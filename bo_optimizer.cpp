#include "bo_optimizer.h"
#include "csv_analyzer.h"
#include "loop_complex.h"
#include "eval_cache.h"
#include "parameter.h"

#include <bayesopt/bayesopt.hpp>        // C++ interface
#include <bayesopt/parameters.hpp>

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Statics
//

BOOptimizer *baysian_opt = nullptr;         // Instance (really just an empty shell at this point)

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// JoSimBO -- ContinuousModel subclass wrapping loop_complex.run_once().
//
// BayesOpt proposes points in [0,1]^n.  Each coordinate is mapped through
// the corresponding parameter's set_normalized_value() before the oracle
// is called.  The EvalCache deduplicates repeated proposals.
//

class JoSimBO : public bayesopt::ContinuousModel {
public:
    JoSimBO(bayesopt::Parameters        params,
            vector<const_parameter*>   &opt_params,
            EvalCache                  &cache,
            FILE                       *sum_fp,
            parameter                  *eval_ptr)
        : ContinuousModel((size_t) opt_params.size(), params)
        , n_evals(0)
        , best_score(DBL_MAX)
        , opt_params(opt_params)
        , cache(cache)
        , sum_fp(sum_fp)
        , eval_ptr(eval_ptr)
        {}

    double evaluateSample(const vectord &query) override
    {
        assert(query.size() == opt_params.size());

        // Build plain double array for the cache (vectord is a uBLAS type)
        vector<double> p(query.size());
        for (size_t i = 0; i < query.size(); i++)
            p[i] = query[i];

        // --- EvalCache lookup ---
        double cached = cache.lookup(p.data());
        if (isfinite(cached))
            return cached;

        // --- Set parameters and call oracle ---
        for (size_t i = 0; i < opt_params.size(); i++)
            opt_params[i]->set_mapped_value(query[i]);

        loop_complex.run_once(sum_fp);
        double score = eval_ptr->get_cur_value();

        // Replace NaN/Inf with DBL_MAX: treats infeasible regions as
        // maximally bad.  BayesOpt minimises so this steers it away.
        if (!isfinite(score))
            score = DBL_MAX;

        cache.store(p.data(), score);

        n_evals++;
        if (score < best_score) {
            best_score = score;
            best_point_found = query;   // store normalized coords
            fprintf(stderr, "BO iter %u: new best = %.6g\n", n_evals, best_score);
        }

        return score;
    }

    bool checkReachability(const vectord &query) override
    {
        // All points in [0,1]^n are structurally reachable.
        // Once the reject expression is wired up, evaluate it here
        // and return false for infeasible points.
        return true;
    }

    unsigned    n_evals;
    double      best_score;
    vectord best_point_found;

private:
    vector<const_parameter*> &opt_params;
    EvalCache          &cache;
    FILE               *sum_fp;
    parameter          *eval_ptr;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer::configure() -- called from define_bo() pragma handler
//


BOOptimizer::BOOptimizer(parameter *obf, unsigned n_it) :
    configured(true),
    n_iterations(n_it),
    obj_funct(obf),
    sum_fp(nullptr)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer::run() -- called from main()
//

void BOOptimizer::run(FILE *result_fp)
{
    assert(configured);

    vector<const_parameter*> opt_params;
    parameter::bo_export(opt_params);                   // get a vector of pointers to the tunable c_params
    const unsigned n = (unsigned) opt_params.size();

    if (n == 0) {
        fprintf(stderr, "BOOptimizer: no optimizable parameters defined\n");
        return;
    }
    
    parameter *ev_ptr = parameter::find_parameter(BAY_OPT_OBJECTIVE);
    if (ev_ptr == nullptr) {
        fprintf(stderr, "BOOptimizer: objective parameter '%s' not defined\n", BAY_OPT_OBJECTIVE);
        return;
    }

    printf("BOOptimizer: %u parameter(s), %u iterations\n", n, n_iterations);

    // --- BayesOpt parameters ---
    bayesopt::Parameters params;
    params.n_iterations  = n_iterations;
//  params.noise         = 1e-10;   // near-deterministic: JoSIM is deterministic
    params.noise         = 1e-4;    // large addition of noise to preven singularity, but also degrades performance
    params.verbose_level = 0;       // suppress BayesOpt's own logging
    params.random_seed = 42;        // fixed seed for reproducibility
    
    // --- EvalCache: capacity = next power of 2 above 4x the budget ---
    size_t cache_cap = 1;
    while (cache_cap < (size_t)(n_iterations * 4))
        cache_cap <<= 1;

    cache_cap *= 16;  // Try larger cache
    
    EvalCache cache(n, cache_cap);

    // --- Build and run optimizer ---
    JoSimBO optimizer(params, opt_params, cache, sum_fp, ev_ptr);

    // BayesOpt default bounds are [0,1]^n -- matches normalised space exactly.
    // No setBoundingBox() call needed.

    vectord best_point((size_t) n);
    try {
        optimizer.optimize(best_point);
    } catch (const std::exception &e) {
        fprintf(stderr, "BOOptimizer: BayesOpt threw: %s\n", e.what());
        fprintf(stderr, "BOOptimizer: iter=%u best so far=%.6g\n",
                optimizer.n_evals, optimizer.best_score);
        // still apply best found so far
        best_point = optimizer.best_point_found;
    } catch (...) {
        fprintf(stderr, "BOOptimizer: BayesOpt threw unknown exception\n");
    }
    
    // --- Apply best found parameters ---
    for (unsigned i = 0; i < n; i++)
        opt_params[i]->set_mapped_value(best_point[i]);

    printf("BOOptimizer: done.  oracle calls=%u  best=%.6g\n",
           optimizer.n_evals, optimizer.best_score);
    cache.print_stats(stdout);

    parameter::save_result(result_fp);
}

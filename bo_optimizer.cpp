#include "bo_optimizer.h"
#include "csv_analyzer.h"
#include "loop_complex.h"
#include "eval_cache.h"
#include "parameter.h"

#include <bayesopt/bayesopt.hpp>
#include <bayesopt/parameters.hpp>
#include <bayesopt/specialtypes.hpp>
#include <bayesopt/prob_distribution.hpp>

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Global instance pointer
//

BOOptimizer *baysian_opt = nullptr;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// JoSimBO -- ContinuousModel subclass wrapping loop_complex.run_once().
//
// In MODE_OPTIMIZE  : evaluateSample() returns the raw cost score (minimised).
// In MODE_ROBUSTNESS: evaluateSample() returns -score so that BayesOpt's
//                     internal machinery stays valid during init, but the real
//                     boundary search is driven by run_robustness() using
//                     getPrediction() directly -- BayesOpt just manages the GP.
//

class JoSimBO : public bayesopt::ContinuousModel {
public:
    enum Mode { OPTIMIZE, ROBUSTNESS };

    JoSimBO(bayesopt::Parameters      params,
            vector<const_parameter*> &opt_params,
            EvalCache                &cache,
            FILE                     *sum_fp,
            parameter                *eval_ptr,
            Mode                      mode)
        : ContinuousModel((size_t) opt_params.size(), params)
        , n_evals(0)
        , best_score(DBL_MAX)
        , opt_params(opt_params)
        , cache(cache)
        , sum_fp(sum_fp)
        , eval_ptr(eval_ptr)
        , mode(mode)
    {}

    //
    // Override generateInitialPoints so the known-good x* is always the
    // first sample evaluated, regardless of n_init_samples.
    // x_star must be set before initializeOptimization() is called.
    //
    void generateInitialPoints(matrixd &xPoints) override
    {
        if (x_star.size() == mDims) {
            // Use x* as the sole initial point -- the GP is anchored there
            // from the very first iteration.
            xPoints.resize(1, mDims);
            for (size_t i = 0; i < mDims; i++)
                xPoints(0, i) = x_star[i];
        } else {
            // Fallback: let BayesOpt choose initial points normally
            ContinuousModel::generateInitialPoints(xPoints);
        }
    }

    double evaluateSample(const vectord &query) override
    {
        assert(query.size() == opt_params.size());

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

        if (!isfinite(score))
            score = (mode == ROBUSTNESS) ? -DBL_MAX : DBL_MAX;

        // In robustness mode BayesOpt sees negated score so its minimiser
        // is consistent: pass (positive margin) -> negative -> BayesOpt happy.
        // The raw score (margin) is what we cache and report.
        double bo_score = (mode == ROBUSTNESS) ? -score : score;

        cache.store(p.data(), bo_score);
        n_evals++;

        if (bo_score < best_score) {
            best_score    = bo_score;
            best_point_found = query;
            if (mode == OPTIMIZE)
                fprintf(stderr, "BO iter %u: new best = %.6g\n", n_evals, best_score);
            else
                fprintf(stderr, "BO iter %u: margin = %.6g\n", n_evals, score);
        }

        return bo_score;
    }

    bool checkReachability(const vectord &query) override
    {
        // Wire up the reject expression here once that pragma is available.
        return true;
    }

    //
    // Surrogate query helpers -- used by the straddle boundary search.
    // These call getPrediction() which is free (no oracle call).
    //

    double surrogate_mean(const vectord &x)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        double m = pd->getMean();
        return isfinite(m) ? m : 0.0;               // treat degenerate prediction as boundary
    }

    double surrogate_std(const vectord &x)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        double s = pd->getStd();
        return (isfinite(s) && s >= 0.0) ? s : 1.0;  // treat as maximally uncertain
    }

    // Straddle acquisition: β·σ(x) - |μ(x)|  (threshold=0, maximise this)
    // Large value = uncertain AND close to boundary = most informative point.
    double straddle(const vectord &x, double beta = 1.96)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        return beta * pd->getStd() - fabs(pd->getMean());
    }

    // Finite-difference gradient of surrogate mean toward zero crossing.
    // Returns the normalised step direction and the gradient magnitude.
    // Step direction points from x toward the boundary (sign chosen so that
    // moving in this direction reduces |μ(x)|).
    // Costs 2*n GP queries, zero oracle calls.
    vectord boundary_gradient(const vectord &x, double &grad_mag, double h = 1e-3)
    {
        size_t n = x.size();
        vectord grad(n), xp(x), xm(x);
        double mu0 = surrogate_mean(x);

        for (size_t i = 0; i < n; i++) {
            xp[i] = min(1.0, x[i] + h);
            xm[i] = max(0.0, x[i] - h);
            double actual_h = xp[i] - xm[i];   // handles boundary clamping
            grad[i] = (surrogate_mean(xp) - surrogate_mean(xm)) / actual_h;
            xp[i] = x[i];
            xm[i] = x[i];
        }

        // Gradient magnitude
        grad_mag = 0.0;
        for (size_t i = 0; i < n; i++)
            grad_mag += grad[i] * grad[i];
        grad_mag = sqrt(grad_mag);

        // Step direction: move in -sign(μ) * ∇μ direction to reduce |μ|
        // (i.e. toward the zero crossing of the mean)
        vectord direction(n);
        double sign_mu = (mu0 >= 0.0) ? 1.0 : -1.0;
        if (grad_mag > 1e-12) {
            for (size_t i = 0; i < n; i++)
                direction[i] = -sign_mu * grad[i] / grad_mag;
        } else {
            // Flat region: pick a random direction
            for (size_t i = 0; i < n; i++)
                direction[i] = (i == (n_evals % n)) ? 1.0 : 0.0;
        }

        return direction;
    }

    // Public state
    unsigned    n_evals;
    double      best_score;
    vectord     best_point_found;
    vectord     x_star;             // known-good starting point, set before optimize()

private:
    vector<const_parameter*> &opt_params;
    EvalCache                &cache;
    FILE                     *sum_fp;
    parameter                *eval_ptr;
    Mode                      mode;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer
//

BOOptimizer::BOOptimizer(parameter *obf, unsigned n_it, Mode m)
    : configured(true)
    , n_iterations(n_it)
    , mode(m)
    , obj_funct(obf)
    , sum_fp(nullptr)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Shared helpers
//

static bayesopt::Parameters make_bo_params(unsigned n_iter, unsigned n_init, double noise = 1e-10)
{
    bayesopt::Parameters p;
    p.n_iterations  = n_iter;
    p.n_init_samples = n_init;
    p.noise         = noise;
    p.verbose_level = 0;
    p.random_seed   = 42;
    return p;
}

static size_t cache_capacity(unsigned n_iter)
{
    size_t cap = 1;
    while (cap < (size_t)(n_iter * 4))
        cap <<= 1;
    return cap * 16;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// MODE_OPTIMIZE -- minimise the objective function
//

void BOOptimizer::run_optimize(FILE *result_fp,
                               vector<const_parameter*> &opt_params,
                               unsigned n)
{
    auto params = make_bo_params(n_iterations, 50);
    EvalCache cache(n, cache_capacity(n_iterations));

    JoSimBO optimizer(params, opt_params, cache, sum_fp,
                      obj_funct, JoSimBO::OPTIMIZE);

    vectord best_point((size_t) n);
    try {
        optimizer.optimize(best_point);
    } catch (const std::exception &e) {
        fprintf(stderr, "BOOptimizer: BayesOpt threw: %s\n", e.what());
        fprintf(stderr, "BOOptimizer: iter=%u best so far=%.6g\n",
                optimizer.n_evals, optimizer.best_score);
        best_point = optimizer.best_point_found;
    } catch (...) {
        fprintf(stderr, "BOOptimizer: BayesOpt threw unknown exception\n");
        best_point = optimizer.best_point_found;
    }

    for (unsigned i = 0; i < n; i++)
        opt_params[i]->set_mapped_value(best_point[i]);

    printf("BOOptimizer[optimize]: done.  oracle calls=%u  best=%.6g\n",
           optimizer.n_evals, optimizer.best_score);
    cache.print_stats(stdout);
    parameter::save_result(result_fp);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// MODE_ROBUSTNESS -- map the feasibility boundary around x*
//
// Strategy:
//   1. Capture x* from get_mapped_value() -- the known-good SA solution.
//   2. Run BayesOpt with x* as the sole initial sample and n_init_samples=1
//      to seed the GP, then drive exploration using the straddle criterion
//      via the manual stepOptimization() loop.  At each step:
//        a. Estimate the surrogate gradient toward the boundary.
//        b. Take a step of size proportional to |μ| / |∇μ| (large gradient
//           = nearly at boundary = small step; flat = large step).
//        c. Add exploration noise scaled by σ to avoid surrogate artifacts.
//        d. Evaluate the oracle at the chosen point, feed back to GP.
//   3. After the budget is exhausted, walk each parameter axis from x* using
//      the GP mean to locate the zero crossing -- per-parameter margins.
//

void BOOptimizer::run_robustness(FILE *result_fp,
                                 vector<const_parameter*> &opt_params,
                                 unsigned n)
{
    // --- Capture x* before any set_mapped_value() calls modify it ---
    vectord x_star((size_t) n);
    for (unsigned i = 0; i < n; i++)
        x_star[i] = opt_params[i]->get_mapped_value();

    // Verify x* is actually a passing point
    {
        for (unsigned i = 0; i < n; i++)
            opt_params[i]->set_mapped_value(x_star[i]);
        loop_complex.run_once(sum_fp);
        double margin_star = obj_funct->get_cur_value();
        if (margin_star <= 0.0) {
            fprintf(stderr,
                "BOOptimizer[robustness]: WARNING: starting point has margin=%.6g (not passing)\n",
                margin_star);
        } else {
            fprintf(stderr,
                "BOOptimizer[robustness]: x* verified, margin=%.6g\n", margin_star);
        }
    }

    // --- BayesOpt params: 1 init sample (x* itself), rest are straddle-guided ---
    auto params = make_bo_params(n_iterations, 1, 1e-4);
    EvalCache cache(n, cache_capacity(n_iterations));

    JoSimBO optimizer(params, opt_params, cache, sum_fp,
                      obj_funct, JoSimBO::ROBUSTNESS);
    optimizer.x_star = x_star;      // generateInitialPoints() will use this

    // --- Seed the GP with x* then run straddle-guided steps ---
    optimizer.initializeOptimization();

    vectord x_cur = x_star;         // current probe point
    const double beta      = 1.96;  // straddle exploration weight
    const double max_step  = 0.15;  // max step in normalised space per iteration
    const double noise_scale = 0.02;// exploration noise amplitude

    for (unsigned iter = 0; iter < n_iterations; iter++) {
        // Estimate gradient of GP mean toward boundary from current point
        double grad_mag;
        vectord direction = optimizer.boundary_gradient(x_cur, grad_mag);

        // Step size: proportional to |μ| / |∇μ|, capped at max_step
        double mu_cur = fabs(optimizer.surrogate_mean(x_cur));
        double step   = (grad_mag > 1e-12)
                        ? min(mu_cur / grad_mag, max_step)
                        : max_step;

        // Straddle-weighted candidate: gradient step + exploration noise
        vectord x_next((size_t) n);
        for (size_t i = 0; i < n; i++) {
            // Small random perturbation scaled by local GP uncertainty
            double sigma_i = optimizer.surrogate_std(x_cur);
            double noise   = noise_scale * sigma_i * (2.0 * ((rand() / (double)RAND_MAX) - 0.5));
            x_next[i] = x_cur[i] + step * direction[i] + noise;
            x_next[i] = fmax(0.0, fmin(1.0, x_next[i]));  // clip to [0,1]
        }

        // Evaluate straddle at candidate vs current -- accept if higher
        double s_next = optimizer.straddle(x_next, beta);
        double s_cur  = optimizer.straddle(x_cur,  beta);
        if (s_next > s_cur)
            x_cur = x_next;

        // Let BayesOpt do one step from the chosen point
        // (feeds the oracle result back into the GP model)
        optimizer.stepOptimization();

        fprintf(stderr, "robustness iter %u: step=%.4f grad=%.4f mu=%.4f\n",
                iter, step, grad_mag, optimizer.surrogate_mean(x_cur));
    }

    // --- Per-parameter margin extraction via GP mean zero-crossing ---
    printf("\nBOOptimizer[robustness]: margin analysis from x*\n");
    printf("%-20s  %8s  %8s  %8s  %8s\n",
           "parameter", "value", "lo_margin", "hi_margin", "margin_%");

    const unsigned N_BISECT = 32;   // bisection steps per axis
    const double   STEP0    = 0.01; // initial step for bracket search

    for (unsigned i = 0; i < n; i++) {
        double xi_star = x_star[i];

        // Walk upward from x* to find upper zero crossing
        vectord x_probe = x_star;
        double hi_cross = xi_star;
        for (double t = xi_star + STEP0; t <= 1.0; t += STEP0) {
            x_probe[i] = t;
            if (optimizer.surrogate_mean(x_probe) <= 0.0) {
                // Found bracket [t-STEP0, t]: bisect
                double lo_b = t - STEP0, hi_b = t;
                for (unsigned k = 0; k < N_BISECT; k++) {
                    double mid = (lo_b + hi_b) / 2.0;
                    x_probe[i] = mid;
                    if (optimizer.surrogate_mean(x_probe) > 0.0)
                        lo_b = mid;
                    else
                        hi_b = mid;
                }
                hi_cross = (lo_b + hi_b) / 2.0;
                break;
            }
            hi_cross = t;           // still passing, extend
        }
        x_probe[i] = x_star[i];

        // Walk downward from x* to find lower zero crossing
        double lo_cross = xi_star;
        for (double t = xi_star - STEP0; t >= 0.0; t -= STEP0) {
            x_probe[i] = t;
            if (optimizer.surrogate_mean(x_probe) <= 0.0) {
                double lo_b = t, hi_b = t + STEP0;
                for (unsigned k = 0; k < N_BISECT; k++) {
                    double mid = (lo_b + hi_b) / 2.0;
                    x_probe[i] = mid;
                    if (optimizer.surrogate_mean(x_probe) > 0.0)
                        hi_b = mid;
                    else
                        lo_b = mid;
                }
                lo_cross = (lo_b + hi_b) / 2.0;
                break;
            }
            lo_cross = t;           // still passing, extend
        }
        x_probe[i] = x_star[i];

        // Report in normalised space; physical conversion via set/get
        double lo_margin_norm = xi_star - lo_cross;
        double hi_margin_norm = hi_cross - xi_star;
        double margin_pct     = 100.0 * (hi_cross - lo_cross) / 2.0;

        printf("%-20s  %8.4f  %8.4f  %8.4f  %7.2f%%\n",
               opt_params[i]->get_name(),
               xi_star, lo_margin_norm, hi_margin_norm, margin_pct);
    }

    printf("\nBOOptimizer[robustness]: oracle calls=%u\n", optimizer.n_evals);
    cache.print_stats(stdout);
    parameter::save_result(result_fp);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer::run() -- dispatch by mode
//

void BOOptimizer::run(FILE *result_fp)
{
    assert(configured);

    vector<const_parameter*> opt_params;
    parameter::bo_export(opt_params);
    const unsigned n = (unsigned) opt_params.size();

    if (n == 0) {
        fprintf(stderr, "BOOptimizer: no optimizable parameters defined\n");
        return;
    }

    if (obj_funct == nullptr) {
        fprintf(stderr, "BOOptimizer: no objective parameter defined\n");
        return;
    }

    printf("BOOptimizer: %u parameter(s), %u iterations, mode=%s\n",
           n, n_iterations,
           (mode == MODE_OPTIMIZE) ? "optimize" : "robustness");

    if (mode == MODE_OPTIMIZE)
        run_optimize(result_fp, opt_params, n);
    else
        run_robustness(result_fp, opt_params, n);
}

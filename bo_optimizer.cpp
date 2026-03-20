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
// In OPTIMIZE   mode: evaluateSample() returns the raw cost score (minimised).
// In ROBUSTNESS mode: evaluateSample() returns -score so that BayesOpt's
//                     internal machinery stays valid during init, but the real
//                     boundary search is driven externally using getPrediction()
//                     directly -- BayesOpt just manages the GP surrogate.
//
// The binary submode (run_binary) uses this class only for the GP and the
// direct oracle call helper (call_oracle), not for BayesOpt's acquisition
// function -- the ray-bisection search drives all evaluation directly.
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
            best_score       = bo_score;
            best_point_found = query;
            if (mode == OPTIMIZE)
                fprintf(stderr, "BO iter %u: new best = %.6g\n", n_evals, best_score);
            else
                fprintf(stderr, "BO iter %u: margin = %.6g\n", n_evals, score);
        }

        return bo_score;
    }

    //
    // call_oracle() -- direct oracle call bypassing BayesOpt's machinery.
    // Used by run_binary() to evaluate arbitrary points during ray-bisection
    // without going through evaluateSample().  Caches the result.
    // Returns the raw score (positive=pass, negative=fail).
    //
    double call_oracle(const vectord &query)
    {
        assert(query.size() == opt_params.size());

        vector<double> p(query.size());
        for (size_t i = 0; i < query.size(); i++)
            p[i] = query[i];

        // Check cache first (avoids re-running JoSIM at already-visited points)
        double cached = cache.lookup(p.data());
        if (isfinite(cached))
            return cached;     // cache stores raw score in binary mode

        for (size_t i = 0; i < opt_params.size(); i++)
            opt_params[i]->set_mapped_value(query[i]);

        loop_complex.run_once(sum_fp);
        double score = eval_ptr->get_cur_value();

        if (!isfinite(score))
            score = -DBL_MAX;   // treat non-finite as deep fail

        cache.store(p.data(), score);
        n_evals++;

        return score;
    }

    bool checkReachability(const vectord &query) override
    {
        // Wire up the reject expression here once that pragma is available.
        return true;
    }

    //
    // Surrogate query helpers -- used by the gradient straddle search.
    // These call getPrediction() which is free (no oracle call).
    //

    double surrogate_mean(const vectord &x)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        double m = pd->getMean();
        return isfinite(m) ? m : 0.0;              // treat degenerate prediction as boundary
    }

    double surrogate_std(const vectord &x)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        double s = pd->getStd();
        return (isfinite(s) && s >= 0.0) ? s : 1.0; // treat as maximally uncertain
    }

    // Straddle acquisition: beta*sigma(x) - |mu(x)|  (threshold=0, maximise this).
    // Large value = uncertain AND close to boundary = most informative point.
    double straddle(const vectord &x, double beta = 1.96)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        return beta * pd->getStd() - fabs(pd->getMean());
    }

    // Finite-difference gradient of surrogate mean toward zero crossing.
    // Returns the normalised step direction and the gradient magnitude.
    // Step direction moves from x toward the boundary (sign chosen to reduce |mu(x)|).
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

        grad_mag = 0.0;
        for (size_t i = 0; i < n; i++)
            grad_mag += grad[i] * grad[i];
        grad_mag = sqrt(grad_mag);

        // Move in -sign(mu)*grad direction to reduce |mu| (toward zero crossing)
        vectord direction(n);
        double sign_mu = (mu0 >= 0.0) ? 1.0 : -1.0;
        if (grad_mag > 1e-12) {
            for (size_t i = 0; i < n; i++)
                direction[i] = -sign_mu * grad[i] / grad_mag;
        } else {
            // Flat region: pick an axis-aligned direction to avoid stalling
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
// BOOptimizer constructor
//

BOOptimizer::BOOptimizer(parameter *obf,
                         unsigned   n_it,
                         Mode       m,
                         SubMode    sm,
                         unsigned   n_r,
                         unsigned   n_brk,
                         unsigned   n_bis,
                         double     thr)
    : configured(true)
    , n_iterations(n_it)
    , mode(m)
    , submode(sm)
    , obj_funct(obf)
    , sum_fp(nullptr)
    , n_rays(n_r)
    , n_bracket(n_brk)
    , n_bisect(n_bis)
    , threshold(thr)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Shared helpers
//

static bayesopt::Parameters make_bo_params(unsigned n_iter, unsigned n_init,
                                           double noise = 1e-10)
{
    bayesopt::Parameters p;
    p.n_iterations   = n_iter;
    p.n_init_samples = n_init;
    p.noise          = noise;
    p.verbose_level  = 0;
    p.random_seed    = 42;
    return p;
}

static size_t cache_capacity(unsigned n_iter)
{
    size_t cap = 1;
    while (cap < (size_t)(n_iter * 4))
        cap <<= 1;
    return cap * 16;
}

// Clip a normalised coordinate to [0,1]
static inline double clip01(double x) { return fmax(0.0, fmin(1.0, x)); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Per-parameter margin extraction -- shared by both robustness submodes.
//
// Walks each parameter axis from x* along the GP surrogate mean to find
// where it crosses threshold.  Cheap (GP queries only, no oracle calls).
// Reports lo_margin and hi_margin in normalised space plus a percentage.
//

static void report_margins(JoSimBO &optimizer,
                            const vectord &x_star,
                            vector<const_parameter*> &opt_params,
                            unsigned n,
                            double threshold,
                            double step0 = 0.01,
                            unsigned n_bisect = 32)
{
    printf("\nBOOptimizer[robustness]: margin analysis from x*\n");
    printf("%-20s  %8s  %8s  %8s  %8s\n",
           "parameter", "value", "lo_margin", "hi_margin", "margin_%");

    for (unsigned i = 0; i < n; i++) {
        double xi_star = x_star[i];
        vectord x_probe = x_star;

        // --- Walk upward to find upper boundary bracket ---
        double hi_cross = 1.0;          // assume passes to the upper bound unless proven otherwise
        bool   hi_found = false;
        for (double t = xi_star + step0; t <= 1.0; t += step0) {
            x_probe[i] = t;
            if (optimizer.surrogate_mean(x_probe) <= threshold) {
                // Bracket found: bisect between t-step0 and t
                double lo_b = t - step0, hi_b = t;
                for (unsigned k = 0; k < n_bisect; k++) {
                    double mid = (lo_b + hi_b) / 2.0;
                    x_probe[i] = mid;
                    if (optimizer.surrogate_mean(x_probe) > threshold)
                        lo_b = mid;
                    else
                        hi_b = mid;
                }
                hi_cross = (lo_b + hi_b) / 2.0;
                hi_found = true;
                break;
            }
        }
        x_probe[i] = xi_star;

        // --- Walk downward to find lower boundary bracket ---
        double lo_cross = 0.0;          // assume passes to the lower bound unless proven otherwise
        bool   lo_found = false;
        for (double t = xi_star - step0; t >= 0.0; t -= step0) {
            x_probe[i] = t;
            if (optimizer.surrogate_mean(x_probe) <= threshold) {
                double lo_b = t, hi_b = t + step0;
                for (unsigned k = 0; k < n_bisect; k++) {
                    double mid = (lo_b + hi_b) / 2.0;
                    x_probe[i] = mid;
                    if (optimizer.surrogate_mean(x_probe) > threshold)
                        hi_b = mid;
                    else
                        lo_b = mid;
                }
                lo_cross = (lo_b + hi_b) / 2.0;
                lo_found = true;
                break;
            }
        }
        x_probe[i] = xi_star;

        double lo_margin_norm = xi_star - lo_cross;
        double hi_margin_norm = hi_cross - xi_star;
        double margin_pct     = 100.0 * (hi_cross - lo_cross) / 2.0;

        printf("%-20s  %8.4f  %8.4f  %8.4f  %7.2f%%%s%s\n",
               opt_params[i]->get_name(),
               xi_star, lo_margin_norm, hi_margin_norm, margin_pct,
               lo_found ? "" : "  (lo bound not found)",
               hi_found ? "" : "  (hi bound not found)");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// MODE_OPTIMIZE -- minimise the objective function via BayesOpt
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
// SUBMODE_BINARY -- ray-bisection boundary search from x*
//
// The oracle returns a signed value with threshold at 0: positive=pass, negative=fail.
// For a digital circuit this is typically +1/-1 exactly.
//
// Strategy:
//   For each of n_rays directions from x* (axis-aligned plus random diagonals):
//     1. Step outward along the ray until the oracle returns fail (bracket found).
//     2. Bisect within the bracket to locate the zero crossing precisely.
//     3. Record the boundary point in the GP via evaluateSample().
//   After all rays, use the GP mean zero-crossing walk for per-parameter margins.
//
// Budget:
//   n_rays * (n_bracket + n_bisect) oracle calls total.
//   If this exceeds n_iterations, n_rays is reduced to fit.
//   Auto n_rays (n_rays==0): defaults to 4*n, covering all axis-aligned directions
//   plus 2*n random diagonal directions.
//

void BOOptimizer::run_binary(FILE *result_fp,
                             vector<const_parameter*> &opt_params,
                             unsigned n)
{
    // --- Capture x* ---
    vectord x_star((size_t) n);
    for (unsigned i = 0; i < n; i++)
        x_star[i] = opt_params[i]->get_mapped_value();

    // Verify x* passes
    {
        for (unsigned i = 0; i < n; i++)
            opt_params[i]->set_mapped_value(x_star[i]);
        loop_complex.run_once(sum_fp);
        double margin_star = obj_funct->get_cur_value();
        if (margin_star <= threshold) {
            fprintf(stderr,
                "BOOptimizer[binary]: WARNING: starting point score=%.6g does not pass threshold=%.6g\n",
                margin_star, threshold);
        } else {
            fprintf(stderr,
                "BOOptimizer[binary]: x* verified, score=%.6g (passes threshold=%.6g)\n",
                margin_star, threshold);
        }
    }

    // --- Build ray directions ---
    // First 2*n rays: axis-aligned +/- for each parameter
    // Remaining rays: random unit-vector diagonals (seeded for reproducibility)
    unsigned actual_n_rays = (n_rays == 0) ? 4 * n : n_rays;

    // Check budget: reduce rays if needed to fit n_iterations
    unsigned calls_per_ray = n_bracket + n_bisect;
    if (calls_per_ray == 0) calls_per_ray = 1;
    if (actual_n_rays * calls_per_ray > n_iterations) {
        actual_n_rays = n_iterations / calls_per_ray;
        if (actual_n_rays == 0) actual_n_rays = 1;
        fprintf(stderr,
            "BOOptimizer[binary]: budget reduced to %u rays (%u calls each, %u total)\n",
            actual_n_rays, calls_per_ray, actual_n_rays * calls_per_ray);
    }

    // Build direction vectors
    vector<vectord> directions;
    directions.reserve(actual_n_rays);

    // Axis-aligned +/- directions
    for (unsigned i = 0; i < n && directions.size() < actual_n_rays; i++) {
        vectord d((size_t) n, 0.0);
        d[i] = 1.0;
        directions.push_back(d);
    }
    for (unsigned i = 0; i < n && directions.size() < actual_n_rays; i++) {
        vectord d((size_t) n, 0.0);
        d[i] = -1.0;
        directions.push_back(d);
    }
    // Fill remaining slots with random unit-vector diagonals
    srand(42);  // fixed seed for reproducibility
    while (directions.size() < actual_n_rays) {
        vectord d((size_t) n);
        double len = 0.0;
        for (size_t i = 0; i < n; i++) {
            d[i] = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
            len += d[i] * d[i];
        }
        len = sqrt(len);
        if (len < 1e-12) continue;  // degenerate sample, retry
        for (size_t i = 0; i < n; i++)
            d[i] /= len;
        directions.push_back(d);
    }

    // --- Seed GP with x* as the one known-good point ---
    auto params = make_bo_params(n_iterations, 1, 1e-4);
    EvalCache cache(n, cache_capacity(n_iterations));

    JoSimBO optimizer(params, opt_params, cache, sum_fp,
                      obj_funct, JoSimBO::ROBUSTNESS);
    optimizer.x_star = x_star;
    optimizer.initializeOptimization();

    // --- Ray-bisection boundary search ---
    const double step0    = 1.0 / (double)(n_bracket + 1); // initial step size along ray
    unsigned     boundary_pts_found = 0;

    for (unsigned ray = 0; ray < actual_n_rays; ray++) {
        const vectord &dir = directions[ray];

        // Step outward from x* along this ray to find a bracket [t_pass, t_fail]
        double t_pass = 0.0;    // last known passing distance along ray
        double t_fail = -1.0;   // first failing distance (-1 = not found yet)

        for (unsigned step = 1; step <= n_bracket; step++) {
            double t = step * step0;

            // Build probe point: x* + t * dir, clipped to [0,1]^n
            vectord x_probe((size_t) n);
            bool    at_boundary = false;    // true if clipping hit [0,1] wall
            for (size_t i = 0; i < n; i++) {
                x_probe[i] = clip01(x_star[i] + t * dir[i]);
                if (x_probe[i] != x_star[i] + t * dir[i])
                    at_boundary = true;
            }

            double score = optimizer.call_oracle(x_probe);
            fprintf(stderr, "ray %u step %u: t=%.4f score=%.4g\n",
                    ray, step, t, score);

            if (score <= threshold) {
                t_fail = t;
                break;
            }
            t_pass = t;

            if (at_boundary) {
                // Hit the [0,1] wall while still passing -- no boundary in this direction
                fprintf(stderr, "ray %u: reached parameter boundary still passing\n", ray);
                break;
            }
        }

        if (t_fail < 0.0) {
            // No failure found along this ray -- boundary is beyond the parameter range
            fprintf(stderr, "ray %u: no boundary found within parameter range\n", ray);
            continue;
        }

        // --- Bisect between t_pass and t_fail to locate the boundary precisely ---
        for (unsigned k = 0; k < n_bisect; k++) {
            double t_mid = (t_pass + t_fail) / 2.0;
            vectord x_probe((size_t) n);
            for (size_t i = 0; i < n; i++)
                x_probe[i] = clip01(x_star[i] + t_mid * dir[i]);

            double score = optimizer.call_oracle(x_probe);
            if (score > threshold)
                t_pass = t_mid;
            else
                t_fail = t_mid;
        }

        // Boundary point is the last failing probe (just inside the fail region)
        vectord x_boundary((size_t) n);
        for (size_t i = 0; i < n; i++)
            x_boundary[i] = clip01(x_star[i] + t_fail * dir[i]);

        fprintf(stderr, "ray %u: boundary at t=%.6f (%.0f oracle calls)\n",
                ray, t_fail, (double) optimizer.n_evals);
        boundary_pts_found++;

        // Feed the boundary point into the GP via stepOptimization so the
        // surrogate learns the boundary shape for the margin extraction step.
        // (The oracle call is cached so this costs nothing extra.)
        optimizer.stepOptimization();
    }

    fprintf(stderr, "BOOptimizer[binary]: %u boundary points found, %u oracle calls\n",
            boundary_pts_found, optimizer.n_evals);

    // --- Extract per-parameter margins from GP mean zero-crossing ---
    report_margins(optimizer, x_star, opt_params, n, threshold,
                   0.01, n_bisect);

    printf("BOOptimizer[binary]: oracle calls=%u\n", optimizer.n_evals);
    cache.print_stats(stdout);
    parameter::save_result(result_fp);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// SUBMODE_GRADIENT -- surrogate-gradient straddle search from x*
//
// Used when the oracle returns a continuous signed margin with a meaningful
// gradient (e.g. the AC synchronizer gate where slope is the figure of merit).
// The GP surrogate gradient guides the search toward the boundary.
//
// This is the original run_robustness() approach, retained here for circuits
// where the gradient is informative.
//

void BOOptimizer::run_gradient(FILE *result_fp,
                               vector<const_parameter*> &opt_params,
                               unsigned n)
{
    // --- Capture x* before any set_mapped_value() calls modify it ---
    vectord x_star((size_t) n);
    for (unsigned i = 0; i < n; i++)
        x_star[i] = opt_params[i]->get_mapped_value();

    // Verify x* passes
    {
        for (unsigned i = 0; i < n; i++)
            opt_params[i]->set_mapped_value(x_star[i]);
        loop_complex.run_once(sum_fp);
        double margin_star = obj_funct->get_cur_value();
        if (margin_star <= threshold) {
            fprintf(stderr,
                "BOOptimizer[gradient]: WARNING: starting point score=%.6g does not pass threshold=%.6g\n",
                margin_star, threshold);
        } else {
            fprintf(stderr,
                "BOOptimizer[gradient]: x* verified, score=%.6g\n", margin_star);
        }
    }

    // Use higher noise for gradient mode: boundary probing clusters evaluations
    // in a small region which stresses the GP kernel matrix at low noise.
    auto params = make_bo_params(n_iterations, 1, 1e-4);
    EvalCache cache(n, cache_capacity(n_iterations));

    JoSimBO optimizer(params, opt_params, cache, sum_fp,
                      obj_funct, JoSimBO::ROBUSTNESS);
    optimizer.x_star = x_star;

    // Seed the GP with x* then run straddle-guided steps
    optimizer.initializeOptimization();

    vectord x_cur = x_star;
    const double beta        = 1.96;  // straddle exploration weight
    const double max_step    = 0.05;  // max step in normalised space (tighter than before)
    const double noise_scale = 0.02;  // exploration noise amplitude

    for (unsigned iter = 0; iter < n_iterations; iter++) {
        double grad_mag;
        vectord direction = optimizer.boundary_gradient(x_cur, grad_mag);

        // Step size: proportional to |mu| / |grad|, capped at max_step
        double mu_cur = fabs(optimizer.surrogate_mean(x_cur) - threshold);
        double step   = (grad_mag > 1e-12)
                        ? min(mu_cur / grad_mag, max_step)
                        : max_step;

        // Gradient step + small exploration noise
        vectord x_next((size_t) n);
        for (size_t i = 0; i < n; i++) {
            double sigma_i = optimizer.surrogate_std(x_cur);
            double noise   = noise_scale * sigma_i
                             * (2.0 * ((rand() / (double)RAND_MAX) - 0.5));
            x_next[i] = clip01(x_cur[i] + step * direction[i] + noise);
        }

        // Accept x_next if its straddle value is higher than x_cur
        if (optimizer.straddle(x_next, beta) > optimizer.straddle(x_cur, beta))
            x_cur = x_next;

        optimizer.stepOptimization();

        fprintf(stderr, "gradient iter %u: step=%.4f grad=%.4f mu=%.4f\n",
                iter, step, grad_mag, optimizer.surrogate_mean(x_cur));
    }

    report_margins(optimizer, x_star, opt_params, n, threshold,
                   0.01, n_bisect);

    printf("BOOptimizer[gradient]: oracle calls=%u\n", optimizer.n_evals);
    cache.print_stats(stdout);
    parameter::save_result(result_fp);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// run_robustness() -- dispatch to binary or gradient submode
//

void BOOptimizer::run_robustness(FILE *result_fp,
                                 vector<const_parameter*> &opt_params,
                                 unsigned n)
{
    switch (submode) {
        case SUBMODE_BINARY:
            run_binary(result_fp, opt_params, n);
            break;
        case SUBMODE_GRADIENT:
            run_gradient(result_fp, opt_params, n);
            break;
        case SUBMODE_PROBABILISTIC:
            fprintf(stderr, "BOOptimizer: probabilistic submode not yet implemented\n");
            break;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer::run() -- top-level dispatch
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

    const char *mode_str = "optimize";
    if (mode == MODE_ROBUSTNESS) {
        switch (submode) {
            case SUBMODE_BINARY:        mode_str = "robustness/binary";        break;
            case SUBMODE_GRADIENT:      mode_str = "robustness/gradient";      break;
            case SUBMODE_PROBABILISTIC: mode_str = "robustness/probabilistic"; break;
        }
    }

    printf("BOOptimizer: %u parameter(s), %u iterations, mode=%s\n",
           n, n_iterations, mode_str);

    if (mode == MODE_OPTIMIZE)
        run_optimize(result_fp, opt_params, n);
    else
        run_robustness(result_fp, opt_params, n);
}

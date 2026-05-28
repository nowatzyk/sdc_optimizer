#include "bo_optimizer.h"
#include "csv_analyzer.h"
#include "loop_complex.h"
#include "eval_cache.h"
#include "parameter.h"
#include "binary_ellipsoid_fit.h"

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
// Used only by MODE_OPTIMIZE and SUBMODE_GRADIENT.
// SUBMODE_BINARY drives oracle calls directly via eval_point() lambda and
// does not use the GP surrogate at all.
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
            xPoints.resize(1, mDims);
            for (size_t i = 0; i < mDims; i++)
                xPoints(0, i) = x_star[i];
        } else {
            ContinuousModel::generateInitialPoints(xPoints);
        }
    }

    double evaluateSample(const vectord &query) override
    {
        assert(query.size() == opt_params.size());

        vector<double> p(query.size());
        for (size_t i = 0; i < query.size(); i++)
            p[i] = query[i];

        double cached = cache.lookup(p.data());
        if (isfinite(cached))
            return cached;

        for (size_t i = 0; i < opt_params.size(); i++)
            opt_params[i]->set_mapped_value(query[i]);

        loop_complex.run_once(sum_fp);
        double score = eval_ptr->get_cur_value();

        if (!isfinite(score))
            score = (mode == ROBUSTNESS) ? -DBL_MAX : DBL_MAX;

        // In robustness mode BayesOpt sees negated score so its minimiser
        // is consistent: pass (positive margin) -> negative -> BayesOpt happy.
        double bo_score = (mode == ROBUSTNESS) ? -score : score;

        cache.store(p.data(), bo_score);
        n_evals++;

        if (bo_score < best_score) {
            best_score       = bo_score;
            best_point_found = query;
            if (mode == OPTIMIZE)
                fprintf(stderr, "BO iter %u: new best = %.6g\n",
                        n_evals, best_score);
            else
                fprintf(stderr, "BO iter %u: margin = %.6g\n",
                        n_evals, score);
        }

        return bo_score;
    }

    bool checkReachability(const vectord &query) override
    {
        // Wire up the reject expression here once that pragma is available.
        return true;
    }

    //
    // Surrogate query helpers -- used by the gradient straddle search.
    //

    double surrogate_mean(const vectord &x)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        double m = pd->getMean();
        return isfinite(m) ? m : 0.0;
    }

    double surrogate_std(const vectord &x)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        double s = pd->getStd();
        return (isfinite(s) && s >= 0.0) ? s : 1.0;
    }

    double straddle(const vectord &x, double beta = 1.96)
    {
        bayesopt::ProbabilityDistribution *pd = getPrediction(x);
        return beta * pd->getStd() - fabs(pd->getMean());
    }

    vectord boundary_gradient(const vectord &x, double &grad_mag,
                              double h = 1e-3)
    {
        size_t n = x.size();
        vectord grad(n), xp(x), xm(x);
        double mu0 = surrogate_mean(x);

        for (size_t i = 0; i < n; i++) {
            xp[i] = min(1.0, x[i] + h);
            xm[i] = max(0.0, x[i] - h);
            double actual_h = xp[i] - xm[i];
            grad[i] = (surrogate_mean(xp) - surrogate_mean(xm)) / actual_h;
            xp[i] = x[i];
            xm[i] = x[i];
        }

        grad_mag = 0.0;
        for (size_t i = 0; i < n; i++)
            grad_mag += grad[i] * grad[i];
        grad_mag = sqrt(grad_mag);

        vectord direction(n);
        double sign_mu = (mu0 >= 0.0) ? 1.0 : -1.0;
        if (grad_mag > 1e-12) {
            for (size_t i = 0; i < n; i++)
                direction[i] = -sign_mu * grad[i] / grad_mag;
        } else {
            for (size_t i = 0; i < n; i++)
                direction[i] = (i == (n_evals % n)) ? 1.0 : 0.0;
        }

        return direction;
    }

    unsigned    n_evals;
    double      best_score;
    vectord     best_point_found;
    vectord     x_star;

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
                         double     thr,
                         unsigned   bef_bgt,
                         unsigned   bef_itr,
                         unsigned   bef_prb)
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
    , bef_budget(bef_bgt)
    , bef_iter(bef_itr)
    , bef_probes(bef_prb)
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

static inline double clip01(double x) { return fmax(0.0, fmin(1.0, x)); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// GP-based per-parameter margin extraction -- used by SUBMODE_GRADIENT only.
//

static void report_margins_gp(JoSimBO &optimizer,
                               const vectord &x_star,
                               vector<const_parameter*> &opt_params,
                               unsigned n,
                               double threshold,
                               double step0    = 0.01,
                               unsigned n_bisect = 32)
{
    printf("\nBOOptimizer[gradient]: margin analysis from x* (GP surrogate)\n");
    printf("%-20s  %8s  %8s  %8s  %8s\n",
           "parameter", "norm_val", "lo_margin", "hi_margin", "margin_%");

    for (unsigned i = 0; i < n; i++) {
        double xi_star = x_star[i];
        vectord x_probe = x_star;

        double hi_cross = 1.0;
        bool   hi_found = false;
        for (double t = xi_star + step0; t <= 1.0; t += step0) {
            x_probe[i] = t;
            if (optimizer.surrogate_mean(x_probe) <= threshold) {
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

        double lo_cross = 0.0;
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

        double lo_margin = xi_star - lo_cross;
        double hi_margin = hi_cross - xi_star;
        double margin_pct = 100.0 * (lo_margin + hi_margin) / 2.0;

        printf("%-20s  %8.4f  %8.4f  %8.4f  %7.2f%%%s%s\n",
               opt_params[i]->get_name(),
               xi_star, lo_margin, hi_margin, margin_pct,
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
// SUBMODE_BINARY -- iterative ray-bisection boundary search with centre relocation
//
// Runs n_rounds = n_relocate+1 rounds of ray-bisection from an evolving centre.
// After each round the centroid of all accumulated boundary points becomes
// the new centre, progressively relocating toward the geometric centre of
// the feasibility region.
//
// Budget notes:
//   The total oracle call budget is n_iterations.  The cache is shared across
//   all rounds, so repeated bisection near the same boundary is cheap.
//   actual_n_rays is set to fit the TOTAL budget (not per-round), since later
//   rounds reuse cached evaluations heavily.
//   n_relocate is currently hardcoded to 3 -- will become a pragma parameter.
//
// Bug fixes vs previous version:
//   - Budget is total not per-round (no longer divides by n_rounds)
//   - t_pass=0 is a valid bracket start (first step fails = boundary is
//     between 0 and step0); bisection now always runs when t_fail >= 0
//   - Boundary point is recorded at t_pass (last confirmed pass), never
//     at x_centre itself unless t_pass genuinely equals 0
//

void BOOptimizer::run_binary(FILE *result_fp,
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
        double score_star = obj_funct->get_cur_value();
        if (score_star <= threshold) {
            fprintf(stderr,
                "BOOptimizer[binary]: WARNING: starting point score=%.6g "
                "does not pass threshold=%.6g\n", score_star, threshold);
        } else {
            fprintf(stderr,
                "BOOptimizer[binary]: x* verified, score=%.6g "
                "(passes threshold=%.6g)\n", score_star, threshold);
        }
    }

    // --- Number of rays and rounds ---
    // n_rays==0 means auto: use 4*n rays per round.
    // Budget is the TOTAL across all rounds -- the cache absorbs repeated
    // bisection near the same boundary, so later rounds cost much less.
    const unsigned n_relocate = 3;          // TBD: make this a pragma parameter
    const unsigned n_rounds   = n_relocate + 1;
    unsigned actual_n_rays    = (n_rays == 0) ? 4 * n : n_rays;

    // Sanity check: warn if total budget seems tight
    unsigned calls_per_ray = n_bracket + n_bisect;
    if (calls_per_ray == 0) calls_per_ray = 1;
    unsigned est_fresh_calls = actual_n_rays * calls_per_ray;  // round 1 only
    if (est_fresh_calls > n_iterations) {
        // Reduce rays so a single round fits in budget
        actual_n_rays = n_iterations / calls_per_ray;
        if (actual_n_rays == 0) actual_n_rays = 1;
        fprintf(stderr,
            "BOOptimizer[binary]: reduced to %u rays to fit single-round "
            "budget of %u calls\n", actual_n_rays, n_iterations);
    }

    fprintf(stderr,
        "BOOptimizer[binary]: %u rounds, %u rays/round, "
        "%u calls/ray (budget %u, later rounds mostly cached)\n",
        n_rounds, actual_n_rays, calls_per_ray, n_iterations);

    if (actual_n_rays < 2 * n)
        fprintf(stderr,
            "BOOptimizer[binary]: WARNING: %u rays < 2*n=%u, "
            "some axis margins will be missing\n", actual_n_rays, 2 * n);

    // --- Direct oracle evaluator shared across all rounds ---
    unsigned n_evals = 0;
    EvalCache cache(n, cache_capacity(n_iterations));

    auto eval_point = [&](const vectord &x) -> double {
        vector<double> p(n);
        for (size_t i = 0; i < n; i++) p[i] = x[i];

        double cached = cache.lookup(p.data());
        if (isfinite(cached)) return cached;

        for (size_t i = 0; i < n; i++)
            opt_params[i]->set_mapped_value(x[i]);
        loop_complex.run_once(sum_fp);
        double score = obj_funct->get_cur_value();
        if (!isfinite(score)) score = -DBL_MAX;  // NaN -> deep fail

        cache.store(p.data(), score);
        n_evals++;
        return score;
    };

    // --- Build ray direction vectors ---
    // +axis (0..n-1), -axis (n..2n-1), random diagonals beyond.
    // Same set reused every round -- direction diversity comes from
    // the changing centre, not from changing directions.
    vector<vectord> directions;
    directions.reserve(actual_n_rays);

    for (unsigned i = 0; i < n && directions.size() < actual_n_rays; i++) {
        vectord d((size_t) n, 0.0); d[i] =  1.0; directions.push_back(d);
    }
    for (unsigned i = 0; i < n && directions.size() < actual_n_rays; i++) {
        vectord d((size_t) n, 0.0); d[i] = -1.0; directions.push_back(d);
    }
    srand(42);  // fixed seed: reproducible diagonal directions
    while (directions.size() < actual_n_rays) {
        vectord d((size_t) n);
        double len = 0.0;
        for (size_t i = 0; i < n; i++) {
            d[i] = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
            len += d[i] * d[i];
        }
        len = sqrt(len);
        if (len < 1e-12) continue;
        for (size_t i = 0; i < n; i++) d[i] /= len;
        directions.push_back(d);
    }

    // Step size for bracket search: divides unit interval into n_bracket steps.
    const double step0 = 1.0 / (double)(n_bracket + 1);

    // Accumulate boundary points across all rounds for centroid computation.
    vector<vectord> all_boundary_pts;

    // --- Single-round ray-bisection lambda ---
    // Runs all rays from x_centre, appends found boundary points to
    // all_boundary_pts, returns number of boundaries found this round.
    auto run_round = [&](const vectord &x_centre, unsigned round) -> unsigned
    {
        fprintf(stderr,
            "\nBOOptimizer[binary]: round %u/%u, centre=(",
            round + 1, n_rounds);
        for (unsigned i = 0; i < n; i++)
            fprintf(stderr, "%s%.4f", i ? ", " : "", x_centre[i]);
        fprintf(stderr, ")\n");

        unsigned found_this_round = 0;

        for (unsigned ray = 0; ray < actual_n_rays; ray++) {
            const vectord &dir = directions[ray];
            double t_pass = 0.0;
            double t_fail = -1.0;   // -1 = not yet found

            // --- Bracket search ---
            for (unsigned step = 1; step <= n_bracket; step++) {
                double t = step * step0;
                vectord xp((size_t) n);
                bool wall = false;
                for (size_t i = 0; i < n; i++) {
                    double xi = x_centre[i] + t * dir[i];
                    xp[i] = clip01(xi);
                    if (xp[i] != xi) wall = true;
                }
                double s = eval_point(xp);
                fprintf(stderr, "  ray %u step %u: t=%.4f score=%.4g\n",
                        ray, step, t, s);

                if (s <= threshold) {
                    t_fail = t;
                    break;          // bracket found: [t_pass, t_fail]
                }
                t_pass = t;         // still passing: advance lower bracket

                if (wall) {
                    // Hit the [0,1] wall still passing -- no boundary here
                    fprintf(stderr,
                        "  ray %u: wall at t=%.4f, still passing\n", ray, t);
                    break;
                }
            }

            if (t_fail < 0.0) {
                // No failure found along this ray within parameter range
                fprintf(stderr, "  ray %u: no boundary found\n", ray);
                continue;
            }

            // --- Bisection: refine bracket [t_pass, t_fail] ---
            // Note: t_pass==0 is valid -- it means the first step already
            // failed, so the bracket is [0, step0].  Bisection still works.
            for (unsigned k = 0; k < n_bisect; k++) {
                double tm = (t_pass + t_fail) / 2.0;
                vectord xp((size_t) n);
                for (size_t i = 0; i < n; i++)
                    xp[i] = clip01(x_centre[i] + tm * dir[i]);
                if (eval_point(xp) > threshold)
                    t_pass = tm;
                else
                    t_fail = tm;
            }

            // Record the last-passing boundary point
            vectord bpt((size_t) n);
            for (size_t i = 0; i < n; i++)
                bpt[i] = clip01(x_centre[i] + t_pass * dir[i]);
            all_boundary_pts.push_back(bpt);
            found_this_round++;

            fprintf(stderr,
                "  ray %u: t_pass=%.6f t_fail=%.6f (%u oracle calls total)\n",
                ray, t_pass, t_fail, n_evals);
        }

        return found_this_round;
    };

    // --- Iterative rounds with centre relocation ---
    vectord x_centre = x_star;

    for (unsigned round = 0; round < n_rounds; round++) {
        run_round(x_centre, round);

        // After all rounds except the last: relocate to centroid
        if (round + 1 >= n_rounds || all_boundary_pts.empty())
            continue;

        // Compute centroid of all boundary points found so far
        vectord centroid((size_t) n, 0.0);
        for (const vectord &bp : all_boundary_pts)
            for (size_t i = 0; i < n; i++)
                centroid[i] += bp[i];
        for (size_t i = 0; i < n; i++)
            centroid[i] /= (double) all_boundary_pts.size();

        // Verify centroid passes; if not, bisect back toward current centre
        double score_c = eval_point(centroid);
        if (score_c > threshold) {
            fprintf(stderr,
                "BOOptimizer[binary]: centroid passes (score=%.4g), "
                "relocating centre\n", score_c);
            x_centre = centroid;
        } else {
            // Centroid fails -- bisect between x_centre (pass) and centroid (fail)
            fprintf(stderr,
                "BOOptimizer[binary]: centroid fails (score=%.4g), "
                "bisecting back toward current centre\n", score_c);
            vectord x_lo = x_centre, x_hi = centroid;
            for (unsigned k = 0; k < 16; k++) {
                vectord xm((size_t) n);
                for (size_t i = 0; i < n; i++)
                    xm[i] = (x_lo[i] + x_hi[i]) / 2.0;
                if (eval_point(xm) > threshold)
                    x_lo = xm;
                else
                    x_hi = xm;
            }
            x_centre = x_lo;
            fprintf(stderr,
                "BOOptimizer[binary]: relocated to (");
            for (unsigned i = 0; i < n; i++)
                fprintf(stderr, "%s%.4f", i ? ", " : "", x_centre[i]);
            fprintf(stderr, ")\n");
        }
    }

    // Helper: convert normalised [0,1] value to physical units via the
    // parameter's own mapping.  set_mapped_value() sets the parameter;
    // get_cur_value() reads back the physical value.
    // Note: this temporarily modifies the parameter value -- it is restored
    // by the caller before any subsequent oracle calls.
    auto to_physical = [&](unsigned idx, double norm) -> double {
        opt_params[idx]->set_mapped_value(norm);
        return opt_params[idx]->get_cur_value();
    };

    // --- Final margin report: axis-aligned rays from final centre ---
    // These rays are cheap: most bisection points will be cache hits.
    printf("\nBOOptimizer[binary]: margin analysis from final centre\n");
    printf("  Centre: (");
    for (unsigned i = 0; i < n; i++)
        printf("%s%s=%.6g", i ? ", " : "",
               opt_params[i]->get_name(), to_physical(i, x_centre[i]));
    printf(")\n\n");

    printf("%-20s  %12s  %12s  %12s  %12s  %12s  %8s\n",
           "parameter", "centre",
           "lo_pass", "lo_fail", "hi_pass", "hi_fail", "margin_%");

    for (unsigned i = 0; i < n; i++) {
        double t_pass_hi = 0.0, t_fail_hi = 0.0;
        double t_pass_lo = 0.0, t_fail_lo = 0.0;
        bool   hi_found = false, lo_found = false;

        // +axis (hi margin)
        {
            double tp = 0.0, tf = -1.0;
            for (unsigned step = 1; step <= n_bracket; step++) {
                double t = step * step0;
                vectord xp((size_t) n);
                bool wall = false;
                for (size_t j = 0; j < n; j++) {
                    double xj = x_centre[j] + t * (j == i ? 1.0 : 0.0);
                    xp[j] = clip01(xj);
                    if (xp[j] != xj) wall = true;
                }
                double s = eval_point(xp);
                if (s <= threshold) { tf = t; break; }
                tp = t;
                if (wall) break;
            }
            if (tf >= 0.0) {
                for (unsigned k = 0; k < n_bisect; k++) {
                    double tm = (tp + tf) / 2.0;
                    vectord xp((size_t) n);
                    for (size_t j = 0; j < n; j++)
                        xp[j] = clip01(x_centre[j] + tm*(j==i?1.0:0.0));
                    if (eval_point(xp) > threshold) tp = tm; else tf = tm;
                }
                t_pass_hi = tp; t_fail_hi = tf; hi_found = true;
            } else {
                t_pass_hi = tp; t_fail_hi = tp;  // reached wall still passing
            }
        }

        // -axis (lo margin)
        {
            double tp = 0.0, tf = -1.0;
            for (unsigned step = 1; step <= n_bracket; step++) {
                double t = step * step0;
                vectord xp((size_t) n);
                bool wall = false;
                for (size_t j = 0; j < n; j++) {
                    double xj = x_centre[j] + t * (j == i ? -1.0 : 0.0);
                    xp[j] = clip01(xj);
                    if (xp[j] != xj) wall = true;
                }
                double s = eval_point(xp);
                if (s <= threshold) { tf = t; break; }
                tp = t;
                if (wall) break;
            }
            if (tf >= 0.0) {
                for (unsigned k = 0; k < n_bisect; k++) {
                    double tm = (tp + tf) / 2.0;
                    vectord xp((size_t) n);
                    for (size_t j = 0; j < n; j++)
                        xp[j] = clip01(x_centre[j] + tm*(j==i?-1.0:0.0));
                    if (eval_point(xp) > threshold) tp = tm; else tf = tm;
                }
                t_pass_lo = tp; t_fail_lo = tf; lo_found = true;
            } else {
                t_pass_lo = tp; t_fail_lo = tp;
            }
        }

        double margin_pct = 100.0 * (t_pass_lo + t_pass_hi) / 2.0;

        // Physical values: set normalised coord, read back physical value.
        // Centre, lo boundary (centre - lo_pass along axis), hi boundary.
        double phys_centre   = to_physical(i, x_centre[i]);
        double phys_lo_pass  = to_physical(i, clip01(x_centre[i] - t_pass_lo));
        double phys_lo_fail  = to_physical(i, clip01(x_centre[i] - t_fail_lo));
        double phys_hi_pass  = to_physical(i, clip01(x_centre[i] + t_pass_hi));
        double phys_hi_fail  = to_physical(i, clip01(x_centre[i] + t_fail_hi));

        printf("%-20s  %12.6g  %12.6g  %12.6g  %12.6g  %12.6g  %7.2f%%%s%s\n",
               opt_params[i]->get_name(),
               phys_centre,
               phys_lo_pass, phys_lo_fail,
               phys_hi_pass, phys_hi_fail,
               margin_pct,
               lo_found ? "" : " (lo:wall)",
               hi_found ? "" : " (hi:wall)");
    }

    // All accumulated boundary points across all rounds (physical units)
    printf("\nAll boundary points (%zu across %u rounds):\n",
           all_boundary_pts.size(), n_rounds);
    for (size_t b = 0; b < all_boundary_pts.size(); b++) {
        printf("  pt %3zu:", b);
        for (size_t i = 0; i < n; i++)
            printf("  %s=%.6g", opt_params[i]->get_name(),
                   to_physical(i, all_boundary_pts[b][i]));
        printf("\n");
    }

    printf("\nBOOptimizer[binary]: oracle calls=%u\n", n_evals);
    cache.print_stats(stdout);
    parameter::save_result(result_fp);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// SUBMODE_GRADIENT -- surrogate-gradient straddle search from x*
//
// Used when the oracle returns a continuous signed margin with a meaningful
// gradient (e.g. AC synchronizer gate).  GP-based margin extraction at end.
//

void BOOptimizer::run_gradient(FILE *result_fp,
                               vector<const_parameter*> &opt_params,
                               unsigned n)
{
    vectord x_star((size_t) n);
    for (unsigned i = 0; i < n; i++)
        x_star[i] = opt_params[i]->get_mapped_value();

    {
        for (unsigned i = 0; i < n; i++)
            opt_params[i]->set_mapped_value(x_star[i]);
        loop_complex.run_once(sum_fp);
        double margin_star = obj_funct->get_cur_value();
        if (margin_star <= threshold)
            fprintf(stderr,
                "BOOptimizer[gradient]: WARNING: starting point score=%.6g "
                "does not pass threshold=%.6g\n", margin_star, threshold);
        else
            fprintf(stderr,
                "BOOptimizer[gradient]: x* verified, score=%.6g\n",
                margin_star);
    }

    // Higher noise: boundary probing clusters evaluations, stressing the
    // GP kernel matrix at low noise.
    auto params = make_bo_params(n_iterations, 1, 1e-4);
    EvalCache cache(n, cache_capacity(n_iterations));

    JoSimBO optimizer(params, opt_params, cache, sum_fp,
                      obj_funct, JoSimBO::ROBUSTNESS);
    optimizer.x_star = x_star;
    optimizer.initializeOptimization();

    vectord x_cur = x_star;
    const double beta        = 1.96;
    const double max_step    = 0.05;
    const double noise_scale = 0.02;

    for (unsigned iter = 0; iter < n_iterations; iter++) {
        double grad_mag;
        vectord direction = optimizer.boundary_gradient(x_cur, grad_mag);

        double mu_cur = fabs(optimizer.surrogate_mean(x_cur) - threshold);
        double step   = (grad_mag > 1e-12)
                        ? min(mu_cur / grad_mag, max_step)
                        : max_step;

        vectord x_next((size_t) n);
        for (size_t i = 0; i < n; i++) {
            double sigma_i = optimizer.surrogate_std(x_cur);
            double noise   = noise_scale * sigma_i
                             * (2.0 * ((rand() / (double)RAND_MAX) - 0.5));
            x_next[i] = clip01(x_cur[i] + step * direction[i] + noise);
        }

        if (optimizer.straddle(x_next, beta) > optimizer.straddle(x_cur, beta))
            x_cur = x_next;

        optimizer.stepOptimization();

        fprintf(stderr, "gradient iter %u: step=%.4f grad=%.4f mu=%.4f\n",
                iter, step, grad_mag, optimizer.surrogate_mean(x_cur));
    }

    report_margins_gp(optimizer, x_star, opt_params, n, threshold,
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
            // run_binary(result_fp, opt_params, n);
            {
                //
                // This integration is just temporary to get things off the ground.
                // Doing all the work in the constructor is not ideal, and cleaning
                // things up after the work is done needs some more care.
                //
                // Also note that there are 4 variables that control how much computation
                // is expended for the optimization and how that CPU time is used.
                // Right now, the default parameters from the constructor are used,
                // but this needs to be controllable from the spice deck.
                // Adding this facility and choosing more sensible defaults is TBD
                //
                auto *bef = new bin_ellipsoid_fit(result_fp, opt_params, obj_funct, sum_fp);
                if (bef_budget > 0) {
                    // Build a plan from the pragma-specified budget and overrides.
                    // n_dim = number of tuneable parameters.
                    unsigned n_dim = (unsigned) opt_params.size();
                    bef_plan plan(n_dim, bef_budget,
                                  (bef_iter   > 0) ? bef_iter   : 5,
                                  (bef_probes > 0) ? bef_probes : 16);
                    bef->run(plan);
                } else {
                    bef->run();     // no budget specified: use bef_plan defaults
                }
            }
            break;
        case SUBMODE_GRADIENT:
            run_gradient(result_fp, opt_params, n);
            break;
        case SUBMODE_PROBABILISTIC:
            fprintf(stderr,
                "BOOptimizer: probabilistic submode not yet implemented\n");
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

//
// bef_plan.h -- Planning and configuration for bin_ellipsoid_fit::run()
//
// All tuning knobs, algorithmic switches, and constraint bounds for the
// binary ellipsoid fit subsystem live here. This is the single authoritative
// place to understand what each parameter does and why it has its default value.
//
// USAGE
// -----
// Simple (budget only -- all derived automatically):
//
//   bef_plan plan(n_dim, 1000);   // 1000 JoSIM calls total
//   bef.run(plan);
//
// With overrides:
//
//   bef_plan plan(n_dim, 1000);
//   plan.n_iter           = 8;    // more outer iterations
//   plan.n_probes_per_ray = 20;   // finer boundary resolution
//   plan.finalise(n_dim);         // recompute derived values
//   bef.run(plan);
//
// From pragma attributes (bo_optimizer calls this path):
//
//   bef_plan plan(n_dim, budget, n_iter, n_probes_per_ray);
//   bef.run(plan);
//
// BUDGET DECOMPOSITION
// --------------------
// The simulation budget is split between the initial exploration (iteration 0)
// and the ellipsoid shell searches (iterations 1..n_iter-1).
//
// Each iteration gets an equal share:
//   budget_per_iter = sim_budget / n_iter
//
// Iteration 0 (initial explore):
//   Axis-aligned rays:  2 * n_dim  (always present, not negotiable)
//   Extra random rays:  n_extra_rays  (derived)
//   Total rays:         (2*n_dim + n_extra_rays)
//   Budget consumed:    (2*n_dim + n_extra_rays) * n_probes_per_ray
//
//   => n_extra_rays = budget_per_iter / n_probes_per_ray - 2*n_dim
//
// Iterations 1..n_iter-1 (shell search):
//   n_shell_probes = budget_per_iter   (one sim per probe point)
//
// Warnings are issued if:
//   - n_extra_rays < 2 * n_dim  (poor angular coverage)
//   - n_shell_probes < 10 * n_dim  (too sparse for meaningful shell sampling)
//
// TIMING ESTIMATE
// ---------------
// After iteration 0 completes, run() measures the wall-clock time and
// prints an estimated total run time to stderr. This estimate assumes
// all simulations take the same time as the average in iteration 0.
//

#pragma once

#include <stdio.h>
#include <assert.h>

struct bef_plan {

    ///////////////////////////////////////////////////////////////////////////
    //
    // User-visible knobs (set via pragma attributes or directly)
    //

    // sim_budget
    //   Total JoSIM simulation budget for the entire run.
    //   This is the primary control: everything else is derived from it
    //   unless explicitly overridden.
    //
    unsigned sim_budget;

    // n_iter
    //   Number of outer iterations.
    //   Iteration 0: initial ray exploration from x_start.
    //   Iterations 1..n_iter-1: shell search + outlier rejection +
    //                            hp_filter (if i>1) + LM solve.
    //   More iterations allow the ellipsoid to relocate and refine.
    //   5 works for 2D; 8-10 is better for 6-8D.
    //   Minimum: 2 (one explore + one shell/solve).
    //
    unsigned n_iter;

    // n_probes_per_ray
    //   Simulation budget per ray, shared between:
    //     mode0 (binary search when endpoint is a fail), and
    //     mode1 (uniform probing when both endpoints pass).
    //   Controls the resolution vs. coverage tradeoff:
    //     Low values (8-12):  fast, broad angular coverage, coarse boundaries.
    //     High values (24+):  precise boundaries but fewer rays, poor coverage.
    //   16 is a good default. Raising it beyond ~24 is usually counterproductive:
    //   you spend a lot of budget nailing one boundary point to high precision
    //   while missing the overall ellipsoid shape.
    //
    unsigned n_probes_per_ray;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Derived values (computed by finalise(), read-only after that)
    //
    // These are computed from the three user-visible knobs above.
    // Override them after calling finalise() only if you have a specific reason.
    //

    // n_extra_rays
    //   Random-direction rays beyond the 2*n_dim axis-aligned rays.
    //   Derived: budget_per_iter / n_probes_per_ray - 2*n_dim
    //   A warning is printed if this is < 2*n_dim (poor angular coverage).
    //
    unsigned n_extra_rays;

    // n_shell_probes
    //   Number of shell probe points per e_shell_search() call.
    //   Derived: budget_per_iter (one sim per probe).
    //   A warning is printed if this is < 10*n_dim (too sparse).
    //
    unsigned n_shell_probes;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Internal knobs (not user-visible via pragma, but can be set in code)
    //

    // n_candidates
    //   For each extra random ray, n_candidates directions are drawn and the
    //   one with the lowest cosine-similarity to all prior directions is chosen.
    //   Higher = more uniform angular coverage but O(n_candidates^2) cost.
    //   32 is a good balance. Rarely needs tuning.
    //
    unsigned n_candidates = 32;

    // shell_a
    //   Sigmoid sharpness used when sampling shell probe directions in
    //   e_shell_search(). Larger = narrower band around the phi=0.5 surface.
    //   Should be >= 50 for a well-converged ellipsoid.
    //   Too small wastes probes deep inside or far outside.
    //
    double shell_a = 50.0;

    // outlier_frac
    //   Fraction of points rejected as outliers each iteration.
    //   Ranked by |fitted - actual|; the worst outlier_frac fraction is excluded
    //   from the next LM fit. Too high risks confirmation bias; too low leaves
    //   noisy boundary points in the data. 0.05 (5%) is conservative.
    //
    double outlier_frac = 0.05;

    ///////////////////////////////////////////////////////////////////////////
    //
    // LM beam search knobs
    //

    // n_lm_iterations
    //   Maximum LM steps per solve_one() call. Each step accumulates all
    //   non-outlier data points and takes one LM parameter update.
    //   200 is conservative; convergence usually happens in 50-100 steps.
    //
    unsigned n_lm_iterations = 200;

    // n_beam_solutions
    //   Maximum number of candidate solutions retained in the beam between
    //   solve() calls. More = more robust against bad starting points but
    //   proportionally more LM work per solve(). 16 is generous; 8 suffices
    //   for most cases.
    //
    unsigned n_beam_solutions = 16;

    // n_beam_derive
    //   Number of top beam solutions from which derive_solution() is called
    //   each round. Each call may add up to 3 derived solutions.
    //   Must be < n_beam_solutions to leave room for derived entries.
    //
    unsigned n_beam_derive = 4;

    // clean_solution_tolerance
    //   When selecting the best solution to carry forward, a solution with
    //   ec >= 0 (LM converged cleanly) is preferred over one with ec < 0,
    //   provided its residual is within this factor of the best residual.
    //   1.1 = accept a clean solution with up to 10% higher residual.
    //
    double clean_solution_tolerance = 1.1;

    // a_init_sans_a
    //   Initial sigmoid sharpness for the sans-a fit variants (where <a> is
    //   held fixed, not subject to LM). Must be in (min_a, max_a).
    //   50 sits in the middle of the useful range.
    //
    double a_init_sans_a = 50.0;

    // a_incr_ramp
    //   Per-LM-step increment applied to <a> in the ramping sans-a variant.
    //   0.5 ramps a from a_init_sans_a toward max_a over n_lm_iterations.
    //   Set to 0 to disable ramping (flat sans-a only).
    //
    double a_incr_ramp = 0.5;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Focal escape recovery
    //

    // n_wall_points
    //   Synthetic fail points added per build_wall() call when a focal point
    //   escapes the parameter space. More = stronger wall but risks data
    //   imbalance (synthetic fails outnumbering real pass points).
    //
    unsigned n_wall_points = 150;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Algorithmic strategy flags
    //
    // These replace the compile-time #defines and allow runtime experimentation
    // without recompilation. Defaults match the currently-active behaviour.
    //

    // center_escape_only  [was: _POK_ELLIPSOID_CENTER_ONLY_ defined]
    //   true  (default): foci may wander outside [0,1] as long as the ellipsoid
    //     centre stays inside. Found to be significantly better on the D2 5D
    //     case: foci naturally return inside after initial excursions.
    //   false: both foci must remain inside [0,1] at all times. More conservative
    //     but prevents LM from exploring useful parameter configurations.
    //
    bool center_escape_only = true;

    // hp_use_ellipsoid_normal  [was: _HP_FILTER_ELLIPSOID_NORMAL_ not defined]
    //   true:  hyperplane normal in hp_filter() = angular bisector of focal
    //     direction vectors at the fail point (approximate ellipsoid normal).
    //   false (default): normal = vector from fail point to ellipsoid centre.
    //     Geometrically less principled but empirically better on D2 latch,
    //     possibly because the centre direction is more stable when the
    //     ellipsoid fit is not yet well-converged.
    //
    bool hp_use_ellipsoid_normal = false;

    // pass_only_outliers  [was: _BIN_EFIT_PASS_ONLY_OUTLIER_ defined]
    //   true  (default): only pass points are candidates for outlier rejection.
    //     Conservative: a fail point marked outlier would encourage the ellipsoid
    //     to extend into a region that genuinely fails.
    //   false: both pass and fail points can be outliers.
    //
    bool pass_only_outliers = true;

    // synthetic_fails_on_rays  [was: _ADD_SYNTHETIC_FAILS_ not defined]
    //   true:  after the first fail on a ray, add synthetic fail points along
    //     the outbound portion without running simulations. Was tried and found
    //     to make results noticeably worse due to over-influence of fail clusters.
    //     May be worth revisiting with reduced weights in nl_lsq_fit.
    //   false (default): no synthetic fails. Only real simulation results used.
    //
    bool synthetic_fails_on_rays = false;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Constraint bounds
    //
    // Hard limits enforced by bef_param_ok() / bef_param_ok_sa().
    // Declared const to make clear they should not change during a run.
    // The defaults match the values historically hardcoded in bef_functions.cpp.
    //
    // Note: once bef_plan is fully wired in, the file-scoped constants in
    // bef_functions.cpp should be removed and bef_param_ok[_sa]() should
    // receive the plan (by const reference or pointer).
    //
    const double min_focal_sum = 0.01;  // fs floor: below this is degenerate
    const double min_a         = 0.001; // sigmoid sharpness floor
    const double max_a         = 100.0; // sigmoid sharpness ceiling
    const double min_foci_d    = 0.01;  // minimum inter-focal distance
    const double min_fs_excess = 0.1;   // fs must exceed foci_d by this fraction
    const double gcs_min_scale = 0.1;   // per-axis scale factor lower bound
    const double gcs_max_scale = 10.0;  // per-axis scale factor upper bound

    ///////////////////////////////////////////////////////////////////////////
    //
    // Constructors
    //

    // Default constructor -- must call finalise(n_dim) before use.
    bef_plan() : sim_budget(0), n_iter(5), n_probes_per_ray(16),
                 n_extra_rays(0), n_shell_probes(0) {}

    // Primary constructor: derive everything from budget + n_dim.
    // n_iter and n_probes_per_ray use their defaults unless overridden.
    bef_plan(unsigned n_dim, unsigned budget,
             unsigned iter = 5, unsigned probes_per_ray = 16)
        : sim_budget(budget), n_iter(iter), n_probes_per_ray(probes_per_ray),
          n_extra_rays(0), n_shell_probes(0)
    {
        finalise(n_dim);
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    // finalise() -- compute derived values and validate.
    //
    // Call this after construction or after overriding any user-visible knob.
    // Prints warnings to stderr if derived values look problematic.
    //
    void finalise(unsigned n_dim)
    {
        assert(n_dim >= 2);
        assert(n_iter >= 2);
        assert(n_probes_per_ray >= 1);
        assert(sim_budget > 0);

        // Equal budget per iteration
        unsigned budget_per_iter = sim_budget / n_iter;

        // Iteration 0: axis rays always get their share first
        unsigned axis_ray_sims = 2 * n_dim * n_probes_per_ray;

        if (budget_per_iter <= axis_ray_sims) {
            // Budget is so tight that even axis rays consume the whole share.
            // Clamp: no extra rays, and shell gets what's left.
            n_extra_rays   = 0;
            n_shell_probes = (sim_budget - axis_ray_sims) / (n_iter - 1);
            fprintf(stderr,
                "bef_plan warning: budget too tight for extra rays "
                "(%u sims/iter, axis rays alone need %u). "
                "n_extra_rays=0, n_shell_probes=%u.\n",
                budget_per_iter, axis_ray_sims, n_shell_probes);
        } else {
            unsigned remaining_iter0 = budget_per_iter - axis_ray_sims;
            n_extra_rays   = remaining_iter0 / n_probes_per_ray;
            n_shell_probes = budget_per_iter;   // remaining iters each get same share
        }

        // Warning: poor angular coverage
        if (n_extra_rays < 2 * n_dim) {
            fprintf(stderr,
                "bef_plan warning: n_extra_rays=%u < 2*n_dim=%u -- "
                "angular coverage will be poor. "
                "Consider increasing sim_budget or decreasing n_probes_per_ray.\n",
                n_extra_rays, 2 * n_dim);
        }

        // Warning: sparse shell sampling
        if (n_shell_probes < 10 * n_dim) {
            fprintf(stderr,
                "bef_plan warning: n_shell_probes=%u < 10*n_dim=%u -- "
                "shell sampling will be too sparse for good boundary coverage. "
                "Consider increasing sim_budget or n_iter.\n",
                n_shell_probes, 10 * n_dim);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    // print() -- human-readable summary of the plan, for logging.
    //
    void print(FILE *fp, unsigned n_dim) const
    {
        fprintf(fp, "bef_plan: sim_budget=%u  n_iter=%u  n_probes_per_ray=%u\n",
                sim_budget, n_iter, n_probes_per_ray);
        fprintf(fp, "  derived: n_extra_rays=%u (axis=%u)  n_shell_probes=%u\n",
                n_extra_rays, 2*n_dim, n_shell_probes);
        fprintf(fp, "  internal: n_candidates=%u  shell_a=%.1f  "
                "outlier_frac=%.3f  n_lm_iter=%u\n",
                n_candidates, shell_a, outlier_frac, n_lm_iterations);
        fprintf(fp, "  beam: n_solutions=%u  n_derive=%u  "
                "clean_tol=%.2f  a_init=%.1f  a_ramp=%.2f\n",
                n_beam_solutions, n_beam_derive,
                clean_solution_tolerance, a_init_sans_a, a_incr_ramp);
        fprintf(fp, "  flags: center_escape_only=%d  hp_ellipsoid_normal=%d  "
                "pass_only_outliers=%d  synthetic_fails=%d\n",
                center_escape_only, hp_use_ellipsoid_normal,
                pass_only_outliers, synthetic_fails_on_rays);
        fprintf(fp, "  bounds: min_a=%.3f  max_a=%.1f  min_foci_d=%.3f  "
                "scale=[%.2f,%.2f]\n",
                min_a, max_a, min_foci_d, gcs_min_scale, gcs_max_scale);
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    // estimate_remaining_time() -- print a timing estimate to stderr.
    //
    // Call this after iteration 0 completes, passing the wall-clock time of
    // iteration 0 in seconds and the number of simulations actually run in it.
    // Remaining iterations each do n_shell_probes simulations; we scale
    // proportionally assuming constant per-simulation time.
    //
    void estimate_remaining_time(double iter0_seconds,
                                 unsigned iter0_sims,
                                 unsigned n_dim) const
    {
        if (iter0_sims == 0) return;

        double t_per_sim    = iter0_seconds / (double) iter0_sims;
        unsigned sims_left  = (unsigned)(n_iter - 1) * n_shell_probes;
        double   t_left_s   = t_per_sim * (double) sims_left;
        double   t_total_s  = iter0_seconds + t_left_s;

        auto fmt_time = [](double s, char *buf) {
            if (s < 60.0)
                snprintf(buf, 32, "%.0f s", s);
            else if (s < 3600.0)
                snprintf(buf, 32, "%.1f min", s / 60.0);
            else
                snprintf(buf, 32, "%.1f h", s / 3600.0);
        };

        char t_iter0[32], t_left[32], t_total[32];
        fmt_time(iter0_seconds, t_iter0);
        fmt_time(t_left_s,      t_left);
        fmt_time(t_total_s,     t_total);

        fprintf(stderr,
            "bef timing: iteration 0 took %s (%u sims, %.3f s/sim)\n"
            "  %u remaining iterations × %u shell probes = %u sims\n"
            "  estimated remaining: %s  |  estimated total: %s\n",
            t_iter0, iter0_sims, t_per_sim,
            n_iter - 1, n_shell_probes, sims_left,
            t_left, t_total);
    }
};

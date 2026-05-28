#pragma once

#include <stdio.h>
#include "eval_cache.h"
#include "parameter.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer -- Bayesian Optimisation via BayesOpt's C++ interface.
//
// Two top-level modes:
//
//   MODE_OPTIMIZE  (default, no "margin" keyword in pragma)
//     Minimises the objective function.  Oracle returns a cost (smaller=better).
//     BayesOpt drives toward the minimum.  Reports the best point found.
//
//   MODE_ROBUSTNESS  ("margin" keyword present)
//     Maps the feasibility boundary around a known-good starting point x*.
//     x* is captured from get_mapped_value() at run() entry -- set it to the
//     SA (or other) optimum before calling run().
//     Three submodes select the boundary-search strategy:
//
//     SUBMODE_BINARY       ("binary" keyword, default)
//       Oracle returns +1 (pass) or -1 (fail).  Ray-bisection search from x*
//       locates the zero crossing along each ray direction.  Cheap and reliable
//       for digital gates where the pass/fail decision is crisp.
//
//     SUBMODE_GRADIENT     ("gradient" keyword)
//       Oracle returns a continuous signed margin with a meaningful gradient.
//       Surrogate-gradient straddle search drives toward the boundary.
//       Suitable for analog/mixed gates like the AC synchronizer.
//
//     SUBMODE_PROBABILISTIC  ("probabilistic" keyword)
//       Oracle returns an error rate or probability of failure.  The boundary
//       is the iso-contour at the user-specified threshold.  The GP models the
//       error rate surface; margin is the region where GP mean < threshold with
//       high confidence.  Most expensive -- each oracle call may require many
//       JoSIM runs.  (Not yet implemented; reserved for future use.)
//
// Normalised space:
//   All search operates in [0,1]^n.  get_mapped_value() / set_mapped_value()
//   on const_parameter handle physical <-> normalised conversion.
//   This class never sees physical component values.
//
// Numeric pragma parameters (all have sensible defaults):
//   n_rays     : number of rays from x* for boundary search (0 = auto: 4*n)
//   n_bracket  : bracket search steps per ray before bisection
//   n_bisect   : bisection steps per ray to locate zero crossing
//   threshold  : pass/fail boundary value (default 0.0)
//

class BOOptimizer {
public:
    enum Mode    { MODE_OPTIMIZE, MODE_ROBUSTNESS };
    enum SubMode { SUBMODE_BINARY, SUBMODE_GRADIENT, SUBMODE_PROBABILISTIC };

    BOOptimizer(parameter *obf,
                unsigned   n_it      = 190,
                Mode       mode      = MODE_OPTIMIZE,
                SubMode    submode   = SUBMODE_BINARY,
                unsigned   n_rays    = 0,       // 0 = auto
                unsigned   n_bracket = 10,
                unsigned   n_bisect  = 20,
                double     threshold = 0.0,
                unsigned   bef_budget  = 0,     // 0 = bef_plan default
                unsigned   bef_iter    = 0,     // 0 = bef_plan default (5)
                unsigned   bef_probes  = 0);    // 0 = bef_plan default (16));

    bool in_use()  { return configured; }
    void run(FILE *result_fp);                  // called from main()
    void specify_summary_file(FILE *sfp) { sum_fp = sfp; }

private:
    bool        configured;
    unsigned    n_iterations;
    Mode        mode;
    SubMode     submode;
    parameter  *obj_funct;                      // objective / margin parameter
    FILE       *sum_fp;

    // Robustness search parameters
    unsigned    n_rays;         // 0 = auto-compute from param count
    unsigned    n_bracket;      // bracket search steps per ray
    unsigned    n_bisect;       // bisection steps per ray
    double      threshold;      // pass/fail boundary value
    // bef_plan overrides (0 = use bef_plan defaults)
    unsigned    bef_budget;     // total JoSIM budget for ellipsoid fit
    unsigned    bef_iter;       // outer iteration count for ellipsoid fit
    unsigned    bef_probes;     // probes per ray for ellipsoid fit
    
    
    void run_optimize    (FILE *result_fp, vector<const_parameter*> &opt_params, unsigned n);
    void run_robustness  (FILE *result_fp, vector<const_parameter*> &opt_params, unsigned n);
    void run_binary      (FILE *result_fp, vector<const_parameter*> &opt_params, unsigned n);
    void run_gradient    (FILE *result_fp, vector<const_parameter*> &opt_params, unsigned n);
};

extern BOOptimizer *baysian_opt;                // pointer to instance when in use

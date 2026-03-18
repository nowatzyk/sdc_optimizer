#pragma once

#include <stdio.h>
#include "eval_cache.h"
#include "parameter.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer -- Bayesian Optimisation via BayesOpt's C++ interface.
//
// Two operating modes, selected at construction time:
//
//   MODE_OPTIMIZE  (default)
//     Minimises the objective function.  Oracle returns a cost (smaller=better).
//     BayesOpt drives toward the minimum.  Reports the best point found.
//
//   MODE_ROBUSTNESS
//     Maps the feasibility boundary around a known-good starting point x*.
//     Oracle returns a signed margin: positive=pass, negative=fail, 0=boundary.
//     A surrogate-gradient straddle search finds the zero-crossing of the GP
//     mean in each direction, reporting per-parameter margins as +/- %.
//     x* is captured from get_mapped_value() at run() entry -- it must be set
//     to the SA (or other) optimum before run() is called.
//
// Normalised space:
//   BayesOpt and the straddle search both work in [0,1]^n.
//   get_mapped_value() / set_mapped_value() on const_parameter handle
//   all physical <-> normalised conversion; this class never sees physical units.
//
// checkReachability():
//   Currently returns true for all points.  Wire up the reject expression here
//   once that pragma is available.
//

class BOOptimizer {
public:
    enum Mode { MODE_OPTIMIZE, MODE_ROBUSTNESS };

    BOOptimizer(parameter *obf, unsigned n_it = 190, Mode mode = MODE_OPTIMIZE);

    bool in_use()  { return configured; }
    void run(FILE *result_fp);                      // called from main()
    void specify_summary_file(FILE *sfp) { sum_fp = sfp; }

private:
    bool        configured;
    unsigned    n_iterations;
    Mode        mode;
    parameter  *obj_funct;                          // objective / margin parameter
    FILE       *sum_fp;

    void run_optimize   (FILE *result_fp, vector<const_parameter*> &opt_params, unsigned n);
    void run_robustness (FILE *result_fp, vector<const_parameter*> &opt_params, unsigned n);
};

extern BOOptimizer *baysian_opt;                    // pointer to instance when in use

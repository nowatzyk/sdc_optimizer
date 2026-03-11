#pragma once

#include <stdio.h>
#include "eval_cache.h"
#include "parameter.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// BOOptimizer -- Bayesian Optimisation via BayesOpt's C++ interface.
//
// Uses ContinuousModel inheritance -- the natural C++ interface.
//
// Configuration:
//   Call BOOptimizer::configure(n_iterations) from the "bo" pragma handler.
//   The optimisable parameters are those in parameter::get_opt_params()
//   (the same set used by simulated annealing, i.e. constant_parameters).
//
// Normalised space:
//   BayesOpt searches in [0,1]^n.  Each parameter maps this through
//   parameter::set_normalized_value() / to_normalized() using its own
//   linear (or log) scale.  The EvalCache also operates in normalised
//   space for scale-independence.
//
// Minimisation:
//   BayesOpt minimises by default.  The oracle returns eval_expr->get_value()
//   which the user writes as a cost (smaller = better).
//   NaN / Inf oracle results are replaced by DBL_MAX.
//
// checkReachability():
//   Currently returns true for all points.  Once the reject pragma is wired
//   up, evaluate it here and return false for infeasible points -- this is
//   more informative to the GP surrogate than NaN -> DBL_MAX substitution.
//

class BOOptimizer {
    //
    // This is a place-holer at this point.
    // Here is the place to store parameters extracted from the spice deck
    //
    bool        configured;
    unsigned    n_iterations;
    parameter   *obj_funct;                         // What to optimiza for (smaller is better)
    FILE        *sum_fp;

public:
    BOOptimizer(parameter *obf, unsigned n_it = 190);

    bool in_use()  { return configured; }
    void run(FILE *result_fp);                      // called from main()

    void specify_summary_file(FILE *sfp) {sum_fp = sfp;};
};

extern BOOptimizer *baysian_opt;                    // Pointer to instance when in use

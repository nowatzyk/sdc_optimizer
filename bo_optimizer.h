#pragma once

#include <stdio.h>
#include "eval_cache.h"

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
public:
    static void configure(unsigned n_iter);         // called from pragma handler
    static bool in_use()  { return configured; }
    static void run(FILE *result_fp, FILE *sum_fp); // called from main()

private:
    static bool     configured;
    static unsigned n_iterations;
};

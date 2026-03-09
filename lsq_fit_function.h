#pragma once

#include "expression.h"

extern "C" {
#include "lsq_fit.h"
}

class LoopComplex;          //Forward declaration

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// lsq_fit_function -- wraps an lsq_fit system as a LoopComplex-managed object.
//
// Lifecycle (called by LoopComplex at the appropriate loop level boundary):
//   initialize() -- resets fit state, ready for new data
//   add_datum()  -- called after every JoSIM run (via LoopComplex::update_all)
//   finalize()   -- solves the fit; results accessible via get_ci() / get_residual()
//
// Self-registers with loop_complex in its constructor.
// The static linked list (root/next) is removed: LoopComplex owns the registry.
//
// find() is retained for the parser's define_function2() lookup.
//

class lsq_fit_function {
public:
    unsigned            level;          // nesting level, set by LoopComplex::register_lsq_fit()

    lsq_fit_function(const char *name, int order,
                     expression *free_var, expression *dep_var,
                     LoopComplex &lc);

    // Lifecycle interface (called by LoopComplex):
    void    initialize();               // reset + init_lsq_fit
    void    add_datum();                // add one data point (if both vars are finite)
    void    finalize();                 // solve fit; populate coefficients / residual

    // Result accessors (used by expression evaluator via function2 callbacks):
    double  fit_ok()     const { return isfinite(residual) ? 1.0 : 0.0; }
    double  get_ci(unsigned i) const { assert(i <= order); return coefficients[i]; }
    double  get_residual()     const { return residual; }

    // Parser lookup -- find by name:
    static lsq_fit_function *find(const char *name);

private:
    const char          *name;
    struct lsq_fit      *lsq_fit_ptr;

    expression          *free_var;
    expression          *dep_var;

    unsigned            order;
    double              residual;
    double              coefficients[5];    // up to 4th order polynomial

    void                clear_results();    // sets residual and coefficients to NAN

    // Registry for find() -- simple vector, not a linked list:
    static vector<lsq_fit_function *> registry;
};


double lsq_fit_get_cx(void *obj_ptr, double x);     // Retrieves the fitted coefficient(s)
double test_lsq_fit(void *obj_ptr, double x);       // Checks if the fit succeeded
double lsq_fit_residual(void *obj_ptr, double x);   // Provides fit residual (quality info)

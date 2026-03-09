#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include "expression.h"

using namespace std;

class lsq_fit_function;     // forward declaration

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// LoopIterator -- abstract base for iteration strategies.
//
// Each iterator owns one parameter and knows how to step through its range.
// Stepping is done in normalised [0,1] space; the parameter's to_physical()
// handles log/linear mapping.  Level 0 = innermost loop.
//
// Concrete subclasses (ScanIterator, BinarySearchIterator) are registered
// with LoopComplex at parse time from their constructors.
//

class LoopIterator {
public:
    unsigned    level;          // 0 = innermost; set by LoopComplex::register_iterator()

    virtual void initialize() = 0;  // set parameter to first value, reset counters
    virtual bool next()       = 0;  // advance to next value; return false when exhausted
    virtual ~LoopIterator()   = default;

protected:
    LoopIterator() : level(0) {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// LoopComplex -- the evaluation oracle.
//
// Built entirely at parse time from pragma definitions.  Owns two separate
// registries for stateful objects so that the correct finalization order is
// always enforced: StatefulExpressions finalize before LsqFits at each level.
//
// run_once() executes the full nested loop, running JoSIM at each innermost
// point, and returns the value of eval_expr (or NaN if none is defined).
//
// The global instance is defined in loop_complex.cpp and declared extern here
// so that parser-side factory functions (define_function1 etc.) can reach it.
//

class LoopComplex {
public:
    LoopComplex();

    //
    // Registration -- called from constructors at parse time
    //
    void register_iterator  (LoopIterator       *it);   // increments current_level
    void register_stateful  (StatefulExpression  *s);   // captures current_level
    void register_lsq_fit   (lsq_fit_function    *f);   // captures current_level; finalized after stateful

    unsigned current_nesting_level() const { return current_level; }

    //
    // eval_expr -- the expression whose value run_once() returns.
    // Set by the parser when it encounters the "eval" pragma.
    // nullptr is valid: run_once() returns NaN in that case (plain sweep mode).
    //
    expression *eval_expr;

    //
    // run_once() -- execute the full nested loop and return eval_expr->get_value().
    //
    // Calls the JoSIM subprocess indirectly via the circuit / fifo machinery
    // that lives in main.cpp.  The boundary between LoopComplex and main is
    // the two callbacks below, which main sets before calling run_once().
    //
    // Returns NaN if no eval_expr is defined.
    //
    double run_once();

    //
    // Callbacks into main for the operations LoopComplex cannot do itself.
    // Must be set before run_once() is called.
    //
    void (*run_josim_cb)()   = nullptr;  // writes cir file, forks JoSIM, reads CSV, waits
    void (*post_run_cb)()    = nullptr;  // any post-run bookkeeping main needs (e.g. list_c_val)

    //
    // n_optimizable_params() -- count of parameters eligible for SA/BO perturbation.
    // Used by Optimizer to size the EvalCache.
    //
    unsigned n_optimizable_params() const;

private:
    vector<LoopIterator*>      iterators;   // in registration order (outermost first)
    vector<StatefulExpression*> stateful;   // expressions: sum, avg, min, max, geomean
    vector<lsq_fit_function*>  lsq_fits;   // finalized after stateful at each level

    unsigned current_level;                 // incremented each time an iterator registers

    // Lifecycle helpers called by run_once()
    void initialize_level   (unsigned lv);
    void finalize_level     (unsigned lv);
    void update_all         ();
};

//
// Global instance -- analogous to 'circuit' in main.cpp.
// Defined in loop_complex.cpp.
//
extern LoopComplex loop_complex;

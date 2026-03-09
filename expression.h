#pragma once

#include <cmath>
#include <vector>
#include <cassert>

using namespace std;

//
// Forward declarations
//
class parameter;
class LoopComplex;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// exp_type -- type of expression node
//
// Note: sum_func is replaced by stateful_func. All stateful subclasses
//       (sum_expression, avg_expression, etc.) use this type. The base
//       class get_value() simply returns the value field, which the
//       subclass manages as its accumulator.
//

enum exp_type {
    constant,
    stateful_func,              // replaces sum_func: covers sum, avg, min, max, geomean
    function1,                  // function with 1 argument
    function2,                  // function with 2 arguments, first being an object pointer
    p_reference,
    addition,
    subtraction,
    multiplication,
    division,
    test_select,
    comp_eq,
    comp_ne,
    comp_gt,
    comp_ge,
    comp_lt,
    comp_le,
    undefined
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// StatefulExpression -- interface for expression nodes that accumulate
// state across JoSIM runs within a loop level.
//
// Subclasses inherit both expression (for pull evaluation via get_value())
// and StatefulExpression (for lifecycle management by LoopComplex).
//
// Registration with LoopComplex happens in the subclass constructor.
// LoopComplex calls initialize/update/finalize at loop level boundaries,
// always before LsqFit finalization at the same level.
//

class StatefulExpression {
public:
    unsigned    level;              // nesting level at registration time
                                    // 0 = innermost loop

    virtual void initialize() = 0; // called before first iteration of level
    virtual void update()     = 0; // called after every JoSIM run
    virtual void finalize()   = 0; // called when level is exhausted

    virtual ~StatefulExpression() = default;

protected:
    StatefulExpression() : level(0) {}
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// expression -- stateless recursive pull-evaluation tree node.
//
// get_value() recursively evaluates the dependency tree and returns
// the current value. NaN propagates naturally through arithmetic;
// isfinite() is used as the NaN test throughout (more robust than
// isnan() under aggressive optimisation flags).
//
// Stateful subclasses (sum_expression etc.) inherit this class and
// StatefulExpression. They use the protected constructor and own
// the value field as their accumulator. The base get_value() returns
// value directly for the stateful_func case, so no override is needed
// except for avg and geomean which compute on demand.
//

class expression {
protected:
    exp_type    type;
    expression  *l_arg, *r_arg;
    expression  *t_arg;
    double      value;                          // accumulator for stateful types; constant value otherwise
    double      (*func1_ptr)(double x);
    double      (*func2_ptr)(void *obj_ptr, double x);
    void        *obj_ptr;

    //
    // Protected constructor for stateful subclasses.
    // Sets type=stateful_func, value=0.0, l_arg=arg.
    // Only callable from subclass constructors.
    //
    expression(exp_type et, expression *arg);

    void        initialize();                   // zeroes all fields to safe defaults

public:
    expression(double x);                       // constant
    expression(double (*f_ptr)(double x),
               expression *arg_ptr);            // function1 (f_ptr must not be nullptr)
    expression(double (*f_ptr)(void *obj_ptr, double x),
               void *o_ptr,
               expression *arg_ptr);            // function2
    expression(expression *la_ptr,
               exp_type et,
               expression *ra_ptr);            // binary op
    expression(parameter *p_ptr);              // parameter reference
    expression(expression *test_ptr,
               expression *true_arg,
               expression *false_arg);         // ternary test-select

    virtual double get_value();                 // recursive pull evaluation
};


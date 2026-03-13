#include "expression.h"
#include "parameter.h"      // for get_cur_value() in p_reference case

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// expression -- private initializer
//

void expression::initialize()
{
    type      = undefined;
    l_arg     = nullptr;
    r_arg     = nullptr;
    t_arg     = nullptr;
    value     = NAN;
    func1_ptr = nullptr;
    func2_ptr = nullptr;
    obj_ptr   = nullptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// expression -- protected constructor for stateful subclasses
//

expression::expression(exp_type et, expression *arg)
{
    initialize();
    assert(et == stateful_func);
    type  = stateful_func;
    l_arg = arg;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// expression -- public constructors
//

expression::expression(double x)
{
    initialize();
    type  = constant;
    value = x;
}

expression::expression(double (*f_ptr)(double), expression *arg_ptr)
{
    initialize();
    assert(f_ptr != nullptr);   // nullptr no longer valid: use sum_expression subclass instead
    type      = function1;
    func1_ptr = f_ptr;
    l_arg     = arg_ptr;
}

expression::expression(double (*f_ptr)(void *, double), void *o_ptr, expression *arg_ptr)
{
    initialize();
    type      = function2;
    func2_ptr = f_ptr;
    obj_ptr   = o_ptr;
    l_arg     = arg_ptr;
}

expression::expression(expression *la_ptr, exp_type et, expression *ra_ptr)
{
    initialize();
    assert(et == addition    || et == subtraction  || et == multiplication ||
           et == division    || et == comp_eq      || et == comp_ne        ||
           et == comp_lt     || et == comp_le      || et == comp_gt        ||
           et == comp_ge);
    type  = et;
    l_arg = la_ptr;
    r_arg = ra_ptr;
}

expression::expression(parameter *param_ptr)
{
    initialize();
    type    = p_reference;
    obj_ptr = param_ptr;
}

expression::expression(expression *tst_ptr, expression *true_arg, expression *false_arg)
{
    initialize();
    type  = test_select;
    value = 0.0;
    t_arg = tst_ptr;
    l_arg = true_arg;
    r_arg = false_arg;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// expression::get_value() -- recursive pull evaluation
//
// NaN propagates naturally through arithmetic operations.
// isfinite() is used as the NaN/Inf test: more robust than isnan()
// under aggressive compiler optimisation flags.
//

double expression::get_value()
{
    double rtn = NAN;
    double l, r;

    switch (type) {

        case constant:
        case stateful_func:
            rtn = value;            // stateful subclass manages value as accumulator
            break;

        case function1:
            l = l_arg->get_value(); // Note: don't filter out NAN : the isfinite() funtion will break
            rtn = func1_ptr(l);
            break;

        case function2:
            l = l_arg->get_value();
            if (isfinite(l)) rtn = func2_ptr(obj_ptr, l);
            break;

        case p_reference:
            rtn = ((parameter *) obj_ptr)->get_cur_value();
            break;

        case addition:
            l = l_arg->get_value();
            if (isfinite(l)) rtn = l + r_arg->get_value();
            break;

        case subtraction:
            l = l_arg->get_value();
            if (isfinite(l)) rtn = l - r_arg->get_value();
            break;

        case multiplication:
            l = l_arg->get_value();
            if (isfinite(l)) rtn = l * r_arg->get_value();
            break;

        case division:
            l = l_arg->get_value();
            if (isfinite(l)) rtn = l / r_arg->get_value();
            break;

        case test_select:
            l = t_arg->get_value();
            if (isfinite(l)) rtn = (l > 0.0) ? l_arg->get_value() : r_arg->get_value();
            break;

        case comp_eq:
            l = l_arg->get_value(); r = r_arg->get_value();
            if (isfinite(l) && isfinite(r)) rtn = (l == r) ? 1.0 : 0.0;
            break;

        case comp_ne:
            l = l_arg->get_value(); r = r_arg->get_value();
            if (isfinite(l) && isfinite(r)) rtn = (l != r) ? 1.0 : 0.0;
            break;

        case comp_gt:
            l = l_arg->get_value(); r = r_arg->get_value();
            if (isfinite(l) && isfinite(r)) rtn = (l > r) ? 1.0 : 0.0;
            break;

        case comp_ge:
            l = l_arg->get_value(); r = r_arg->get_value();
            if (isfinite(l) && isfinite(r)) rtn = (l >= r) ? 1.0 : 0.0;
            break;

        case comp_lt:
            l = l_arg->get_value(); r = r_arg->get_value();
            if (isfinite(l) && isfinite(r)) rtn = (l < r) ? 1.0 : 0.0;
            break;

        case comp_le:
            l = l_arg->get_value(); r = r_arg->get_value();
            if (isfinite(l) && isfinite(r)) rtn = (l <= r) ? 1.0 : 0.0;
            break;

        default:
            assert(0);
    }

    return rtn;
}

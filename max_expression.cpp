#include "max_expression.h"
#include "loop_complex.h"

max_expression::max_expression(expression *arg, LoopComplex &lc)
    : expression(stateful_func, arg)
{
    value = NAN;
    lc.register_stateful(this);
}

void max_expression::update()
{
    double v = l_arg->get_value();
    if (isfinite(v) && (!isfinite(value) || v > value)) value = v;
}

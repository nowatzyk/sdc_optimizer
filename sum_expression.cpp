#include "sum_expression.h"
#include "loop_complex.h"

sum_expression::sum_expression(expression *arg, LoopComplex &lc)
    : expression(stateful_func, arg)
{
    lc.register_stateful(this);
}

void sum_expression::update()
{
    double v = l_arg->get_value();
    if (isfinite(v)) value += v;
}

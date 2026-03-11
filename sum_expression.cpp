#include "sum_expression.h"
#include "loop_complex.h"

sum_expression::sum_expression(expression *arg)
    : expression(stateful_func, arg)
{
    loop_complex.register_stateful(this);
}

void sum_expression::update()
{
    double v = l_arg->get_value();
    if (isfinite(v)) value += v;
}

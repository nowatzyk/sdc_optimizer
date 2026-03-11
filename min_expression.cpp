#include "min_expression.h"
#include "loop_complex.h"

min_expression::min_expression(expression *arg)
    : expression(stateful_func, arg)
{
    value = NAN;
    loop_complex.register_stateful(this);
}

void min_expression::update()
{
    double v = l_arg->get_value();
    if (isfinite(v) && (!isfinite(value) || v < value)) value = v;
}

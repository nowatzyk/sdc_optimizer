#include "geomean_expression.h"
#include "loop_complex.h"

geomean_expression::geomean_expression(expression *arg, LoopComplex &lc)
    : expression(stateful_func, arg), count(0)
{
    lc.register_stateful(this);
}

void geomean_expression::update()
{
    double v = l_arg->get_value();
    if (isfinite(v) && v > 0.0) { value += log(v); count++; }
}

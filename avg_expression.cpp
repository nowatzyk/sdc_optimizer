#include "avg_expression.h"
#include "loop_complex.h"

avg_expression::avg_expression(expression *arg, LoopComplex &lc)
    : expression(stateful_func, arg), count(0)
{
    lc.register_stateful(this);
}

void avg_expression::update()
{
    double v = l_arg->get_value();
    if (isfinite(v)) { value += v; count++; }
}

#pragma once
#include "expression.h"

//
// geomean_expression -- geometric mean of its argument across JoSIM runs.
// value field accumulates log(arg); get_value() returns exp(value/count).
// Non-positive or NaN arguments are skipped.
//
class geomean_expression : public expression, public StatefulExpression {
    unsigned count;
public:
    geomean_expression(expression *arg);
    void   initialize() override { value = 0.0; count = 0; }
    void   update()     override;
    void   finalize()   override { }
    double get_value()  override { return (count > 0) ? exp(value / count) : NAN; }
};

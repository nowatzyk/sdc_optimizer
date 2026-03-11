#pragma once
#include "expression.h"

//
// avg_expression -- computes average of its argument across JoSIM runs.
// value field holds the running sum; get_value() returns sum/count.
// Returns NaN if no finite values were seen.
//
class avg_expression : public expression, public StatefulExpression {
    unsigned count;
public:
    avg_expression(expression *arg);
    void   initialize() override { value = 0.0; count = 0; }
    void   update()     override;
    void   finalize()   override { }
    double get_value()  override { return (count > 0) ? value / (double) count : NAN; }
};

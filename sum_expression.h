#pragma once
#include "expression.h"

//
// sum_expression -- accumulates sum of its argument across JoSIM runs.
// get_value() returns the running sum (value field, managed as accumulator).
// NaN arguments are skipped.
//
class sum_expression : public expression, public StatefulExpression {
public:
    sum_expression(expression *arg);
    void initialize() override { value = 0.0; }
    void update()     override;
    void finalize()   override { }
};

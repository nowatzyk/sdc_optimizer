#pragma once
#include "expression.h"

//
// min_expression -- tracks minimum of its argument across JoSIM runs.
// Initialises to NaN (no samples yet). get_value() inherited: returns value.
//
class min_expression : public expression, public StatefulExpression {
public:
    min_expression(expression *arg);
    void initialize() override { value = NAN; }
    void update()     override;
    void finalize()   override { }
};

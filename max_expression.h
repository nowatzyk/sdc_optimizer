#pragma once
#include "expression.h"

//
// max_expression -- tracks maximum of its argument across JoSIM runs.
// Initialises to NaN (no samples yet). get_value() inherited: returns value.
//
class max_expression : public expression, public StatefulExpression {
public:
    max_expression(expression *arg, LoopComplex &lc);
    void initialize() override { value = NAN; }
    void update()     override;
    void finalize()   override { }
};

//
// Data structures for the CSV analyzer
//
#pragma once

const unsigned max_ln_length = 100000;      // We expect some very long lines...
                                            
extern const char* EVAL_PARAMETER;          // Eval parameter keyword
extern const char* REJECT_PARAMETER;        // Keyword fo a reject functon (name of a parameter)
extern const char* BAY_OPT_OBJECTIVE;       // Keyword to define the objective function for Baysian opt

void run_josim();                           // Execute one simulation run 


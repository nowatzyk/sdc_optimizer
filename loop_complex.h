#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include "expression.h"
#include "parameter.h"

using namespace std;

class lsq_fit_function;                     // forward declaration

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Run level: This class runs the actual simulation
//

class run_level {
    unsigned        level;                  // advisory / debugging
    parameter       *looping_p_ptr;         // Optional looping parameter

    vector <StatefulExpression *> s_expr_ptrs; // Stateful expressions that belong to this level
    vector <parameter *> p_updt_ptrs;       // Parameter update pointers
    
    run_level       *below;                 // The next run level below this one or nullptr
    
public:
    run_level(run_level *rl_ptr = nullptr); // constructor: meant to recursively contruct the
                                            // the loop control structure
    
    unsigned        is_looping ()           // Just a test if the current top level is executed more than once
                    {return (looping_p_ptr != nullptr);};
    
    void            add_looping_param(parameter *p_ptr);
    void            add_s_expr(StatefulExpression *se_ptr);
    void            add_param(parameter *p_ptr);
    
    void            run_this_level(FILE *sum_fp);  // Execute this level
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// LoopComplex -- the evaluation oracle.
//
// Built entirely at parse time from pragma definitions.  Owns two separate
// registries for stateful objects so that the correct finalization order is
// always enforced: StatefulExpressions finalize before LsqFits at each level.
//
// run_once() executes the full nested loop, running JoSIM at each innermost
// point. The results of this computation is avalable as updated values in
// the parameters: constant parameters may be updated by outer optimizers
// and expression-type parameters will be set to the result of the smulation.
//

class LoopComplex {
    
    run_level       *run_level_ptr;         // Points to the top level      

public:
    LoopComplex();

    //
    // Registration -- called from constructors at parse time
    //

    void register_stateful  (StatefulExpression  *se);   // captures current_level
    void register_parameter (parameter *pp, unsigned is_looping = 0);   // Registers a parameter

    void run_once(FILE *sum_fp);            // Perform a simulation run (incl. all looping and the kitchen sink
};

//
// Global instance -- analogous to 'circuit' in main.cpp.
// Defined in loop_complex.cpp.
//
extern LoopComplex loop_complex;

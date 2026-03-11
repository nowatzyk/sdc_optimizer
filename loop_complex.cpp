#include "loop_complex.h"
#include "csv_analyzer.h"       // for lsq_fit_function
#include "parameter.h"          // for update_assignments()
#include "nodes_of_interest.h"
#include "lsq_fit_function.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Global instance
//
LoopComplex loop_complex;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Run level
//

run_level::run_level(run_level* rl_ptr)
    : looping_p_ptr(nullptr)
{
    if (rl_ptr == nullptr) {
        // This is the bottom
        level = 0;
        below = nullptr;
    } else {
        level = rl_ptr->level + 1;
        below = rl_ptr;
    }
}

void run_level::run_this_level(FILE *sum_fp)
    //
    // This is primarily a recursive execution of a set of netsted loops. Each instance of
    // run_level will execute one loop, which is defined by a looping parameter instance.
    // The outermost/top-level loop is called, which in turn for each iteration will call the
    // instance of the next level down until the iterator (parameter.next()) returns a 1
    // which signals that the loop is complete.
    // At the bottom, the actual JoSIM simulation is called, which terminates the recursion.
    // Ideally, a level would be added only once a looping parameter is being defined. In this
    // case a level is always associated with a looping parameter. However, that is problematic
    // at the lowest level: for example parameters may be used to just configure a circuit
    // whithout any looping. Running a single simulation should be practical. So instead of
    // having specialized 0-level, a level may have no looping at all. This allows the construction
    // of an empty run_level instance that just runs the simulation and does nothing else.
    // Everything else is added later.
    //
{
        for (auto *lp : lsq_f_ptrs )  lp->initialize();
        for (auto *ep : s_expr_ptrs)  ep->initialize();
        
        if (looping_p_ptr != nullptr)
            looping_p_ptr->initialize();    // If there is a looping contruct, initialize it
            
        do {
            for (auto *pp : p_updt_ptrs)  pp->update();
            
            if (below != nullptr)
                below->run_this_level(sum_fp);    // recursively execute all looping levels
            else
                run_josim();                // Run the JoSIM simulator
                
            for (auto *ep : s_expr_ptrs)  ep->update();
            for (auto *lp : lsq_f_ptrs )  lp->add_datum();
            
            if (below == nullptr) {
                //
                // This is the bottom/innnermost/lowest level that does the reporting
                // for each JoSIM run
                //
                parameter::list_c_val(sum_fp);
                fprintf(sum_fp, " 0\n");
                n_josim_runs++;             // Done with this run
            }
            
        } while ((looping_p_ptr != nullptr) && (looping_p_ptr->next() == 0));
            
        for (auto *ep : s_expr_ptrs)  ep->finalize();
        for (auto *lp : lsq_f_ptrs )  lp->finalize();
        
        parameter::list_c_val(sum_fp);
        fprintf(sum_fp, " %u\n", level + 1);
}

void run_level::add_looping_param(parameter* p_ptr)
{
    assert(looping_p_ptr == nullptr);       // Make sure that this is done only once!
    
    looping_p_ptr = p_ptr;
    p_updt_ptrs.push_back(p_ptr);           // Note: looping parameters usually need updating too,
                                            //       so this takes care of it. No seperate registration
                                            //       required.
}

void run_level::add_s_expr(StatefulExpression* se_ptr)
{
    s_expr_ptrs.push_back(se_ptr);
}

void run_level::add_param(parameter* p_ptr)
{
    p_updt_ptrs.push_back(p_ptr);
}

void run_level::add_lsq_fit(lsq_fit_function* lf_ptr)
{
    lsq_f_ptrs.push_back(lf_ptr);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// LoopComplex
//
// This is basically just a warpper arround the run level core
//

LoopComplex::LoopComplex()
{
    run_level_ptr = new run_level(nullptr);     // Creates the bottom level / innermost loop
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Registration
//

void LoopComplex::register_stateful(StatefulExpression *se)
{
    assert(se);
    run_level_ptr->add_s_expr(se);
}

void LoopComplex::register_lsq_fit(lsq_fit_function *lf)
{
    assert(lf);
    run_level_ptr->add_lsq_fit(lf);
}

void LoopComplex::register_parameter(parameter* pp, unsigned int is_looping)
{
    assert((pp != nullptr) && (is_looping <= 1));
    
    if (is_looping == 0)
        run_level_ptr->add_param(pp);
    else {
        // Note: each level can have only one looping parameter
        if (run_level_ptr->is_looping() != 0)
            // Already has a looping parameter, so a new level is needed
            run_level_ptr = new run_level(run_level_ptr);
        run_level_ptr->add_looping_param(pp);
    }
}

void LoopComplex::run_once(FILE *sum_fp)
{
    run_level_ptr->run_this_level(sum_fp);
}

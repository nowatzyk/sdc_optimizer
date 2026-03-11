//
// A few functions to use simulated annealing to optimmize a function
//

#ifndef ANNEAL_DEFINED_
#define ANNEAL_DEFINED_

#include <vector>
using namespace std;

#include "expression.h"
#include "parameter.h"
#include "xrand.h"

#define _ADAPTIVE_RANGEING_         // If defined, the change scale is adjusted
                                    // individually for each parameter


///////////////////////////////////////////////////////////////////////////////////////////////////

struct an_sched {                   // anneal schedule
    double      temp;               // Starting temp for this phase
    unsigned    n_steps;            // #of steps
};

class sa_parameter;                 // forward declaration, just for cosmetics sake

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  The simulated annealer
//

class sim_anneal {
    vector<sa_parameter*> sa_p_ptrs;// parameters subject to optimization via sim-anealing
    vector<an_sched> sa_schedule;   // How to go about it

    parameter  *eval_ptr;           // The objective function
    parameter  *reject_ptr;         // an optional early reject function
    
    sa_parameter *changed_param;    // Which was changed last
    double       changed_value;     // What was its pre-change value
    
    double      best_eval;          // The best evaluation
    unsigned    have_best:1;        // if set, the best value is defined
    unsigned    changed_flag:1;     // Change-machinery: indicate validity
    
    FILE        *lfp;               // Log file pointer
    FILE        *sum_fp;            // summary file

    double      comp_eval();        // Objective function to tune!
                                    // Lower values are better, 0 is the best possible

public:
    sim_anneal(char *log_file_name, parameter *ev_ptr, parameter *rj_ptr);

    void add_sa_sched(double temp, unsigned n);  // Used to add to the SA schedule
    void add_parameter(const_parameter *p_ptr);  // adds an parameter

    int write(const char *fn);      // Write all parameters to a file
    int read(const char *fn);       // Read all parameters from a file

    void save_best(double best_ev); // Save the best
    double restore_best();          // restore the best
                                    
    void change(double cr);         // select a paramter and change it
    void unchange();                // un-do the last change

    unsigned n_tune() {return sa_p_ptrs.size();}; // Returns number of tunable parameters
    unsigned long n_steps();        // Returns the total number of annealing steps to be performed

    void optimize();                // The actual optimization function
    void save_result(FILE *fp);        // Save the results
    void specify_summary_file(FILE *sfp) {sum_fp = sfp;};
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  the sa_paramter class adds info to parameters subject to simulated annealing
//

class sa_parameter {

    const_parameter *p_ptr;         // Name of this paramter (needed to access it at run-time)

    double      best_val;           // Best value encountered so far

#ifdef _ADAPTIVE_RANGEING_
    int         n_reject;           // #of of rejects in a row
    int         n_accept;           // #of accepts in a row
    int         acc_pos;            // last accept was a positive change
    double      range;              // Change-range
#endif

public:
    sa_parameter (const_parameter *p_ptr);

    void save_best() {best_val = p_ptr->get_mapped_value();}
                                    // Shall save the best value encountered so far
    void restore_best() {p_ptr->set_mapped_value(best_val);}
                                    // Shall restore this paramter to the best value encountered so far
                                    
    friend class sim_anneal;
};

extern sim_anneal *sim_anneal_ptr;  // The SA machinery (if used, otherwise nullptr)

#endif

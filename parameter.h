//
// Declaration of the paparemer class
//
#pragma once

#define _ADAPTIVE_RANGEING_                 // If defined, the change scale is adjusted
                                            // individually for each parameter
const double MIN_RANGE = 1.0e-6;            // Min. parameter change range
extern const char*  EVAL_PARAMETER;         // The evaluation function parameter

#define REJ_MAX     20                      // #of rejects before reducing change scale
                                            //   in simulated annealing

#define MAX_RESET   1000                    // patience factor: if opt wasn't improved in this
                                            // many steps, reset to optimum

#define MAX_CHANGE_TRIES 50                 // #of times a parameter change can be rejected

#define EPS         1.0e-10                 // Greedy search cut-off

///////////////////////////////////////////////////////////////////////////////////////////////////

struct zel_element {                        // Zero-expression list element
    expression  *e_ptr;                     // Pointer to expression
    zel_element *next;                      // Next or nullpter
};

struct fit_element {                        // lsq_fit list element
    lsq_fit_function  *fit_ptr;             // Pointer to lsq fit object
    fit_element *next;                      // Next or nullpter
};

struct param_iterator {                     // A parameter iterator
    zel_element *zel_ptr;                   // stuff that needs to be zeroed
    fit_element *fit_ptr;                   // lsq_fits that need to be computed
    parameter   *param_ptr;                 // Pointer to the parameter
};

struct an_sched {                           // annealing schedule
    double          temp;                   // Starting temp for this phase
    unsigned        n_steps;                // #of steps
};

class parameter {
    //
    //  Common elements (housekeeping, etc)
    //
    char        *name;
    enum {scan, assignment, bin_search, sim_anneal, undefined} type;   // Parameter type
   
    parameter      *next;                   // List pointer
    static parameter *root;                 // Anchor
    static unsigned nesting_level;          // How many iterator type parameters are nested
    static param_iterator iterator[max_param_iterators]; // Note: innermost loop is at level 0
                                            
    //
    // Variables used by the scan function:
    //
    double          min_value, max_value;   // Range of values to step through
                                            // Note: these are also used by bin search and simulated annealing
    unsigned        n_steps;                // #of n_steps
    double          d_value;                // Step size
    unsigned        step_cntr;              // Step counter
    
    //
    // assignment type parameter
    //
    expression      *expr;                  // Pointer to expression
     
    //
    // Simulated annealing state:
    //
    static unsigned enable_sim_anneal;      // If set != 0, simulated annealing is in use
    double          cur_value;              // Current value
    double          best_val;               // Best value for this parameter encountered so far
    static unsigned	have_best;              // is !0 if best value is defined
#ifdef _ADAPTIVE_RANGEING_
    int             n_reject;               // #of of rejects in a row
    int             n_accept;               // #of accepts in a row
    int             acc_pos;                // last accept was a positive change
    double          range;                  // Change-range
#endif
    static vector<parameter*> p_ptr;        // Collection of parameters Subject to anealing
    static unsigned changed_flag;           // Change-machinery: indicate validity    
    static parameter *changed_param;        // Which was changed last
    static double	changed_value;          // What was its pre-change value
    static parameter *eval_ptr;             // The evaluation function
    static parameter *reject_ptr;           // Optional reject function
    
    //
    // The working variables of the simulated annealing optimization function
    //
    static double   t0, t1, dt;             // annealing temperature
    static double   old_eval, new_eval, min_eval;
    static unsigned n_sim_anneal;           // Counts the smulated annealing cycles
    static int      n_rej, t_max;
    static double   d_scale;                // scales the range of changes
    static unsigned out_cnt;                // Used to control the volume of the log-output
    static vector<an_sched> as;             // The anealing schedule
    static unsigned cur_sched;              // Which schedule is current
    static unsigned sched_stp_cntr;         // Current step
    static FILE     *sa_log;                // Optinional log file
    
    //
    // Internal functions
    //
    void            i_reset() {step_cntr = 0;}; // Restarts an iterator type
    int             i_next(); 
    void            default_init(char *nm); // Initialize a generic parameter class member

    void            print_name(FILE *fp);   // Output functions
    void            print_value(FILE *fp);
    
    static unsigned sim_anneal_cycle();     // Perform one simulated annealing update cycle
    static void     restore_best();
    static void     save_best();
    static void     change_param(double c_range);
    static void     unchange_param();

public:
    parameter(char *nm, double v_min, double v_max, unsigned n_steps);   // Creates a scan-type parameter
    
    parameter(char *nm, double v_min, double v_init, double v_max, parameter *eval_ptr, parameter *reject_ptr);
                                            // creates a simulated annealing type
    
    parameter(char *nm, class expression *expr); // creates an assignment type parameter
    
    static void reset();                    // Go back to initial state
    static int  advance(unsigned &level);   // Advance to next value combination, returns != 0 when exhausted
    
    static unsigned list_names(FILE *fp);   // adds names to file (preceeded by space, followed by nothing)
                                            // returns number of names listed
    static void list_c_val(FILE *fp);       // list the current values, like above

    static parameter *find_parameter(const char *nm);   // find parameter by its name (symbol)
    double      get_cur_value();
    
    static void enqueue_zero_op(expression *e_ptr); // Summation op to be zeroed
    static void enqueue_fit_function(lsq_fit_function *f_ptr); // Fits to be updated
    
    static unsigned open_sa_log_file(char *nm); // Opens a log file for the simulated annealing feature
    static unsigned read_sa_schedule(char *nm); // Read the simulated annealing schedule
    
    static unsigned sim_anneal_in_use() {return enable_sim_anneal;}
    static void print_sa_results(FILE *sar_fp);
};

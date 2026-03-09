//
// Declaration of the paparemer class
//
#pragma once

#include <math.h>
#include <stdio.h>
#include <vector>

class expression;                           // Forward declaration (instead of pulling in the expression.h file)

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

struct an_sched {                           // annealing schedule
    double          temp;                   // Starting temp for this phase
    unsigned        n_steps;                // #of steps
};

class parameter {
protected:
    //
    //  Common elements (housekeeping, etc)
    //
    char        *name;
    enum {constant, assignment, scan, bin_search, undefined} type;   // Parameter type
   
    static std::vector<parameter *> parameters;  // Collection of all parameters
    static unsigned nesting_level;          // How many iterator type parameters are nested
    
    // attributes:
    unsigned        no_print:1;             // If set, do not include in summary
    unsigned        tunable:1;              // If set, will be subject to simulated annealing, BO, etc.
    unsigned        log_map:1;              // Uses logarithmic mapping of tuning range

    // Common:
    double          min_value, max_value;   // Value range
    double          value;                  // The current value
    
    // Internal functions:
    void            print_name(FILE *fp);   // Output functions
    void            print_value(FILE *fp);

    parameter(char *nm,
              double v_min = __DBL_MAX__, double v_max = -__DBL_MAX__,
              unsigned npr = 0, unsigned tun = 0, unsigned l_map = 1);

public:
    virtual ~parameter() = default;         // Recommended C++ practice, alas not needed here
    
    double get_cur_value() {return value;}; // returns the current value of the parameter
    
    virtual void    update() {};            // This must be called *before* the simulation step
                                            // By default do nothing
                                            
    virtual void    initialize() {};        // Initializer: for looping parameters, like scan, etc.
                                            // This will be called before the loop is executed for 
                                            // the firest time
                                            
    virtual unsigned next() {return 1; };   // iterator: must be called after each loop iteration
                                            // Returns 0 when the loop is not complete and 1 otherwise
                                            
    static parameter *find_parameter(const char *nm);   // find parameter by its name (symbol)
    
    static unsigned list_names(FILE *fp);   // adds names to file (preceeded by space, followed by nothing)
                                            // returns number of names listed
    static void list_c_val(FILE *fp);       // list the current values, like above
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Derived classes 
//

class const_parameter : public parameter {
public:
    const_parameter (char *nm,
                    double val,
                    double v_min = __DBL_MAX__, double v_max = -__DBL_MAX__,
                    unsigned npr = 0, unsigned tun = 0, unsigned l_map = 1);
    
    double get_mapped_value();              // Returns the value in [0,1] provided that a range is defined
    void   set_mapped_value(double m_val);  // Does the reverse
};

class expr_parameter : public parameter {
    expression      *expr;                  // Pointer to expression
public:
    expr_parameter (char *nm,
                    expression  *expr,
                    double v_min = __DBL_MAX__, double v_max = -__DBL_MAX__,
                    unsigned npr = 0);
    
    void            update() override;
};
    
class scan_parameter : public parameter {

    unsigned        n_steps;                // #of n_steps
    double          d_value;                // Step size
    unsigned        step_cntr;              // Step counter

public:
    scan_parameter (char *nm,
                    double v_min, double v_max, unsigned n_steps,
                    unsigned npr = 0, unsigned l_map = 0);   // Creates a scan-type parameter

    void initialize() override {step_cntr = 0;};
    void update() override;
    unsigned next() override {return ++step_cntr >= n_steps; };
};

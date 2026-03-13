//
// Declaration of the paparemer class
//
#pragma once

#include <math.h>
#include <stdio.h>
#include <vector>
using namespace std;

class expression;                           // Forward declaration (instead of pulling in the expression.h file)
class sim_anneal;                           // dito
class const_parameter;                      //  ..

///////////////////////////////////////////////////////////////////////////////////////////////////

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
    unsigned        locked:1;               // Used to detect dependency cycles

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
    
    virtual double get_cur_value() {return NAN;}; // value depends on parameter type
                                            
    virtual void    initialize() {};        // Initializer: for looping parameters, like scan, etc.
                                            // This will be called before the loop is executed for 
                                            // the firest time. Non-looping parameters simply do nothing here.
                                            
    virtual unsigned next() {return 1; };   // iterator: must be called after each loop iteration
                                            // Returns 0 when the loop is not complete and 1 otherwise
                                            // For a non-looping parameter, this is a No-op.
                                            
    static parameter *find_parameter(const char *nm);   // find parameter by its name (symbol)
    
    static unsigned list_names(FILE *fp);   // adds names to file (preceeded by space, followed by nothing)
                                            // returns number of names listed
    static void list_c_val(FILE *fp);       // list the current values, like above
    
    static void sa_p_export(sim_anneal *sa_ptr); // Export all tunable paramters to the
                                            // simulated annealing subsystem
    static void bo_export(vector<const_parameter *> &tunable_params); // same for the BO subsystem
    
    static void save_result(FILE *of);      // For optimizer runs, save the results

    char *get_name() {return name;};
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
    
    void   print_self(FILE *fp);            // Prints a version of itself

    double get_cur_value() override {return value;}; // just returns the value
};

class expr_parameter : public parameter {
    expression      *expr;                  // Pointer to expression
public:
    expr_parameter (char *nm,
                    expression  *expr,
                    double v_min = __DBL_MAX__, double v_max = -__DBL_MAX__,
                    unsigned npr = 0);
    
    double get_cur_value() override;        // Returns the current value
};
    
class scan_parameter : public parameter {

    unsigned        n_steps;                // #of n_steps
    double          d_value;                // Step size
    unsigned        step_cntr;              // Step counter

public:
    scan_parameter (char *nm,
                    double v_min, double v_max, unsigned n_steps,
                    unsigned npr = 0, unsigned l_map = 0);   // Creates a scan-type parameter

    void            initialize()    override {step_cntr = 0;};
    double          get_cur_value() override;
    unsigned        next()          override {return ++step_cntr >= n_steps; };
};

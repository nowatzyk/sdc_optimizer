//
// Declaration of the parameter class
//
#pragma once

#include <math.h>
#include <stdio.h>
#include <vector>
using namespace std;

class expression;                           // Forward declaration
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
    unsigned        locked:1;              // Used to detect dependency cycles
    unsigned        warning_issued:1;       // Used to prevent excessive warnings

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
    virtual ~parameter() = default;         // Recommended C++ practice
    
    virtual double get_cur_value() {return NAN;}; // value depends on parameter type
                                            
    virtual void    initialize() {};        // Initializer: for looping parameters, like scan, etc.
                                            // This will be called before the loop is executed for 
                                            // the firest time. Non-looping parameters simply do nothing here.
                                            
    virtual unsigned next() {return 1; };   // iterator: must be called after each loop iteration
                                            // Returns 0 when the loop is not complete and 1 otherwise
                                            // For a non-looping parameter, this is a No-op.

    //
    // print_to_file() -- emit the parameter's current value to a file stream.
    //
    // The default implementation formats get_cur_value() as "%.12lg", which
    // is the correct behaviour for all numeric parameters.
    //
    // pwl_parameter overrides this to emit a multi-token PWL string instead of
    // a single number. spice_elements::print() calls this method so that the
    // SPICE deck gets either a number or a PWL string transparently.
    //
    virtual void print_to_file(FILE *fp);

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
    
    double map_01_to_parm(double val);      // Map a value in [0,1] space to the perameter value
    double map_parm_to_01(double val);      // Do the reverse of above
    
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
    unsigned        next()          override // Keep the state of the final iteration
                    {   if ((step_cntr+1) >= n_steps) return 1;
                        else {step_cntr++;            return 0; };
                    };
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// pwl_parameter -- emits a SPICE PWL waveform string from a named time pattern.
//
// The pattern defines a sequence of pulse centre times (evaluated lazily from
// expressions at deck-assembly time). For each centre time t_c the parameter
// emits four PWL points that form one rectangular pulse:
//
//   t_c - t_width/2 - t_rise   v_low    <- start of rise edge
//   t_c - t_width/2            v_high   <- end of rise edge / start of flat top
//   t_c + t_width/2            v_high   <- end of flat top / start of fall edge
//   t_c + t_width/2 + t_fall   v_low    <- end of fall edge
//
// The full output string begins with "0 <v_low>" (baseline at t=0) and ends
// with the v_low point after the last pulse, making a valid SPICE PWL vector.
//
// All five control values (t_rise, t_fall, t_width, v_high, v_low) are
// expression* pointers evaluated at deck-assembly time, consistent with the
// existing parameter expression system.  This allows them to reference tunable
// parameters so that waveform timing can be explored by SA/BO.
//
// Usage in circuit file:
//   *Pragma parameter inp_a pwl out_a rise 0.2p width 1p fall 0.2p high 1m low 0
//   Vinput_a N1a gnd pwl(~inp_a~)
//
// The five keyword arguments may appear in any order.
//

class time_pattern;                         // forward declaration (defined in nodes_of_interest.h)

class pwl_parameter : public parameter {
    time_pattern    *pattern;               // The source time pattern (pulse centres)
    expression      *t_rise;                // Rise time (seconds)
    expression      *t_fall;                // Fall time (seconds)
    expression      *t_width;               // Pulse width (seconds, flat-top duration)
    expression      *v_high;                // High voltage level
    expression      *v_low;                 // Low voltage (baseline) level

public:
    pwl_parameter(char *nm,
                  time_pattern *pat,
                  expression   *t_rise,
                  expression   *t_fall,
                  expression   *t_width,
                  expression   *v_high,
                  expression   *v_low);

    //
    // print_to_file() -- generates and emits the full PWL string.
    // Called by spice_elements::print() during deck assembly.
    //
    void print_to_file(FILE *fp) override;

    //
    // get_cur_value() is not meaningful for a PWL parameter since the output
    // is a multi-point string, not a scalar. Returns NAN to make misuse visible.
    //
    double get_cur_value() override { return NAN; }
};

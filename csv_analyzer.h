//
// Data structures for the CSV analyzer
//
#pragma once

const unsigned max_ln_length = 100000;      // We expect some very long lines...
const unsigned max_ts = 128;                // Max number of time series

const double threshold_frac = 0.6;          // Peak search threshold, 75% of max 
const double threshold_hyst = 0.01;         // Hysteresis to avoid noise induced false peak locations

const unsigned max_pks = 1024;              // needed to allocate arrays for peaks:
                                            // Peaks are located first and then analyzed

const unsigned parameter_name_lim = 64;     // Max symbol name length for parameters
const unsigned spice_src_max_char = 1024;   // Max line length in spice source file
const char spice_escape = '~';              // Escape character in spice decks for parameter substitutions

const unsigned n_slope_pks = 4;             // #of data points to be used in slope determination

const unsigned ts_name_lim =64;             // Dito for TS names

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Spice deck parser stuff
//



struct t_of_peaks {
    unsigned    n_pks;                      // #of peaks
    double      *t_pks;                     // Time of peaks
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Class definitions
//

class nodes_of_interest {
    static nodes_of_interest *root;         // Start of NOI list
    nodes_of_interest *next;                // Pointer to the next list element
    static unsigned n_noi;                  // #of nodes of interest
    
    char        *name;                      // Name of this NOI
    unsigned    va_index;                   // Index into value array
public:
    static int  find(const char *name);     // Find the named NOI and return its VA
    nodes_of_interest(const char *name);    //create a new node of interest (Note: doesn't check for duplicate names: you need to do that before creaing a new NOI!)
    static unsigned get_n_noi() {return n_noi;}; // Just returns the current number of NOI's
};

class parameter {
    char        *name;
    enum {scan, expression, bin_search, sim_anneal} type;   // Parameter type
                                            
    //
    // Variables used by the scan function:
    //
    double      min_value, max_value;       // Range of values to step through
    unsigned    n_steps;                    // #of n_steps
    double      d_value;                    // Step size
    unsigned    step_cntr;                  // Step counter
    
    parameter      *next;                   // List pointer
    static parameter *root;                 // Anchor
    
    void         i_reset() {step_cntr = 0;};    // Duh!
    int          i_next(); 

    void         print_name(FILE *fp);
    void         print_value(FILE *fp);
    
public:
    parameter(char *nm, double v_min, double v_max, unsigned );   // Creates a scan-type parameter
    
    static void reset();                    // Go back to initial state
    static int  advance();                  // Advance to next value combination, returns != 0 when exhausted
    static void list_names(FILE *fp);       // adds names to file (preceeded by space, followed by nothing)
    static void list_c_val(FILE *fp);       // list the current values, like above

    static parameter *find_parameter(const char *nm);   // find parameter by its name (symbol)
    double      get_cur_value() {return (min_value + (double) step_cntr * d_value);};
};

struct spe_text {                           // Spice element: Text
    char        *text;                      // Pointer to the 0-termninated text (belt&suspenders...)
};

struct spe_param {                          // Spice element: Parameter
    parameter   *param;
};

class spice_elements  {
    enum {text, parameter, new_line} se_type;         // Type of element
    spice_elements   *next;                 // Pointer to the next element
    union {
        spe_text    txt;                    // This is text element
        spe_param   par;                    // This is a parameter element
    };
    
public:
    spice_elements(char *txt);              // Construct a new text element
    spice_elements(class parameter *par);   // Construct a new parameter element
    spice_elements();                       // Construct a new line element
    void        add_next(spice_elements *nxt) {next = nxt;};
    spice_elements *get_next() {return next;};
    
    void print(FILE *of);
    
    int is_nl() {return (se_type == new_line) ? 1 : 0;};
};

class spice_deck {
    spice_elements  *first, *last;          // Pointer to the actual code

public:
    spice_deck();                           // constructor
    
    void        add_line(char *txt);        // Adds a line of text
    void        add_param_ref(parameter *par); // adds a reference to a parameter
    void        add_nl();                   // adds a new-line character to spice deck
    
    int read_cir_file(const char *fn);      // Reads a spice deck from a file into memory
    int write_cir_file(const char *fn);     // writes the spice deck to a file
};

class time_series {
    char        *name;                      // Name of this column
    double      *data;
    unsigned    n_data;
    unsigned    max_data;
    double      v_min, v_max;
    
    static unsigned n_init;                 // #of values for initial allocation
    
    static time_series *root;               // Root of list
    time_series *next;
    
public:
    time_series(char *name);
    ~time_series();
    
    void        reset();                    // Reset the time series
    void        add_datum(double);
    double      get_max() {return v_max;};
    double      get_min() {return v_min;};
    char        *get_name() {return name;};
    double      get_val(unsigned i) {assert(i < n_data); return data[i];};
    
    unsigned    get_n() {return n_data;};
    
    double      peak_search(unsigned from, unsigned to, double thr, double eps, unsigned &end);
    
    static void set_default_length(unsigned n);  // Just set a default length
};

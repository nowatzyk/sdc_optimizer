//
// Data structures for the CSV analyzer
//
#pragma once

const unsigned max_ln_length = 100000;      // We expect some very long lines...
const unsigned max_ts = 256;                // Max number of time series

const double threshold_frac = 0.6;          // Peak search threshold, 75% of max 
const double threshold_hyst = 0.01;         // Hysteresis to avoid noise induced false peak locations

const unsigned max_pks = 1024;              // needed to allocate arrays for peaks:
                                            // Peaks are located first and then analyzed

const unsigned max_param_iterators = 16;    // Should be plenty
const unsigned spice_src_max_char = 1024;   // Max line length in spice source file
const char spice_escape = '~';              // Escape character in spice decks for parameter substitutions

const double edge_search_min_chg = 0.6667;  // For an edge search on a phase, this is the min change required
                                            // to be considered an edge (in fractions of a turn, = 2PI phase)
const double edge_search_t_win = 10.0e-12;  // max edge transition time 10 pico-seconds
const double edge_search_slope_frac = 1.0; // Defines the fraction of the min_chg (see above) that is used to
                                            // to determine the transition point

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Spice deck parser stuff
//

struct t_of_peaks {
    unsigned    n_pks;                      // #of peaks
    double      *t_pks;                     // Time of peaks
};

struct f1_table {                           // Function pointer f1_table
    const char  *name;                      // Name of this function (or nullptr)
    double      (*func1_ptr) (double x);    // Pointer to the funcion
};

struct f2_table {                           // Function pointer f2_table
    const char  *name;                      // Name of this function (or nullptr)
    double      (*func2_ptr) (class nodes_of_interest *noi_ptr, double x);    // Pointer to the funcion
};


///////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Class definitions
//

class nodes_of_interest {
    static nodes_of_interest *root;         // Start of NOI list
    nodes_of_interest *next;                // Pointer to the next list element
    static unsigned n_noi;                  // #of nodes of interest
    
    const char  *name;                      // Name of this NOI
    int         col_index;                  // Which column of the Josim output (or -1)

public:
    static nodes_of_interest *find(const char *name);     // Find the named NOI and return its VA
    nodes_of_interest(const char *name);    // create a new node of interest
                                            // Note: doesn't check for duplicate names:
                                            // you need to do that before creaing a new NOI!)
    
    static unsigned get_n_noi() {return n_noi;}; // Just returns the current number of NOI's
    int         get_col_index() {return col_index;};
    unsigned    set_col_index(unsigned n);  // Set the column index, return != 0 if there is a problem
    const char *get_name() {return name;};
    static nodes_of_interest *find_undef(); // Returns point to the first NOI with col_index < 0
    
    int print(char *f_name);
    static void print_all();                // Saves all NOI's to files
};

// Some related function definitione:
double locate_peak(nodes_of_interest *noi_ptr, double x);
double count_peaks(nodes_of_interest *noi_ptr, double x);
double locate_rise(nodes_of_interest *noi_ptr, double x);
double locate_fall(nodes_of_interest *noi_ptr, double x);
double locate_edge(nodes_of_interest *noi_ptr, double x);
double count_edges(nodes_of_interest *noi_ptr, double x);

class parameter;                            // Forward declaration so that it can be used in the expression class
    
enum    exp_type {                          // Type of expression
                constant,
                sum_func,
                function1,                  // Function with 1 argument
                function2,                  // function with 2 arguments
                p_reference,
                addition,
                subtraction,
                multiplication,
                division,
                test_select
                };

class expression {
            
    exp_type    type;                       // Type of this expression element
            
    expression  *l_arg, *r_arg;             // Pointers to the left/right arguments
    expression  *t_arg;                     // A 3rd argument, needed for thest-select op
    double      value;                      // the value (if this is a constant or a sum)
    double      (*func1_ptr)(double x);     // Function1 pointer
    double      (*func2_ptr)(class nodes_of_interest *noi_ptr, double x);     // Function2 pointer
    parameter   *p_ptr;                     // Parameter pointer
    nodes_of_interest *noi_ptr;             // NOI pointer
    // Note: TBD save some space and make this a union
    
public:
    expression  (double x);                 // creates a constant expression
    expression  (double (*f_ptr)(double x), expression *arg_ptr);   // creates a function expression with 1 arg 
    expression  (double (*f_ptr)(class nodes_of_interest *noi_ptr, double x), 
                 class nodes_of_interest *noi_ptr,
                 expression *arg_ptr);
    expression  (expression *la_ptr, exp_type et, expression *ra_ptr); // creates a prmitive op
    expression  (parameter *p_ptr);         // Creates a parameter reference
    expression  (expression *test_ptr, expression *true_arg, expression *false_arg);
    
    double      get_value();                // Returns the current value of this expression
    void        zero_sum() {assert(type == sum_func); value = 0.0;}   // Should only be used on sums
};

struct zel_element {                        // Zero-expression list element
    expression  *e_ptr;                     // Pointer to expression
    zel_element *next;                      // Next or nullpter
};

struct param_iterator {                     // A parameter iterator
    zel_element *zel_ptr;                   // stuff that needs to be aeroed
    parameter   *param_ptr;                 // Pointer to the parameter
};

class parameter {
    char        *name;
    enum {scan, assignment, bin_search, sim_anneal} type;   // Parameter type
                                            
    //
    // Variables used by the scan function:
    //
    double      min_value, max_value;       // Range of values to step through
    unsigned    n_steps;                    // #of n_steps
    double      d_value;                    // Step size
    unsigned    step_cntr;                  // Step counter
    
    //
    // assignment type parameter
    //
    expression  *expr;                      // Pointer to expression
    
    parameter      *next;                   // List pointer
    static parameter *root;                 // Anchor
    static unsigned nesting_level;          // How many iterator type parameters are nested
    static param_iterator iterator[max_param_iterators]; // Note: innermost loop is at level 0
    
    void         i_reset() {step_cntr = 0;};    // Duh!
    int          i_next(); 

    void         print_name(FILE *fp);
    void         print_value(FILE *fp);
    
public:
    parameter(char *nm, double v_min, double v_max, unsigned );   // Creates a scan-type parameter
    parameter(char *nm, class expression *expr); // creates an assignment type parameter
    
    static void reset();                    // Go back to initial state
    static int  advance();                  // Advance to next value combination, returns != 0 when exhausted
    
    static void list_names(FILE *fp);       // adds names to file (preceeded by space, followed by nothing)
    static void list_c_val(FILE *fp);       // list the current values, like above

    static parameter *find_parameter(const char *nm);   // find parameter by its name (symbol)
    double      get_cur_value();
    
    static void enqueue_zero_op(expression *e_ptr); // Summation op to be zeroed
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
    double      *data;
    unsigned    n_data;
    unsigned    max_data;
    double      v_min, v_max;
    
    static unsigned n_init;                 // #of values for initial allocation
    
    static time_series *root;               // Root of list
    time_series *next;
    
public:
    time_series();
    ~time_series();
    
    void        reset();                    // Reset the time series
    void        add_datum(double);
    double      get_max() {return v_max;};
    double      get_min() {return v_min;};
    double      get_val(unsigned i) {assert(i < n_data); return data[i];};
    
    unsigned    get_n() {return n_data;};
    
    double      peak_search(unsigned from, unsigned to, double thr, double eps, unsigned &end);
    double      edge_search(unsigned from, unsigned to, double min_chg, double t_window,
                            unsigned &e_type, unsigned &end);
    
    void        print_all(FILE *fp);        // Dump all data to a file
    
    static void set_default_length(unsigned n);  // Just set a default length
};

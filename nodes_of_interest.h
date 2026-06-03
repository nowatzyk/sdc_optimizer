#pragma once

#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

#include "expression.h"

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// time_series -- stores one column of JoSIM output and supports peak/edge searching.
//

class time_series {
    double      *data;
    unsigned    n_data;
    unsigned    max_data;
    double      v_min, v_max;

    static unsigned     n_init;         // initial allocation size
    static time_series *root;
    time_series        *next;

public:
    time_series();
    ~time_series();

    void        reset();
    void        add_datum(double v);

    double      get_max()           { return v_max; }
    double      get_min()           { return v_min; }
    double      get_val(unsigned i) { assert(i < n_data); return data[i]; }
    unsigned    get_n()             { return n_data; }

    double      peak_search(unsigned from, unsigned to,
                            double thr, double eps, unsigned &end);
    double      edge_search(unsigned from, unsigned to,
                            double min_chg, double t_window,
                            unsigned &e_type, unsigned &end);

    void        print_all(FILE *fp);

    static void set_default_length(unsigned n);
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// nodes_of_interest -- tracks a named JoSIM output column and provides
// expression-callable analysis functions.
//

class nodes_of_interest {
    static vector<nodes_of_interest *> all_noi_s;

    const char *name;
    int         col_index;                  // column index in JoSIM output, or -1 if not yet matched
    
    static unsigned is_valid_flag;          // A flag that indicates if a result is available:
                                            // It is cleared prior to a simulation run and set 
                                            // once the JoSIM completed successfully

public:
    nodes_of_interest(const char *name);

    static nodes_of_interest *find(const char *name);
    static nodes_of_interest *find_undef(); // first NOI with col_index < 0
    static unsigned           get_n_noi()       { return all_noi_s.size(); }

    int         get_col_index()                 { return col_index; }
    unsigned    set_col_index(unsigned n);       // returns != 0 if there is a conflict
    const char *get_name()                      { return name; }

    int         print(char *f_name);
    static void print_all();
    
    static unsigned is_valid() {return is_valid_flag; };
    static void set_validity(unsigned isv) {is_valid_flag = (isv != 0); };
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Time Pattern object
//
//  This is used to define time intervals for matching the output of a circuit against an expected
//  behavior:
//
//  ---t0----t1-------t2------t3------...-->  time-line
//
//  The Pattern object is a vector of times [t0, t1,... tn]. A good match consists of no event
//  in [0,t0), followed by one event in [t0,t1), no event in [t1,t2), ....
//  Events can be a peak, a {rising,falling,aby}edge. Peaks are usually used on voltage nodes,
//  while edges make more sense on phase.
//  t0 can be 0.0, so that the first interval is empty, which essentially inverts the sense.
//
class time_pattern {
    const char          *name;              // Name of this pattern
    nodes_of_interest   *noi_ptr;           // NOI to which this pattern applies
    static vector<time_pattern *> all_tps;  // The collection of time patterns

    // Expression-valued time storage (replaces the pre-evaluated double vector
    // for patterns defined via the new expression syntax).
    //
    // Each entry is (expression*, is_relative):
    //   is_relative=0: expression evaluates to an absolute time
    //   is_relative=1: expression evaluates to a delta added to the previous time
    //
    // add_time(double) wraps the value in a constant expression for backward compat.
    // get_times() evaluates all expressions, resolves relative deltas, validates
    // monotonicity, and returns the resulting absolute times.
    
    struct time_expr_entry {
        expression  *expr;
        unsigned     is_relative;
    };
    vector<time_expr_entry>  time_exprs;    // expression-valued time list

public:
    void   add_time_expr(expression *e, unsigned is_relative); // add expression entry
    vector<double> get_times() const;       // evaluate and return absolute times
    
    const char *get_name() const {return name; };

    time_pattern(char *nm, nodes_of_interest *np);
    void add_time(double tt);               // Adds a transition time
    static time_pattern *find(const char *name); // Finds a pattern via name loop-up (linear search)
    nodes_of_interest   *get_noi_ptr() {return noi_ptr;};
    
    unsigned cnt_missmatches(vector<double> &t_peaks);
};


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Analysis functions -- callable as function2 expression callbacks.
// obj_ptr is a nodes_of_interest*.
//

double locate_peak  (void *obj_ptr, double x);
double count_peaks  (void *obj_ptr, double x);
double locate_rise  (void *obj_ptr, double x);
double locate_fall  (void *obj_ptr, double x);
double locate_edge  (void *obj_ptr, double x);
double count_edges  (void *obj_ptr, double x);

double match_peaks  (void *obj_ptr, double x);
double match_r_edg  (void *obj_ptr, double x);
double match_f_edg  (void *obj_ptr, double x);
double match_a_edg  (void *obj_ptr, double x);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Globals owned by nodes_of_interest.cpp but shared with main / loop_complex.
//

extern vector<time_series *> josim_out_columns;
extern unsigned     n_josim_runs;
extern double       sim_time_start;
extern double       sim_time_incr;

double sim_time(double n);      // returns simulation time for fractional row n

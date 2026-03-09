#pragma once

#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>
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

    const char  *name;
    int          col_index;     // column index in JoSIM output, or -1 if not yet matched

public:
    nodes_of_interest(const char *name);

    static nodes_of_interest *find(const char *name);
    static nodes_of_interest *find_undef();     // first NOI with col_index < 0
    static unsigned           get_n_noi()       { return all_noi_s.size(); }

    int         get_col_index()                 { return col_index; }
    unsigned    set_col_index(unsigned n);       // returns != 0 if there is a conflict
    const char *get_name()                      { return name; }

    int         print(char *f_name);
    static void print_all();
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

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Globals owned by nodes_of_interest.cpp but shared with main / loop_complex.
//

extern vector<time_series *> josim_out_columns;
extern unsigned     n_josim_runs;
extern double       sim_time_start;
extern double       sim_time_incr;

double sim_time(double n);      // returns simulation time for fractional row n

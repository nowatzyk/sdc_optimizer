//
// This is the interface to the parser, which is using plain C
//
#ifndef _PARSER_INTERFACE_
#define _PARSER_INTERFACE_

#include <stdio.h>

extern FILE *yyin;                      // fd for spice source file
int yyparse();                          // The Bison parser (in plain C)

struct units {
    char        c;                      // character representation
    double      val;                    // Value
};

struct range_pair {
    double      from;
    double      to;
};

struct f1_table {                           // Function pointer f1_table
    const char          *name;              // Name of this function (or nullptr)
    double              (*func1_ptr) (double x); // Pointer to the funcion
};

enum function_type {noi_type, lsq_fit_type, none_of_the_above};

struct f2_table {                           // Function pointer f2_table
    const char          *name;              // Name of this function (or nullptr)
    double              (*func2_ptr) (void *obj_ptr, double x);
                                            // Pointer to the associated object and argument
    enum function_type  f_type;             // used to find the correct object pointer
};

extern struct units unit_table[];       // Common SI scale identifiers

extern unsigned  yy_n_parse_err;        // #of errors during parsing the circuit file

///////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Parser interface functions
//

void add_line_to_spice_deck(char *string);
void add_subst_to_spice_deck(char *p_name);
void add_new_line_to_spice_deck();

void define_tran(double t_incr, double t_stop, double t_start, double dT_max);

void define_param_scan(char *name, void *rng, double n_steps, unsigned flags);
void define_param_constant(char *name, double value, void *rng, unsigned flags);
void define_param_expression(char *name, void *expr, void *rng, unsigned flags);

void define_sim_anneal(char *log_file_nm);
void define_add2SA_sched(double temp, double n_iter);
void define_monitor(char *name);
void define_snapshot(char *name, double start, double frequency);
void define_lsq_fit(char *name, double order, void *x, void *y);
void define_bo(double n_itr);

void *define_range(double from, double to);

void *define_add(void *x, void *y);
void *define_sub(void *x, void *y);
void *define_mul(void *x, void *y);
void *define_div(void *x, void *y);
void *define_const(double x);
void *define_ref(char *name);
void *define_function1(char *name, void *x);
void *define_function2(char *name, char *ts_name, void *y);
void *define_test_sel(void *t, void *a, void *b);
void *define_eq(void *x, void *y);
void *define_ne(void *x, void *y);
void *define_gt(void *x, void *y);
void *define_ge(void *x, void *y);
void *define_le(void *x, void *y);
void *define_lt(void *x, void *y);
void *define_not(void *x);

#endif

//
// This is the interface to the parser, which is using plain C
//
#ifndef _PARSER_INTERFACE_
#define _PARSER_INTERFACE_

#include <stdio.h>

extern FILE *yyin;                          // fd for spice source file
int yyparse();                              // The Bison parser (in plain C)

struct units {
    char        c;                          // character representation
    double      val;                        // Value
};

struct range_pair {
    double      from;
    double      to;
};

struct f1_table {                           // Function pointer f1_table
    const char          *name;              // Name of this function (or nullptr)
    double              (*func1_ptr) (double x); // Pointer to the funcion
};

enum function_type {noi_type, lsq_fit_type, tm_pattern, none_of_the_above};

struct f2_table {                           // Function pointer f2_table
    const char          *name;              // Name of this function (or nullptr)
    double              (*func2_ptr) (void *obj_ptr, double x);
                                            // Pointer to the associated object and argument
    enum function_type  f_type;             // used to find the correct object pointer
};

//
// p_time_element -- one element of a time pattern list.
//
// The time is now stored as an expression pointer so that pattern times can be
// parameterised and evaluated at deck-assembly time.
//
// is_relative == 0: time_expr evaluates to an absolute time
// is_relative == 1: time_expr evaluates to a delta added to the previous time
//
// Note: the old representation used a negative double to signal relative times.
//       This is replaced by an explicit flag so that expressions can produce
//       any numeric value without ambiguity.
//
struct p_time_element {
    void               *time_expr;          // expression* -- evaluated at deck-assembly time
    unsigned            is_relative;        // 0 = absolute, 1 = relative (+delta)
    struct p_time_element *next;
};

extern struct units unit_table[];           // Common SI scale identifiers

extern unsigned  yy_n_parse_err;            // #of errors during parsing the circuit file

//
// bo_attr -- BO pragma attribute node, built into a linked list by the parser.
// One node per keyword seen in the baysian_opt pragma line.
// Passed as void* through the Bison grammar (C-compatible) to define_bo().
//

typedef enum {
    BO_ATTR_MARGIN,         // enable robustness / margin analysis mode
    BO_ATTR_BINARY,         // submode: binary pass/fail oracle (+1/-1)
    BO_ATTR_GRADIENT,       // submode: continuous oracle with gradient
    BO_ATTR_PROBABILISTIC,  // submode: stochastic oracle (error rate)
    BO_ATTR_N_RAYS,         // number of rays from x* for boundary search
    BO_ATTR_N_BRACKET,      // bracket search steps per ray
    BO_ATTR_N_BISECT,       // bisection steps per ray
    BO_ATTR_THRESHOLD,      // pass/fail boundary value (default 0.0)
    // bef_plan controls for the binary ellipsoid fit subsystem:
    BO_ATTR_BEF_BUDGET,     // total JoSIM simulation budget for the ellipsoid fit
    BO_ATTR_BEF_ITER,       // number of outer fit iterations (default 5)
    BO_ATTR_BEF_PROBES      // probes per ray in the initial exploration (default 16)
} bo_attr_type;

struct bo_attr {
    bo_attr_type    type;
    double          value;
    struct bo_attr *next;
};

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
void define_param_pwl(char *name, char *pattern_name,
                      void *t_rise, void *t_fall, void *t_width,
                      void *v_high, void *v_low);
void *define_range(double from, double to);

void define_sim_anneal(char *log_file_nm);
void define_add2SA_sched(double temp, double n_iter);

void define_monitor(char *name);
void define_snapshot(char *name, double start, double frequency);
void define_lsq_fit(char *name, double order, void *x, void *y);

void define_bo(double n_itr, void *attr_list);
void *bo_attr_cat(void *list, void *node);
void *define_bo_attr(int type, double value);

void define_pattern(char *name, char *noi, void *p);
void *define_p_cat(void *p1, void *p2);
void *define_p_term(unsigned rel, void *expr);  // expr is expression*
void *define_p_rep(void *t, double n);

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

void start_include(char *file_name);
#endif

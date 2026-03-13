//
// parser_interface.cpp -- implementations of all parser callback functions.
//
// These are called by the Bison/Flex parser as it processes the annotated
// spice deck.  They are declared in parser_interf.h (hand-written, C-callable).
//
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>

#include "csv_analyzer.h"       // for units, f1_table, f2_table, function_type
#include "expression.h"
#include "sum_expression.h"
#include "avg_expression.h"
#include "min_expression.h"
#include "max_expression.h"
#include "geomean_expression.h"
#include "loop_complex.h"
#include "lsq_fit_function.h"
#include "nodes_of_interest.h"
#include "spice_deck.h"
#include "parameter.h"
#include "anneal.h"
#include "bo_optimizer.h"

extern "C" {
#include "parser_interf.h"
#include "lex.yy.h"                 // for yylineno
}


unsigned yy_n_parse_err = 0;        // Counts # of parse errors

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function tables -- looked up by define_function1 / define_function2
//

static double isfinite_op(double x) { return isfinite(x) ? 1.0 : 0.0; }

static double not_op(double x)
{
    // Note: what to do if x == NAN is up for debate. Here it is treated as a not-zero, thus a 0 is returned
    return (x == 0.0) ? 1.0 : 0.0;
}

f1_table func1_tab[] = {
    {"Sum",      nullptr},          // nullptr signals: use sum_expression subclass
    {"Avg",      nullptr},          // avg_expression
    {"Min",      nullptr},          // min_expression
    {"Max",      nullptr},          // max_expression
    {"Geomean",  nullptr},          // geomean_expression
    {"sqrt",     sqrt},
    {"abs",      fabs},
    {"ln",       log},
    {"exp",      exp},
    {"isfinite", isfinite_op},
    {nullptr,    nullptr}
};

f2_table func2_tab[] = {
    {"peak",     locate_peak,       noi_type},
    {"n_peaks",  count_peaks,       noi_type},
    {"t_rise",   locate_rise,       noi_type},
    {"t_fall",   locate_fall,       noi_type},
    {"t_edge",   locate_edge,       noi_type},
    {"n_edges",  count_edges,       noi_type},
    {"fit_ok",   test_lsq_fit,      lsq_fit_type},
    {"fit_coef", lsq_fit_get_cx,    lsq_fit_type},
    {"fit_resi", lsq_fit_residual,  lsq_fit_type},
    {"match_peaks",    match_peaks, tm_pattern},
    {"match_rise_edge",match_r_edg, tm_pattern},
    {"match_fall_edge",match_f_edg, tm_pattern},
    {"match_any_edge", match_a_edg, tm_pattern},
    {nullptr,    nullptr,      none_of_the_above}
};

units unit_table[] = {
    {' ', 1.0},                             // Nothing: unity
    {'G', 1.0e9},                           // Giga
    {'M', 1.0e6},                           // Mega
    {'K', 1.0e3},                           // Kilo
    {'m', 1.0e-3},                          // milli
    {'u', 1.0e-6},                          // micro
    {'n', 1.0e-9},                          // nano
    {'p', 1.0e-12},                         // pico
    {'L', 125.0e-6},                        // Liks, unit of 125 micro-Ampere
    {'O', 2.632e-12},                       // oHenry
    {'T', 2.0 * M_PI},                      // Turns
    {0, 1.0}
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Spice deck callbacks
//


void define_tran(double t_incr, double t_stop, double t_start, double dT_max)
//
// This is intercepted only to get a conservative estimate on the #of output line to
// allocate the storage for the time series.
// Note: the parser does not cover all flavors of the .tran statement, just what I need right now
//       TBD: This needs to be fixed. Robust workarround: if not parsed correctly, just
//       pass it through and leave the allocation as is.
//
{
    char buf[128];                  // just a sufficiently large buffer
    sprintf(buf, ".tran %.12lg  %.12lg %.12lg %.12lg\n", t_incr, t_stop, t_start, dT_max);
    add_line_to_spice_deck(strdup(buf));    // Just add the original statement
    
    unsigned n = ceil(t_stop / t_incr);
    if (n < 100)
        fprintf(stderr, "Warning: the .tran statement in the spice deck produces less than 100 rows of output\n");
    
    time_series::set_default_length(n + 16); // adds a little slack
    
    sim_time_start = t_start;
    sim_time_incr = t_incr;
}

void add_line_to_spice_deck(char *string)
// Just a wraper to add a text line to the spice void add_line_to_spice_deck(char* string)
{
    assert(string != nullptr);
    circuit.add_line(string);
}

void add_subst_to_spice_deck(char *p_name)
// Add a parameter reference
{
    parameter *p_ptr = parameter::find_parameter(p_name);
    if (p_ptr == nullptr) {
        fprintf(stderr, "Line %d: reference to undefined parameter '%s'\n", yylineno, p_name);
        yy_n_parse_err++;
    } else
        circuit.add_param_ref(p_ptr);
    
    free(p_name);                   // This was allocated via strdup() in the lexer and is no longer needed
                                    // Yeah, this isn't good practice, so sue me
}

void add_new_line_to_spice_deck()
{
    circuit.add_nl();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parameter definition callbacks
//

void * define_range(double f, double t)
{
    if (f >= t) {
        fprintf(stderr, "Line %d: [from < to] required\n", yylineno);
        yy_n_parse_err += 1;
        return nullptr;
    }

    struct range_pair *rp = (struct range_pair *) malloc(sizeof(struct range_pair));
    rp->from = f;
    rp->to = t;
    
    return rp;
}

void define_param_expression(char *name, void* expr, void *rng, unsigned flags)
// Defines a new expression type parameter
{
    if (parameter::find_parameter(name) != nullptr) {
        fprintf(stderr, "Line %d: redefinition of '%s'\n", yylineno, name);
        yy_n_parse_err += 1;
        return;
    }
    
    double v_min = -__DBL_MAX__;
    double v_max =  __DBL_MAX__;
    if (rng != nullptr) {
        v_min = ((range_pair*) rng)->from;
        v_max = ((range_pair*) rng)->to;
        free(rng);
    }
    
    if (flags & 6) {
        fprintf(stderr, "Line %d: Cannot tune/log-map expressions\n", yylineno);
        yy_n_parse_err += 1;
        return;
    }
    
    new expr_parameter(name, (expression *) expr, v_min, v_max, flags & 1);
}

void define_param_constant(char *name, double value, void *rng, unsigned flags)
// Defines a new constant type parameter that might be tunable
{
    if (parameter::find_parameter(name) != nullptr) {
        fprintf(stderr, "Line %d: redefinition of '%s'\n", yylineno, name);
        yy_n_parse_err += 1;
        return;
    }
    
    double v_min = -__DBL_MAX__;
    double v_max =  __DBL_MAX__;
    if (rng != nullptr) {
        v_min = ((range_pair*) rng)->from;
        v_max = ((range_pair*) rng)->to;
        free(rng);
    } else if (flags & 6) {
        fprintf(stderr, "Line %d: tuning/mapping requires a range specification\n", yylineno);
        yy_n_parse_err += 1;
        return;
    }
    
    if ((v_min > value) || (value > v_max)) {
        fprintf(stderr, "Line %d: value out of range\n", yylineno);
        yy_n_parse_err += 1;
        return;
    }
    
    if ((flags & 2) && (v_min <= 0.0)) {
        fprintf(stderr, "Line %d: log-map needs positive a lower bound\n", yylineno);
        yy_n_parse_err += 1;
        return;
    }
    
    new const_parameter(name, value, v_min, v_max, flags & 1, (flags >> 2) & 1, (flags >> 1) & 1);
}

void define_param_scan(char *name, void *rng, double n_steps, unsigned flags)
{
    if (parameter::find_parameter(name)) {
        fprintf(stderr, "Line %d: duplicate parameter definition for '%s'\n",
                yylineno, name);
        yy_n_parse_err++;
        return;
    }
    
    if (flags & 4) {
        fprintf(stderr, "Line %d: can't tune a scan\n", yylineno);
        yy_n_parse_err++;
        return;
    }
   
    if ((n_steps < 0.0) || (n_steps > 1.0e5)) {
        fprintf(stderr, "Line %d: questionable number of steps\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    unsigned n = (unsigned) nearbyint(n_steps);
    
    double v_min = -__DBL_MAX__;
    double v_max =  __DBL_MAX__;
    if (rng != nullptr) {
        v_min = ((range_pair*) rng)->from;
        v_max = ((range_pair*) rng)->to;
        free(rng);
    } else if (flags & 2) {
        fprintf(stderr, "Line %d: log scanning needs a positive lower bound\n", yylineno);
        yy_n_parse_err += 1;
        return;
    }
    
    new scan_parameter(name, v_min, v_max, n, flags & 1, (flags >> 1) & 1);
}

void define_sim_anneal(char *log_file_nm)
{
   if (sim_anneal_ptr != nullptr) {
        fprintf(stderr, "Line %d: Multiple simulated annealing directives\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    
    if (baysian_opt != nullptr) {
        fprintf(stderr, "Line %d: Conflicts with Baysian optimization\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    
    parameter *ev_ptr = parameter::find_parameter(EVAL_PARAMETER);
    if (nullptr == ev_ptr) {
        fprintf(stderr, "Line %d: '%s' paramter missing\n", yylineno, EVAL_PARAMETER);
        yy_n_parse_err++;
        return;        
    }

    sim_anneal_ptr = new sim_anneal(log_file_nm, ev_ptr, parameter::find_parameter(REJECT_PARAMETER));
}

void define_bo(double n_itr)
{
    if (baysian_opt != nullptr) {
        fprintf(stderr, "Line %d: Multiple simulated baysian opt directives\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    
    if (baysian_opt != nullptr) {
        fprintf(stderr, "Line %d: Conflicts with simulated annealing\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    
    if (n_itr < 1.0) {
        fprintf(stderr, "Line %d: Bad number of iterations\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    unsigned n = (unsigned) nearbyint(n_itr);
    
    parameter *of_ptr = parameter::find_parameter(BAY_OPT_OBJECTIVE);
    if (nullptr == of_ptr) {
        fprintf(stderr, "Line %d: '%s' paramter missing\n", yylineno, BAY_OPT_OBJECTIVE);
        yy_n_parse_err++;
        return;        
    }
    
    baysian_opt = new BOOptimizer(of_ptr, n);
}

void define_add2SA_sched(double temp, double n_iter)
{
    if (sim_anneal_ptr == nullptr) {
        fprintf(stderr, "Line %d: Define simulated annealing first\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    
    if ((temp < 0.0) || (n_iter < 1.0)) {
        fprintf(stderr, "Line %d: temerature < 0 or #of steps < 1\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    
    unsigned n = (unsigned) nearbyint(n_iter);
    sim_anneal_ptr->add_sa_sched(temp, n);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Misc:
//

void define_monitor(char *name)
{
    if (nodes_of_interest::find(name) != nullptr) {
        fprintf(stderr, "Line %d: '%s' is already monitored\n", yylineno, name);
        yy_n_parse_err++;
        return;
    }
    new nodes_of_interest(name);
}

void define_snapshot(char *name, double start, double freq)
{
    extern char    *snapshot_file_name;
    extern unsigned snapshot_first;
    extern unsigned snapshot_frequency;

    int s = (int) nearbyint(start);
    int f = (int) nearbyint(freq);
    if ((s < 0) || (f < 1)) {
        fprintf(stderr,
            "Bad snapshot pragma: start >= 0 and frequency >= 1 required\n");
        yy_n_parse_err++;
        return;
    }
    snapshot_file_name  = name;
    snapshot_first      = (unsigned) s;
    snapshot_frequency  = (unsigned) f;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Expression construction callbacks
//

void *define_const(double x)        { return new expression(x); }

void *define_ref(char *name)
{
    parameter *p_ptr = parameter::find_parameter(name);
    if (p_ptr == nullptr) {
        fprintf(stderr, "Line %d: reference to undefined parameter '%s' - needs to be defined first\n",
                yylineno, name);
        yy_n_parse_err += 1;
        return nullptr;
    }
    
    return new expression(p_ptr);
}

void *define_add(void *x, void *y)
{ return new expression((expression*)x, addition,       (expression*)y); }
void *define_sub(void *x, void *y)
{ return new expression((expression*)x, subtraction,    (expression*)y); }
void *define_mul(void *x, void *y)
{ return new expression((expression*)x, multiplication, (expression*)y); }
void *define_div(void *x, void *y)
{ return new expression((expression*)x, division,       (expression*)y); }
void *define_eq (void *x, void *y)
{ return new expression((expression*)x, comp_eq,        (expression*)y); }
void *define_ne (void *x, void *y)
{ return new expression((expression*)x, comp_ne,        (expression*)y); }
void *define_gt (void *x, void *y)
{ return new expression((expression*)x, comp_gt,        (expression*)y); }
void *define_ge (void *x, void *y)
{ return new expression((expression*)x, comp_ge,        (expression*)y); }
void *define_lt (void *x, void *y)
{ return new expression((expression*)x, comp_lt,        (expression*)y); }
void *define_le (void *x, void *y)
{ return new expression((expression*)x, comp_le,        (expression*)y); }

void *define_not(void *x)
{ return new expression(not_op, (expression*)x); }

void *define_test_sel(void *t, void *a, void *b)
{ return new expression((expression*)t, (expression*)a, (expression*)b); }

void *define_function1(char *name, void *x)
{
    for (unsigned i = 0; func1_tab[i].name != nullptr; i++) {
        if (strcmp(name, func1_tab[i].name)) continue;

        if (func1_tab[i].func1_ptr == nullptr) {
            // Stateful expression: dispatch by name
            expression *arg = (expression *) x;
            if      (!strcmp(name, "Sum"))     return new sum_expression    (arg);
            else if (!strcmp(name, "Avg"))     return new avg_expression    (arg);
            else if (!strcmp(name, "Min"))     return new min_expression    (arg);
            else if (!strcmp(name, "Max"))     return new max_expression    (arg);
            else if (!strcmp(name, "Geomean")) return new geomean_expression(arg);
            else { assert(0); return nullptr; }
        }

        return new expression(func1_tab[i].func1_ptr, (expression *) x);
    }

    fprintf(stderr, "Line %d: undefined function '%s'\n", yylineno, name);
    yy_n_parse_err++;
    return nullptr;
}

void *define_function2(char *name, char *obj_name, void *y)
{
    double (*f_ptr)(void *, double) = nullptr;
    function_type f_type = none_of_the_above;

    for (unsigned i = 0; func2_tab[i].name != nullptr; i++) {
        if (!strcmp(name, func2_tab[i].name)) {
            f_ptr   = func2_tab[i].func2_ptr;
            f_type  = func2_tab[i].f_type;
            break;
        }
    }
    if (!f_ptr) {
        fprintf(stderr, "Line %d: undefined function '%s'\n", yylineno, name);
        yy_n_parse_err++;
        return nullptr;
    }

    void *obj_ptr = nullptr;
    switch (f_type) {
        case noi_type:      obj_ptr = nodes_of_interest::find(obj_name);  break;
        case lsq_fit_type:  obj_ptr = lsq_fit_function::find(obj_name);   break;
        case tm_pattern:    obj_ptr = time_pattern::find(obj_name);       break;
        default: assert(0);
    }
    if (!obj_ptr) {
        fprintf(stderr, "Line %d: object '%s' is not defined\n",
                yylineno, obj_name);
        yy_n_parse_err++;
        return nullptr;
    }

    return new expression(f_ptr, obj_ptr, (expression *) y);
}

void define_lsq_fit(char *name, double order, void *x, void *y)
// creates an LSQ fit operation
//
// Adding an LSQ fit via the parameter class is a bit awkward and probably not a good
// idea. It is awkward, because a LSQ fit does not have an obviuos value associated with
// it.
{
    if (lsq_fit_function::find(name)) {
        fprintf(stderr, "Line %d: duplicate name '%s'\n", yylineno, name);
        yy_n_parse_err++;
        return;
    }
    int i_order = (int) nearbyint(order);
    if ((i_order < 1) || (i_order > 4)) {
        fprintf(stderr,
            "Line %d: only 1st to 4th order polynomials supported\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    // Constructor self-registers with loop_complex
    new lsq_fit_function(name, i_order,
                         (expression *) x, (expression *) y,
                         loop_complex);
}

void *define_p_cat(void *p1, void *p2)
// concatenate two pattern time elelemts
{
    struct p_time_element *t_ptr = (struct p_time_element *) p1;
    while (t_ptr->next != nullptr)          // finds the end of the list
        t_ptr->next = t_ptr;
    t_ptr = (struct p_time_element *) p2;
    return p1;
}

void *define_p_term(unsigned rel, double t)
{
    struct p_time_element *t_ptr = (struct p_time_element *) malloc(sizeof(struct p_time_element));
    t_ptr->next = nullptr;
    t_ptr->time = NAN;
    
    if (t < 0.0) {
        fprintf(stderr, "Line %d: time step must be >= 0\n", yylineno);
        yy_n_parse_err++;
        return t_ptr;
    }
    
    if (rel) t_ptr->time = -t;
    else     t_ptr->time = -t;
    
    return t_ptr;
}

void *define_p_rep(void *t, double n)
{
    struct p_time_element *t_ptr = (struct p_time_element *) t;
    if (n < 1.0) {
        fprintf(stderr, "Line %d: repetition count must >= 1\n", yylineno);
        yy_n_parse_err++;
        return t_ptr;
    }
    unsigned i = (unsigned) nearbyint(n);
    
    struct p_time_element *l_ptr = t_ptr;  // Pointer to the last element
    while (l_ptr->next != nullptr)         // finds the end of the list
        l_ptr->next = l_ptr;
    
    for (struct p_time_element *m_ptr = l_ptr; i > 0; i--) { // replicate n times 
        struct p_time_element *r_ptr = t_ptr;
        do {
            struct p_time_element *n_ptr = (struct p_time_element *) malloc(sizeof(struct p_time_element));
            n_ptr->time = r_ptr->time;
            n_ptr->next = nullptr;
            m_ptr->next = n_ptr;
            m_ptr = n_ptr;
            
            if (r_ptr == l_ptr)
                break;
            r_ptr = r_ptr->next;
        } while (1);
    }
    
    return t_ptr;
}

void define_pattern(char *name, char *noi, void *p)
{
    nodes_of_interest *noi_ptr = nodes_of_interest::find(noi);
    if (noi_ptr == nullptr) {
        fprintf(stderr, "Line %d: '%s' is not defined\n", yylineno, noi);
        yy_n_parse_err++;
        return;
    }
    
    time_pattern *t_ptr = time_pattern::find(name);
    if (t_ptr != nullptr) {
        fprintf(stderr, "Line %d: '%s' already in use\n", yylineno, name);
        yy_n_parse_err++;
        return;
    }
    
    t_ptr = new time_pattern(name, noi_ptr);
    
    double t = -__DBL_MAX__;
    struct p_time_element *te_ptr = (struct p_time_element *) p;
    while (te_ptr != nullptr) {
        if (te_ptr->time < 0.0) t -= te_ptr->time;
        else if (te_ptr->time <= t) {
            fprintf(stderr, "Line %d: time must be monotonically increasing (use +dt for relative time)\n",
                    yylineno);
            yy_n_parse_err++;
            return;             // Yeah, leaks memory, but will be aborting, so who cares.
        } else
            t = te_ptr->time;
        t_ptr->add_time(t);
        
        struct p_time_element *n_ptr = te_ptr->next;
        free(te_ptr);
        te_ptr = n_ptr;
    }
}

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
void define_param_scan(char *name, double v_start, double v_stop, double n_steps);
void define_monitor(char *name);

void define_para_expression(char *name, void *expr);
void *define_add(void *x, void *y);
void *define_sub(void *x, void *y);
void *define_mul(void *x, void *y);
void *define_div(void *x, void *y);
void *define_const(double x);
void *define_ref(char *name);
void *define_function1(char *name, void *x);
void *define_function2(char *name, char *ts_name, void *y);

#endif

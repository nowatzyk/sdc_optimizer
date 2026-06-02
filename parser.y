
%{
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "../parser_interf.h"

void yyerror(const char *s);
extern int yylineno;
int yylex(void);

%}

%union {
    char   *text;
    double val;
    void   *pointer;
    unsigned flags;
    struct {
        void *t_rise, *t_fall, *t_width, *v_high, *v_low;
    } pwl_ctrl;
}
    
%token <text> END LINE SYMBOL PSUBSTITUTION
%token <val> NUMBER
%type <pointer> expr terminal range opt_range pattern p_time bo_attr_list bo_attr
%type <flags> attribute attributes
%type <val> bo_n_iter
%type <pwl_ctrl> pwl_attrs
%token TRAN EOL PRAGMA SCAN MONITOR COMMENT PARAMETER EQUAL COMMA SET OSQB CSQB INCLUDE
%token SIM_ANNEAL SA_SCHEDULE SNAPSHOT LSQ_FIT NO_PRINT LOG_MAP TUNEABLE BAYSIAN_OPT PATTERN FOR
%token MARGIN BO_BINARY BO_GRADIENT BO_PROBABILISTIC
%token PWL PWL_RISE PWL_FALL PWL_WIDTH PWL_HIGH PWL_LOW
%token BO_N_RAYS BO_N_BRACKET BO_N_BISECT BO_THRESHOLD
%token BO_BEF_BUDGET BO_BEF_ITER BO_BEF_PROBES

%left TEST_OP OTHERWISE
%left COMP_GT COMP_GE COMP_EQ COMP_LT COMP_LE COMP_NE
%left PLUS MINUS
%left MULT DIVIDE
%left OPAR CPAR
%left NOT

%%
input:
        lines END EOL                       { add_line_to_spice_deck($2); add_new_line_to_spice_deck();}
    ;

lines:
    |   lines line
    ;
    
line:
        LINE                                { add_line_to_spice_deck($1); }
    |   PSUBSTITUTION                       { add_subst_to_spice_deck($1); }
    |   COMMENT EOL
    |   TRAN NUMBER NUMBER NUMBER NUMBER EOL { define_tran($2, $3, $4, $5); }
    |   PRAGMA pragma_body EOL
    |   PRAGMA INCLUDE SYMBOL EOL           { start_include($3); }
    |   EOL                                 { add_new_line_to_spice_deck(); }
    ;

pragma_body:
        PARAMETER SYMBOL SCAN range NUMBER attributes
                                            { define_param_scan($2, $4, $5, $6); }
    |   PARAMETER SYMBOL EQUAL expr opt_range attributes
                                            { define_param_expression($2, $4, $5, $6); }
    |   PARAMETER SYMBOL SET NUMBER opt_range attributes
                                            { define_param_constant($2, $4, $5, $6); }
    |   PARAMETER SYMBOL PWL SYMBOL pwl_attrs  { define_param_pwl($2, $4,
                                                    $5.t_rise, $5.t_fall,
                                                    $5.t_width, $5.v_high,
                                                    $5.v_low); }
    |   MONITOR SYMBOL                      { define_monitor($2); }
    |   SNAPSHOT SYMBOL NUMBER COMMA NUMBER { define_snapshot($2, $3, $5); }
    |   LSQ_FIT SYMBOL NUMBER COMMA expr COMMA expr  { define_lsq_fit($2, $3, $5, $7); }
    |   SIM_ANNEAL SYMBOL                   { define_sim_anneal($2); }
    |   SA_SCHEDULE NUMBER COMMA NUMBER     { define_add2SA_sched($2, $4); }
    |   BAYSIAN_OPT bo_n_iter bo_attr_list  { define_bo($2, $3); }
    |   PATTERN SYMBOL FOR SYMBOL EQUAL pattern { define_pattern($2, $4, $6); }
    ;

//
// BO attribute list -- a linked list of bo_attr nodes built left-to-right.
// Passed as a single void* to define_bo(), which walks the list to extract
// all settings, applying defaults for anything not specified.
//
bo_attr_list :                              { $$ = NULL; }
    |   bo_attr_list bo_attr                { $$ = bo_attr_cat($1, $2); }
    ;

//
// Individual BO attributes.  Mode keywords carry no value (0.0 placeholder).
// Numeric keywords carry their value directly.
//
bo_attr :
        MARGIN                             { $$ = define_bo_attr(BO_ATTR_MARGIN,       0.0); }
    |   BO_BINARY                          { $$ = define_bo_attr(BO_ATTR_BINARY,       0.0); }
    |   BO_GRADIENT                        { $$ = define_bo_attr(BO_ATTR_GRADIENT,     0.0); }
    |   BO_PROBABILISTIC                   { $$ = define_bo_attr(BO_ATTR_PROBABILISTIC,0.0); }
    |   BO_N_RAYS     NUMBER               { $$ = define_bo_attr(BO_ATTR_N_RAYS,       $2); }
    |   BO_N_BRACKET  NUMBER               { $$ = define_bo_attr(BO_ATTR_N_BRACKET,    $2); }
    |   BO_N_BISECT   NUMBER               { $$ = define_bo_attr(BO_ATTR_N_BISECT,     $2); }
    |   BO_THRESHOLD  NUMBER               { $$ = define_bo_attr(BO_ATTR_THRESHOLD,    $2); }
    |   BO_BEF_BUDGET NUMBER               { $$ = define_bo_attr(BO_ATTR_BEF_BUDGET,   $2); }
    |   BO_BEF_ITER   NUMBER               { $$ = define_bo_attr(BO_ATTR_BEF_ITER,     $2); }
    |   BO_BEF_PROBES NUMBER               { $$ = define_bo_attr(BO_ATTR_BEF_PROBES,   $2); }
    ;

bo_n_iter :                                 { $$ = 190.0; }
    |   NUMBER                              { $$ = $1; }
    ;
    
pattern:
        p_time                              { $$ = $1; }
    |   pattern COMMA p_time                { $$ = define_p_cat($1, $3); }
    ;

p_time:
        expr                                { $$ = define_p_term(0, $1); }
    |   PLUS expr                           { $$ = define_p_term(1, $2); }
    |   OSQB pattern CSQB NUMBER            { $$ = define_p_rep($2, $4); }
    ;

opt_range:                                  { $$ = NULL; }
    |   range                               { $$ = $1; }
    ;
    
range:
        OSQB NUMBER COMMA NUMBER CSQB       { $$ = define_range($2, $4); }
    ;
    
attributes:                                 { $$ = 0; }
    |   attributes attribute                { $$ = $1 | $2; }
    ;
    
attribute :
        NO_PRINT                            { $$ = 1; }
    |   LOG_MAP                             { $$ = 2; }
    |   TUNEABLE                            { $$ = 4; }
    ;
    
expr:
        terminal                            { $$ = $1; }
    |   expr PLUS expr                      { $$ = define_add($1, $3); }
    |   expr MINUS expr                     { $$ = define_sub($1, $3); }
    |   expr MULT expr                      { $$ = define_mul($1, $3); }
    |   expr DIVIDE expr                    { $$ = define_div($1, $3); }
    |   expr COMP_EQ expr                   { $$ = define_eq($1, $3); }
    |   expr COMP_NE expr                   { $$ = define_ne($1, $3); }
    |   expr COMP_GT expr                   { $$ = define_gt($1, $3); }
    |   expr COMP_GE expr                   { $$ = define_ge($1, $3); }
    |   expr COMP_LT expr                   { $$ = define_lt($1, $3); }
    |   expr COMP_LE expr                   { $$ = define_le($1, $3); }
    |   OPAR expr CPAR                      { $$ = $2; }
    |   expr TEST_OP expr OTHERWISE expr    { $$ = define_test_sel($1, $3, $5); }
    |   NOT expr                            { $$ = define_not($2); }
    ;
    
terminal:
        NUMBER                              { $$ = define_const($1); }
    |   SYMBOL OPAR expr CPAR               { $$ = define_function1($1, $3); }
    |   SYMBOL OPAR SYMBOL COMMA expr CPAR  { $$ = define_function2($1, $3, $5); }
    |   SYMBOL                              { $$ = define_ref($1); }
    ;

//
// pwl_attrs -- the five control value expressions for a PWL parameter.
// Keywords may appear in any order; each sets one field in the pwl_ctrl struct.
// A missing keyword leaves the field as nullptr, which define_param_pwl()
// will catch and report as an error.
//
pwl_attrs :                             { $$.t_rise=$$.t_fall=$$.t_width=
                                          $$.v_high=$$.v_low=NULL; }
    |   pwl_attrs PWL_RISE  expr        { $$ = $1; $$.t_rise  = $3; }
    |   pwl_attrs PWL_FALL  expr        { $$ = $1; $$.t_fall  = $3; }
    |   pwl_attrs PWL_WIDTH expr        { $$ = $1; $$.t_width = $3; }
    |   pwl_attrs PWL_HIGH  expr        { $$ = $1; $$.v_high  = $3; }
    |   pwl_attrs PWL_LOW   expr        { $$ = $1; $$.v_low   = $3; }
    ;

%%

void yyerror(const char *s) {
    fprintf(stderr, "Error in line %d : %s\n", yylineno, s);
}

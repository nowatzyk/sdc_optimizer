
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
}
    
%token <text> END LINE SYMBOL PSUBSTITUTION
%token <val> NUMBER
%type <pointer> expr terminal range opt_range
%type <flags> attribute attributes
%token TRAN EOL PRAGMA SCAN MONITOR COMMENT PARAMETER EQUAL COMMA SET OSQB CSQB
%token SIM_ANNEAL SA_SCHEDULE SNAPSHOT LSQ_FIT NO_PRINT LOG_MAP TUNEABLE

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
    |   EOL                                 { add_new_line_to_spice_deck(); }
    ;

pragma_body:
        PARAMETER SYMBOL SCAN range NUMBER attributes
                                            { define_param_scan($2, $4, $5, $6); }
    |   PARAMETER SYMBOL EQUAL expr opt_range attributes
                                            { define_param_expression($2, $4, $5, $6); }
    |   PARAMETER SYMBOL SET NUMBER opt_range attributes
                                            { define_param_constant($2, $4, $5, $6); }
    |   MONITOR SYMBOL                      { define_monitor($2); }
    |   SNAPSHOT SYMBOL NUMBER COMMA NUMBER { define_snapshot($2, $3, $5); }
    |   LSQ_FIT SYMBOL NUMBER COMMA expr COMMA expr  { define_lsq_fit($2, $3, $5, $7); }
    |   SIM_ANNEAL SYMBOL                   { define_sim_anneal($2); }
    |   SA_SCHEDULE NUMBER COMMA NUMBER     { define_add2SA_sched($2, $4); }
    ;
    
opt_range:                                  { $$ = NULL; }
    |   range                               { $$ = $1; }
    ;
    
attributes:                                 { $$ = 0; }
    |   attributes attribute                { $$ = $1 | $2; }
    ;
    
attribute :
        NO_PRINT                            { $$ = 1; }
    |   LOG_MAP                             { $$ = 2; }
    |   TUNEABLE                            { $$ = 4; }
    ;
    
range:
        OSQB NUMBER COMMA NUMBER CSQB       { $$ = define_range($2, $4); }
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
        
%%

void yyerror(const char *s) {
    fprintf(stderr, "Error in line %d : %s\n", yylineno, s);
}

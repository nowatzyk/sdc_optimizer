
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
}
    
%token <text> END LINE SYMBOL PSUBSTITUTION
%token <val> NUMBER
%type <pointer> expr terminal
%token TRAN EOL PRAGMA SCAN MONITOR COMMENT PARAMETER ASSIGN COMMA
%token SIM_ANNEAL

%left TEST_OP OTHERWISE
%left COMP_GT COMP_GE COMP_EQ COMP_LT COMP_LE COMP_NE
%left PLUS MINUS
%left MULT DIVIDE
%left OPAR CPAR

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
        PARAMETER SYMBOL SCAN NUMBER NUMBER NUMBER  { define_param_scan($2, $4, $5, $6);}
    |   PARAMETER SYMBOL ASSIGN expr        { define_para_expression($2, $4); }
    |   PARAMETER SYMBOL SIM_ANNEAL NUMBER NUMBER NUMBER  { define_sim_anneal($2, $4, $5, $6);}
    |   MONITOR SYMBOL                      { define_monitor($2);}
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

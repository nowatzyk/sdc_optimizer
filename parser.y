
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
%token TRAN EOL PRAGMA SCAN MONITOR COMMENT PARAMETER OPAR CPAR ASSIGN COMMA

%left PLUS MINUS
%left MULT DIVIDE

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
    |   MONITOR SYMBOL                      { define_monitor($2);}
    ;
    
expr:
        terminal                            { $1; }
    |   expr PLUS terminal                  { define_add($1, $3); }
    |   expr MINUS terminal                 { define_sub($1, $3); }
    |   expr MULT terminal                  { define_mul($1, $3); }
    |   expr DIVIDE terminal                { define_div($1, $3); }
    |   OPAR expr CPAR                      { $2; }
    ;
    
terminal:
        NUMBER                              { define_const($1); }
    |   SYMBOL OPAR expr CPAR               { define_function1($1, $3); }
    |   SYMBOL OPAR SYMBOL COMMA expr CPAR  { define_function2($1, $3, $5); }
    |   SYMBOL                              { define_ref($1); }
    ;
        
%%

void yyerror(const char *s) {
    fprintf(stderr, "Error in line %d : %s\n", yylineno, s);
}


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
%token TRAN EOL PRAGMA SCAN MONITOR COMMENT PARAMETER

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
    |   MONITOR SYMBOL                      { define_monitor($2);}
    ;
        
%%

void yyerror(const char *s) {
    fprintf(stderr, "Error in line %d : %s\n", yylineno, s);
}

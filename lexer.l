%option noyywrap nodefault yylineno

%{
#include "parser.h"   // Bison-generated header

#include <math.h>
#include "../parser_interf.h"
%}

%x IN_PRAGMA

%%

^\*Pragma/[ \t].*    { BEGIN(IN_PRAGMA); return PRAGMA; }

^\.tran/[ \t].*      { BEGIN(IN_PRAGMA); return TRAN; }

^\.end/[ \t\n].*     { yylval.text = strdup(yytext); return END; }

<*>\n                { BEGIN(INITIAL); return EOL; }

<IN_PRAGMA>[-0-9.eE]+[GMKmunpLOT]? { double val;
                        if (1 != sscanf(yytext, "%lf", &val)) {
                            fprintf(stderr, "Bad number in circuit file on line %i\n", yylineno);
                            val = NAN;
                        };
                        char postfix = yytext[strlen(yytext)-1];
                        for (unsigned i = 0; unit_table[i].c != 0; i++)
                            if (postfix == unit_table[i].c) {
                                val *= unit_table[i].val;
                                break;
                            };
                        yylval.val = val;
                        return NUMBER;
                    }
                    
<IN_PRAGMA>"scan"    { return SCAN; }

<IN_PRAGMA>"parameter" { return PARAMETER; }

<IN_PRAGMA>"monitor" { return MONITOR; }

<IN_PRAGMA>[ \t]*   ; // ignore white space
                    
<IN_PRAGMA>[a-zA-Z][a-zA-Z0-9_|]+ { yylval.text = strdup(yytext); return SYMBOL; }

<IN_PRAGMA>.        { fprintf(stderr, "Line %d: Unrecognized pragma character in *.cir file: '%s'\n", yylineno, yytext); } 


^"*".*              { return COMMENT; }

~[a-zA-Z][a-zA-Z0-9_|]+~ {
                        yytext[strlen(yytext) - 1] = 0;         // strip trailing "~"
                        yylval.text = strdup(yytext+1);
                        return PSUBSTITUTION;
                    }
                    
[^~\n]+             { yylval.text = strdup(yytext); return LINE; }

.                   { fprintf(stderr, "Line %d: Unmatched character in *.cir file: '%s'\n", yylineno, yytext); } 

%%

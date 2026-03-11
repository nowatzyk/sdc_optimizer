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

<IN_PRAGMA>[-]?[0-9.]+[eE]?[-]?[0-9]*[GMKmunpLOT]? { double val;
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
                    
<IN_PRAGMA>"NAN"     { yylval.val = NAN; return NUMBER; }
                    
<IN_PRAGMA>"scan"    { return SCAN; }

<IN_PRAGMA>"sim_anneal" { return SIM_ANNEAL; }
<IN_PRAGMA>"sa_schedule" { return SA_SCHEDULE; }

<IN_PRAGMA>"baysian_opt" { return BAYSIAN_OPT; }

<IN_PRAGMA>"parameter" { return PARAMETER; }

<IN_PRAGMA>"monitor" { return MONITOR; }

<IN_PRAGMA>"snapshot" { return SNAPSHOT; }

<IN_PRAGMA>"lsq_fit" { return LSQ_FIT; }

<IN_PRAGMA>"no_print" { return NO_PRINT; }

<IN_PRAGMA>"log_map" { return LOG_MAP; }

<IN_PRAGMA>"tune"    { return TUNEABLE; }

<IN_PRAGMA>[ \t]*   ; // ignore white space
                    
<IN_PRAGMA>[a-zA-Z][a-zA-Z0-9_]* { yylval.text = strdup(yytext); return SYMBOL; }
<IN_PRAGMA>\"[^\n\t\" ]+\" { // Symbol can be arbitrary text that is in double quotes
                        yytext[strlen(yytext) - 1] = 0;
                        yylval.text = strdup(yytext+1);
                        return SYMBOL;
                    }

<IN_PRAGMA>"="      { return EQUAL; }
<IN_PRAGMA>":="     { return SET; }
<IN_PRAGMA>"+"      { return PLUS; }
<IN_PRAGMA>"-"      { return MINUS; }
<IN_PRAGMA>"("      { return OPAR; }
<IN_PRAGMA>")"      { return CPAR; }
<IN_PRAGMA>"*"      { return MULT; }
<IN_PRAGMA>"/"      { return DIVIDE; }
<IN_PRAGMA>","      { return COMMA; }
<IN_PRAGMA>"?"      { return TEST_OP; }
<IN_PRAGMA>":"      { return OTHERWISE; }
<IN_PRAGMA>"=="     { return COMP_EQ; }
<IN_PRAGMA>"!="     { return COMP_NE; }
<IN_PRAGMA>">"      { return COMP_GT; }
<IN_PRAGMA>">="     { return COMP_GE; }
<IN_PRAGMA>"<="     { return COMP_LE; }
<IN_PRAGMA>"<"      { return COMP_LT; }
<IN_PRAGMA>"!"      { return NOT; }
<IN_PRAGMA>"["      { return OSQB; }
<IN_PRAGMA>"]"      { return CSQB; }

<IN_PRAGMA>.        { fprintf(stderr, "Line %d: Unrecognized pragma character in *.cir file: '%s'\n", yylineno, yytext); } 


^"*".*              { return COMMENT; }

~[a-zA-Z][a-zA-Z0-9_]+~ {
                        yytext[strlen(yytext) - 1] = 0;         // strip trailing "~"
                        yylval.text = strdup(yytext+1);
                        return PSUBSTITUTION;
                    }
                    
[^~\n]+             { yylval.text = strdup(yytext); return LINE; }

.                   { fprintf(stderr, "Line %d: Unmatched character in *.cir file: '%s'\n", yylineno, yytext); } 

%%

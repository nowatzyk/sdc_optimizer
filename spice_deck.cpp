#include "spice_deck.h"
#include "parameter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

extern "C" {
#include "parser.h"
#include "lex.yy.h"
}

extern unsigned yy_n_parse_err;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Global instance
//

spice_deck circuit;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// spice_elements
//

spice_elements::spice_elements(char *t)
{
    se_type  = text;
    next     = nullptr;
    txt.text = t;
}

spice_elements::spice_elements(class parameter *par_ptr)
{
    se_type   = parameter;
    next      = nullptr;
    par.param = par_ptr;
}

spice_elements::spice_elements()
{
    se_type = new_line;
    next    = nullptr;
}

void spice_elements::print(FILE *of)
{
    switch (se_type) {
        case text:
            fputs(txt.text, of);
            break;
        case parameter:
            fprintf(of, "%.12lg", par.param->get_cur_value());
            break;
        case new_line:
            fputc('\n', of);
            break;
        default:
            assert(0);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// spice_deck
//

spice_deck::spice_deck()
    : first(nullptr), last(nullptr)
{
}

void spice_deck::add_line(char *txt)
{
    spice_elements *sep = new spice_elements(txt);
    if (first == nullptr) first = sep;
    else                  last->add_next(sep);
    last = sep;
}

void spice_deck::add_param_ref(parameter *par)
{
    spice_elements *sep = new spice_elements(par);
    if (first == nullptr) first = sep;
    else                  last->add_next(sep);
    last = sep;
}

void spice_deck::add_nl()
{
    if ((last != nullptr) && (last->is_nl() == 1))
        return;     // suppress consecutive newlines
    spice_elements *sep = new spice_elements();
    if (first == nullptr) first = sep;
    else                  last->add_next(sep);
    last = sep;
}

int spice_deck::read_cir_file(const char* fn)
//
// Reads the spice deck into memory
// The intent is to supply it to JoSim multiple times, but with the ability to
// change variable, parameters etc.
// I need a tool for the design space exploration, so running JoSim manually gets old fast.
//
// Return codes:
//   >= 0 : number of lines successfully read
//     -1 : failed to open the file
//
{
    yyin = fopen(fn, "r");
    if (yyin == nullptr)
        return -1;
    
    if ((yyparse() != 0) || (yy_n_parse_err > 0)) {
        fprintf(stderr, "Parsing of the spice circuit input failed.\n");
        exit(1);
    }
    
    fclose(yyin);
    
    return yylineno;
}

int spice_deck::write_cir_file(const char* fn)
//
// Write the spice deck to a file
//
{
    FILE *out = fopen(fn, "w");
    if (out == nullptr)
        return -1;
    
    for(spice_elements *spe_ptr = first; spe_ptr != nullptr; spe_ptr = spe_ptr->get_next())
        spe_ptr->print(out);
    
    fclose(out);
    return 0;
}

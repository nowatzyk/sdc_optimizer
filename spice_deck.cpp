#include "spice_deck.h"
#include "parameter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

extern unsigned yy_n_parse_err;
extern unsigned n_josim_runs;

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
    double v;
    
    switch (se_type) {
        case text:
            fputs(txt.text, of);
            break;

        case parameter:
            v = par.param->get_cur_value();
            if (!isfinite(v)) {
                fprintf(stderr, "Reference to '%s' parameter yieded NAN in run %u\n",
                        par.param->get_name(), n_josim_runs);
                exit(1);
            }
            fprintf(of, "%.12lg", v);
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

int spice_deck::open_include_file(char* fn)
//
// Facility to alow include files. Recursive include files are OK
//
{
    FILE *incl_fp = fopen(fn, "r");
    free (fn);
    if (incl_fp == nullptr)
        return -1;

    include_file_stack.push_back(get_current_buffer());
    include_line_no_stk.push_back(yylineno);
    yylineno = 1;
    yyin = incl_fp;
    yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));
    return 0;
}

int spice_deck::yywrap()
{
    if (include_file_stack.size() > 0) {
        fclose(yyin);
        
        // Restore previous buffer
        yy_switch_to_buffer(include_file_stack.back());
        yylineno = include_line_no_stk.back();
        include_file_stack.pop_back();
        include_line_no_stk.pop_back();
        return 0;
    }
    return 1;
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

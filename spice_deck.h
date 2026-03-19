#pragma once

#include <cstdio>
#include <cassert>
#include <vector>
using namespace std;

extern "C" {
    #include "parser.h"
    #include "lex.yy.h"
    YY_BUFFER_STATE get_current_buffer(void);
}

class parameter;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// spice_elements -- one token in the spice deck (text, parameter reference, or newline).
//

extern const char spice_escape;             // Escape character in spice decks for parameter substitutions
extern const unsigned spice_src_max_char;   // Max line length in the spice source file

struct spe_text  { char      *text;  };
struct spe_param { parameter *param; };

class spice_elements {
    enum { text, parameter, new_line } se_type;
    spice_elements *next;
    union {
        spe_text  txt;
        spe_param par;
    };

public:
    spice_elements(char *txt);              // text element
    spice_elements(class parameter *par);   // parameter reference
    spice_elements();                       // newline element

    void            add_next(spice_elements *nxt) { next = nxt; }
    spice_elements *get_next()                    { return next; }
    int             is_nl()  { return (se_type == new_line) ? 1 : 0; }

    void print(FILE *of);
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// spice_deck -- in-memory representation of the circuit file.
//
// Read once at startup, written to the JoSIM input FIFO before each run
// with current parameter values substituted.
//

class spice_deck {
    spice_elements *first, *last;
    vector<YY_BUFFER_STATE>  include_file_stack;
    vector<int>              include_line_no_stk;
    
public:
    spice_deck();

    void add_line     (char *txt);
    void add_param_ref(parameter *par);
    void add_nl       ();
    
    int yywrap();                           // Function to teake care of EOF in include files

    int  read_cir_file (const char *fn);    // parses and loads; returns line count or -1
    int  open_include_file(char *fn);       // starts to read from the include file
    int  write_cir_file(const char *fn);    // writes with substituted values; returns 0 or -1
};

extern spice_deck circuit;

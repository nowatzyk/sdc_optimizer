//
//  This file contails the code for the parameter class.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <vector>
using namespace std;

#include "parameter.h"
#include "expression.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Statics:
//

vector<parameter *> parameter::parameters;      // The collection of all parameters
unsigned parameter::nesting_level = 0;          // = #of looping constructs


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The parameter base class
//

parameter::parameter(char* nm, double v_min, double v_max, unsigned npr, unsigned tun, unsigned l_map) : 
    name(nm), type(undefined),
    no_print(npr), tunable(tun), log_map(l_map), min_value(v_min), max_value(v_max), value(NAN)
{
    assert(v_min < v_max);
    parameters.push_back(this);         // TBD this needs to be verified
}

void parameter::print_name(FILE* fp)
{
    fprintf(fp, "%s", name);
}

void parameter::print_value(FILE* fp)
{
    fprintf(fp, " %.15lg", value);
}

parameter *parameter::find_parameter(const char* nm)
{
    for (unsigned i = 0; i < parameters.size(); i++)
        if (!strcmp(nm, parameters[i]->name))
            return parameters[i];

    return nullptr;
}

unsigned parameter::list_names(FILE* fp)
{
    unsigned cnt = 1;
    for (unsigned i = 0; i < parameters.size(); i++)
        if (parameters[i]->no_print == 0) {
            fprintf(fp, " (%u):", cnt++);
            parameters[i]->print_name(fp);
        }
    
    return cnt;
}

void parameter::list_c_val(FILE* fp)
{
    for (unsigned i = 0; i < parameters.size(); i++)
        if (parameters[i]->no_print == 0)
            parameters[i]->print_value(fp);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  The constant constant parameter class 
//
//  Note: this is the only class that can be made tunable. Tuning an extression, or loop
//        operator does not make sense
//

const_parameter::const_parameter(char* nm, double val, double v_min, double v_max,
                                 unsigned int npr, unsigned int tun, unsigned int l_map) :
                                 parameter(nm, v_min, v_max, npr, tun, l_map)
{
    // If tuning is enabled: a sensible range must be given
    assert((tun == 0) || (v_min != -__DBL_MAX__ && v_max != __DBL_MAX__));
    
    // if logarithmic mapping is enabled, the min-range must be > 0
    assert((l_map == 0) || (v_min > 0));
    
    value = val;
}

double const_parameter::get_mapped_value()
{
    double m_val = NAN;             // Default to when the valuse is out of bound
    
    if ((min_value <= value) && (value <= max_value)) {
        if (log_map)
            // Logarithmic mapping:
            m_val = (log(value) - log(min_value)) / (log(max_value) - log(min_value));
        else
            // Linear mapping:
            m_val = (value - min_value) / (max_value - min_value);
    }
    
    return m_val;
}

void const_parameter::set_mapped_value(double m_val)
{
    assert((0.0 <= m_val) && (m_val <= 1.0));
    if (log_map)
        value = exp(m_val * log(max_value) + (1.0 - m_val) * log(min_value));
    else
        value = min_value + m_val * (max_value - min_value);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// the assignment parameter class
//

expr_parameter::expr_parameter(char* nm, expression* e_ptr, double v_min, double v_max, unsigned npr) :
    parameter(nm, v_min, v_max, npr, 0, 0)
{
    expr = e_ptr;
}

void expr_parameter::update()
{
    double new_value = expr->get_value();
    
    if (isfinite(new_value)) {
        //
        // If a range is specified then the value will be clamped to this range.
        //
        // Other semantics are possible and may be specified via an attribute.
        // For example set to NAN if the expression is out of bound, or issue a warning, or abort
        //
        if (new_value < min_value) new_value = min_value;
        if (new_value > max_value) new_value = max_value;
    }
    
    value = new_value;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// the scan paramter class
//

scan_parameter::scan_parameter(char* nm, double v_min, double v_max, unsigned n_s,
                               unsigned npr, unsigned l_map) :
    parameter(nm, v_min, v_max, npr, 0, l_map)
{
    assert((v_min < v_max) && (-__DBL_MAX__ < v_min) && (v_max < __DBL_MAX__) && (n_s > 1));
    n_steps = n_s;
    step_cntr = 0;
    if (log_map) {
        assert(v_min > 0.0);
        d_value = (log(v_max) - log(v_min)) / (double) (n_s - 1);
    } else
        d_value = (v_max - v_min) / (double) (n_s - 1);
}

void scan_parameter::update()
{
    if (log_map)
        value = min_value * exp(d_value * (double) step_cntr);
    else
        value = min_value + d_value * (double) step_cntr;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

extern "C" {
#include "parser.h"                         // Bison generated headers (in build dir)
#include "lex.yy.h"                         // Flex-generated header (in build dir)
#include "parser_interf.h"                  // needed to integrate the parser
#include "lsq_fit.h"
#include "fit_functions.h"
}

#include "csv_analyzer.h"

const char* csv_out_path = "/tmp/JoSim2csv_analyzer.csv";      // Output fifo from Josim to this program
const char* cir_inp_path = "/tmp/JoSim2csv_analyzer.cir";      // Input fifo to Josim


////////////////////////////////////////////////////////////////////
//
// Static/gobal things:
//

nodes_of_interest* nodes_of_interest::root = nullptr; // Root of the NOI list
unsigned nodes_of_interest::n_noi = 0;      // Counts the number of NOI's

parameter* parameter::root = nullptr;       // start of the parameter list
param_iterator parameter::iterator[max_param_iterators] = {{nullptr, nullptr}}; // mechanism to allow nested loops
unsigned parameter::nesting_level = 0;      // Nesting level
 
units unit_table[] = {
    {' ', 1.0},                             // Nothing: unity
    {'G', 1.0e9},                           // Giga
    {'M', 1.0e6},                           // Mega
    {'K', 1.0e3},                           // Kilo
    {'m', 1.0e-3},                          // milli
    {'u', 1.0e-6},                          // micro
    {'n', 1.0e-9},                          // nano
    {'p', 1.0e-12},                         // pico
    {'L', 125.0e-6},                        // Liks, unit of 125 micro-Ampere
    {'O', 2.632e-12},                       // oHenry
    {'T', 2.0 * M_PI},                      // Turns
    {0, 1.0}
};

unsigned yy_n_parse_err = 0;                // Counts # of parse errors

unsigned time_series::n_init = 1024;        // First allocation default
time_series* time_series::root = nullptr;   // Root of the time series

spice_deck circuit;                         // circuit under test

f1_table func1_tab[] = {                    // Table of functions with one argument
    {"Sum", nullptr},                       // The summartion function is special and it is
                                            // implemented in the expression class
    {"sqrt", sqrt},                         // Square root
    {"abs", fabs},                          // Absolute
    {"ln", log},                            // Natural log
    {"exp", exp},                           // e^x
    {nullptr, nullptr}
};

f2_table func2_tab[] = {                    // Functions with 2 arguments
    {"peak", locate_peak},
    {"n_peaks", count_peaks},
    {"t_rise", locate_rise},
    {"t_fall", locate_fall},
    {"t_edge", locate_edge},
    {"n_edges", count_edges},
    {nullptr, nullptr}
};
    
time_series *josim_out_columns[max_ts] = {nullptr}; // Storage for Josim output
unsigned n_josim_out_columns = 0;           // #of these columns that are in use
unsigned n_josim_runs = 0;                  // # of josim runs executed
                                            // Note: the first one is used to set up things
                                            
double sim_time_start = -__DBL_MAX__;       // Simulation time start (in Seconds)
double sim_time_incr = 0.0;                 // Simulation time increment
                                            
char josim_output_buf[max_ln_length];       // Where to put the output from Josim (large, thus not put on heap)

double sim_time(double n)
// Returns the simulation time for the n-th row
// Note: row numbers may be fractional
{
    return sim_time_start + sim_time_incr * n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class functions
//

nodes_of_interest *nodes_of_interest::find(const char* name)
// Finds the NOI by name or return nullptr
{
    for (nodes_of_interest *noi_ptr = root; noi_ptr != nullptr; noi_ptr = noi_ptr->next)
        if (!strcmp(name, noi_ptr->name))
            return noi_ptr;
    return nullptr;
}

nodes_of_interest::nodes_of_interest(const char* nm)
// Constructor
{
    next = root;
    root = this;
    col_index = -1;
    name = nm;                              // was allocated via strdup in the lexer
    n_noi++;
}

unsigned nodes_of_interest::set_col_index(unsigned ind)
{
    if ((col_index < 0) || (col_index == ind)) {
        col_index = ind;
        return 0;
    }
    return 1;
}

nodes_of_interest *nodes_of_interest::find_undef()
{
    for (nodes_of_interest *n_ptr = root; n_ptr != nullptr; n_ptr = n_ptr->next)
        if (n_ptr->col_index < 0)
            return n_ptr;
    return nullptr;
}

int nodes_of_interest::print(char* f_name)
{
    FILE *of = fopen(f_name, "w");
    if (of == nullptr)
        return 1;
    fprintf(of, "# NOI name: %s\n", name);
    if ((col_index < 0) || (col_index >= max_ts) || (josim_out_columns[col_index] == nullptr))
        fprintf(of, "# has no data\n");
    else
        josim_out_columns[col_index]->print_all(of);
    fclose(of);
    
    return 0;       // success
}

void nodes_of_interest::print_all()
{
    unsigned i = 0;
    for (nodes_of_interest *n_ptr = root; n_ptr != nullptr; n_ptr = n_ptr->next) {
        char buf[64];
        sprintf(buf, "NOI_%u.dat", i);
        n_ptr->print(buf);
        i++;
    }
}


//////

expression::expression(double x)
{
    type = constant;
    value = x;
    l_arg = nullptr;
    r_arg = nullptr;
    t_arg = nullptr;
    func1_ptr = nullptr;
    p_ptr = nullptr;
    func2_ptr = nullptr;
    noi_ptr = nullptr;
}

expression::expression(double (*f_ptr) (double), expression* arg_ptr)
{
    type = (f_ptr == nullptr) ? sum_func : function1;
    value = 0.0;
    l_arg = arg_ptr;
    r_arg = nullptr;
    t_arg = nullptr;
    func1_ptr = f_ptr;
    p_ptr = nullptr;
    func2_ptr = nullptr;
    noi_ptr = nullptr;
}

expression::expression(expression* la_ptr, exp_type et, expression* ra_ptr)
{
    assert((et == addition) || (et == subtraction) || (et == multiplication) || (et == division) ||
           (et == comp_eq) || (et == comp_ne)  || (et == comp_lt) || (et == comp_le) || (et == comp_gt) || (et == comp_ge) );
    type = et;
    value = 0.0;
    l_arg = la_ptr;
    r_arg = ra_ptr;
    t_arg = nullptr;
    func1_ptr = nullptr;    
    p_ptr = nullptr;
    func2_ptr = nullptr;
    noi_ptr = nullptr;
}

expression::expression(parameter* param_ptr)
{
    type = p_reference;
    value = 0.0;
    l_arg = nullptr;
    r_arg = nullptr;
    t_arg = nullptr;
    func1_ptr = nullptr;
    p_ptr = param_ptr;
    func2_ptr = nullptr;
    noi_ptr = nullptr;
}

expression::expression(double (*f_ptr) (class nodes_of_interest *n_ptr, double x),
                       class nodes_of_interest* n_ptr, expression* arg_ptr)
{
    type = function2;
    value = 0.0;
    l_arg = arg_ptr;
    r_arg = nullptr;
    t_arg = nullptr;
    func1_ptr = nullptr;
    p_ptr = nullptr;
    func2_ptr = f_ptr;
    noi_ptr = n_ptr;
}

expression::expression(expression* tst_ptr, expression* true_arg, expression* false_arg)
{
    type = test_select;
    value = 0.0;
    l_arg = true_arg;
    r_arg = false_arg;
    t_arg = tst_ptr;
    func1_ptr = nullptr;    
    p_ptr = nullptr;
    func2_ptr = nullptr;
    noi_ptr = nullptr;
}

double expression::get_value()
{
    switch (type) {
        case constant:
            return value;
            
        case sum_func:
            value += l_arg->get_value();        // Note: this one is different - it has a side-effect
            return value;
            
        case function1:
            return func1_ptr(l_arg->get_value());
            
        case function2:
            return func2_ptr(noi_ptr, l_arg->get_value());
            
        case p_reference:
            return p_ptr->get_cur_value();
            
        case addition:
            return l_arg->get_value() + r_arg->get_value();
            
        case subtraction:
            return l_arg->get_value() - r_arg->get_value();
             
        case multiplication:
            return l_arg->get_value() * r_arg->get_value();
            
        case division:
            return l_arg->get_value() / r_arg->get_value();
            
        case test_select:
            return (t_arg->get_value() > 0.0) ? l_arg->get_value() : r_arg->get_value();
            
        case comp_eq:
            return (l_arg->get_value() == r_arg->get_value()) ? 1.0 : 0.0;
            
        case comp_ne:
            return (l_arg->get_value() != r_arg->get_value()) ? 1.0 : 0.0;
            
        case comp_gt:
            return (l_arg->get_value() > r_arg->get_value()) ? 1.0 : 0.0;
            
        case comp_ge:
            return (l_arg->get_value() >= r_arg->get_value()) ? 1.0 : 0.0;
            
        case comp_lt:
            return (l_arg->get_value() < r_arg->get_value()) ? 1.0 : 0.0;
            
        case comp_le:
            return (l_arg->get_value() <= r_arg->get_value()) ? 1.0 : 0.0;
            
        default: assert(0);
    }
}


//////

double parameter::get_cur_value()
{
    switch (type) {
        case scan:
            return (min_value + (double) step_cntr * d_value);
            
        case assignment:
            return expr->get_value();
            
        default:
            assert (0);
    }
}


int parameter::i_next()
{
    step_cntr++;
    
    if (step_cntr >= n_steps) {
        step_cntr = 0;
        return 1;
    } 
    
    return 0;
}

void parameter::print_name(FILE* fp)
{
    fprintf(fp, " %s", name);
}

void parameter::print_value(FILE* fp)
{
    fprintf(fp, " %.15lg", get_cur_value());
}

void parameter::reset()
{
    parameter *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->i_reset();
        p_ptr = p_ptr->next;
    }
}

void parameter::list_names(FILE* fp)
{
    parameter *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->print_name(fp);
        p_ptr = p_ptr->next;
    }
}

void parameter::list_c_val(FILE* fp)
{
    parameter *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->print_value(fp);
        p_ptr = p_ptr->next;
    }
}

int parameter::advance(unsigned &lvl)
//
// Advance the parameter values. Returns 1 when all combinations are done, 0 otherwise
//
// <lvl> is set to highest nesting level that was advanced. This is used in the summary report
// so that a plot function can use this info to only use data points that are a function of
// the innermost loop(s) having been completed.
//
{
    for (unsigned level = 0; level < nesting_level; level++) {

        if (iterator[level].param_ptr->i_next() == 0) {
            lvl = level;
            return 0;
        }
        
        // Level <level> completed. Need to zero any sums over this level:
        for (zel_element *z_ptr = iterator[level].zel_ptr; z_ptr != nullptr; z_ptr = z_ptr->next)
            z_ptr->e_ptr->zero_sum();
    }

    lvl = nesting_level - 1;
    return 1;
}

parameter::parameter(char* nm, double v_min, double v_max, unsigned n)
// Constructor for a scan-type parameter
{ 
    type = scan;            //This is a scan type parameter
    
    name = nm;
    min_value = v_min;
    max_value = v_max;
    n_steps = n;
    
    expr = nullptr;
    
    next = root;
    root = this;
    
    step_cntr = 0;
    if (n <= 1) d_value = 1.0;          // Does not matter, never used (only 1 step)
    else        d_value = (v_max - v_min) / (double) (n - 1);
    
    iterator[nesting_level++].param_ptr = this;  // This parameter is a looping one
}

parameter::parameter(char* nm, class expression* expr_ptr)
// assignment type parameter
{
    type = assignment;
    
    name = nm;

    min_value = 0.0;        // does not matter...
    max_value = 1.0;
    n_steps = 0;
    step_cntr = 0;
    d_value = 1.0;
    
    expr = expr_ptr;        // This is the only thing that matters!
    
    next = root;
    root = this;
}

parameter *parameter::find_parameter(const char* nm)
{
    parameter *pptr = root;
    while (pptr != nullptr) {
        if (!strcmp(nm, pptr->name))
            return pptr;
        pptr = pptr->next;
    }
    return nullptr;
}

void parameter::enqueue_zero_op(expression* expr_ptr)
{
    zel_element *z_ptr = new zel_element;
    z_ptr->e_ptr = expr_ptr;
    z_ptr-> next = iterator[nesting_level].zel_ptr;
    iterator[nesting_level].zel_ptr = z_ptr;
}


/////////

spice_elements::spice_elements(char* line)
{
    se_type = text;
    next = nullptr;
    txt.text = line;
}

spice_elements::spice_elements(class parameter* par_ptr)
{
    se_type = parameter;
    next = nullptr;
    par.param = par_ptr;
}

spice_elements::spice_elements()
{
    se_type = new_line;
    next = nullptr;
}

void spice_elements::print(FILE* of)
{
    switch(se_type) {
        case text:
            fputs(txt.text, of);
            break;
            
        case parameter:
            fprintf(of, "%.12lg", par.param->get_cur_value());
            break;
            
        case new_line:
            fputc('\n', of);
            break;
        
        default: assert(0);
    }
}


/////////

spice_deck::spice_deck()
//
// Constructor
//
{
    first = nullptr;
    last = nullptr;
}

void spice_deck::add_line(char* txt)
{
    spice_elements *sep = new spice_elements(txt);
    if (first == nullptr)
        first = sep;
    else
        last->add_next(sep);
    last = sep;
}

void spice_deck::add_param_ref(parameter *par)
// adds a parameter reference
{
    spice_elements *sep = new spice_elements(par);
    if (first == nullptr)
        first = sep;
    else
        last->add_next(sep);
    last = sep;
}

void spice_deck::add_nl()
{
    if ((last != nullptr) && (last->is_nl() == 1))
        return;             // Ignore multiple new-lines
    
    spice_elements *sep = new spice_elements();
    if (first == nullptr)
        first = sep;
    else
        last->add_next(sep);
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


////////////

double time_series::peak_search(unsigned int from, unsigned int to, double thr, double eps, unsigned &end)
//
// Find the next peak staring at <from>
//
// Returns: the simulation time of the peak
//
// Failure returns: -1.0 = starting peak is above threshold
//                  -2.0 = ending peak is above threshold
//                  -3.0 = no peak found in interval
//
{
    assert((from < to) && (to < n_data));
    
    if (data[from] >= thr) return -1.0;
    if (data[to] >= thr) return -2.0;
    
    while ((from < to) && (data[from] < (thr + eps))) from++;
    if (from >= to) return -3.0;
    
    unsigned up = from + 1;
    while ((up < to) && (data[up] > (thr - eps))) up++;
    if (up >= to) return -3.0;
    
    end = up;       // Hint for the next peak_search
    
    //
    // The data points under considerations are data[from],....,data[up]
    //
    unsigned n = 1 + up - from;     // #of data points
    if (!(n & 1)) {                 // make it an odd number of data points to simplify fitting
        n++;
        from--;
    }
    unsigned center = (from + up) / 2;

    if (n < 3)
        return (double) center;      // too few data points to attempt peak fit
        
    //
    // Lets fit a parabola
    //
    double S0 = 0.0;
    double S1 = 0.0;
    double S2 = 0.0;
    double S3 = 0.0;
    double S4 = 0.0;
    double T0 = 0.0;
    double T1 = 0.0;
    double T2 = 0.0;
    for (unsigned ind = from; ind <= up; ind++) {
        double i = (double) ind - (double) center;
        S0 += 1.0;
        S1 += i;
        S2 += i*i;
        S3 += i*i*i;
        S4 += i*i*i*i;
        T0 += data[ind];
        T1 += data[ind] * i;
        T2 += data[ind] * i * i;
    }
    
    // Just sanity checking my math - can take these terms out later
    assert(S0 == (double) n);
    assert(S1 == 0.0);
    assert(S3 == 0.0);
    
    double tmp = 2.0 * S2 * (S0 * T2 - S2 * T0);
    if (fabs(tmp) < 1.0e-10)
        return -4.0;                    // fit failed
        
    double imax = - T1*(S0*S4 - S2*S2) / tmp;
    
    return sim_time((double) center + imax);
}

double time_series::edge_search(unsigned int from, unsigned int to, double min_chg,
                                double t_window, unsigned int& e_type, unsigned int& end)
//
// Search for an edge in the time series. An edge is defined as a change in the value of the
// time series that exceeds <min_chg>  (>0) over a time period of <t_window> seconds.
// The search begins at <from> and ends at <to>. The location of the edge is reported in fractional
// time steps, basically a floating point typed index into the array.
//
// The edge is determined by fitting an 3-order polynomical to the edge data and determining the
// inflection point, that is when the curvature changes, which is the 0-crossing point of the second
// derivative.
//
// Upon success, the reported edge type <e_type> is set to 1 for a rising edge and 0 for a falling edge,
// and <end> is set to a datum index that can be used as a starting point to locate the next edge.
//
// Returns the simulation time of the edge (in seconds)
// Failures are reported by a negative result:
// -3.0 : no edge found in interval.
//
{
    static lsq_fit *lsq_fit_sys_ptr = nullptr;      // The LSQ fit
    
    assert((min_chg > 0.0) && (t_window > 0.0) && (sim_time_incr > 1.0e-15));
    unsigned d_win = (unsigned) nearbyint(t_window / sim_time_incr);
    assert((d_win > 0) && (from < to) && (to <=  n_data));
    
    //
    // Search for an edge:
    //
    unsigned e_pos = from;
    do {
        if ((e_pos + d_win) >= to)
            return -3.0;                // Window exceeds data range
        if (fabs(data[e_pos] - data[e_pos + d_win]) >= min_chg)
            break;                      // Found an edge
        e_pos++;
    } while(1);
    
    //
    // To center the search window, a preliminary threshold is determined
    // and where the time series cross this threshold:
    //
    double val_bgn = data[e_pos];
    double val_end = data[e_pos + d_win];
    unsigned is_rising = val_end > val_bgn;       // Determine type of edge
    double threshold = 0.5 * (val_bgn + val_end); // THR = 1/2 way between endpoints
    while (((data[e_pos] < threshold) == is_rising) && ((e_pos + 1) < to))
            e_pos++;
    assert(e_pos < to);
    
    //
    // Now data[e_pos] is the datum just below the threshold and data[e_pos + 1] is
    // the datum just above the threshold. Switch above/below if this is a falling edge.
    // The main point is that the threshold is between the data [e_pos] and [e_pos + 1].
    //
    // Important nitpicking: if the threshold is equal to either datum, then the value intevals
    // differ:
    // rising edge:         [val_begin, thr)[thr, val_end]
    // falling edge:        [val_begin, thr](thr, val_end]
    //
    // Now the window od tata points will be widened symmetrically about the e_pos location.
    // The threshold is refined after each step. The window is widened one data point at a time.
    // The data must remain monotonically ordered. Once the windo cannot be extended without
    // breaking monotonicity or reaching the data range limits, this widening process is complte.
    //
    unsigned win_bgn = e_pos;             // The first datum index within the current window
    unsigned win_end = e_pos + 1;         // The last datum index within the current window
    val_bgn = data[win_bgn];
    val_end = data[win_end];
    
    while (((win_bgn - 1) >= from) && ((win_end + 1) < to)) {  // The window extension loop
        double d_bgn = val_bgn - data[win_bgn - 1];
        double d_end = data[win_end + 1] - val_end;
        if ( ((d_bgn < 0) == is_rising) ||          // end of monotonicity reached ?
             ((d_end < 0) == is_rising)    ) break; // Yes!

        // extend along the steepest gradient first
        if (fabs(d_bgn) > fabs(d_end)) {
            win_bgn--;
            val_bgn = data[win_bgn];
        } else {
            win_end++;
            val_end = data[win_end];
        }
        
        if (fabs(val_bgn - val_end) > edge_search_slope_frac * min_chg)
            break;                       // Done: once enough slope is covered
    }
      
    //
    // Now locate the edge within this search window
    //
    double edge_time = -3.0;
    if ((win_end - win_bgn) >= 3) {
        // we have 4 or more data points. The plan is to fit a 3rd order polynomial to this data
        // and look for the inflection point (2nd derivative == 0).
        
        if (lsq_fit_sys_ptr == nullptr)
            lsq_fit_sys_ptr = new_lsq_fit (1, 4, polynomial_o3);
        int ec = init_lsq_fit (lsq_fit_sys_ptr);            // Get ready for the fit
        assert(ec == 0);
        
        for (unsigned i = win_bgn; i <= win_end; i++) {     // add the data points
            // Note: the dependent variable is time: x=Phase, y=time (or data index)
            //       this rotation of the coordinate system results in a better fit of a 3rd order polynomial
            ec = add_lsq_fit  (lsq_fit_sys_ptr, &(data[i]), (double) i);
            assert(ec == 0);
        }
        
        ec = solve_lsq_fit(lsq_fit_sys_ptr);
        if (ec == 0) {
            // got a solution:
            double x = -coeff_lsq_fit(lsq_fit_sys_ptr, 2) / (3.0 * coeff_lsq_fit(lsq_fit_sys_ptr, 3));
            // x: value at the inflection point
            double edge_pos = eval_lsq_fit (lsq_fit_sys_ptr, &x);
            // edge_pos: fractional position of the edge (should be near e_pos)
            if ((edge_pos > (double) win_bgn) && (edge_pos < (double) win_end)) {
                // Fit result passes sanity check
                edge_time = sim_time(edge_pos);
            }
        }
        if (0) { // Debugging aid, remove once deemed stable: this data is nice to verify the fit
            FILE *of = fopen("q.dat", "w");
            assert(of);
            double dy = (data[win_end] - data[win_bgn]) / 200.0;
            for (double y = data[win_bgn]; y <= data[win_end]; y += dy) {
                double ep = eval_lsq_fit (lsq_fit_sys_ptr, &y);
                fprintf(of, "%.6lg %.6lg %.6lg\n", sim_time(ep), y, (edge_time < sim_time(ep)) ?
                    data[win_end] : data[win_bgn]);
            }
            fclose(of);
        }
    }
    
    if (edge_time < 0.0) {
        // Fit failed, or too few data points. Use the threshold crossing
        // as an approximation to the edge position:
        val_bgn = data[win_bgn];
        val_end = data[win_end];
        threshold = 0.5 * (val_bgn + val_end); // updated THR = 1/2 way between endpoints
        
        double edge_pos = (double) win_bgn + (double) (win_end - win_bgn) * 
                          (threshold - data[win_bgn]) / (data[win_end] - data[win_bgn]);
        edge_time = sim_time(edge_pos);
    }
    
    end = win_end;
    e_type = is_rising;
    return edge_time;
}

time_series::time_series()
{
    max_data = n_init;
    data = (double *) malloc(sizeof(double) * max_data);
    reset();
}

void time_series::set_default_length(unsigned int n)
{
    assert(n > 0);
    n_init = n;
}

void time_series::reset()
{
    n_data = 0;
    v_min = __DBL_MAX__;
    v_max = - __DBL_MAX__;
}


void time_series::add_datum(double value)
{
    if (n_data >= max_data) {
        max_data += max_data / 2;   // add 50%
        data = (double *) realloc(data, sizeof(double) * max_data);
    }
    data[n_data++] = value;
    if (v_max < value) v_max = value;
    if (v_min > value) v_min = value;
}

void time_series::print_all(FILE* fp)
{
    for (unsigned i = 0; i < n_data; i++)
        fprintf(fp, "%.6lg %.12lg\n", sim_time((double) i), data[i]);
}


int find_peaks(time_series *time, time_series *ts, double *pk_loc, unsigned max_pk)
//
// Find the peaks (up to max_pk)
//
// Returns #of peaks found, -1 if too many
//
{
    unsigned i = 0;
    unsigned imax = ts->get_n() - 1;
    
    int n_pks = 0;
    
    do {
        unsigned next;
        double ip = ts->peak_search(i, imax, threshold_frac * ts->get_max(), threshold_hyst * ts->get_max(), next);
        if (ip < 0)
            break;
        unsigned peak = floor(ip);
        assert(peak <= imax);
        
        double f = ip - floor(ip);
        double t = time->get_val(peak) * (1.0 - f) + time->get_val(peak + 1) * f;
        
//    printf(" P@ %.6lg\n", t);
        
        if (n_pks >= max_pk)
            return -1;
        pk_loc[n_pks++] = t;
        
        i = next;
        
    } while (1);
    
    return n_pks;
}

/*
int locate_peaks(time_series *ts[], unsigned n_ts, const char *nm, double* &pka)
{
    int i = find_ts(ts, n_ts, nm);
    if (i < 0)
        return -1;    // no such ts
        
    char buf[128];
    sprintf(buf, "%s.dat", nm);
    FILE *of = fopen(buf, "w");
    for (unsigned j = 0; j < ts[0]->get_n(); j++)
        fprintf(of, "%.9lg %.9lg\n", ts[0]->get_val(j) * 1.0e12, ts[i]->get_val(j));
    fclose(of);
    
    double tmp[max_pks];
    i = find_peaks(ts[0], ts[i], tmp,  max_pks);
    if (i <= 0)
        return -2;    // too many peaks / no peaks
        
    pka = new double[i];
    for(int j = 0; j < i; j++) pka[j] = tmp[j];
    
    sprintf(buf, "%s_pk.dat", nm);
    of = fopen(buf, "w");
    for (unsigned j = 0; j < i; j++)
        fprintf(of, "%.9lg 0.2e-3\n", pka[j] * 1.0e12);
    fclose(of);
    
    return i;
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parser interface functions
//

void define_tran(double t_incr, double t_stop, double t_start, double dT_max)
//
// This is intercepted only to get a conservative estimate on the #of output line to
// allocate the storage for the time series.
// Note: the parser does not cover all flavors of the .tran statement, just what I need right now
//
{
    char buf[128];                  // just a sufficiently large buffer
    sprintf(buf, ".tran %.12lg  %.12lg %.12lg %.12lg\n", t_incr, t_stop, t_start, dT_max);
    add_line_to_spice_deck(strdup(buf));    // Just add the original statement
    
    unsigned n = ceil(t_stop / t_incr);
    if (n < 100)
        fprintf(stderr, "Warning: the .tran statement in the spice deck produces less than 100 rows of output\n");
    
    time_series::set_default_length(n + 16); // adds a little slack
    
    sim_time_start = t_start;
    sim_time_incr = t_incr;
}

void add_line_to_spice_deck(char *string)
// Just a wraper to add a text line to the spice void add_line_to_spice_deck(char* string)
{
    assert(string != nullptr);
    circuit.add_line(string);
}

void add_subst_to_spice_deck(char *p_name)
// Add a parameter reference
{
    parameter *p_ptr = parameter::find_parameter(p_name);
    if (p_ptr == nullptr) {
        fprintf(stderr, "Line %d: reference to undefined parameter '%s'\n", yylineno, p_name);
        yy_n_parse_err++;
    } else
        circuit.add_param_ref(p_ptr);
    
    free(p_name);                   // This was allocated via strdup() in the lexer and is no longer needed
                                    // Yeah, this isn't good practice, so sue me
}

void define_param_scan(char *name, double v_start, double v_stop, double n_steps)
// Defines a parameter scan
{
    if (nullptr != parameter::find_parameter(name)) {
        fprintf(stderr, "Line %d; Duplicate parameter definition for '%s'\n", yylineno, name);
        yy_n_parse_err += 1;
        return;
    };
    
    unsigned n = nearbyint(n_steps);
    new parameter(name, v_start, v_stop, n);
}

void define_monitor(char *name)
{
    nodes_of_interest *n_ptr = nodes_of_interest::find(name);
    if (n_ptr != nullptr) {
        fprintf(stderr, "Line %d: '%s' is already monitored\n", yylineno, name);
        yy_n_parse_err += 1;
        return;
    }
    
    n_ptr = new nodes_of_interest(name);
}

void add_new_line_to_spice_deck()
{
    circuit.add_nl();
}

// expression related functions:

void define_para_expression(char *name, void *expr)
// Defines a new expression type parameter
{
    parameter *p_ptr = parameter::find_parameter(name);
    if (p_ptr != nullptr) {
        fprintf(stderr, "Line %d: redefinition of '%s'\n", yylineno, name);
        yy_n_parse_err += 1;
        return;
    }
    
    p_ptr = new parameter(name, (expression *) expr);
}

void *define_add(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, addition, (expression *) y);
}

void *define_sub(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, subtraction, (expression *) y);
}

void *define_mul(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, multiplication, (expression *) y);
}

void *define_div(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, division, (expression *) y);
}

void *define_eq(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, comp_eq, (expression *) y);
}

void *define_ne(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, comp_ne, (expression *) y);
}

void *define_gt(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, comp_gt, (expression *) y);
}

void *define_ge(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, comp_ge, (expression *) y);
}

void *define_lt(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, comp_lt, (expression *) y);
}

void *define_le(void *x, void *y)
// Just an add op
{
    return new expression((expression *) x, comp_le, (expression *) y);
}

void *define_const(double x)
// Just a constant
{
    expression *e_ptr = new expression(x);
    return e_ptr;
}

void *define_ref(char *name)
{
    parameter *p_ptr = parameter::find_parameter(name);
    if (p_ptr == nullptr) {
        fprintf(stderr, "Line %d: reference to undefined parameter '%s' - needs to be defined first\n", yylineno, name);
        yy_n_parse_err += 1;
        return nullptr;
    }
    
    return new expression(p_ptr);
}

void *define_function1(char *name, void *x)
// Produce a void * define_function1(char* name, void* x)
{
    for (unsigned i = 0; func1_tab[i].name != nullptr; i++)
        if (!strcmp(name, func1_tab[i].name)) {
            expression *ep = new expression(func1_tab[i].func1_ptr, (expression *) x);
            if (func1_tab[i].func1_ptr == nullptr)
                // This is a Sum-operation, that needs to be zero-ed if the current
                // parameter looping level completes
                parameter::enqueue_zero_op(ep);
            return ep;
        }
    fprintf(stderr, "Line %d: undefined function '%s'\n", yylineno, name);
    yy_n_parse_err += 1;
    return nullptr;
}

void *define_function2(char *name, char *ts_name, void *y)
// defines function on time time_series::~time_series()
{
    double      (*f_ptr) (class nodes_of_interest *noi_ptr, double x) = nullptr;
    for (unsigned i = 0; func2_tab[i].name != nullptr; i++)
        if (!strcmp(name, func2_tab[i].name)) {
            f_ptr = func2_tab[i].func2_ptr;
            break;
        }
            
    nodes_of_interest *n_ptr = nodes_of_interest::find(ts_name);
    
    if ((f_ptr == nullptr) || (n_ptr == nullptr)) {
        fprintf(stderr, "Line %d: function or monitored column undefined\n", yylineno);
        yy_n_parse_err += 1;
        return nullptr;
    }

    return new expression(f_ptr, n_ptr, (expression *) y);
}

void *define_test_sel(void *t, void *a, void *b)
// This defines the usual C-stype conditional expression (test) ? a : b
{
    return new expression((expression *) t, (expression *) a, (expression *) b);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Analysis functions
//

double locate_peak(nodes_of_interest *noi_ptr, double x)
//
// The node of interest is expected to be a voltage time series and this fuction will determine the
// time of the x-th peak. <x> is expected to be an integer and peaks are numbered 0,1,2,...
//
// A NAN is returned upon failure
//
{
    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < max_ts));
    time_series *ts_ptr = josim_out_columns[ind];
    assert(ts_ptr != nullptr);
    
    unsigned i = 0;
    unsigned imax = ts_ptr->get_n() - 1;
    do {
        unsigned next;
        double tp = ts_ptr->peak_search(i, imax, threshold_frac * ts_ptr->get_max(),
                                        threshold_hyst * ts_ptr->get_max(), next);
        if (tp < 0)
            return NAN;                 // Fail: no peak found
        
        if (n_peak == 0)
            return tp;        // Success: that is the desired peak

        n_peak--;                       // Not this one ...
        i = next;
        
    } while (1);
}

double count_peaks(nodes_of_interest *noi_ptr, double x)
//
// The node of interest is expected to be a voltage time series and this fuction counts the number
// of peaks.
//
// The argment <x> isn't used for now. May be used to supply threshold fractions later
//
{
    unsigned n_peaks = 0;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < max_ts));
    time_series *ts_ptr = josim_out_columns[ind];
    assert(ts_ptr != nullptr);
    
    unsigned i = 0;
    unsigned imax = ts_ptr->get_n() - 1;
    do {
        unsigned next;
        double tp = ts_ptr->peak_search(i, imax, threshold_frac * ts_ptr->get_max(),
                                        threshold_hyst * ts_ptr->get_max(), next);
        if (tp < 0)
            return (double) n_peaks;     // Fail: no more peaks

        n_peaks++;                       // Got one
        i = next;
        
    } while (1);
    
    return NAN;                         // Should not get here
}

double locate_rise(nodes_of_interest *noi_ptr, double x)
//
// Similar to the peak seach, but locate the rising edge of a phase transition.
// The node of interest ought to be a phase column and a rising transition mean that the phase
// increases by about 1 turn (2PI). Rising transitions are relative, so a stair-case function
// will have multiple rising transitions. The time of the transition is defined to be the
// Inflection point, where the slope stops rising and starts falling: zero cossing of the
// second derivative.
//
{
    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < max_ts));
    time_series *ts_ptr = josim_out_columns[ind];
    assert(ts_ptr != nullptr);
   
    unsigned imin = 0;
    unsigned imax = ts_ptr->get_n() - 1;
    do {
        unsigned next;
        unsigned edge_type;
        double tx = ts_ptr->edge_search(imin, imax, edge_search_min_chg * 2.0 * M_PI, edge_search_t_win,
                                        edge_type, next);
        if (tx < 0.0)
            break;                      // Fail: no edge found
        imin = next;                    // Prepare next sarch
        
        if (edge_type != 1)
            continue;                   // Wrong kind f edge
        
        if (n_peak == 0)
            return tx;                  // Success: that is the desired peak

        n_peak--;                       // Not this one ...
        
    } while (1);   
    
    return NAN;
}

double locate_fall(nodes_of_interest *noi_ptr, double x)
{
    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < max_ts));
    time_series *ts_ptr = josim_out_columns[ind];
    assert(ts_ptr != nullptr);
   
    unsigned imin = 0;
    unsigned imax = ts_ptr->get_n() - 1;
    do {
        unsigned next;
        unsigned edge_type;
        double tx = ts_ptr->edge_search(imin, imax, edge_search_min_chg * 2.0 * M_PI, edge_search_t_win,
                                        edge_type, next);
        if (tx < 0.0)
            break;                      // Fail: no edge found
        imin = next;                    // Prepare next sarch
        
        if (edge_type != 0)
            continue;                   // Wrong kind f edge
        
        if (n_peak == 0)
            return tx;                  // Success: that is the desired peak

        n_peak--;                       // Not this one ...
        
    } while (1);   
    
    return NAN;
}

double locate_edge(nodes_of_interest *noi_ptr, double x)
// Dito, but locates any edge (rising or falling
{
    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < max_ts));
    time_series *ts_ptr = josim_out_columns[ind];
    assert(ts_ptr != nullptr);
   
    unsigned imin = 0;
    unsigned imax = ts_ptr->get_n() - 1;
    do {
        unsigned next;
        unsigned edge_type;
        double tx = ts_ptr->edge_search(imin, imax, edge_search_min_chg * 2.0 * M_PI, edge_search_t_win,
                                        edge_type, next);
        if (tx < 0.0)
            break;                      // Fail: no edge found
        imin = next;                    // Prepare next sarch
        
        if (n_peak == 0)
            return tx;                  // Success: that is the desired peak

        n_peak--;                       // Not this one ...
        
    } while (1);   
    
    return NAN;
}

double count_edges(nodes_of_interest *noi_ptr, double x)
// Counts the number of edges (regardless of direction)
//
// <x> isn't used now, but may be used later to supply edge search parameters
{
    unsigned n_edges = 0;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < max_ts));
    time_series *ts_ptr = josim_out_columns[ind];
    assert(ts_ptr != nullptr);
   
    unsigned imin = 0;
    unsigned imax = ts_ptr->get_n() - 1;
    do {
        unsigned next;
        unsigned edge_type;
        double tx = ts_ptr->edge_search(imin, imax, edge_search_min_chg * 2.0 * M_PI, edge_search_t_win,
                                        edge_type, next);
        // Note: it is possible to save some time by writing a simpler edge finder that doesn't bother
        //       with fitting, interpolating and timing...
        
        if (tx < 0.0)
            return (double) n_edges;    // Done: no more edges
        imin = next;                    // Prepare next sarch

        n_edges++;                      // Got one more edge
        
    } while (1);   
    
    return NAN;                         // unreachable
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The main program
//

int main(int argc, char *argv[])
{
    // Create the named pipes (FIFO) (unless already present)
    if ((access(csv_out_path, F_OK) != 0) && (mkfifo(csv_out_path, 0666) == -1)) {
        perror("mkfifo failed");
        exit(1);
    }
    
    if ((access(cir_inp_path, F_OK) != 0) && (mkfifo(cir_inp_path, 0666) == -1)) {
        perror("mkfifo failed");
        exit(1);
    }
    
    if (argc != 2) {
        fprintf(stderr, "csv_analyzer <spice source file>\n");
        exit(1);
    }
    
    if (1 > circuit.read_cir_file(argv[1])) {
        fprintf(stderr, "Failed to read spice deck from '%s'\n", argv[1]);
        exit(1);
    }
        
    FILE *sum_fp = fopen("Summary.dat", "w");
    assert(sum_fp != nullptr);
    fprintf(sum_fp, "#");
    parameter::list_names(sum_fp);
    fprintf(sum_fp, " level\n");

    //
    // Initialization and setup done.
    //
    while (1) {
        //
        // Do the actual work:
        //   1. Run JoSIM with the edited circuit file
        //   2. Digest the JoSim output
        //   3. Analyze the result
        //   4. Repeat with different parameter values
        //
        
        ///////////////////////////////////////////////////////////////////////////////////////////
        //  1   Run JoSim                                                                            //
        ///////////////////////////////////////////////////////////////////////////////////////////
        pid_t child_pid = fork();
        if (child_pid == 0) {
            //
            // This is the child process
            //
            freopen("run.log", "w", stdout);    // Keep this file for debugging purposes
        
            int ie = execl("/usr/local/bin/josim-cli", "-a", "1", "-o", csv_out_path, cir_inp_path, nullptr);
            perror("Failed to stat Josim");
            exit(1);
        }
    
        if (circuit.write_cir_file(cir_inp_path)) {
            fprintf(stderr, "failed to write the spice deck to JoSim\n");
            exit(1);
        }
    
        ///////////////////////////////////////////////////////////////////////////////////////////
        //  2   Parse the output from JoSim                                                          //
        ///////////////////////////////////////////////////////////////////////////////////////////
        FILE *inpf = fopen(csv_out_path, "r");
        assert(inpf != nullptr);

        //
        // Process first line of CVS file
        //
        if (!fgets(josim_output_buf, max_ln_length - 1, inpf) ||
            strlen(josim_output_buf) > (max_ln_length -4)) {            // Meant to detect buffer overflow
            fprintf(stderr, "Failed to read first line\n");
            exit(1);
        }

        char *buf_p;
        {
            //
            // Discard white space and doublle quotes
            //
            char *tmp = josim_output_buf;
            for (buf_p = josim_output_buf; *buf_p; buf_p++)
                if ((*buf_p != '"') && (*buf_p != ' ') && (*buf_p != '\t') && (*buf_p != '\n'))
                    *tmp++ = *buf_p;
            *tmp = 0;
        }
    
        buf_p = josim_output_buf;
        for (unsigned ic = 0; ic < max_ts; ic++) {
            char *tmp = buf_p;
            while (*tmp && (*tmp != ',')) tmp++;    // Find end of string
            int last = !*tmp;       // Set if this is the last colun
            *tmp = 0;               // Replance comma with 0
            //
            // buf_p is now pointing to the name string for this column
            //
            nodes_of_interest *noi_ptr = nodes_of_interest::find(buf_p);
            if (n_josim_runs == 0) {
                // First Josim run: allocate the ts for the NOI's
                if (noi_ptr == nullptr) {
                    // Don't care about this column
                    josim_out_columns[ic] = nullptr;
                } else {
                    josim_out_columns[ic] = new time_series();
                    if (noi_ptr->set_col_index(ic)) {
                        fprintf(stderr, "Josim output has multiple columns named '%s'\n", tmp);
                        exit(1);
                    }
                }
            } else {
                // Verify that the column order is the same as the first running
                if (noi_ptr != nullptr) {
                    if (ic != noi_ptr->get_col_index()) {
                        fprintf(stderr, "Josim run %u: column order has changed\n", n_josim_runs + 1);
                        exit(1);
                    }
                    josim_out_columns[ic]->reset();         // Reset the TS to accept new data
                }
            }
            
            if (last) {
                if (n_josim_runs == 0) {
                    n_josim_out_columns = ic + 1;
                    printf("Detected %u columns\n", n_josim_out_columns);
                }
                break;
            }
            
            if ((ic + 1) >= max_ts) {
                fprintf(stderr, "Josim output: too many columns - increase max_ts\n");
                exit(1);
            }
            
            buf_p = tmp + 1;
        }
        
        if (n_josim_runs == 0) {
            // Verify that all NOI's are defined
            nodes_of_interest *n_ptr = nodes_of_interest::find_undef();
            if (n_ptr != nullptr) {
                fprintf(stderr, "Node of interest '%s' not found in Josim output\n", n_ptr->get_name());
                exit(1);
            }
        }

        //
        // Now read the data (in CSV format)
        //
        int line_no = 1;
        do {
            line_no++;
            if (!fgets(josim_output_buf, max_ln_length - 1, inpf))
                break;                             // Input file exhausetd
            buf_p = josim_output_buf;
            double val;
            int nc;
            for (int i = 0; i < n_josim_out_columns; i++) {
                if (1 > sscanf(buf_p,"%lf%n", &val, &nc)) {
                    fprintf(stderr, "Josim run %u: Failed to digest line %d\n", n_josim_runs + 1, line_no);
                    exit(1);
                }
                if (i == 0) {
                    // Verify time column
                    if (fabs(val - sim_time((double)(line_no - 2))) > 0.1e-12) {
                        fprintf(stderr, "Josim output on line %d: time expected %.10lg time found %.10lg\n",
                                line_no, sim_time((double)(line_no - 2)), val);
                        exit(1);
                    }
                }
                if (josim_out_columns[i] != nullptr)
                    josim_out_columns[i]->add_datum(val);   // Selectively add data
                buf_p += nc +1;
            }
        } while (1);
        
        fclose(inpf);
        
        int status;
        waitpid(child_pid, &status, 0);         // reap the child process
    
        if (n_josim_runs == 0) printf("Read %d rows of data, exit status=%d\n", line_no, status);
    
        ///////////////////////////////////////////////////////////////////////////////////////////
        // 3 Analyze the simulation results                                                      // 
        ///////////////////////////////////////////////////////////////////////////////////////////
        //
        // Analysis is primarily done via parameter expressions
        //
//        nodes_of_interest::print_all();
        parameter::list_c_val(sum_fp);

        unsigned level;
        n_josim_runs++;                         // Done with this run
        unsigned all_done = parameter::advance(level);
        fprintf(sum_fp, " %u\n", level);
        
        if (all_done)
            break;
    }
    
    fclose(sum_fp);

    exit(0);
}

void misc()
{
    //
    // Just some code from previous iterations that may be useful
    //
/*
        double slope = 0.0;
        unsigned n_sp = (min_n_pks < n_slope_pks) ? min_n_pks : n_slope_pks;
        if (n_sp >= 2) {
            //
            // If there are more than 2 points, we can attempt to compute a slope
            //
            double Sy = 0.0;
            double Sx = 0.0;
            double Sxx = 0.0;
            double Sxy = 0.0;
            for (unsigned i = 0; i < n_sp; i++) {
                double x = top[0].t_pks[i] - top[1].t_pks[i];   // Time difference N5 vs. N5A
                double y = top[8].t_pks[i] - top[9].t_pks[i];   // Time difference N9 vs. N9A
                Sx += x;
                Sy += y;
                Sxx += x*x;
                Sxy += x*y;
            }
            //
            // We are fitting a linear expression y(i) = a + b * x(i)
            //
            double t = (double) n_sp * Sxx - Sx * Sx;
            if (fabs(t) > 1.0e-30) {
                //
                // There seems to be a solution:
                // Otherwise, we juset leave the slope at 0
                //
                double a = (Sy * Sxx - Sx * Sxy) / t;
                double b = ((double) n_sp * Sxy - Sx * Sy) / t;
                //
                // of interest is only the slope <b> here. <a> shouble be near 0,
                // it is computed only to aid debugging
                //
                slope = b;
            }
        }
        
        double   min_dly = __DBL_MAX__,   avg_dly = 0.0,   max_dly = - __DBL_MAX__;
        double A_min_dly = __DBL_MAX__, A_avg_dly = 0.0, A_max_dly = - __DBL_MAX__;
        for (unsigned i = 0; i < min_n_pks; i++) {
            double t = top[8].t_pks[i] - top[0].t_pks[i];
            min_dly = fmin(t, min_dly);
            avg_dly += t;
            max_dly = fmax(t, max_dly);
            
            t = top[9].t_pks[i] - top[1].t_pks[i];
            A_min_dly = fmin(t, A_min_dly);
            A_avg_dly += t;
            A_max_dly = fmax(t, A_max_dly);
        }
        avg_dly   =   avg_dly / (double) min_n_pks;
        A_avg_dly = A_avg_dly / (double) min_n_pks;
        
        parameter::list_c_val(sum_fp);
        fprintf(sum_fp, " %.9lg %.9lg %.9lg %.9lg %.9lg %.9lg %.9lg\n", slope,
                  min_dly * 1.0e12,   avg_dly * 1.0e12,   max_dly * 1.0e12,
                A_min_dly * 1.0e12, A_avg_dly * 1.0e12, A_max_dly * 1.0e12);
*/
}

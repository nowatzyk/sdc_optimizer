#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

extern "C" {
#include "parser.h"                         // Bison generated headers (in build dir)
#include "lex.yy.h"                         // Flex-generated header (in build dir)
#include "parser_interf.h"                  // needed to integrate the parser
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
    
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class functions
//

int nodes_of_interest::find(const char* name)
// Finds the NOI by name or returns -1
{
    for (nodes_of_interest *noi_ptr = root; noi_ptr != nullptr; noi_ptr = noi_ptr->next)
        if (!strcmp(name, noi_ptr->name))
            return noi_ptr->va_index;
    return -1;
}

nodes_of_interest::nodes_of_interest(const char* name)
// Constructor
{
    next = root;
    root = this;
    va_index = n_noi++;
    name = strdup(name);                    // Yeah, I should not mix malloc() and new, but I don't care right now.
}

//////

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
    fprintf(fp, " %.15lg", min_value + (double) step_cntr * d_value);
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

int parameter::advance()
//
// Advance the parameter values. Returns 1 when all combinations are done, 0 otherwise
//
{
    parameter *p_ptr = root;
    while(p_ptr != nullptr) {
        if (p_ptr->i_next() == 0)
            return 0;
        p_ptr = p_ptr->next;
    }
    
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
    
    next = root;
    root = this;
    
    step_cntr = 0;
    if (n <= 1) d_value = 1.0;          // Does not matter, never used (only 1 step)
    else        d_value = (v_max - v_min) / (double) (n - 1);
}



parameter * parameter::find_parameter(const char* nm)
{
    parameter *pptr = root;
    while (pptr != nullptr) {
        if (!strcmp(nm, pptr->name))
            return pptr;
        pptr = pptr->next;
    }
    return nullptr;
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
// Returns: the frational peak location. data points start at location 0
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
    
    return (double) center + imax;
}


time_series::time_series(char* nm)
{
    name = nm;
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

int find_ts(time_series *ts[], unsigned n_ts, const char *s)
{
    for (int i = 0; i < n_ts; i++)
        if (!strcmp(ts[i]->get_name(), s))
            return i;
    return -1;
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
    new time_series(name);
}

void add_new_line_to_spice_deck()
{
    circuit.add_nl();
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
    
    circuit.write_cir_file("q.txt");
    
/*
    char *buf = new char[max_ln_length];
    time_series *ts[max_ts];
    for (unsigned i = 0; i < max_ts; i++) ts[i] = nullptr;
    unsigned n_ts = 0;

    t_of_peaks top[n_noi];
    for (int i = 0; i < n_noi; i++)
        top[i].t_pks = nullptr;
    
    FILE *sum_fp = fopen("Summary.dat", "w");
    assert(sum_fp != nullptr);
    fprintf(sum_fp, "#");
    parameter::list_names(sum_fp);
//    fprintf(sum_fp, " slop(N5/N9) min_delay avg_delay max_delay A:min_delay A:avg_delay A:max_delay\n");
    fprintf(sum_fp, " peak-times for");
    for (unsigned i; i < n_noi; i++)
        fprintf(sum_fp, " %s", nodes_of_interest[i]);
    fprintf(sum_fp, "\n");

    //
    // Initialization and setup done.
    //
    for (unsigned n_sim = 0; 1; n_sim++) {
        n_ts = 0;
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
    
        if (sd.write_cir_file(cir_inp_path)) {
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
        if (!fgets(buf, max_ln_length - 1, inpf) || strlen(buf) > (max_ln_length -4)) {
            fprintf(stderr, "Failed to read first line\n");
            exit(1);
        }

        char *buf_p;
        {
            //
            // Discard white space and doublle quotes
            //
            char *tmp = buf;
            for (buf_p = buf; *buf_p; buf_p++)
                if (*buf_p != '"' && *buf_p != ' ' && *buf_p != '\t' && *buf_p != '\n')
                    *tmp++ = *buf_p;
            *tmp = 0;
        }
    
        buf_p = buf;
        int last;
        do {
            char *tmp = buf_p;
            while (*tmp && (*tmp != ',')) tmp++;    // Find end of string
            last = !*tmp;
            *tmp = 0;
            if (ts[n_ts] != nullptr) ts[n_ts]->reset(buf_p);
            else                     ts[n_ts] = new time_series(buf_p);
            n_ts++;
            buf_p = tmp + 1;
        } while (last == 0);
    
        // printf("Detected %u columns\n", n_ts);
    
        //
        // Now read the data (in CSV format)
        //
        int line_no = 1;
        do {
            line_no++;
            if (!fgets(buf, max_ln_length - 1, inpf))
                break;                             // Input file exhausetd
            buf_p = buf;
            double val;
            int nc;
            for (int i = 0; i < n_ts; i++) {
                if (1 > sscanf(buf_p,"%lf%n", &val, &nc)) {
                    fprintf(stderr, "Failed to digest line %d\n", line_no);
                    exit(1);
                }
                ts[i]->add_datum(val);
                buf_p += nc +1;
            }
        } while (1);
        
        fclose(inpf);
    
        // printf("Read %d rows of data\n", line_no);
    
        ///////////////////////////////////////////////////////////////////////////////////////////
        // 3 Analyze the simulation results                                                      // 
        ///////////////////////////////////////////////////////////////////////////////////////////
        //
        // Save the time delta...
        //
        unsigned min_n_pks = ~0u, max_n_pks = 0;
        for (unsigned i = 0; i < n_noi; i++) {
            //
            // For each node of interest, locate all the peaks in the time series
            //
            if (top[i].t_pks != nullptr)
                delete[] top[i].t_pks;
            top[i].t_pks = nullptr;
            top[i].n_pks = 0;
            int j = locate_peaks(ts, n_ts, nodes_of_interest[i], top[i].t_pks);
            if (j <= 0) {
                // Failed to locate peaks:
                min_n_pks = 0;
                top[i].n_pks = 0;
                continue;
            }
            top[i].n_pks = j;
            if (min_n_pks > j)
                min_n_pks = j;
            if (max_n_pks < j)
                max_n_pks = j;
        }
        // Note: if the min/max #of peaks are not the same, the simulation is suspect!
        //       this usually happens when some pules do not traverse the jtl until the end
        //       or if the delay exceeds the period so that at the end of the simulation
        //       some peaks are beyond the simulation time limit.
        
        
        for (unsigned i = 0; i < n_noi; i++) {
            //
            // Print the first peack > 50ps
            //
            unsigned pk_found = 0;
            for (unsigned j = 0; j < top[i].n_pks; j++) {
                double t = top[i].t_pks[j] * 1.0e12;
                if (t > 50.0) {
                    // peak is past 50ps (initialization phase)
                    fprintf(sum_fp, " %.15lg", t);
                    pk_found = 1;
                    break;
                }
            }
            
            if (!pk_found) fprintf(sum_fp, " 0");
        }
        fprintf(sum_fp, "\n");
*/
/*
        sprintf(buf, "delay_%u.dat", n_sim);
        FILE *of = fopen(buf, "w");
        assert(of);
        fprintf(of, "#");
        parameter::list_names(of);
        fprintf(of, "\n#");
        parameter::list_c_val(of);
        fprintf(of, "\n");
    
        for(unsigned i = 0; i < min_n_pks; i++) {
            fprintf(of, "%.9lg %.9lg %.9lg %.9lg %.9lg %.9lg %.9lg\n", 
                    (top[0].t_pks[i] - top[1].t_pks[i]) * 1.0e12,
                    (top[2].t_pks[i] - top[3].t_pks[i]) * 1.0e12,
                    (top[4].t_pks[i] - top[5].t_pks[i]) * 1.0e12,
                    (top[6].t_pks[i] - top[7].t_pks[i]) * 1.0e12,
                    (top[8].t_pks[i] - top[9].t_pks[i]) * 1.0e12,
                    (top[8].t_pks[i] - top[0].t_pks[i]) * 1.0e12,
                    (top[9].t_pks[i] - top[1].t_pks[i]) * 1.0e12
            );
        }
   
        fclose(of);
        
        //
        // Overall analysis:
        //
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
/*
        
        if (parameter::advance())
            break;
    }
    
    fclose(sum_fp);
*/
    exit(0);
}

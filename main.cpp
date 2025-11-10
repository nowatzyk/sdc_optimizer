#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

const unsigned max_ln_length = 100000;      // We expect some very long lines...
const unsigned max_ts = 128;                // Max number of time series

const double threshold_frac = 0.6;         // Peak search threshold, 75% of max 
const double threshold_hyst = 0.01;         // Hysteresis to avoid noise induced false peak locations

const unsigned max_pks = 1024;              // needed to allocate arrays for peaks:
                                            // Peaks are located first and then analyzed
                                            
const char* csv_out_path = "/tmp/JoSim2csv_analyzer.csv";      // Output fifo from Josim to this program
const char* cir_inp_path = "/tmp/JoSim2csv_analyzer.cir";      // Input fifo to Josim

const unsigned spice_src_max_char = 1024;   // Max line length in spice source file
const char spice_escape = '~';              // Escape character in spice decks for pragma substitutions

const unsigned n_slope_pks = 4;             // #of data points to be used in slope determination

struct units {
    char        c;                          // character representation
    double      val;                        // Value
} unit_table[] = {
    {' ', 1.0},                             // Nothing: unity
    {'m', 1.0e-3},                          // milli
    {'u', 1.0e-6},                          // micro
    {'n', 1.0e-9},                          // nano
    {'p', 1.0e-9},                          // pico
    {'L', 125.0e-6},                        // Liks, unit of 125 micro-Ampere
    {'O', 2.632e-12},                       // oHenry
    {'T', 2.0 * M_PI},                      // Turns
    {0, 1.0}
};

const unsigned n_noi = 10;                  // #of nodes of interest
const char *nodes_of_interest[n_noi] = {
    "V(N5)", "V(N5A)", "V(N6)", "V(N6A)", "V(N7)", "V(N7A)", "V(N8)", "V(N8A)", "V(N9)", "V(N9A)"
};

struct t_of_peaks {
    unsigned    n_pks;                      // %of peaks
    double      *t_pks;                     // Time of peaks
};

class pragma {
    char        *name;
    double      min_value, max_value;       // Range of values to step through
    unsigned    n_steps;                    // #of n_steps
    double      d_value;                    // Step size
    unsigned    unit;                       // Multiplyer to be applied to value before sending to JoSim
    
    unsigned    step_cntr;                  // Step counter
    pragma      *next;                      // List pointer
    static pragma *root;                    // Anchor
    void         i_reset() {step_cntr = 0;};    // Duh!
    int          i_next()  {if (++step_cntr >= n_steps) {n_steps = 0; return 1;}; return 0;};
    void         print_name(FILE *fp);
    void         print_value(FILE *fp);
    
    pragma(char *nm, double v_min, double v_max, unsigned n, unsigned u);   // private constructor
public:
    
    static void reset();                    // Go back to initial state
    static int  advance();                  // Advance to next value combination, returns != 0 when exhausted
    static void list_names(FILE *fp);       // adds names to file (preceeded by space, followed by nothing)
    static void list_c_val(FILE *fp);       // list the current values, like above
    
    static int  new_pragma(char *def);      // Creates a pragma
    static pragma *find_pragma(char *nm);   // find pragma by its name (symbol)
    double      get_cur_value() {return (min_value + (double) step_cntr * d_value) * unit_table[unit].val;};
};

pragma* pragma::root = nullptr;             // start of the list

void pragma::print_name(FILE* fp)
{
    if (unit > 0) fprintf(fp, " %s", name);
    else          fprintf(fp, " %s[%c]", name, unit_table[unit].c);
}

void pragma::print_value(FILE* fp)
{
    fprintf(fp, " %.9lf", min_value + (double) step_cntr * d_value);
}

void pragma::reset()
{
    pragma *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->i_reset();
        p_ptr = p_ptr->next;
    }
}

void pragma::list_names(FILE* fp)
{
    pragma *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->print_name(fp);
        p_ptr = p_ptr->next;
    }
}

void pragma::list_c_val(FILE* fp)
{
    pragma *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->print_value(fp);
        p_ptr = p_ptr->next;
    }
}

int pragma::advance()
//
// Advance the pragma values. Returns 1 when all combinations are done, 0 otherwise
//
{
    pragma *p_ptr = root;
    while(p_ptr != nullptr) {
        if (p_ptr->i_next() == 0)
            return 0;
    }
    
    return 1;
}

pragma::pragma(char* nm, double v_min, double v_max, unsigned n, unsigned u)
// Constructor
{
    name = strdup(nm);
    min_value = v_min;
    max_value = v_max;
    n_steps = n;
    unit = u;
    
    next = root;
    root = this;
    
    step_cntr = 0;
    d_value = (v_max - v_min) / (double) n;
}

int pragma::new_pragma(char* def)
//
// Creates a new pragma
//
{
    char buf[spice_src_max_char];
    double v_min, v_max;
    unsigned n;
    
    char u1, u2;
    unsigned unit;
    if (6 != sscanf(def, "%s%lf%c%lf%c%u", buf, &v_min, &u1, &v_max, &u2, &n))
        return -1;                          // missing or wrong parameters
    if (u1 != u2)
        return -2;                          // Units must be the same 
    for (unit = 0; unit_table[unit].c != 0; unit++)
        if (unit_table[unit].c == u1)
            break;
    if (unit_table[unit].c == 0)
        return -3;                          // Undefined unit character
    if ((v_min >= v_max) || (n < 1) || (strlen(buf) < 1))
        return -4;                         // bad values
    if (find_pragma(buf) != nullptr)
        return -5;                          // already defined
        
    new pragma(buf, v_min, v_max, n, unit);
    return 0;
}

pragma * pragma::find_pragma(char* nm)
{
    pragma *pptr = root;
    while (pptr != nullptr) {
        if (!strcmp(nm, pptr->name))
            return pptr;
        pptr = pptr->next;
    }
    return nullptr;
}




class spice_deck {
   unsigned     n_lines;                    // #of lines
   unsigned     max_lines;                  // just for allocation purposes
   char         **src;                      // The source file

public:
    spice_deck();                           // constructor
    
    int read_cir_file(const char *fn);      // Reads a spice deck from a file into memory
    int write_cir_file(const char *fn);     // writes the spice deck to a file
};

spice_deck::spice_deck()
//
// Constructor
//
{
    n_lines = 0;
    max_lines = 1024;
    src = (char **) malloc(sizeof(char *) * max_lines);
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
    FILE *inp = fopen(fn, "r");
    if (inp == nullptr)
        return -1;
    
    char buf[spice_src_max_char];
    unsigned n_pragma = 1;              // Just to keep the lines numbers straight
    
    while(fgets(buf, spice_src_max_char - 1, inp)) {

        if (strlen(buf) > (spice_src_max_char - 2))  {
            fprintf(stderr, "Line %u in %s is too long (increase spice_src_max_char)\n", n_lines + n_pragma, fn);
            exit(1);
        }

        if (!strncmp(buf, "*Pragma", 7)) {
            if (pragma::new_pragma(buf + 7)) {
                fprintf(stderr, "Bad pragma def in line %u of %s\n", n_lines + n_pragma, fn);
            }
            n_pragma++;
            continue;
        }
        
        if (n_lines >= max_lines) {
            max_lines += max_lines / 2;         // Add 50%
            src = (char **) realloc(src, sizeof(char *) * max_lines);
        }
        src[n_lines++] = strdup(buf);
    }
    
    fclose(inp);
    
    return n_lines;
}

int spice_deck::write_cir_file(const char* fn)
//
// Write the spice deck to a file
//
{
    FILE *out = fopen(fn, "w");
    if (out == nullptr)
        return -1;
    
    char buf[2*spice_src_max_char];     // NOTE - bug lurking here:
    // pragma substitutions can make the line longer and may cause an buffer overflow.
    // Not likely because few substitution are made, but I'm too lazy right now to make this bullet proof.
    
    for (unsigned i = 0; i < n_lines; i++) {
        char *cs_ptr = src[i];
        char *cd_ptr = buf;
        while (*cs_ptr) {
            if (*cs_ptr == spice_escape) {
                cs_ptr++;
                unsigned i = 0;
                while (cs_ptr[i] != spice_escape) {
                    if (cs_ptr[i] == 0) {
                        fprintf(stderr, "Unmatched pragma references: ~some_name~ , trailing ~ missing in circuit file\n");
                        exit(1);
                    }
                    cd_ptr[i] = cs_ptr[i];
                    i++;
                }
                cd_ptr[i] = 0;
                pragma *p_ptr = pragma::find_pragma(cd_ptr);
                if (p_ptr == nullptr) {
                    fprintf(stderr, "Undefined pragma reference in *.cir file: '%s'\n", cd_ptr);
                    exit(1);
                }
                int n;
                sprintf(cd_ptr, " %.9lg %n", p_ptr->get_cur_value(), &n);
                
                cs_ptr += i + 1;
                cd_ptr += n;
            } else
                *cd_ptr++ = *cs_ptr++;
        }
        *cd_ptr = 0;
        
        if (0 > fputs(buf, out)) {    // some problem
            fclose(out);
            return -2;
        }
    }
    
    fclose(out);
    return 0;
}




class time_series {
    char    *name;          // Name of this coumn
    double  *data;
    unsigned    n_data;
    unsigned    max_data;
    double      v_min, v_max;
public:
    time_series(const char *name);
    ~time_series();
    
    void        reset(const char *nm);          // Reset the time series and give it a new (same) name
    void add_datum(double);
    double get_max() {return v_max;};
    double get_min() {return v_min;};
    char *get_name() {return name;};
    double get_val(unsigned i) {assert(i < n_data); return data[i];};
    
    unsigned get_n() {return n_data;};
    
    double peak_search(unsigned from, unsigned to, double thr, double eps, unsigned &end);
};

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


time_series::time_series(const char* nm)
{
    name = strdup(nm);
//    printf(">>%s<<\n", name); 
    max_data = 128;
    n_data = 0;
    v_min = __DBL_MAX__;
    v_max = - __DBL_MAX__;
    data = (double *) malloc(sizeof(double) * max_data);
}

void time_series::reset(const char* nm)
{
    free(name);
    name = strdup(nm);
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
    
    spice_deck sd;
    if (!sd.read_cir_file(argv[1])) {
        fprintf(stderr, "Failed to read spice deck from '%s'\n", argv[1]);
        exit(1);
    }
    
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
    pragma::list_names(sum_fp);
    fprintf(sum_fp, " slop(N5/N9) min_delay avg_delay max_delay A:min_delay A:avg_delay A:max_delay\n");

    //
    // Initialization and setup done.
    //
    for (unsigned n_sim = 0; 1; n_sim++) {
        //
        // Do the actual work:
        //   1. Run JoSIM with the edited circuit file
        //   2. Digest the JoSim output
        //   3. Analyze the result
        //   4. Repeat with different pragma values
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
            top[i].n_pks = locate_peaks(ts, n_ts, nodes_of_interest[i], top[i].t_pks);
            if (min_n_pks > top[i].n_pks)
                min_n_pks = top[i].n_pks;
            if (max_n_pks < top[i].n_pks)
                max_n_pks = top[i].n_pks;
        }
        // Note: if the min/max #of peaks are not the same, the simulation is suspect!
        //       this usually happens when some pules do not traverse the jtl until the end
        //       or if the delay exceeds the period so that at the end of the simulation
        //       some peaks are beyond the simulation time limit.

        sprintf(buf, "delay_%u.dat", n_sim);
        FILE *of = fopen(buf, "w");
        assert(of);
        fprintf(of, "#");
        pragma::list_names(of);
        fprintf(of, "\n#");
        pragma::list_c_val(of);
        fprintf(of, "\n");
    
        for(unsigned i = 0; i < min_n_pks; i++) {
            fprintf(of, "%.9lg %.9lg %.9lg %.9lg %.9lg\n", 
                    (top[0].t_pks[i] - top[1].t_pks[i]) * 1.0e12,
                    (top[2].t_pks[i] - top[3].t_pks[i]) * 1.0e12,
                    (top[4].t_pks[i] - top[5].t_pks[i]) * 1.0e12,
                    (top[6].t_pks[i] - top[7].t_pks[i]) * 1.0e12,
                    (top[8].t_pks[i] - top[9].t_pks[i]) * 1.0e12);
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
        
        pragma::list_c_val(sum_fp);
        fprintf(sum_fp, " %.9lg %.9lg %.9lg %.9lg %.9lg %.9lg %.9lg\n", slope,
                  min_dly * 1.0e12,   avg_dly * 1.0e12,   max_dly * 1.0e12,
                A_min_dly * 1.0e12, A_avg_dly * 1.0e12, A_max_dly * 1.0e12);
        
        
        if (pragma::advance())
            break;
    }
    
    fclose(sum_fp);

    exit(0);
}

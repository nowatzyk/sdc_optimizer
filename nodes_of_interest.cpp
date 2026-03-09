///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Nodes of interests referres to circuit nodes that are defined in the spice circuit description
// file. The spice deck will have a number of ".print" directives, each of which produces one 
// column in the csv output. The names are taken from the first line of the csv output from JoSIM.
// 
// Note: JoSIM folds all names to upper-case. This framework does not. (maybe it should?)
//

#include "csv_analyzer.h"           // for threshold_frac, threshold_hyst, edge_search constants
#include "nodes_of_interest.h"

extern "C" {
#include "lsq_fit.h"
#include "fit_functions.h"
}

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <cassert>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Globals
//

vector <time_series *> josim_out_columns;
unsigned     n_josim_runs        = 0;
double       sim_time_start      = -__DBL_MAX__;
double       sim_time_incr       = 0.0;

const double edge_search_min_chg = 0.6667;  // For an edge search on a phase, this is the min change required
                                            // to be considered an edge (in fractions of a turn, = 2PI phase)
const double edge_search_t_win = 10.0e-12;  // max edge transition time 10 pico-seconds
const double edge_search_slope_frac = 1.0;  // Defines the fraction of the min_chg (see above) that is used to
                                            // to determine the transition point

const double threshold_frac = 0.6;          // Peak search threshold, 75% of max 
const double threshold_hyst = 0.01;         // Hysteresis to avoid noise induced false peak locations

std::vector<nodes_of_interest *> nodes_of_interest::all_noi_s;  //The collection of NOI's 

double sim_time(double n)
// Returns the simulation time for the n-th row
// Note: row numbers may be fractional
{
    return sim_time_start + sim_time_incr * n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// time_series -- statics
//

unsigned     time_series::n_init = 1024;
time_series *time_series::root   = nullptr;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// time_series -- construction / reset
//

time_series::time_series()
{
    max_data = n_init;
    data     = (double *) malloc(sizeof(double) * max_data);
    reset();
}

time_series::~time_series()
{
    free(data);
}

void time_series::reset()
{
    n_data = 0;
    v_min  =  __DBL_MAX__;
    v_max  = -__DBL_MAX__;
}

void time_series::set_default_length(unsigned n)
{
    assert(n > 0);
    n_init = n;
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

void time_series::print_all(FILE *fp)
{
    for (unsigned i = 0; i < n_data; i++)
        fprintf(fp, "%.6lg %.12lg\n", sim_time((double) i), data[i]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// time_series::peak_search
//

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


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// time_series::edge_search
//

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
            lsq_fit_sys_ptr = new_lsq_fit (1, 4, polynomial_o4); // Note: only 3rd order used!
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

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// nodes_of_interest -- construction / lookup
//

nodes_of_interest::nodes_of_interest(const char *nm)
// Constructor
{
    all_noi_s.push_back(this);
    col_index = -1;
    name      = nm;                              // was allocated via strdup in the lexer
}

nodes_of_interest *nodes_of_interest::find(const char* name)
// Finds the NOI by name or return nullptr
{
    for (unsigned i = 0; i < all_noi_s.size(); i++)
        if (!strcmp(name, all_noi_s[i]->name))
            return all_noi_s[i];
    return nullptr;
}

nodes_of_interest *nodes_of_interest::find_undef()
{
    for (unsigned i = 0; i < all_noi_s.size(); i++)
        if (all_noi_s[i]->col_index < 0)
            return all_noi_s[i];
    return nullptr;
}

unsigned nodes_of_interest::set_col_index(unsigned ind)
{
    if ((col_index < 0) || ((unsigned)col_index == ind)) {
        col_index = (int) ind;
        return 0;
    }
    return 1;
}

int nodes_of_interest::print(char* f_name)
{
    FILE *of = fopen(f_name, "w");
    if (of == nullptr)
        return 1;
    fprintf(of, "# NOI name: %s\n", name);
    if ((col_index < 0) || (col_index >= josim_out_columns.size()) || (josim_out_columns[col_index] == nullptr))
        fprintf(of, "# has no data\n");
    else
        josim_out_columns[col_index]->print_all(of);
    fclose(of);
    
    return 0;       // success
}

void nodes_of_interest::print_all()
{
    for (unsigned i = 0; i < all_noi_s.size(); i++) {
        char buf[64];
        sprintf(buf, "NOI_%u.dat", i);
        all_noi_s[i]->print(buf);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Analysis callbacks (function2 expressions)
//

double locate_peak(void *obj_ptr, double x)
//
// The node of interest is expected to be a voltage time series and this fuction will determine the
// time of the x-th peak. <x> is expected to be an integer and peaks are numbered 0,1,2,...
//
// A NAN is returned upon failure
//
{
    nodes_of_interest *noi_ptr = (nodes_of_interest *) obj_ptr;

    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < josim_out_columns.size()));
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

double count_peaks(void *obj_ptr, double x)
//
// The node of interest is expected to be a voltage time series and this fuction counts the number
// of peaks.
//
// The argment <x> isn't used for now. May be used to supply threshold fractions later
//
{
    nodes_of_interest *noi_ptr = (nodes_of_interest *) obj_ptr;

    unsigned n_peaks = 0;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < josim_out_columns.size()));
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

double locate_rise(void *obj_ptr, double x)
//
// Similar to the peak seach, but locate the rising edge of a phase transition.
// The node of interest ought to be a phase column and a rising transition mean that the phase
// increases by about 1 turn (2PI). Rising transitions are relative, so a stair-case function
// will have multiple rising transitions. The time of the transition is defined to be the
// Inflection point, where the slope stops rising and starts falling: zero cossing of the
// second derivative.
//
{
    nodes_of_interest *noi_ptr = (nodes_of_interest *) obj_ptr;

    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < josim_out_columns.size()));
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

double locate_fall(void *obj_ptr, double x)
{
    nodes_of_interest *noi_ptr = (nodes_of_interest *) obj_ptr;
    
    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < josim_out_columns.size()));
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

double locate_edge(void *obj_ptr, double x)
// Dito, but locates any edge (rising or falling
{
    nodes_of_interest *noi_ptr = (nodes_of_interest *) obj_ptr;
    
    int n_peak = (int) nearbyint(x);
    if (n_peak < 0)
        return NAN;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < josim_out_columns.size()));
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

double count_edges(void *obj_ptr, double x)
// Counts the number of edges (regardless of direction)
//
// <x> isn't used now, but may be used later to supply edge search parameters
{
    nodes_of_interest *noi_ptr = (nodes_of_interest *) obj_ptr;
    
    unsigned n_edges = 0;
    
    int ind = noi_ptr->get_col_index();
    assert((ind >= 0) && (ind < josim_out_columns.size()));
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

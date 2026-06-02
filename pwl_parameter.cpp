//
// pwl_parameter.cpp -- implementation of the PWL waveform parameter.
//
// See parameter.h for the full description.
//

#include "parameter.h"
#include "expression.h"
#include "nodes_of_interest.h"  // for time_pattern

#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// parameter::print_to_file() -- default implementation for numeric parameters.
//
// Replaces the inline fprintf in spice_elements::print() so that pwl_parameter
// can override it without changing the spice_elements machinery.
//

void parameter::print_to_file(FILE *fp)
{
    double v = get_cur_value();
    if (!isfinite(v)) {
        fprintf(stderr, "Reference to '%s' parameter yielded NAN\n", get_name());
        fprintf(fp, "NAN");
        return;
    }
    fprintf(fp, "%.12lg", v);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// pwl_parameter constructor
//

pwl_parameter::pwl_parameter(char *nm,
                             time_pattern *pat,
                             expression   *tr,
                             expression   *tf,
                             expression   *tw,
                             expression   *vh,
                             expression   *vl)
    : parameter(nm,
                -__DBL_MAX__, __DBL_MAX__,
                1,              // no_print: exclude from scalar summaries
                0,              // not tunable
                0)              // no log mapping
    , pattern(pat)
    , t_rise(tr), t_fall(tf), t_width(tw)
    , v_high(vh), v_low(vl)
{
    assert(pattern && t_rise && t_fall && t_width && v_high && v_low);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// pwl_parameter::print_to_file()
//
// Each pulse centre time t_c from the pattern produces exactly 4 PWL points:
//
//   A: (t_c - half_w - t_r,  v_lo)   <- signal still low, rise about to start
//   B: (t_c - half_w,        v_hi)   <- rise complete, flat top starts
//   C: (t_c + half_w,        v_hi)   <- flat top ends, fall starts
//   D: (t_c + half_w + t_f,  v_lo)   <- fall complete, signal low again
//
// The complete output is:
//   "0 v_lo  A B C D  [A B C D ...]"
//
// where the initial "0 v_lo" establishes the baseline and each subsequent
// pulse is immediately preceded by its own point A (which also serves as the
// low-level dwell between pulses).
//
// Times are emitted in seconds (bare floating-point, no unit suffix).
//

void pwl_parameter::print_to_file(FILE *fp)
{
    // Evaluate control expressions at current parameter values
    double tr = t_rise ->get_value();
    double tf = t_fall ->get_value();
    double tw = t_width->get_value();
    double vh = v_high ->get_value();
    double vl = v_low  ->get_value();

    if (!isfinite(tr) || !isfinite(tf) || !isfinite(tw) ||
        !isfinite(vh) || !isfinite(vl)) {
        fprintf(stderr,
            "pwl_parameter '%s': control expression(s) evaluated to NAN -- "
            "emitting zero waveform\n", get_name());
        fprintf(fp, "0 0 1 0");
        return;
    }
    if (tr < 0.0 || tf < 0.0 || tw < 0.0) {
        fprintf(stderr,
            "pwl_parameter '%s': rise/fall/width must be >= 0\n", get_name());
        fprintf(fp, "0 0 1 0");
        return;
    }

    double half_w = tw * 0.5;

    // Evaluate all pattern times.
    // time_pattern::get_times() walks the p_time_element list, evaluating each
    // expression and resolving relative deltas, and returns absolute times.
    vector<double> times = pattern->get_times();

    if (times.empty()) {
        fprintf(stderr, "pwl_parameter '%s': pattern '%s' has no time points\n",
                get_name(), pattern->get_name());
        fprintf(fp, "0 %.12lg", vl);
        return;
    }

    // Validate that pulses don't overlap (each pulse spans
    // [t_c - half_w - tr, t_c + half_w + tf])
    for (unsigned i = 0; i + 1 < times.size(); i++) {
        double end_i   = times[i]   + half_w + tf;
        double start_i1 = times[i+1] - half_w - tr;
        if (start_i1 < end_i) {
            fprintf(stderr,
                "pwl_parameter '%s': pulse %u (centre %.4gp) overlaps "
                "pulse %u (centre %.4gp) -- waveform may be incorrect\n",
                get_name(),
                i,   times[i]   * 1e12,
                i+1, times[i+1] * 1e12);
        }
    }

    // Initial baseline
    fprintf(fp, "0 %.12lg", vl);

    for (double t_c : times) {
        double tA = t_c - half_w - tr;   // rise starts
        double tB = t_c - half_w;        // rise ends
        double tC = t_c + half_w;        // fall starts
        double tD = t_c + half_w + tf;   // fall ends

        fprintf(fp, " %.12lg %.12lg", tA, vl);  // point A: low, rise about to start
        fprintf(fp, " %.12lg %.12lg", tB, vh);  // point B: high, flat top starts
        fprintf(fp, " %.12lg %.12lg", tC, vh);  // point C: high, fall about to start
        fprintf(fp, " %.12lg %.12lg", tD, vl);  // point D: low, pulse complete
    }
}

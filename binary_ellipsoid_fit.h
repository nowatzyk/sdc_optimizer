//
// This is a sub-system to the BO functions. Specifically, it approxiamtes a multi-dimensional
// Shmoo plot with a n-dimensional hyper-ellipsoid. It tries to find the largest (most volume)
// such ellipsoid, reports iis center, and the marging for each parameter
//
#pragma once

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

extern "C" {                                    // The fit-component functions are plain C-code:
    #include "lsq_fit.h"
};

#include "nl_lsq_fit.h"
#include "bef_functions.h"
#include "bef_plan.h"

#include "eval_cache.h"
#include "loop_complex.h"
#include "parameter.h"
#include "bo_optimizer.h"
#include "xrand.h"

const uint64_t  bef_pnt_value    = 0x0001ull;   // If bit set: value = 1 (pass), fail otherwise
const uint64_t  bef_pnt_synth    = 0x0002ull;   // Synthetic point (not result of a Josim run)
const uint64_t  bef_pnt_squashed = 0x0004ull;   // This was a pass, but was turned to faill because of a
                                                // fail from a closer evaluation. This can happen because
                                                // the pass region may not be convex. Squashing points 
                                                // is a way to discared non-convex points
const uint64_t  bef_pnt_outlier  = 0x0008ull;   // This point was a pss, but is considered to be an outlier
                                                // because it is too far from the consensus ellipsoid
const uint64_t  bef_convexified  = 0x0010ull;   // Filtered out via the convexify heuristic

class explored_points {
    unsigned        n_dim;                      // # of dimensions (= #of FP numbers for coordinate)
    unsigned        n_alloc;                    // allocation size
    unsigned        n_used;                     // Number of points valid

    unsigned        item_size;                  // #of doubles per item
    double          *data;
    
    void resize(unsigned new_m);                // does the obvious, but does not need to be public
    
public: 
    explored_points (unsigned n_dim, unsigned size = 8192);
    
    ~explored_points() { std::free(data); };
    
    uint64_t        &meta(unsigned i) {
        return *reinterpret_cast<uint64_t*>(data + i * item_size);
    }
    
    double* value(unsigned i) { return data + (1 + i * item_size); };
                                                // Note: there is no range check!!!!
                                                // access value[i] only for i in [0,n_dim - 1]
   
    unsigned add_pnt(double *pnt);              // Adds one point (with all attributes ste to 0

    unsigned size() { return n_used; };
    
    void print_stat(FILE *fp);                  // Diagnostic tool: print stats
};

struct lm_solution {                            // Captures one LM LSQ fit result
    double              begin_res;              // Residual at the start
    double              end_res;                // Residual at the end
    unsigned            n_iter;                 // #of LM iterations performed
    int                 ec;                     // Error condition returned from the last step 
    unsigned            pok_fail;               // Reason for POK fail (or 0 if successfull)
    double              *param_init;            // Starting parameter set (POK passes)
    double              *param_final;           // Parameter set at completeion
    double              a_init;                 // initial sigmoid shape factor: if > 0.0, then a sans a fit
    double              a_incr;                 // change of a per iteration
};

class bin_ellipsoid_fit {
    unsigned            n_dim;                  // #of dimensions/parameters to consider
    double              *x_start;               // Interior start point
    
    explored_points     *exp_pnts;              // Points visited
                                                // needed for using nl_lsq_fit() in external data mode
    
    nl_lsq_fit          *e_fit;                 // the NL lsq fit machinery
    nl_lsq_fit          *e_fit_sa;              // dito, but without the sigmoid shape factor <a> (sans a)
    unsigned            n_e_params;             // #of parameters for the ellipsoid (a function of n_dim)
    double              *e_params;              // Parameters of the fitted ellipsoid
    vector<lm_solution> solutions;              // LM solutions
    
    EvalCache           *ev_cache;              // Perhaps saves some JoSIM calls (doubtful)
    
    unsigned            warn_high, warn_low;    // Bit-vectors to suppress range limit warnings
    
    FILE                *result_fp;             // Reporting output file
    
    vector<const_parameter *>& opt_params;      // The actual parametes to be expolred
    parameter           *of_ptr;                // The objective function: >0 is pass, fail otherwise (inc. NAN)
    FILE                *sum_fp;                // Logs all evalualtion

    //
    // Timing members for run-time estimate (set in constructor, used in run())
    //
    struct timespec     t_ctor_start_;          // Wall-clock time at start of constructor
    struct timespec     t_ctor_end_;            // Wall-clock time at end of constructor (after initial explore)
    unsigned            n_ctor_sims_;           // #of simulations run during initial exploration

    //
    // Active plan pointer -- set at the start of run(), cleared at the end.
    // nullptr outside of run().
    //
    const bef_plan      *plan_;

    unsigned            eval_pnt(double *pnt, int *ind = nullptr);  // Evaluate one point
    void                explore(double *pnt);   // Explore the passing region from <pnt>
    void                ray_search(double *pnt, double *dir);
                                                // probe points along the ray from point <pnt> in direction <dir>
    void                ray_search_mode0(double ds, double de, double *pnt, double *dir, unsigned n_probes);
    void                ray_search_mode1(double ds, double de, double *pnt, double *dir, unsigned n_probes);
                                                // recursive probe functions
    void                estimate_initial_e_params(); // Make up an estimate for the fit to start from
    void                set_up_proto_solution(double ai, double da);    // creates a solution entity
    void                derive_solution(lm_solution &sol);              // derives a solution entity
    unsigned            param_recovery(lm_solution &sol);               // Fix failing parameter set (with margin)
    unsigned            min_param_recovery(double *p_ptr, unsigned sa); // dito, but with minimal changes to parameters
    int                 solve();                // Perform the one fit iteration
    void                solve_one(lm_solution &sol); // Perform one LM LSQ fit
    void                de_duplicate();         // removes duplicates

    void                reject_outliers();      // Filter out outliers
    void                e_shell_search(unsigned n, double a = 50.0); // explore points on the ellipsoid shell
    void                hp_filter();            // Hyper-plane convexifier
    int                 build_wall(unsigned n); // Adds synthetic points to discurage ellipsoid escape
    
    unsigned            n_rays;                 // #of extra rays to cast
    unsigned            n_candidates;           // #of candidates to consider per non-axis ray
    unsigned            n_probes_p_ray;         // #of probes per ray (search budget)
    double              *ray_dirs;              // Storage area for the ray directions
    
    double              outlier_frac;           // The fraction of data points to be considered outliers
                                                // This needs to be in [0,1]. Say 0.01 for 1%
    void                print_results();        // Output the print_results

    void                print_elliosoid(char *fn, unsigned x, unsigned y);  // a debugging aid
                                                // Outputs the points of the ellipsoid intersection
                                                // with the hyperplane defined by the dimensions x and y
    void                print_e_major_axis(char *fn);   // a debugging aid
                                                // Outputs a file with the point along the major
                                                // ellipsoid axis

    void                print_all(char *bfn, double z_cut_off = 1.0);   // Yet another debugging aid:
                                                // Draws all ellipsoids in unscaled [0,1] space for all dimension pairs
    void                print_one(unsigned ix, unsigned iy, FILE *of, double z_cut_off);  // used by above
    
    void                print_a_point(FILE *of, const double *pnt, const double *off = nullptr, double d = 0.0);
    void                print_a_point_ec(FILE *of, const double *pnt, const double *off = nullptr, double d = 0.0);
                                                // same, but in ellipsoid coordinates: scalling is applied

public:
    //
    // Constructor: allocate data structures, validate the start point, and
    // run the initial axis-aligned exploration. Does NOT run the fit loop.
    // Call run() after construction (with an optional plan) to do the work.
    //
    // The old constructor parameters (n_iter, n_ray_mul, n_can, n_p_p_ray, outl_frac)
    // are now supplied via bef_plan. Default values are preserved as plan defaults
    // so existing call sites can migrate incrementally:
    //
    //   Old:  new bin_ellipsoid_fit(fp, params, of, sfp, 5, 3, 32, 16, 0.05);
    //   New:  auto *b = new bin_ellipsoid_fit(fp, params, of, sfp);
    //         b->run();                             // uses all defaults
    //   Or:   bef_plan plan(n_dim, budget);         // budget-controlled plan
    //         b->run(plan);
    //
    bin_ellipsoid_fit (FILE *result_fp,         // where to place the results
                       vector<const_parameter*> &opt_params,  // the parameters that need to optimized
                       parameter *of_ptr,       // pointer to the objective function
                       FILE *sum_fp             // Passed on to the loop complex to log the smulation outputs
                       );

    //
    // run() -- execute the fit loop.
    //
    // Uses the supplied plan for all tuning parameters. If no plan is given,
    // a default plan is used that matches the previously hardcoded values:
    //   n_iter=5, n_ray_mul=3, n_candidates=32, n_probes_per_ray=16, outlier_frac=0.05
    //
    void run(const bef_plan &plan = bef_plan());
};


double vect_scalar_product(const double *va, const double *vb, unsigned n_dim);
void vect_normalize(double *vect, unsigned n_dim);

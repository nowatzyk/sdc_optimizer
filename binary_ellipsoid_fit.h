//
// This is a sub-system to the BO functions. Specifically, it approxiamtes a multi-dimensional
// Shmoo plot with a n-dimensional hyper-ellipsoid. It tries to find the largest (most volume)
// such ellipsoid, reports iis center, and the marging for each parameter
//
#pragma once

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>

extern "C" {                                    // The fit-component functions are plain C-code:
    #include "lsq_fit.h"
};

#include "nl_lsq_fit.h"
#include "bef_functions.h"

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
};


class bin_ellipsoid_fit {
    unsigned            n_dim;                  // #of dimensions/parameters to consider
    double              *x_start;               // Interior start point
    
    explored_points     *exp_pnts;              // Points visited
                                                // needed for using nl_lsq_fit() in external data mode
    
    nl_lsq_fit          *e_fit;                 // the NL lsq fit machinery
    unsigned            n_e_params;             // #of parameters for the ellipsoid (a function of n_dim)
    double              *e_params;              // Parameters of the fitted ellipsoid
    
    EvalCache           *ev_cache;              // Perhaps saves some JoSIM calls (doubtful)
    
    unsigned            warn_high, warn_low;    // Bit-vectors to suppress range limit warnings
    
    FILE                *result_fp;             // Reporting output file
    
    vector<const_parameter *>& opt_params;      // The actual parametes to be expolred
    parameter           *of_ptr;                // The objective function: >0 is pass, fail otherwise (inc. NAN)
    FILE                *sum_fp;                // Logs all evalualtion
    
    unsigned            eval_pnt(double *pnt, int *ind = nullptr);  // Evaluate one point
    void                explore(double *pnt);   // Explore the passing region from <pnt>
    void                ray_search(double *pnt, double *dir);
                                                // probe points along the ray from point <pnt> in direction <dir>
    void                ray_search_mode0(double ds, double de, double *pnt, double *dir, unsigned n_probes);
    void                ray_search_mode1(double ds, double de, double *pnt, double *dir, unsigned n_probes);
                                                // recursive probe functions
    void                estimate_initial_e_params(); // Make up an estimate for the fit to start from
    int                 solve();                // perform one iteration of the LSQ fit 
    void                reject_outliers();      // Filter out outliers
    void                e_shell_search(unsigned n); // explore points on the ellipsoid shell
    void                hp_filter();            // Hyper-plane convexifier
    
    unsigned            n_iterations;           // #of LM-fitting steps to be performed
    unsigned            n_rays;                 // #of extra rays to cast
    unsigned            n_candidates;           // #of candidates to consider per non-axis ray
    unsigned            n_probes_p_ray;         // #of probes per ray (search budget)
    double              *ray_dirs;              // Storage area for the ray directions
    
    double              outlier_frac;           // The fraction of data points to be considered outliers
                                                // This needs to be in [0,1]. Say 0.01 for 1%
    void                print_results();        // Output the print_results

    void                print_elliosoid(char *fn, unsigned x, unsigned y);
                                                // Outputs the points of the ellipsoid intersection
                                                // with the hyperplane defined by the dimensions x and y
public:
    bin_ellipsoid_fit (FILE *result_fp,         // where to place the results
                       vector<const_parameter*> &opt_params,  // the parameters that need to optimized
                       parameter *of_ptr,       // pointer to the objective function
                       FILE *sum_fp,            // Passed on to the loop complex to log the smulation outputs
                       unsigned n_iter = 4,     // #of relocations of the ellipsoid / lsq fitting steps
                       unsigned n_ray_mul = 3,  // multiply the #of parameter by this number to get the 
                                                // number of extra (beyond 2 per axis) rays to cast
                       unsigned n_can = 32,     // #of randomly choosen candidates from which to choose
                                                // the best direction for an extra ray to cast
                       unsigned n_p_p_ray = 16, // #of points to probe for each ray
                       double outl_frac = 0.05  // Outlier fraction (default is 5%, or 1 in 20) 
                       );

};


double vect_scalar_product(const double *va, const double *vb, unsigned n_dim);
void vect_normalize(double *vect, unsigned n_dim);

//
// The non-linear leaset-square fit should be used only if the linear version
// is not sufficient. It is considerably more cumbersome, slower and prone
// to failure (non-convergence) than the linear version.
//

#ifndef _NL_LSQ_FIT_HEADER_
#define _NL_LSQ_FIT_HEADER_

#define _USE_MATH_DEFINES 
#include <math.h>
#include <float.h>
#include <assert.h>

//#define _PARALLEL_NL_LSQ_FIT_                 // When defined allows concurrent solvers
//
// The lsq_fit - functions are included because this code uses its linear equation solver
// (the differentiation functions can be used to verify the derivatives)
//
extern "C" {
#include "lsq_fit.h"                            // least square fit functions
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Configuration and control:
//

//#define _NL_LSQ_FIT_USE_KAHAN_SUMMATION_      // If defined, use Kahan summation
#define _NL_LSQ_FIT_INIT_LAM 1000.0             // Initial Levenberg damping parameter
#define _NL_LSQ_FIT_LINE_LSRM 50.0             // If defined, line search is used and this (was 200)
                                                // constant limits the delta-parameter scale
                                                // 50 is a good start
#define _NL_LSQ_FIT_LINE_RETRY_LIMIT 20         // Re-try limit (before giving up), 20 worked OK in some cases

//#define _NL_LSQ_FIT_CONV_CRIT 0.0001          // Convergence creteria:
#define _NL_LSQ_FIT_CONV_CRIT 0.001             // 1 in 1000 is good enough for the ellipsoid fit use
                                                // The solver declares victory once the incremental change <dp>
                                                // of the parameter vector is less than this value, i.e.
                                                //     |dp|/|p| < _NL_LSQ_FIT_CONV_CRIT
                                                // is true, where |x| = sqrt(sum x(i)*x(i)) is the Eucilian norm.

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// <ext_double> is a replacement type for <double> to implement Kahan summation:
//

class ext_double {                              // Extended precision double
                                                // to be used for Kahan summation
    double          s;
    double          c;

public:
    ext_double() {;};                           // Constructor does nothing
    ext_double(double x) {s = x; c = 0.0;};     // Initializing constructor

    inline double operator = (double x)         // Assignement operator
                {c = 0.0; return s = x;};

    inline double operator += (double x)        // Add a <double>
                {   volatile double y = x - c;
                    volatile double t = s + y;
                    c = t - s;
                    c -= y;
                    return s = t;
                };

    inline double operator += (ext_double x)    // Add an <ext_double>
                {   *this += x.c;
                    return *this += x.s;
                };

    inline double operator *= (double x)       // Multiplication with a <double>
                {c *= x; s *= x; return s;};

    inline operator double() const             // De-reference operator
                { return s;};
};

#ifdef _NL_LSQ_FIT_USE_KAHAN_SUMMATION_
typedef ext_double accum_ty;
#else
typedef double accum_ty;
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The non-linear fit class:
//

class nl_lsq_fit {                              // Least square fit object
    int     n_params;                           // #of of parameters
    int     n_vars;                             // #of variables
    int     n_funcs;                            // #of functions

    double  (**Fp)(const double *x, const double *p, int ip);
                                                // Function pointer array
    void    (**Dp)(double *d, const double *x, const double *p, int ip);
                                                // Pointer to the derivative computing functions

    int     (*P_ok)(double *p, int ip);         // Parameter check function

    accum_ty    *JtJ;                           // lower triangular part of the Normal matrix
    accum_ty    *JtR;                           // J^t * [yi - f(xi,beta)]
    accum_ty    res;                            // Residual accumulator (for S(p) )

    double  *parameters;                        // Parameter array

    int     ipb;                                // ip-base (see desription in constructor)

    enum {  ALLOCATED,                          // Data structures are allocated, but NOT initialized
            INITIALIZED,
            DATA_ADDED,
            SOLVED
         }  state;                              // State of this fit

    double  *data_points;                       // Internal copy of the data to be fitted (may be 0)
    int     max_dp;                             // Current allocation (in #of data points)
    int     n_dp;                               // #of data points

    //
    // Iteration control variables:
    //
    double  lam;                                // The Levenberg's damping factor
    double  res_last;                           // residual from the last iteration
    double *para_last;                          // Last parameter set
    int     n_iteration;                        // Iteration number = #of of iterations completed
    int     n_retry;                            // Retry counter

    // Common temporary storage
    // Note: this class is NOT thread safe! To make it so would require to allocate/deallocate
    //       temporary storage (the stuff below this line) dynamically.

#ifndef _PARALLEL_NL_LSQ_FIT_
    static double   *t1, *t2;                   // two arrays of size <t_size>
    static double   *tA;                        // an array of size <t_size>^2
    static int      t_size;                     // how large they are
#else
    double  *t1, *t2;                           // two arrays of size <t_size>
    double  *tA;                                // an array of size <t_size>^2
    int     t_size;                             // how large they are
#endif

    // private functions
    int     sub_solve(double lam, double *d);   // The main function of the LM algorithm
    double  res_incr(const double *d, double &a, double &b);
                                                // compute the new residual (incremental change wrt. to <res>)
    void    update();                           // replay the internal data
    void    clear();                            // Prepare for data accumulation

public:
    nl_lsq_fit (                                // Constructor
        int n_param,                            // #of parameters
        int n_var,                              // #of variables
        int n_func,                             // #of functions
        double (**F)(const double *x, const double *p, int ip),
                                                // Function pointer
                                                // Note: all functions share all parameters
        void (**D)(double *d, const double *x, const double *p, int ip),
                                                // Pointer to derivatives:
                                                // For each function there must be an array of [n_param]
                                                // function pointers that point to the derivative
                                                // of the function wrt. each parameter
        int (*P_ok)(double*, int) = nullptr,    // Parameter check function (to express constraints)
        int ip_base = 0,                        // info-parameter base
        int int_storage = 128                   // If set to > 0, internal (to this class) storage is used
                                                // (this is the default !)
    );

    ~nl_lsq_fit ();                             // Destructor

    int init    (const double *pa = nullptr);   // Initialize the system, get ready for a fit:
                                                // <p> points to an array of initial parameter value
                                                // or is nullptr, in which case the current set of parameters
                                                // will be used.

    //
    // There are two ways to use this fit-function: with or without internal storage:
    //
    // WITH internal storage:
    // ======================
    // * the <int_sorage> paramter is a initial guess at the number of data points. Can be 1, storage
    //   is added as needed, but this costs some overhead. A good guess can save multiple allocations.
    // * data-points are only added once: they are essentially copied to an array within this class.
    // * the slover iterates over this data
    // * the intended calling sequence is:
    //
    //  NEW nl_lsq_fit -> init() -> [ add_datumc() ]+ -> [ solve_1s() ]+ -> get_params()
    //
    //  [...]+ denotes calls more than once
    //  The number of solver-calls is up to the user. The solver will report via retrurn code 1
    //  that the solution has converged. The user may use fewer or more iterations and/or may
    //  use other convergence criteria based on the reported residual.
    //
    // WITHOUT internal storage:
    // =========================
    // * this class does not remember the data points, instead the solver only performs one update
    //   of the parameter set. If this solution is deemed insufficient, the calling program needs to
    //   add all datapoints again. The order does not matter, but it better be the same data, otherwise
    //   your milage varies :-( (GIGO principle applies). This mode is usefull in cases where the
    //   data is stored somewhere in the calling program, or where the calling program has a way
    //   to recreate the data-set at a cost that is less that just storing it in this class.
    // * the intended calleing sequence is:
    //
    //  NEW nl_lsq_fit -> init() -> [ [ add_datum() ]+ -> solve_1s() ]+ -> get_params()
    //
    //
    // NOTES on the calling sequences:
    // * get_params() and eval() may be called at any time and will use the current parameters
    // * not all impropper calling sequences are checked. This would cause unnecessary overhead.
    //   for example, you should not use "add_datum()" when internal storage is used. However
    //   nothing will prevent this because "add_datum()" is called from "add_datumc()" for obvious
    //   reasons.
    //
    int add_datum   (const double *x, double f, double w = 1.0); // Add a point to the fit (external storage)
    int add_datumc  (const double *x, double f, double w = 1.0); // Dito, but for use with internal sorage


    int solve_1s    (double &cur_res, double &new_res); // Compute a new solution for the parameters (1 step)

    void get_params (double *pa);               // Retrieve the new parameter
    
    void set_params (const double *pa);         // Set parameters:
                                                // This is the proverbial rope to hang yourself:
                                                // this function allows to change the parameter set
                                                // during the solver iterations. This can be used
                                                // to change direction if the solver is doing something stupid.
                                                // But recognizing this and doing somthing about it is
                                                // difficult and requires you to have an idea what is
                                                // going on and why
    
    double eval (const double *x);              // Compute the fitted value for <x>
};

////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Debugging aid:
//
void check_diffs(int np, int ip, double *x, double *pa,
         double (*F)(const double *x, const double *p, int),
         void (*D)(double *d, const double *x, const double *p, int));

#endif

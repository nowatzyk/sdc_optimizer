//
// Non-linear least square fit functions
//
// This set of functions uses the Levenberg-Marquardt algorithm to iteratively determine
// the parameters of a non-linear function so that the sum of the squared differences
// between the function and a set of data points is minimized.
//
// There are a ton of caveats on when this will work: for example the initial guess for
// the solution must be reasonably good for the iteration to converge. If the solution
// has multiple minima, and if the initial guess is too far from a minima, the iteration
// will fail to converge on a solution. In any event, these function will converge on
// the nearest minima not on a global minima. Given the iterative nature of this procedure,
// these functions are a lot slower than the linear LSQ fit functions which are far
// more robust and should be used whenever possible.
//
// For more details, look up the Wikipedia description of this algorithm.
//
// Note: these functions do not depend on any external libraries other than sharing
//       the linear equation solver from lsq_fit.c
//
// Note: Kahan summation is an option to deal with ill-condition problems or fits that
//       are very good and that involve a large number of data points.
//

// Transtrum's PhD distertaion in 2011 claims that delayed gratification (= use of
// asymmetric changes to lambda) improves convergence. He suggest a small up and a large down.
// For the OCS mirror fit-problem, that doesn't work well: being biased towards
// quadratic interpolation has problems: converges more often and can have some better
// results, but also fails badly often (tried changing lam_down to 0.25). This set 
// works OK:
static const double lam_up = 2.0;       // Change lambda by this factor for up-hill steps
                        // (= when the residual got worse)
static const double lam_down = 0.8;     // Change lambda by this factor for down-hill steps

#include <stdio.h>
#include <stdlib.h>

#include "nl_lsq_fit.h"

template<class T> void ignore( const T& ) { }

//
// Note: nl_lsq_fit was written for an application where it was called a large number of times.
//       thus temporary storage is permamently allocated to reduce overhead. This class was
//       not used in a concurrent (parallel processing, multi-threaded) fashion, thus this
//       temporary storage is allocated only once and shared. It is never used to store stuff
//       beyond the scope/duration of the function where it is referenced. To enable parallel
//       execution, this is the only part that needs to be changed.
//
#ifndef _PARALLEL_NL_LSQ_FIT_
double* nl_lsq_fit::t1 = 0;                 // Private, temporary storage
double* nl_lsq_fit::t2 = 0;                 // same size as t1!
double* nl_lsq_fit::tA = 0;                 // use of the matrix of the linear system
int     nl_lsq_fit::t_size = 0;
#endif

nl_lsq_fit::nl_lsq_fit( int n_param, int n_var, int n_func,
                        double (**F)(const double *x, const double *p, int ip),
                        void (**D)(double *d, const double *x, const double *p, int ip),
                        int (*Pok)(double *p, int ip),
                        int ip_base,
                        int int_storage)
//
// Creates an instance of a non-linear least-square fit system:
//
// <n_param> : is the number of parameters to be determined
//             the function(s) to be fitted depend on these parameters is some fashion
//             there may be more than one function, but they all share one linear array
//             of parameters
// <n_var>   : is the number of variables that the function(s) depend on. Variables are
//             simply a linear vector of double's
// <n_func>  : is the number of functions. The result from each function is added together
//             to generate the actual function to be fitted. Unlike the <lsq_fit.c> package,
//             functions are not weighted with a linear coefficient. If a coefficient is
//             needed, the value that is returned by a function should be multiplied by a
//             parameter.
//             NOTE: this facility is intended to simplify cases where the function to be
//             fitted is assembled from multiple, structurally identical components that
//             are parameterized differently. If this is not needed, just set <n_func> = 1
//             and supply only one function (at the expense of wasteing one level of indirection)
//
// <F>       : the function(s) to be fitted. Each function shall depend on the variable array
//             <x> and the parameter array <p>. Neither of which should be modified by the
//             function. <ip> is an optional parameter that is passed to the function to
//             facilitate parameter selection in case of multiple functions and/or multiple
//             fit-systems. It is the sum of <ip_base> + <function-number [0,1,...]>
//             <ip> has no other side-effects and may be ignored safely for ordinary use of
//             this fit-system.
//
// <D>       : is an array of pointers to functions that compute the derivatives of F(x,p)
//             with respect to <p>. Given that there are <n_param> parameters, there must be
//             <n_param> detrivaties returned in the <d> array, which has been allocated
//             by this class. If a function does not depend on some parameter, the corresponding
//             derivative must be set to 0. This fit-system depends on derivatives that can
//             be computed from a mathematical expression of the function <F>. If there is no
//             closed-form expression for the derivative, it is possible to compute the
//             derivative numerically using the functions supplied in <lsq_fit>. These functions
//             should also be used to verify that the computed derivatives are correct (debugging).
//
// <Pok>     : is an optional function pointer (can be nullptr) that is used to test if a parameter change
//             is valid. This function will be called at the end of each iteration with the new
//             parameters. If it returns a value != 0, then the parameters are accepted and
//             the iteration proceeds. If the function returns a 0, then the magnitude of
//             paramter changes is scaled back successively by a factor of 0.5 until either
//             an acceptable parameter change is achieved, or 10 scale back attempts failed, it which
//             case the solver returns a failure and the fit should be restarted with a better
//             initial estimate for the parameters.
//
//             Note: PoK is intended not to touch the parameteres which were declared as const-pointer,
//                   but I wanted to try some interference to see what happens...
//
// <ip_base> : This is just an integer that is added to the function number (starting from 0)
//             and passed to the function and their derivative. It serves no other purpose and
//             may be ignored. The intended use it to distinguish functions by their parameters
//             in the case that there are multiple fit-systems that share functions and/or there
//             are multiple instances of the same function that compose the function to be
//             fitted. In these cases, the function/derivativ-computation can use <ip_base> to
//             compute an index to the parameters or to locate auxilary paramters. It is also
//             passed to the P_ok() function as the 2nd argument.
//             Here is some rope you can hang yourself with. 
//
// <int_storage> : If <int_storage> is > 0, the data points for the fit will be copied to an
//             array that is allocated within this class. Thus data points (f, x[]) tuples need
//             to be added only once: subsequent iteration will use the data points from the private
//             array. This is the default mode of operation. If <int_storage> is set to 0, the fit
//             will rely on external storage: the same set of data points must be added in each iteration,
//             which is more cumbersome but can save memory when the data is easily avalaible to
//             the calling program. This can also be used to deal with cases where the data points
//             are compressed and/or computed on the fly (in a reproducible fashion).
//
{
    assert(n_param > 0 && n_func > 0 && n_var > 0);     // Mild sanity check
    assert(int_storage >= 0);
    
#ifdef _PARALLEL_NL_LSQ_FIT_
    t_size = 0;                             // nothing allocated yet
#endif
    
    n_params = n_param;                     // Copy vital statistics
    n_funcs = n_func;
    n_vars = n_var;
    ipb = ip_base;

    if (int_storage > 0) {                  // Will use internal storage
        max_dp = int_storage;               // initail allocation (guess, could be wrong)
        data_points = new double[(n_vars + 1) * max_dp];
        //
        // Note: this is awkward. The politically correct way would be to use the vector-template. But
        // that would require another level of indirection (because of the unknown #of variables <n_var>)
        // or two separate allocations. There is no realloc() in C++, so extending the array involves
        // explicit copying. This solution scarifices readablility for speed and low overhead by sticking
        // to the C-way of doing business.
        //
    } else {                                // Use external storage
        max_dp = 0;
        data_points = nullptr;
    }

    Fp = F;                                 // Remember the functions
    Dp = D;
    P_ok = Pok;

    // storage alloaction:
    JtJ = new accum_ty[(n_params * (n_params + 1)) / 2];    // Normal matrix
    JtR = new accum_ty[n_params];           // Residual vector
    parameters = new double[n_params];      // Parameter array
    para_last  = new double[n_params];

    if (t_size < n_params) {                // Allocate, shared temp. storage
        if (t_size != 0) {
            delete[] t1;
            delete[] t2;
            delete[] tA;
        }
        t_size = n_params;
        t1 = new double[t_size];
        t2 = new double[t_size];
        tA = new double[t_size * t_size];
    }

    state = ALLOCATED;                      // Ready to roll
}

nl_lsq_fit::~nl_lsq_fit()
    //
    // Destructor: free allocation
    //
{
    if (data_points != nullptr) {
        delete[] data_points;
        data_points = nullptr;              // Just being paranoid
    }
    delete[] JtJ;
    JtJ = nullptr;
    delete[] JtR;
    JtR = nullptr;
    delete[] parameters;
    delete[] para_last;
    parameters = nullptr;
    para_last = nullptr;
}

void nl_lsq_fit::clear()
    //
    // Prepare for a round of data collection (part of init)
    //
{
    //
    // 1. Reset accumulators and data counters
    //
    for (int i = 0; i < n_params; i++)
        JtR[i] = 0.0;
    for (int i = 0; i < ((n_params * (n_params + 1)) / 2); i++)
        JtJ[i] = 0.0;
    res  = 0.0;
    n_dp = 0;
    
    n_retry = 0;                            // Reset retry counter

    //
    // 2. Switch state
    //
    state = INITIALIZED;
}

int nl_lsq_fit::init(const double *pa)
    //
    // Initialize the system for a fit
    //
    // Returns != 0 if the parameters do not pass the parameter check function.
    //
{
    //
    // 1. Copy the initial parameter estimate
    //
    if (pa != nullptr) {
        for (int i = 0; i < n_params; i++)
            parameters[i] = pa[i];
    } else {
        if (state == ALLOCATED) {
            //
            // This is a really bad idea to start a fit without a parameter estimate!
            //
            return 1;                       // Don't bother trying
        }
    }

    if (P_ok) {                             // Parameter check function present:
        if (!(P_ok)(parameters, ipb)) {
            state = ALLOCATED;
            return 1;                       // Invalid parameters
        }
    }

    //
    // 2. Reset the system
    //
    clear();

    lam = _NL_LSQ_FIT_INIT_LAM;             // Initial damping factor
    n_iteration = 0;

    return 0;
}

void nl_lsq_fit::get_params(double *pa)
    //
    // Retrieve the fit result
    //
{
    for (int i = 0; i < n_params; i++)
        pa[i] = parameters[i];
}

void nl_lsq_fit::set_params(const double *pa)
    //
    // Retrieve the fit result
    //
{
    for (int i = 0; i < n_params; i++)
        parameters[i] = pa[i];
}

double nl_lsq_fit::eval(const double *x)
    //
    // Compute the fitted value for <x>
    //
{
    double value = 0;

    for (int i = 0; i < n_funcs; i++)
        value += (Fp[i])(x, parameters, i + ipb);

    return value;
}

int nl_lsq_fit::add_datumc(const double *x, double f, double w)
    //
    // Add one data point
    //
    // Data points are rejected (= function return != 0) if the
    // function or its derivative produces an un-normal FP number
    // (such as NaN, INF, etc.)
    //
{
    assert (data_points != 0);              // internal storage: copy data points

    double *tdp = data_points;
    if (n_dp >= max_dp) {                   // Need to allocate some more space:
        max_dp += 128 + max_dp / 5;         // Grow by 128 plus 20%
        data_points = new double[max_dp * (1 + n_vars)];
        for (int i = 0; i < (n_dp * (1 + n_vars)); i++)
            data_points[i] = tdp[i];        // Copy the data points stored so far
        delete[] tdp;
        tdp = data_points;
    }
    tdp += n_dp * (1 + n_vars);
    *tdp++ = f;                             // copy data: f,x0,..,xn-1
    for (int i = 0; i < n_vars; i++)
        *tdp++ = x[i];

    // Note: n_dp is incremented below!
    return add_datum(x, f, w);              // Do the actual work
}

int nl_lsq_fit::add_datum(const double *x, double y, double w)
    //
    // Add one data point:
    // <x> : function argument
    // <y> : desired function value
    // <w> : optional weight of this sample: default = 1
    //
    // The computational meat
    //
    // Data points are rejected (= function return != 0) if the
    // function or its derivative produces an un-normal FP number
    // (such as NaN, INF, etc.)
    //
{
    if (state != DATA_ADDED) {              // if not in the expected (= most common state) ...
        if (state == INITIALIZED)
            state = DATA_ADDED;
        else {
            if (state == SOLVED && (data_points == 0)) {
                clear();                    // This is a subsequent data scan with external storage
                state = DATA_ADDED;
            } else
                return 2;                   // This is a mistake.
        }
    }

    //
    // 1. Compute the function and the derivatives at <x>
    //
    double f = 0.0;
    for (int i = 0; i < n_params; i++)
        t1[i] = 0.0;
    for (int i = 0; i < n_funcs; i++) {
        double t = (Fp[i])(x, parameters, ipb + i); // Compute F(x,p)

        if (isnan(t) || isinf(t))
            return 1;                       // Skip ill-formed numbers

        f += t; 

        for (int j = 0; j < n_params; j++)
            t2[j] = 0.0;                    // zero the diff-vector:
                                            // This allows users to only assign the non-zero
                                            // parts of the diff-vector, which is useful
                                            // if there are multiple functions.

        (Dp[i])(t2, x, parameters, ipb + i);// Compute d/dp(i) F(x,p)
        for (int j = 0; j < n_params; j++) {
            t = t2[j];

            if (isnan(t) || isinf(t))
                return 1;                   // Skip ill-formed numbers

            t1[j] += t * w;
        }
    }

    //
    // 2. Update residue = Sum (y(i) - f(x(i),p))^2
    //                      i
    double r = (y - f) * w;
    res += r * r;

    //
    // 3. Update the normal matrix J^t*J
    //
    int k = 0;
    for (int i = 0; i < n_params; i++)
        for (int j = 0; j <= i; j++)
            JtJ[k++] += t1[i] * t1[j];

    //
    // 4. Update the residual vector J^t * [yi - f(xi,p)]
    //
    for (int i = 0; i < n_params; i++)
        JtR[i] += t1[i] * r;


    n_dp++;                                 // data point added

    return 0;
}

void nl_lsq_fit::update()
    //
    // This function requires the internal data mode:
    // It initializes the fit system and replay the data-points
    // The intended use is to perform one iteration of the LM
    // solver useing the cached data.
    //
{
    assert(state == SOLVED && data_points != 0 && n_dp > 0);

    int n_data_points = n_dp;               // Remember how many dp's there were

    clear();                                // Reset the system

    for (int i = 0; i < n_data_points; i++) // Add all data points
        add_datum(data_points + (1 + i * (n_vars + 1)), data_points[i * (n_vars + 1)]);
}

int nl_lsq_fit::sub_solve(double lam, double *d)
    //
    // This is the main function of the LM-algorithm: Solve
    //
    // (JtJ + lam*diag(JtJ)) * d = JtR
    //
    // Note: this function will destroy <t1> and <tA>. <t2> may be used
    //       to pass the delta <d>
    //
{
    //
    // 1. Assemble the matrix in <tA>
    //
    int k = 0;
    for (int i = 0; i < n_params; i++) {
        for (int j = 0; j <= i; j++) {
            double t = JtJ[k++];
            if (i == j) {                   // Diagonal element:
                t *= 1.0 + lam;             // Apply the Levenberg damping factor lamda <lam>
                tA[i * n_params + j] = t;
            } else {                        // Not diagonal element:
                tA[i * n_params + j] = t;   // complete the
                tA[j * n_params + i] = t;   //  matrix
            }
        }

        t1[i] = JtR[i];                     // A copy is made because the Ax = b solver clobbers <b> = <t1>
    }

    //for (int i = 0; i < n_params; i++) {
    //  for (int j = 0; j < n_params; j++)
    //      printf(" %12.4g", tA[i * n_params + j]);
    //  printf(" |  x%d  |  %12.4g\n", i, t1[i]);
    //}

    //
    // 2. Solve the system:
    //
    int ec = lin_equ (n_params, tA, d, t1);
    //for (int i = 0; i < n_params; i++)
    //  printf("  x[%d]= %14.6e\n", i, d[i]);

    return ec;
}

static int nl_key (const void *a, const void *b)
    // Sort residual accumulator by size
{
    double aa = *((accum_ty *) a), bb = *((accum_ty *) b);
    aa = fabs(aa);
    bb = fabs(bb);

    if (aa > bb) return  1;
    if (aa < bb) return -1;
    return 0;
}

double nl_lsq_fit::res_incr(const double *d, double &a, double &b)
    //
    // Computes the estimate for the new residual if <d> were to be added to the parameters
    //
    // <a>, <b> will be set to the polynomial coefficient that can be used to
    // compute the incremental residual for other points: supposed that the parameter
    // vectctor <d> were scaled by a value x (each element of <d> is multiplied by x)
    // then the incremental residue could be expressed as:
    //
    //  res_inc = a * x^2 + b * x
    //
    // There is no constant term.
    //
{
    //
    // Note: this function is notoriously ill-conditioned. The change of the residual
    //       tends to be the difference of numbers with a large range. The difference also
    //       can be small when said numbers are large. Thus this function stores the
    //       contributing terms in an array first, sorts the array and then adds them
    //       up. This looks like a Rube-Goldberg contraption, but is quite necessary for
    //       some of the problems that this was designed to deal with.
    //
    int n = 2 * n_params + (n_params * (n_params - 1)) / 2;

#ifndef _PARALLEL_NL_LSQ_FIT_
    static accum_ty *d_res = 0;             // Residual accumulator
    static int max_d_res = 0;

    if (max_d_res < n) {                    // Allocate residual accumulator (once)
        max_d_res = n;
        if (d_res != 0)
            delete[] d_res;
        d_res = new accum_ty[max_d_res];
    }
#else
    accum_ty *d_res = new accum_ty[n];
#endif

    int k = 0;
    accum_ty *dp_a = d_res;         
    accum_ty *dp_b = d_res + n;
    for (int i = 0; i < n_params; i++) {
        accum_ty t;                         // To make things work with either double and Kahan summation

        for (int j = 0; j < i; j++) {
            t = JtJ[k++];
            t *= 2.0 * d[i] * d[j];
            *dp_a++ = t;                    // quadratic terms are stored from the beginning
        }
        t = JtJ[k++];
        t *= d[i] * d[i];
        *dp_a++ =  t;                       // Diagonal elements to be added only once!

        t = JtR[i];
        t *= -2.0 * d[i];
        *(--dp_b) = t;                      // linear terms are saved from the end
    }
 
    assert(dp_a == dp_b);                   // Sanity check (comment out later)

    //
    // Compute the quadratic term:
    //
    accum_ty delta_res = 0.0;
    qsort(d_res, n - n_params, sizeof(accum_ty), nl_key);   // Sort numbers, smallest first
    delta_res = 0.0;
    for (int i = 0; i < (n - n_params); i++)
        delta_res += d_res[i];              // Add them up
    a = delta_res;

    //
    // Compute the linear term:
    //
    delta_res = 0.0;
    qsort(dp_b, n_params, sizeof(accum_ty), nl_key);    // Sort numbers, smallest first
    delta_res = 0.0;
    for (int i = 0; i < n_params; i++)
        delta_res += dp_b[i];               // Add them up
    b = delta_res;

    //
    // Compute the total, incremental residue (x = 1):
    //
    qsort(d_res, n, sizeof(accum_ty), nl_key);// Sort numbers, smallest first
    delta_res = 0.0;
    for (int i = 0; i < n; i++)
        delta_res += d_res[i];              // Add them up
        
#ifdef _PARALLEL_NL_LSQ_FIT_
    delete[] d_res;
#endif
    return delta_res;
}

int nl_lsq_fit::solve_1s(double &cur_res, double &new_res)
    //
    // Solve one iteration of the LM algorithm
    //
    // Return codes:
    //     0 = Success, needs more iterations
    //     1 = Success, convergence creteria met
    //    -1 = #of of retries exausted due to lack of convergence
    //    -2 = Failed to solve equation system (for example due to linear dependency)
    //    -3 = #of retries exhausted due to parameter constraints
    //    -4 = sequence error: this function was called out of order
    //    -5 = too few data points
    //
    // If the new residual <new_res> exceed the current residual <cur_res>, this
    // iteration resulted in a retry step.
    //
{
    if (state != DATA_ADDED) {
        if (state == SOLVED && data_points != 0 && n_dp > 0)
            update();                       // Perform a data scan in internal storage mode
        else
            return -4;
    }

    if (n_dp <= n_params)
        return -5;                          // Note: it is pointless to try to solve this system
                                            // if there are fewer data points than
                                            // there are parameters to fit.

    if (n_iteration > 0) {                  // Check for forward progress

        if (res_last < res) {               // Trouble: things got worse
            if (n_retry > _NL_LSQ_FIT_LINE_RETRY_LIMIT)
                return -1;                  // Retry counter exhauseted

            n_retry++;                      // retrying
            for (int i = 0; i < n_params; i++)
                parameters[i] = para_last[i];   // Restore last successful parameter set

            lam *= lam_up;                  // Increase lambda

            cur_res = res_last;
            new_res = res;                  // Note: the caller can figure out from this that a
                                            //       retry is happening (new is worse than current)
            state = SOLVED;
            return 0;
        }

        n_retry = 0;                        // Progress: reset retry counter
    }
    res_last = res;

    int ec = sub_solve(lam, t2);            // Solve the system
    if (ec != 0)
        return -2;

    double a, b;                            // Compute residual
    double ri = res_incr(t2, a, b);             
    ignore(ri);                             // Supress compiler warning about ri: it is used when the printout
                                            // below is uncommented. res_incr() chnages <a> and <b>, so this
                                            // is not redundant.
    
    double x = 1.0;             // Line search parameter
#ifdef _NL_LSQ_FIT_LINE_LSRM
    if (fabs(a) > 1e-20) {
        x = -b / (2.0 * a);
        if (x <= 0.0)
            x = 1.0;        // Doesn't make sense
        else if (x > _NL_LSQ_FIT_LINE_LSRM)
            x = _NL_LSQ_FIT_LINE_LSRM;      // Let's not get carried away..
    }
#endif

    double S_p = 0.0, S_dp = 0.0;           // Used to compute convergence creterion
    for (int i = 0; i < n_params; i++) {
        double t = parameters[i];
        para_last[i] = t;                   // Save current parameters
        S_p += t*t;

        t = t2[i] * x;
        parameters[i] += t;                 // Update parameters
        S_dp += t*t;
    }
    ec = ((S_dp / S_p) < (_NL_LSQ_FIT_CONV_CRIT * _NL_LSQ_FIT_CONV_CRIT));

    int sc = 0;                             // shift-cutting flag
    if (P_ok) {                             // There are constraints on the parameters
        for (int i = 0; !(P_ok)(parameters, ipb); i++) {
            if (i > _NL_LSQ_FIT_LINE_RETRY_LIMIT)
                return -3;
            x *= 0.5;       // Use shift-cutting
            sc = 1;
            for (int j = 0; j < n_params; j++)
                parameters[j] = para_last[j] + t2[j] * x;
        }
    }

    if (!sc)
        lam *= lam_down;                    // Unless things go wrong, reduce lam for next round

    cur_res = res;
    new_res = res + (a * x + b) * x;

    //printf(">> r(n)= %.4e  ri= %.4e  r(n+1)= %.4e a= %.4e  b= %.4e x= %7.2f ->%.4e e= %.4e\n", (double) res, ri, res + ri, a, b, x, new_res, ri - (a + b));

    n_iteration++;                          // Count iterations
    state = SOLVED;                         // Change state

    return ec;                              // Success
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Debugging stuff:
//

static double *DBG_x_ptr;
static int DBG_ip;
static double (*DBG_F_ptr)(const double *x, const double *p, int ip);

double DBG_wrapper(double *p)
    //
    // Just a wrapper to interface with the diff() function
    //
{
    return (*DBG_F_ptr)(DBG_x_ptr, p, DBG_ip);
}

void check_diffs(int np, int ip, double *x, double *pa,
                double (*F)(const double *x, const double *p, int),
                void (*D)(double *d, const double *x, const double *p, int))
    //
    // Debugging aid: used to check if the diff-functions are correct
    //
    // <np> : #of parameters
    // <ip> : just passed to the function
    // <x>  : The point where the test is takeing place
    // <pa> : Parameter array
    // <F>  : The function
    // <D>  : Computes a vector of the derivatives
    //
{
    assert (np > 0);

    static double *dpa = nullptr;           // Diff-array
    static int max_dpa = 0;

    if (max_dpa < np) {                     // Allocate diff-array
        if (dpa)
            delete[] dpa;
        max_dpa = np;
        dpa = new double[max_dpa];
    }

    // Set-up the wrapper function
    DBG_x_ptr = x;
    DBG_F_ptr = F;
    DBG_ip = ip;

    // Compute diff's:
    for (int i = 0; i < np; i++)
        dpa[i] = 0.0;                       // Initialize to 0
    (*D)(dpa, x, pa, ip);

    printf(">>> Check_diffs:\n                     Num. Diff             using D\n");
    for (int i = 0; i < np; i++) {
        double d_num = diff(pa, DBG_wrapper, i);
        double d_d = dpa[i];
        double t = d_d + d_num;
        int chk = 1;
        if (fabs(t) > 1.0e-12) {
            t = (d_d - d_num) / t;
            if (fabs(t) < 1.0e-5)
            chk = 0;
        }
        printf("d/dp[%2d] F(x,p) = %20.10e  %20.10e%s\n", i, d_num, d_d, (chk) ? "  check" : "");
    }
}


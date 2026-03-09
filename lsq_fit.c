/*  Customized f_fit for M-data
 *
 *  f_fit <data file name>
 *
 *  F(x1,...,xn) = SUM Ai * Fi (x1,...,xn)
 *
 *  Note: Fi() must be orthogonal (otherwise: no solutions are possible)
 *
 */

#define _ISOC99_SOURCE
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

#include "lsq_fit.h"

#define EPS  1e-60		    /* for near-zero tests */

int lin_equ (int n, double *A, double *X, double *D)
/******************************************************************\
* 								                                   *
*  Linear equation solver (Gaussian elimiation with pivot search)  *
* 								                                   *
*  solves    SUM(i=0..n)  A[j][i] * X[i] = D[j]  		           *
* 								                                   *
\******************************************************************/
{
    int i, j, k;
    double  t, p;
    long double lp;

    for (i = 0; i < n; ++i) {
		k = i;
		p = fabs (A[k*n + i]);
		for (j = i + 1; j < n; ++j) {	// find pivot
			if (fabs (A[j*n + i]) > p) {
				k = j;
				p = fabs (A[k*n + i]);
			}
		}

		if (p < EPS)					// check for singularity
			return 1;				// no solution

		if (k != i) {
			for (j = i; j < n; ++j) {
				t = A[k*n + j];
				A[k*n + j] = A[i*n + j];
				A[i*n + j] = t;
			}
			t = D[k];
			D[k] = D[i];
			D[i] = t;
		}

		lp = A[i*n + i];
		for (j = i + 1; j < n; ++j) {	// elimination
			long double q;
			q = (long double) A[j*n + i] / lp;
			for (k = i + 1; k < n; ++k)
				A[j*n + k] = (long double) A[j*n + k] - q * (long double) A[i*n + k];
			D[j] = (long double) D[j] - q * (long double) D[i];
		}
	}

    for (i = n - 1; i >= 0; --i) {	// roll back
        long double q;
		q = D[i];
		for (j = i + 1; j < n; ++j)
			q -= (long double) X[j] * (long double) A[i*n + j];
		X[i] = q / (long double) A[i*n + i];
	}

    return 0;						// successfull run
}

struct lsq_fit *new_lsq_fit (int n_var, int n_func, double (**F)(double *))
     //
     // allocates memory for a new fit system
     //
{
  struct lsq_fit *l;

  if (n_var < 1 || n_func < 1)
    return 0;

  l = (struct lsq_fit *) malloc (sizeof(struct lsq_fit));
  if (!l)
    return 0;		// oops! no more memory

  l->n_var = n_var;
  l->n_func = n_func;
  l->F = F;
  l->FP = 0;
  l->state = ALLOCATED;

  l->D = (double *) malloc (n_func * sizeof(double));
  if (l->D == 0)
    return 0;

  l->C = (double *) malloc (n_func * sizeof(double));
  if (l->C == 0)
    return 0;

  // Note: only the lower triangular part of A (incl. diagonal)
  //       is allocated
  l->A = (double *) malloc (((n_func * (n_func + 1)) / 2) * sizeof(double));
  if (l->A == 0)
    return 0;

#ifdef _LSQ_FIT_EXTENDED_PRECISION_
  l->Dc = (double *) malloc (n_func * sizeof(double));
  if (l->Dc == 0)
    return 0;

  l->Ac = (double *) malloc (((n_func * (n_func + 1)) / 2) * sizeof(double));
  if (l->Ac == 0)
    return 0;
#endif

  return l;
}

struct lsq_fit *new_lsq_fit_p (int n_var, int n_func, double (*FP)(double *, int))
     //
     // allocates memory for a new fit system
     // Same as above, but uses a parameterized function
     //
{
  struct lsq_fit *l;

  if (n_var < 1 || n_func < 1)
    return 0;

  l = (struct lsq_fit *) malloc (sizeof(struct lsq_fit));
  if (!l)
    return 0;		// oops! no more memory

  l->n_var = n_var;
  l->n_func = n_func;
  l->F = 0;
  l->FP = FP;
  l->state = ALLOCATED;
  l->sum_f2 = 0.0;

  l->D = (double *) malloc (n_func * sizeof(double));
  if (l->D == 0)
    return 0;

  l->C = (double *) malloc (n_func * sizeof(double));
  if (l->C == 0)
    return 0;

  // Note: only the lower triangular part of A (incl. diagonal)
  //       is allocated
  l->A = (double *) malloc (((n_func * (n_func + 1)) / 2) * sizeof(double *));
  if (l->A == 0)
    return 0;

#ifdef _LSQ_FIT_EXTENDED_PRECISION_
  l->Dc = (double *) malloc (n_func * sizeof(double));
  if (l->Dc == 0)
    return 0;

  l->Ac = (double *) malloc (((n_func * (n_func + 1)) / 2) * sizeof(double *));
  if (l->Ac == 0)
    return 0;
#endif

  return l;
}

int init_lsq_fit (struct lsq_fit *l)
     /*
      * Reset & initialize the fit-structure
      *
      * A non-0 return value indicates trouble!
      */
{
  int i;
  double *p;

  if (!l || l->n_var < 1 || l->n_func < 1 ||
      (l->state != ALLOCATED && l->state != INITIALIZED &&
       l->state != DATA_ADDED && l->state != SOLVED) ||
      !(l->A) || !(l->D))
    return 1;	// sanity check

  for (i = 0; i < l->n_func; i++)
    l->D[i] = 0.0;

  for (i = (l->n_func * (l->n_func + 1)) / 2, p = l->A; i--;)
      *p++ = 0.0;

#ifdef _LSQ_FIT_EXTENDED_PRECISION_
  for (i = 0; i < l->n_func; i++)
    l->Dc[i] = 0.0;

  for (i = (l->n_func * (l->n_func + 1)) / 2, p = l->Ac; i--;)
      *p++ = 0.0;

#endif

  l->sum_f2 = 0.0;
  l->c_f2 = 0.0;

  l->state = INITIALIZED;

  return 0;
}

double coeff_lsq_fit(struct lsq_fit *l, int i)
     //
     // Returns the coefficient for function <i>
     // in the solution
     //
{
      assert(0 <= i && i < l->n_func);	// Bound check
      assert(l && l->state == SOLVED);	// Is the a solution?

      return l->C[i];					// Yes, all OK!
}

int add_lsq_fit(struct lsq_fit *l, double *x, double f)
     //
     // Add a data-point to the lsq fit system
     //
     // non 0 return indicates trouble
     //
{
    double *p, *q, t;
    int i, j;
    static int x_size = 0;
    static double *X = 0;
	volatile double y, tt, sum;
#ifdef _LSQ_FIT_EXTENDED_PRECISION_
	double *pc;
#endif

    if (!l || !x ||
	(l-> state != INITIALIZED && l->state != DATA_ADDED))
      return 1;

    if (x_size < l->n_func) {
	x_size = l->n_func;
	// X = (double *) realloc (X, sizeof(double) * x_size);
	if (X != 0) free(X);
	    X = (double *) malloc(sizeof(double) * x_size); // Microsoft's realloc doesn't gork a 0-pointer!
	if (!X)
	    return 1;
    }

	//
	// Evaluate each function only once:
	//
    if (l->F) {
	for (i = l->n_func; i--;)
	    X[i] = (*(l->F[i]))(x);
    } else {
	for (i = 0; i < l->n_func; i++)
	    //
	    // Note: it is important that f(x,0) is called first!
	    //       the interface is calling functions in ascending order. This
	    //       allows the user to compute complicated functions in the first call
	    //       and then use partial result from the first call.
	    //
	    X[i] = (*(l->FP))(x, i);  		// just one function, but with a parameter
    }

    //
    // Accumulate stuff:
    //
    p = l->A;
#ifdef _LSQ_FIT_EXTENDED_PRECISION_
    pc = l->Ac;
#endif
    for (i = 0; i < l->n_func; i++) {
	t = X[i];
#ifndef _LSQ_FIT_EXTENDED_PRECISION_
	l->D[i] += f * t;

	// Note: only the lower triangular matrix is populated (incl. diagonal)
	for (j = 0, q = X; j <= i; j++)
	    *p++ += t * *q++;
#else
	sum = l->D[i];
	y = f * t - l->Dc[i];
	tt = sum + y;
	l->Dc[i] = (tt - sum) - y;
	l->D[i] = tt;

	// Note: only the lower triangular matrix is populated (incl. diagonal)
	for (j = 0, q = X; j <= i; j++) {
	    sum = *p;
	    y = t * *q++ - *pc;
	    tt = sum + y;
	    *pc++ = (tt - sum) - y;
	    *p++  = tt;
	}
#endif
    }

    {
	// This summation is problematic (allways use Kahan summation: low cost)
	sum = l->sum_f2;
	y = f * f - l->c_f2;
	tt = sum + y;
	l->c_f2 = tt - sum;
	l->c_f2 -= y;
	l->sum_f2 = tt;
    }

    l->state = DATA_ADDED;
    return 0;
}

int solve_lsq_fit (struct lsq_fit *l)
    //
    // Performs the least square fit
    //
    // non-0 return means trouble
    //
{
    int i, j, n;
    double *p;

    // The working storage is kept because usually, this solver
    // is called a lot!
    static int x_size = 0;
    static double *DD = 0;
    static double *AA = 0;

    if (!l || l-> state != DATA_ADDED || l->n_func < 1)
	return 1;
    n = l->n_func;				// saves some typing (mem-refs ?)

    if (x_size < n) {
	x_size = n;
	//DD = (double *) realloc (DD, sizeof(double) * n);
	    if (DD != 0) free(DD);
	    DD = (double *) malloc (sizeof(double) * n);
	if (!DD)
	    return 2;
	
	//AA = (double *) realloc (AA, sizeof(double) * n * n);
	if (AA != 0) free(AA);
	AA = (double *) malloc (sizeof(double) * n * n);
	if (!AA)
	    return 2;
    }

    p = l->A;
    for (i = 0; i < n; i++) {			// Copy the data (solver is destructive)
	DD[i] = l->D[i];

	for (j = 0; j <= i; j++)		// Complete the matrix
	    AA[j*n + i] = AA[i*n + j] = *p++;
    }

    if(lin_equ (n, AA, l->C, DD)) {
	l->state = ALLOCATED;
	return 3;				// we are hoosed!
    }

    l->state = SOLVED;

    return 0;
}

double eval_lsq_fit (struct lsq_fit *l, double *x)
    //
    // evaluate the lsq_fit for some argument <x>
    //
{
    double t = 0.0;
    int i;

    assert (l && l->state == SOLVED && x);	// Buggy host program

    if (l->F) {
	for (i = 0; i < l->n_func; i++)
	    t += l->C[i] * (*(l->F[i]))(x);
    } else {
	for (i = 0; i < l->n_func; i++)
	    t += l->C[i] * (*(l->FP))(x, i);
    }

    return t;
}

double eval_function (double *x, double (**F)(double *), const double *C, int n_func)
    //
    // evaluate the functions <F> for some argument <x> with the coeficients <C>
    //
    // This does the same as above, but does not require the overhead of keeping the
    // lsq_fit structure.
    //
    // Danger: there are no checks that the argument array match the functions or that
    //         the number of functions makes sense...
    //
{
    double t = 0.0;
    int i;

    for (i = 0; i < n_func; i++)
	t += C[i] * (*(F[i]))(x);
 
    return t;
}

double eval_function_p (double *x, double (*F)(double *, int i), const double *C, int n_func)
    //
    // evaluate the functions <F> for some argument <x> with the coeficients <C>
    //
    // This does the same as above, but uses one parameterized function
    //
    // Danger: there are no checks that the argument array match the functions or that
    //         the number of functions makes sense...
    //
{
    double t = 0.0;
    int i;

    for (i = 0; i < n_func; i++)
	t += C[i] * (*(F))(x, i);
 
    return t;
}

//#define SORT_OF_WORKS_BUF_HAS_NUMERICAL_STABILITY_PROBLEMS
#ifdef SORT_OF_WORKS_BUF_HAS_NUMERICAL_STABILITY_PROBLEMS

double rsquare_lsq_fit(struct lsq_fit *l)
     //
     // Returns R^2 = Sum (f(x) - y)^2
     //
     // where f(x) are the fitted values and y are the actual ones
     //
{
      double *p, R2, t;
      int i, j;

      assert(l && l->state == SOLVED);

      R2 = l->sum_f2;		// sum( f^2 )

      p =l->A;
      for (i = 0; i < l->n_func; i++) {
	    t = l->C[i];

	    for (j = 0; j < i; j++)
		  R2 += 2.0 * (l->C[j] * t * *p++);

	    R2 += t * t * *p++;	// The diagonal elements are added only once!

	    R2 -= 2.0 * l->D[i] * t;
      }

	  //
	  // Note: the above computation can go negative if this is a very
	  //       good fit: R2 is the result of subtracting two large numbers
	  //       to result a small number, which can cause FP rounding problems.
	  //       the fabs() below helps a little by preving a negative number
	  //       to be returned.
	  //
      return fabs(R2);
}
#else

static int lsq_rs_key(const void *a, const void *b)
	// Sort FP numbers by size
{
	double aa = *((double *) a), bb = *((double *) b);
	aa = fabs(aa);
	bb = fabs(bb);

	if (aa > bb) return  1;
	if (aa < bb) return -1;
	return 0;
}

double rsquare_lsq_fit(struct lsq_fit *l)
     //
     // Returns R^2 = Sum (f(x) - y)^2
     //
     // where f(x) are the fitted values and y are the actual ones
     //
{
      double *p, *rp;
      int i, j, n;
	  static double *R = 0;
	  static int n_R = 0;

	  volatile double sum, c, y, t;

      assert(l && l->state == SOLVED);

#ifdef _LSQ_FIT_EXTENDED_PRECISION_
	  n = 1 + 2 * (2 * l->n_func + ((l->n_func - 1) * l->n_func) /2);
#else
	  n = 1 + 2 * l->n_func + ((l->n_func - 1) * l->n_func) /2;
#endif

	  // n = #of elements needed to add up

	  if (n_R < n) {			// Allocate space
		  n_R = n;
		  if (R != 0)
			  free(R);
		  R = (double *) malloc(sizeof(double) * n_R);
	  }

	  rp = R;					// Point to summation array

      *rp++ = l->sum_f2;		// sum( f^2 )

      p =l->A;
      for (i = 0; i < l->n_func; i++) {
	    t = l->C[i];

	    for (j = 0; j < i; j++)
		  *rp++ = 2.0 * (l->C[j] * t * *p++);

	    *rp++ = t * t * *p++;	// The diagonal elements are added only once!

	    *rp++ = -2.0 * l->D[i] * t;
		//printf("D[%2d]= %20.8e  Dc= %20.8e\n", i, l->D[i], l->Dc[i]);
      }

#ifdef _LSQ_FIT_EXTENDED_PRECISION_
      p =l->Ac;
      for (i = 0; i < l->n_func; i++) {
	    t = l->C[i];

	    for (j = 0; j < i; j++)
		  *rp++ = -2.0 * (l->C[j] * t * *p++);

	    *rp++ = -t * t * *p++;	// The diagonal elements are added only once!

	    *rp++ = 2.0 * l->Dc[i] * t;
      }
#endif

	  assert(n == rp - R);		// Sanity check (proper allocation)

	  qsort(R, n, sizeof(double), lsq_rs_key);		// Sort numbers, smallest first

	//
	// Kahan summartion
	//
	sum = 0.0;
	c = l->c_f2;
	rp = R;
	for (i = 0; i < n; i++) {
		//printf(">> %2d %20.8e\n", i, *rp);
		y = *rp++ - c;
		t = sum + y;
		c = t - sum;
		c -= y;
		sum = t;
	}

	return fabs(sum);
}
#endif

void free_lsq_fit(struct lsq_fit *l)
     //
     // Free that data structures
     //
{
      assert(l && l->C && l->A && l->D);

      free(l->C);
      free(l->A);
      free(l->D);
#ifdef _LSQ_FIT_EXTENDED_PRECISION_
      free(l->Ac);
      free(l->Dc);
#endif

      l->C = 0;				// school of parnoid programming
      l->A = 0;
      l->D = 0;

      free(l);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Newton-Raphson stuff (see explanation in the *.h file)
//
// Note: lines comented out with //O1>> ...
//       are caused by a simple optimization that deviates from the textbook implementation
//       by computing -F instead of F to save a function call to negate the F-vector
//
// Note: Optimizations that were not implemented:
//       1. inline the vector functions (copy, scale, add)
//       2. Compact allocation of the nri-structure: just allocate it with one
//          large malloc-call ot the right size and then adjust the pointers accordingly.
//          that can save quite a bit given the fragmented nature of it.
//       3. Get rid of the factor of 2 (this doesn't do much)
//

#define NRI_MAX_RETRY 8				// NRI can over-estimate the change to the current solution,
									// in which case the next solution is worse. In case this
									// happens, the change is scaled back by 1/2 and a new attempt is
									// made. If there is still an overshoot, the delta is scaled back
									// again. This constant limits how often this may happen.
									// Once the delta is 1/2^NRI_MAX_RETRY things look pretty grim and
									// convergence is rather unlikely, so the NRI is aborted.
									// Alternates to this include:
									// * scaling by other factors than 1/2
									// * allowing a sub-optimal step in the hope it may dig the procedure out of
									//   a local minima.

struct function_set *new_func_set(int n)
	//
	// Just creates an empty function set (constructor)
	//
{
	struct function_set *fs;
	assert(n > 0);

	fs = (struct function_set *) malloc(sizeof(struct function_set));
	assert(fs);

	fs->nf_max = n;
	fs->nf = 0;
	fs->val = (double *) malloc(sizeof(double) * n);
	assert(fs->val);
	fs->gp = (double(**) (double *)) malloc(sizeof(func_ptr) * n); // Allocates an array of function pointers

	return fs;
}

int add_function(struct function_set *fs, func_ptr fp)
	//
	// Adds a function to the function set and returns its index
	//
{
	int i;

	assert(fs);
	assert(fp);

	for (i = 0; i < fs->nf; i++) {
		if (fs->gp[i] == fp)
			return i;						// This function is already in the set
	}

	if (fs->nf >= fs->nf_max) {				// Enough space left?
		assert(0);
		exit(1);							// Make sure over-runs do not go further
	}

	fs->gp[fs->nf] = fp;					// Save the new function pointer
	i = fs->nf;								// Actually redundant
	fs->nf = i + 1;

	return i;
}

void nri_eval_functions(struct function_set *fs, double *x)
	//
	// Evaluates all functions and saves the results
	//
{
	int i;
	for (i = 0; i < fs->nf; i++)
		fs->val[i] = (*(fs->gp[i]))(x);
}

double comp_single_func(struct function_1 *f1, struct function_set *fs)
	//
	// Compute one function
{
	int i;
	double t = 0.0;							// accumulator

	for (i = 0; i < f1->n; i++)
		t += f1->coeff[i] * fs->val[f1->index[i]];

	return t;
}

double comp_double_func(struct function_2 *f2, struct function_set *fs)
	//
	// Compute one function
{
	int i;
	double t = 0.0;							// accumulator

	for (i = 0; i < f2->n; i++)
		t += f2->coeff[i] * fs->val[f2->index1[i]] * fs->val[f2->index2[i]];

	return t;
}

void NRI_compute_F (struct nr_invert *nri, double *F, const double *s)
	//
	// Part of the Newton-Raphson iteration:
	// Compute the function F(s,x)
	//
	// <s> is a vector of <n_fits> elements
	// <F> is the result-vector of <n_var> elements
	// <x> is not needed by this function because
	// it relies on the evaluated generating functions
	// that are cached in the function set
	//
	//O1>> This function computes -F, not F : this saves flipping the sign later
{
	int i, j;

	for (i = 0; i < nri->n_var; i++) {
		// For each f(i) of F do:
		double t = 0.0;

		for (j = 0; j < nri->n_fits; j++) {

			// Sum over the s-dependent terms
			if (nri->Fs_comp[i][j])
				//O1>> t += s[j] * comp_single_func(nri->Fs_comp[i][j], nri->fs);
				t -= s[j] * comp_single_func(nri->Fs_comp[i][j], nri->fs);
		}

		// Sum over the s-independent, drivertive terms
		if (nri->Fd_comp[i])
			//O1>> t += comp_double_func(nri->Fd_comp[i], nri->fs);
			t -= comp_double_func(nri->Fd_comp[i], nri->fs);

		F[i] = t;							// Done with f(i)
	}
}

void NRI_compute_jacobian (struct nr_invert *nri, double *JF, const double *s)
	//
	// Part of the Newton-Raphson iteration:
	// Compute the Jacobian Matrix of F
	//
	// <JF> is an <n_var>x<n_var> matrix that is stored as one linear
	// vector. Indexing is JF[i][j] = JF[i * n_var + j], which is compatible
	// with the linear equation solver in this package.
	//
{
	int i, j, k, l;

	l = 0;									// <l> is used to enumerate the unique elements of the JK-matrix
	for (i = 0; i < nri->n_var; i++) {
		// For each f(i) do:

		for (j = i; j < nri->n_var; j++, l++) {
			// Compute d/dx(j) f(i)

			double t = 0.0;
			for (k = 0; k < nri->n_fits; k++) {
				// For each fit system do:

				// Sum over the s-dependent terms
				if (nri->JFs_comp[l][k])
					t += s[k] * comp_single_func(nri->JFs_comp[l][k], nri->fs);
			}

			// Sum over the s-independent, drivertive terms
			if (nri->JFd_comp[l])
				t += comp_double_func(nri->JFd_comp[l], nri->fs);

			JF[i * nri->n_var + j] = t;
			JF[j * nri->n_var + i] = t;		// Note: JK is symmetric, diagonals are set twice
		}
	}
}

double norm (double *x, int n)
	//
	// Just compute the Euclidian norm of the vector <x>
	//
{
	double t = 0.0;
	int i;

	for (i = 0; i < n; i++)
		t += x[i] * x[i];

	return sqrt(t);
}

void vs_mpy(double *x, const double s, int n)
	//
	// multiply each element of the vector <x> by <s>
	//
{
	int i;

	for (i = 0; i < n; i++)
		x[i] *= s;
}

void v_add(double *x, const double *y, const double *z, int n)
	//
	// add vectors <y> and <z> and assign the result to <x>
	//
{
	int i;

	for (i = 0; i < n; i++)
		x[i] = y[i] + z[i];
}

void v_copy(double *x, const double *y, int n)
	//
	// copy vector <y> to <x>
	//
{
	int i;

	for (i = 0; i < n; i++)
		x[i] = y[i];
}

int NRI_iterate (struct nr_invert *nri, double *x, const double *s, double eps, int n)
	//
	// Perform the Newton-Raphson iteration:
	//
	// <nri> is the system, which must be set-up and defined
	// <x>   is the solution vector (output), which must be initialized to a good
	//       starting point
	// <s>   is the vector of the desired value (input)
	// <n>   is the iteration limit
	//
{
	int i, i_retry;
	double r_n, r_np1 = 0.0;				// Residuals 

	static int size = 0;					// Size (= n_var) of the problem
	static double *F = 0;					// The Function value vector
	static double *x_last = 0;				// a copy to the last solution estimate for undo in case of
											//   an overshoot
	static double *h = 0;					// the delta to the current solution estimate
	static double *JF = 0;					// the Jacobian value matrix

	assert(nri);							// Just make sure that the nri data structure is set up
	assert(eps >= 0.0);						// 0 is included because that gives the ability to open-loop
											// this code: set eps=0 and just do a fixed number of iterations.

	if (!F || size < (nri->n_var)) {		// Static allocation of storage for the intermediate vars
		size = nri->n_var;
		if (F != 0) {						// Enlarging previous allocation (MS's realloc() can'r deal with 0-pointers !!)
			free(F);
			free(h);
			free(x_last);
			free(JF);
		}
		F = (double *) malloc(sizeof(double) * size);
		h = (double *) malloc(sizeof(double) * size);
		x_last = (double *) malloc(sizeof(double) * size);
		JF = (double *) malloc(sizeof(double) * size * size);
		assert(F && h && JF);
	}

	i_retry = 0;							// Retry limit counter

	for (i = 0; i < n; i++) {				// the big loop!
		if (nri->update)					// If there is an update function, call it first!
			(nri->update)(x);
		nri_eval_functions(nri->fs, x);		// Update the functions with the curent values

		NRI_compute_F (nri, F, s);			// compute F
		r_n = r_np1;
		r_np1 = norm(F, nri->n_var);		// Get new norm

		if (i > 0) {						// not on the first iteration, do
			//
			// This is the "modified" part of NR
			//
			if (r_np1 >= r_n) {				// The norm got worse:

				i_retry++;					// Make sure that this does not get stuck here
				if (i_retry > NRI_MAX_RETRY)// too many retries?
					return 1;				// Yes -> out of here

				vs_mpy(h, 0.5, nri->n_var);	// Reduce the step size
				v_add(x, x_last, h, nri->n_var);		// get a new x
				i--;						// This iteration does not count
				r_np1 = r_n;				// Undo damage from failed update
				continue;
			} else
				i_retry = 0;				// Made forward progress
		}

		NRI_compute_jacobian (nri, JF, s);	// Compute JF

		//O1>> vs_mpy(F, -1.0, nri->n_var);		// Clunky (easy to get rid off)

		// Now solve: JF * h = -F
		if (lin_equ (nri->n_var, JF, h, F))
			return 2;						// No solution

		v_copy(x_last, x, nri->n_var);		// Save last x (for the modified NR)
		v_add(x, x_last, h, nri->n_var);	// Update x

		if (norm(h, nri->n_var) <= eps)		// Are we done?
			return 0;						// Success!
	}

	return 3;								// iterations exhausted
}

void gen_sd(struct function_1 *cfa, struct function_set *fs,
			double *C, double (**G_1d)(double *), int nf, int max_func)
	//
	// Generate the s-dependent part
	//
	// <cfa>	is the coefficient/function array that will accumulate the solution.
	//          It adds the coeficient of functions that are used more than once.
	// <fs>		is the set of functions. Each distinct function (as identified by
	//			its pointer, is evaluated only once. That doesn't save anything if
	//			the functions are just polynomials, by this also simplifies expressions
	//			(common-subexpression reduction) and pays of for more interesting
	//			functions (say sin())
	// <C>		is the coefficient array. It is tested for 0, but this code may be
	//			extended to dro very small coefficients too...
	// <G_1d>	points to the first order derivatives of the functions
	// <nf>		is the number of fujctions
	// <max_func> is the size of <cfa>
	//
{
	int i, j, ip;

	for (i = 0; i < nf; i++) {
		double (*p)(double *);

		p = G_1d[i];						// <p> points to d/dx f(i)
		if (!p || C[i] == 0.0)
			continue;						// p is 0, meaning this derivative is always 0
											// or its coefficient is always 0

		ip = add_function(fs, p);			// Get index to function

		for (j = 0; j < cfa->n; j++)		// Is this a new function?
			if (ip == cfa->index[j])
				break;

		if (j < cfa->n)						// No: combine coefficients
			cfa->coeff[j] += -2.0 * C[i];
		else {								// Yes: allocate new term
			if (cfa->n >= max_func) {		// check for over-run
				assert(0);
				exit(1);					// make sure this is allways checked!
			}
			cfa->index[cfa->n] = ip;
			cfa->coeff[cfa->n] = -2.0 * C[i];
			cfa->n += 1;
		}
	}
}

void gen_id(struct function_2 *cfa, struct function_set *fs,
			double *C, double (**G_0d)(double *), double (**G_1d)(double *), int nf, int max_func)
	//
	// Generate the s-independent terms
	//
	// <cfa>	is the coefficient/function array that will accumolate the solution.
	//          It adds the coeficient of functions that are used more than once.
	// <fs>		is the set of functions. Each distinct function (as identified by
	//			its pointer, is evaluated only once. That doesn't save anything if
	//			the functions are just polynomials, by this also simplifies expressions
	//			(common-subexpression reduction) and pays of for more interesting
	//			functions (say sin())
	// <C>		is the coefficient array. It is tested for 0, but this code may be
	//			extended to dro very small coefficients too...
	// <G_0d>   is the original function set
	// <G_1d>	points to the first order derivatives of the functions
	// <nf>		is the number of fujctions
	// <max_func> is the size of <cfa>
	//
{
	int i, j, k, ip0, jp0, jp1;

	for (i = 0; i < nf; i++) {
		double (*p0)(double *);

		p0 = G_0d[i];						// <p0> points to f(i)
		if (!p0 || C[i] == 0.0)
			continue;						// p0 is 0, meaning this function is always 0
											// (not likely: it is part of the fit-functions)
											// or its coefficient is always 0

		ip0 = add_function(fs, p0);			// Get index to this function

		for (j = 0; j < nf; j++) {			// enumerate all the derivatives
			double (*p1)(double *);
			jp0 = ip0;						// Make a copy in case it needs to be swapped

			p1 = G_1d[j];					// <p1> points to d/dx f(i)
			if (!p1 || C[j] == 0.0)
				continue;					// p1 is 0, meaning this function is always 0
											// or its coefficient is always 0

			jp1 = add_function(fs, p1);		// Get index to this function
			if (jp0 > jp1) {				// Canonical representation
				k = jp0; jp0 = jp1; jp1 = k;
			}

			for (k = 0; k < cfa->n; k++)	// Is this a new function pair?
				if (jp0 == cfa->index1[k] && jp1 == cfa->index2[k])
					break;

			if (k < cfa->n)	{				// No: combine coefficients
				cfa->coeff[k] += 2.0 * C[i] * C[j];
			} else {						// Yes: allocate new term
				if (cfa->n >= max_func) {	// check for over-run
					assert(0);
					exit(1);				// Make sure this is caught (even if assert's are compiled out)
				}
				cfa->index1[cfa->n] = jp0;
				cfa->index2[cfa->n] = jp1;
				cfa->coeff[cfa->n] = 2.0 * C[i] * C[j];
				cfa->n += 1;
			}
		}
	}
}

struct function_1 *clone_f1(struct function_1 *f1p)
	//
	// Make a copy of <f1p>
	//
{
	int i;
	struct function_1 *p;

	if (f1p->n <= 0)
		return 0;					// empty set

	p = (struct function_1 *) malloc(sizeof(struct function_1));
	assert(p);
	p->n = f1p->n;
	p->coeff = (double *) malloc(sizeof(double) * f1p->n);
	p->index = (short *) malloc(sizeof(short) * f1p->n);
	assert(p->coeff && p->index);

	for (i = 0; i < f1p->n; i++) {
		p->coeff[i] = f1p->coeff[i];
		p->index[i] = f1p->index[i];
	}

	return p;
}

struct function_2 *clone_f2(struct function_2 *f2p)
	//
	// Make a copy of <f2p>
	//
{
	int i;
	struct function_2 *p;

	if (f2p->n <= 0)
		return 0;					// empty set

	p = (struct function_2 *) malloc(sizeof(struct function_2));
	assert(p);
	p->n = f2p->n;
	p->coeff = (double *) malloc(sizeof(double) * f2p->n);
	p->index1 = (short *) malloc(sizeof(short) * f2p->n);
	p->index2 = (short *) malloc(sizeof(short) * f2p->n);
	assert(p->coeff && p->index1 && p->index2);

	for (i = 0; i < f2p->n; i++) {
		p->coeff[i] = f2p->coeff[i];
		p->index1[i] = f2p->index1[i];
		p->index2[i] = f2p->index2[i];
	}

	return p;
}

struct nr_invert *NRI_create(int n_sys, struct lsq_fit **lsf_ptr, struct function_set **fsp,
							double (***G_1d)(double *), double (****G_2d)(double *), void (*uf)(double *))
	//
	// create the machinery to invert one or more least-square fits:
	//
	// A LSF computes a parameterized function f(x) so that SUM (f(xi) - si)^2 is minimied.
	// Thus f(xi) is approximately si for the <i> training points.
	//
	// This facility is meant to compute <x> for a given <s>:  f^-1(s) = x
	//
	// Lots of caveats apply...
	//
	// <n_sys> is the number of least square fit systems. This can be 1, but it also
	//         supports multiple (related, similar) systems that are meant to be solved
	//         simultaniously:
	//         let: f(i)(x) = s(i) be one LSF system, then <x> is computed so that
	//              SUM (s(i) - f(i)(x))^2 is minimized for all systems.
	// <lsf_ptr> points to one or more SOLVED least square system(s)
	// <fsp>   is optional and points to a function set pointer. Multiple nri's may
	//         share the same function set to save some space.
	// <G_1d>  is an array of <n_var> pointers to arrays of <n_func> function pointers
	//         which are initialized with the functions that produce the first partial
	//         derivatives wrt. the corrsponding variable. A zero-function pointer
	//         shall be used if the corresponding derivative is allways 0. This means
	//         do NOT point to a fuction that allways returns 0 (that would work, but
	//         would result in lots of overhead for unnecessary function calls).
	// <G_2d>  is an array of <n_var> pointers to arrays of <n_var> pointers to
	//         to arrays of <n_func> function pointers, which return the 2nd order
	//         partial derivaties wrt. to two variables. The 2D matrix (first 2 levels
	//         of indirection) shall be completely defined, even though this matrix is
	//         symmetric. But it may point to the same utlimate array. For example, for
	//         a 3x3 system there are only 6 unique partial derivatives, not 9, because
	//         of commutativity.
	//
	// Color me astounded to write code for a function-pointer array with 4 levels
	// of indirection. It does actually look sensible :-)
	//
	// a 0 is returned if something does not work out
	//
	// Note: this function is mostly overhead like memory alloaction and various sanity checks
	//       the actual math is in the gen_* functions
	//
{
	static struct function_1 f1_proto;		// Prototypes for the function descriptor sets
	static struct function_2 f2_proto;		// This is a staging area so that the actual
	static int f_proto_size = 0;			// allocation can be teight

	double (**G_0d)(double *);				// The actual functions that make up the fit system

	struct nr_invert *nri;
	struct function_set *func_set;
	int i, j, k;

	int nv;									// #of variables (n_var)
	int n_func = 0;							// max #of functions for the fit systems
	int max_func;							// max. number of functions needed in set (worst case)

	//
	// 0. Prepare for work
	//
	assert(n_sys > 0);						// Need at least one system

	nv = lsf_ptr[0]->n_var;
	n_func = lsf_ptr[0]->n_func;
	G_0d = lsf_ptr[0]->F;
	assert(G_0d && n_func > 0);
	for (i = 0; i < n_sys; i++) {			// Sanity check
		if (lsf_ptr[i]->FP != 0)
			return 0;						// Does not work with parametrized functions

		if (lsf_ptr[i]->state != SOLVED)
			return 0;						// System must be solved

		if (lsf_ptr[i]->n_var != nv || n_func != lsf_ptr[i]->n_func)
			return 0;						// Incompatible systems (different # of variables or functions)

		for (j = 0; j < n_func; j++)
			if (G_0d[j] != lsf_ptr[i]->F[j])
				return 0;					// Incompatible function set
	}

	max_func = n_func +						// the plain functions
			   nv * n_func +				// first-order derivatives
			   nv * n_func * n_func * 2;	// second order derivatives
	// Note: I'm not sure if this is right. If it is too small, the <gen_*> funtions will
	//       find out and abort.

	if (f_proto_size < max_func) {			// need to allocate space for the prototypes
		if (f_proto_size > 0) {				// Free previous allocation
			free(f1_proto.coeff);
			free(f1_proto.index);
			free(f2_proto.coeff);
			free(f2_proto.index1);
			free(f2_proto.index2);
		}
		f_proto_size = max_func;
		f1_proto.coeff = (double *) malloc(sizeof(double) * max_func);
		f1_proto.index = (short *) malloc(sizeof(short) * max_func);
		f2_proto.coeff = (double *) malloc(sizeof(double) * max_func);
		f2_proto.index1 = (short *) malloc(sizeof(short) * max_func);
		f2_proto.index2 = (short *) malloc(sizeof(short) * max_func);
	}


	//
	// 1. Allocate nri and fill in the simple stuff
	//
	nri = (struct nr_invert *) malloc(sizeof(struct nr_invert));
	assert(nri);
	nri->n_fits = n_sys;
	nri->n_var = nv;
	nri->update = uf;

	if (fsp && (*fsp)) {					// A function set is provided
		assert((*fsp)->n_var == nv);
		func_set = *fsp;
	} else {								// Need to allocate a fresh function set
		func_set = new_func_set(max_func);	// worst case allocation
		func_set->n_var = nv;
	}
	nri->fs = func_set;
	if (fsp)
		*fsp = func_set;					// Save it for re-use

	// Allocate various pointer arrays (just tedious - nothing deep)
	nri->Fs_comp  = (struct function_1 ***) malloc(sizeof(struct function_1 **) * nv);
	nri->Fd_comp  = (struct function_2 **)  malloc(sizeof(struct function_1 *)  * nv);
	nri->JFs_comp = (struct function_1 ***) malloc(sizeof(struct function_1 **) * ((nv * (nv + 1)) / 2));
	nri->JFd_comp = (struct function_2 **)  malloc(sizeof(struct function_1 *)  * ((nv * (nv + 1)) / 2));
	assert(nri->Fs_comp && nri->Fd_comp && nri->JFs_comp && nri->JFd_comp);
	for (i = 0; i < nv; i++) {
		nri->Fs_comp[i] = (struct function_1 **) malloc(sizeof(struct function_1 *) * n_sys);
		assert(nri->Fs_comp[i]);
	}
	for (i = 0; i < ((nv * (nv + 1)) / 2); i++) {
		nri->JFs_comp[i] = (struct function_1 **) malloc(sizeof(struct function_1 *) * n_sys);
		assert(nri->JFs_comp[i]);
	}

	//
	// 2. Generate the s-dependent part of F
	//
	for (i = 0; i < nv; i++) {				// For each variable generate f(i)
		for (j = 0; j < n_sys; j++) {		// For each system do
			f1_proto.n = 0;
			gen_sd(&f1_proto, func_set, lsf_ptr[j]->C, G_1d[i], n_func, max_func);
			nri->Fs_comp[i][j] = clone_f1(&f1_proto);
		}
	}

	//
	// 3. Generate the s-independent part of F
	//
	for (i = 0; i < nv; i++) {				// For each variable generate f(i)
		f2_proto.n = 0;

		for (j = 0; j < n_sys; j++)			// For each system do
			gen_id(&f2_proto, func_set, lsf_ptr[j]->C, G_0d, G_1d[i], n_func, max_func);

		nri->Fd_comp[i] = clone_f2(&f2_proto);
	}

	//
	// 4. Generate the s-dependent part of JF
	//
	k = 0;									// Triangular enumeration
	for (i = 0; i < nv; i++) {				// For each variable do a partial derivative
		for (j = i; j < nv; j++) {			// Dito for the second variable
			int l;
			for (l = 0; l < n_sys; l++) {	// For each system do
				f1_proto.n = 0;
				gen_sd(&f1_proto, func_set, lsf_ptr[l]->C, G_2d[i][j], n_func, max_func);
				nri->JFs_comp[k][l] = clone_f1(&f1_proto);
			}
			k++;
		}
	}

	//
	// 5. Generate the s-independent part of JF
	//
	k = 0;									// Triangular enumeration
	for (i = 0; i < nv; i++) {				// For each variable do a partial derivative
		for (j = i; j < nv; j++) {			// Dito for the second variable
			int l;
			f2_proto.n = 0;

			for (l = 0; l < n_sys; l++) {	// For each system do
				gen_id(&f2_proto, func_set, lsf_ptr[l]->C, G_0d, G_2d[i][j], n_func, max_func);
				gen_id(&f2_proto, func_set, lsf_ptr[l]->C, G_1d[i], G_1d[j], n_func, max_func);
			}

			nri->JFd_comp[k] = clone_f2(&f2_proto);
			k++;
		}
	}

	return nri;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Numerical differentiation functions to debug the diff-functions supplied to NRI
//

//
// diff() can be configured to use 8, 6 or 4 points:
// 
#define D1_8POINTS							// Use 8 points for diff2()
//#define D1_6POINTS						// Use 6 points for diff2()
// If neither are defined, 5 points will be used
//
// More points is computationally more expensive but can yield better results.
// Errors in the order of 1 in 1000 are common for small <h>. Simple functions
// may yield better results with fewer points. Selecting a larger <h> and more
// points may work better for more complex functions. Accuracy is fundametally
// dicy. Tune by choosing the #of points and by setting D1_H_MIN
#define D1_H_MIN 1024.0						// Something like 1 to 1e4

double diff(double *x, double (*F)(double *), int n)
	//
	// Numerical 1st order differentiation of <F> along dimension <n>
	//
	// This uses the finite difference method with 4 to 8 points
	//
{
	//
	// Weight and offset tables differ for the #of points used
	//
#ifdef D1_8POINTS
#define D1_N_POINTS 8
#define D1_SCALE 840.0
	static int idx[8] = { -4, 4, 3, -3, -2, 2, 1, -1};
	static int iw[8] = {3, -3, 32, -32, 168, -168, 672, -672};
#else
#ifdef D1_6POINTS
#define D1_N_POINTS 6
#define D1_SCALE 60.0
	static int idx[6] = {-3, 3, 2, -2, -1, 1};
	static int iw[6] = {-1, 1, -9, 9, -45, 45};
#else
#define D1_N_POINTS 4
#define D1_SCALE 12.0
	static int idx[4] = {-2, 2, 1, -1};
	static int iw[4] = {1, -1, 8, -8};
#endif
#endif

	volatile double h;						// to prevent compiler optimizations
	double xc;								// Center value (to be saved/restored)
	volatile double sum, c, t, y;			// Used for Kahan summation
	int i;

	assert(n >= 0);

	xc = x[n];								// Save x[n], we will restore it later

	//
	// 1. Determine h and make sure that it is a FP model number
	//
	h = fabs(xc);
	if (h < D1_H_MIN)
		h = D1_H_MIN;						// Deal with small x problem
	h *= sqrt(DBL_EPSILON);					// Initial guess at h
	t = xc;
	c = t - h;
	h = t - c;								// Ensure that h is a model number

	//
	// 2. Use Kahan summation to add the weighted f(x)
	//
	sum = 0.0;
	c = 0.0;
	for (i = 0; i < D1_N_POINTS; i++) {
		x[n] = xc + h * (double) idx[i];
		y = (double) (iw[i]) * (*F)(x);

		y -= c;								// Kahan summation starts here
		t = sum + y;
		c = t - sum;
		c -= y;
		sum = t;
	}

	x[n] = xc;								// Restore x

	return sum / (D1_SCALE * h);
}

//
// diff2() can be configured to use 9, 7 or 5 points:
// 
#define D2_9POINTS							// Use 9 points for diff2()
//#define D2_7POINTS						// Use 7 points for diff2()
// If neither are defined, 5 points will be used
//
// More points is computationally more expensive but can yield better results.
// Errors in the order of 1 in 1000 are common for small <h>. Simple functions
// may yeild better results with fewer points. Selecting a larger <h> and more
// points may work better for more complex functions. Accuracy is fundametally
// dicy. Tune by choosing the #of points and by setting D2_H_MIN
//
#define D2_H_MIN 256.0						// Something like 1 to 1e4

double diff2(double *x, double (*F)(double *), int n)
	//
	// Numerical 2nd order differentiation of <F> along dimension <n>
	//
	// This uses the finite difference method with 9 points
	//
{
	//
	// Weight and offset tables differ for the #of points used
	//
#ifdef D2_9POINTS
#define D2_N_POINTS 9
#define D2_SCALE 5040.0
	static int idx[9] = { -4, 4, -3, 3, -2, 2, -1, 1, 0};
	static int iw[9] = {-9, -9, 128, 128, -1008, -1008, 8064, 8064, -14350};
#else
#ifdef D2_7POINTS
#define D2_N_POINTS 7
#define D2_SCALE 180.0
	static int idx[7] = {-3, 3, -2, 2, -1, 1, 0};
	static int iw[7] = {2, 2, -27, -27, 270, 270, -490};
#else
#define D2_N_POINTS 5
#define D2_SCALE 12.0
	static int idx[5] = {-2, 2, -1, 1, 0};
	static int iw[5] = {-1, -1, 16, 16, -30};
#endif
#endif

	volatile double h;						// to prevent compiler optimizations
	double xc;								// Center value (to be saved/restored)
	volatile double sum, c, t, y;			// Used for Kahan summation
	int i;

	assert(n >= 0);

	xc = x[n];								// Save x[n], we will restore it later

	//
	// 1. Determine h and make sure that it is a FP model number
	//
	h = fabs(xc);
	if (h < D2_H_MIN)
		h = D2_H_MIN;						// Deal with small x problem
	h *= sqrt(DBL_EPSILON);					// Initial guess at h
	t = xc;
	c = t - h;
	h = t - c;								// Ensure that h is a model number

	//
	// 2. Use Kahan summation to add the weighted f(x)
	//
	sum = 0.0;
	c = 0.0;
	for (i = 0; i < D2_N_POINTS; i++) {
		x[n] = xc + h * (double) idx[i];
		y = (double) (iw[i]) * (*F)(x);

		y -= c;								// Kahan summation starts here
		t = sum + y;
		c = t - sum;
		c -= y;
		sum = t;
	}

	x[n] = xc;								// Restore x

	return sum / (D2_SCALE * h * h);
}

//
// diff_2d() can be configured to use 16, 36 or 64 points:
// 
#define D11_64POINTS							// Use 64 points for diff_2d()
//#define D11_36POINTS							// Use 36 points for diff_2d()
// If neither are defined, 16 points will be used
//
// More points is computationally more expensive but can yield better results.
// Errors in the order of 1 in 1000 are common for small <h>. Simple functions
// may yeild better results with fewer points. Selecting a larger <h> and more
// points may work better for more complex functions. Accuracy is fundametally
// dicy. Tune by choosing the #of points and by setting D11_H_MIN
#define D11_H_MIN 256.0							// Something like 1 to 1e4

double diff_2d(double *x, double (*F)(double *), int n1, int n2)
	//
	// Computes 2nd order derivatives: d/dx[n1] ( d/dx[n2] F(x) )
	//
{
	//
	// Weight and offset tables differ for the #of points used
	//
#ifdef D11_64POINTS
#define D11_N_POINTS 64
#define D11_SCALE 705600.0
	static int idx[64] = {
           -4,   -4,    4,    4,    3,   -4,    4,    3,
           -4,    4,   -3,   -3,    4,    4,   -2,   -2,
            2,   -4,   -4,    2,    3,   -3,   -3,    3,
           -1,   -1,   -4,    1,    4,    4,    1,   -4,
           -2,   -3,    2,   -2,    3,    3,   -3,    2,
            1,   -1,    3,    3,   -1,    1,   -3,   -3,
           -2,   -2,    2,    2,   -1,    1,    1,   -1,
		    2,    2,   -2,   -2,    1,   -1,   -1,    1  };
	static int idy[64] = {
           -4,    4,    4,   -4,   -4,   -3,   -3,    4,
            3,    3,    4,   -4,    2,   -2,   -4,    4,
            4,    2,   -2,   -4,    3,    3,   -3,   -3,
            4,   -4,    1,    4,   -1,    1,   -4,   -1,
            3,   -2,   -3,   -3,   -2,    2,    2,    3,
            3,    3,    1,   -1,   -3,   -3,   -1,    1,
           -2,    2,    2,   -2,    2,    2,   -2,   -2,
		   -1,    1,    1,   -1,    1,    1,   -1,   -1 };
	static int iw[64] = {
               9,      -9,       9,      -9,      96,     -96,      96,     -96,
              96,     -96,      96,     -96,     504,    -504,     504,    -504,
             504,    -504,     504,    -504,    1024,   -1024,    1024,   -1024,
            2016,   -2016,    2016,   -2016,    2016,   -2016,    2016,   -2016,
            5376,   -5376,    5376,   -5376,    5376,   -5376,    5376,   -5376,
           21504,  -21504,   21504,  -21504,   21504,  -21504,   21504,  -21504,
           28224,  -28224,   28224,  -28224,  112896, -112896,  112896, -112896,
		   112896,-112896,  112896, -112896,  451584, -451584,  451584, -451584  };
#else
#ifdef D11_36POINTS
#define D11_N_POINTS 36
#define D11_SCALE 3600.0
	static int idx[36] = {
	       3,   -3,   -3,    3,   -2,   -2,
          -3,   -3,    3,    3,    2,    2,
           3,    3,    1,    1,   -1,   -1,
          -3,   -3,   -2,   -2,    2,    2,
          -2,   -1,   -1,    1,    2,   -2,
		   1,    2,   -1,   -1,    1,    1  };
	static int idy[36] = {
           3,    3,   -3,   -3,    3,   -3,
           2,   -2,   -2,    2,   -3,    3,
           1,   -1,    3,   -3,   -3,    3,
          -1,    1,   -2,    2,    2,   -2,
           1,   -2,    2,    2,   -1,   -1,
		  -2,    1,   -1,    1,    1,   -1  };
	static int iw[36] = {
           1,   -1,    1,   -1,    9,   -9,
           9,   -9,    9,   -9,    9,   -9,
          45,  -45,   45,  -45,   45,  -45,
          45,  -45,   81,  -81,   81,  -81,
         405, -405,  405, -405,  405, -405,
		 405, -405, 2025,-2025, 2025,-2025  };
#else
#define D11_N_POINTS 16
#define D11_SCALE 144.0
	static int idx[16] = {-2,  2,  2, -2,    2,  1, -2, -1,    1,  2, -1, -2,    -1,   1,   1,  -1};
	static int idy[16] = {-2, -2,  2,  2,   -1,  2,  1, -2,   -2,  1,  2, -1,    -1,  -1,   1,   1};
	static int  iw[16] = { 1, -1,  1, -1,    8, -8,  8, -8,    8, -8,  8, -8,    64, -64,  64, -64};
#endif
#endif

	volatile double h;						// to prevent compiler optimizations
	double xc, yc;							// Center values (to be saved/restored)
	volatile double sum, c, t, y;			// Used for Kahan summation
	int i;

	assert(n1 >= 0 && n2 >= 0);

	if (n1 == n2)							// Special case if both dimensions are same
		return diff2(x, F, n1);

	xc = x[n1];								// save center values
	yc = x[n2];

	//
	// 1. Determine h and make sure that it is a FP model number
	//
	h = fabs(xc);
	if (h < D11_H_MIN)
		h = D11_H_MIN;						// Deal with small x problem
	h *= sqrt(DBL_EPSILON);					// Initial guess at h
	t = xc;
	c = t - h;
	h = t - c;								// Ensure that h is a model number

	//
	// 2. Use Kahan summation to add the weighted f(x,y)
	//
	sum = 0.0;
	c = 0.0;
	for (i = 0; i < D11_N_POINTS; i++) {
		x[n1] = xc + h * (double) idx[i];
		x[n2] = yc + h * (double) idy[i];
		y = (double) (iw[i]) * (*F)(x);

		y -= c;								// Kahan summation starts here
		t = sum + y;
		c = t - sum;
		c -= y;
		sum = t;
	}

	x[n1] = xc;								// Restore x/y
	x[n2] = yc;

	return sum / (D11_SCALE * h * h);
}

void verify_1st_order_diffs (int n_var, int n_func, double *x, double (**F)(double *), double (***G_1d)(double *), void	(*update_f)(double *))
	//
	// Verify the first order derivatives by
	// comparing them to results from numerical differentiation:
	//
	// <n_var> : #number of variables
	// <n_func> : #of functions
	// <x>      : the point where the verification is supposed to be done
	// <F>      : the functions
	// <G_1d>   : the first order derivatives
	// <update_f>: optional update function for the derivatives
	//
{
	int i, j;

	assert(n_var > 0 && n_func > 0);

	if (update_f)
		(*update_f)(x);				// Call update function, if there is one

	printf("Verify first order derivatives at:\n");
	for (i = 0; i < n_var; i++)
		printf("    x[%d] = %20.5g\n", i, x[i]);

	for (i = 0; i < n_func; i++) {
		double num_diff, comp_diff, err;

		for (j = 0; j < n_var; j++) {
			num_diff = diff(x, F[i], j);
			if (G_1d[j][i])
				comp_diff = (*G_1d[j][i])(x);
			else
				comp_diff = 0.0;

			err = comp_diff - num_diff;
			printf("  d/dx[%d] F[%2d] = n: %20.10g c: %20.10g err= %.3e\n", j, i,
				num_diff, comp_diff, err);
		}
	}
}

void verify_2nd_order_diffs (int n_var, int n_func, double *x, double (**F)(double *), double (****G_2d)(double *), void	(*update_f)(double *))
	//
	// Verify the first order derivatives by
	// comparing them to results from numerical differentiation:
	//
	// <n_var> : #number of variables
	// <n_func> : #of functions
	// <x>      : the point where the verification is supposed to be done
	// <F>      : the functions
	// <G_2d>   : the second order derivatives
	// <update_f>: optional update function for the derivatives
	//
{
	int i, j, k;

	assert(n_var > 0 && n_func > 0);

	if (update_f)
		(*update_f)(x);				// Call update function, if there is one

	printf("Verify second order derivatives at:\n");
	for (i = 0; i < n_var; i++)
		printf("    x[%d] = %20.5g\n", i, x[i]);

	for (i = 0; i < n_func; i++) {
		double num_diff, comp_diff, err;

		for (j = 0; j < n_var; j++) {
			for (k = 0; k < n_var; k++) {
				num_diff = diff_2d(x, F[i], j, k);
				if (G_2d[j][k][i])
					comp_diff = (*G_2d[j][k][i])(x);
				else
					comp_diff = 0.0;

				err = comp_diff - num_diff;
				printf("  d/dx[%d]dx[%d] F[%2d] = n: %20.10g c: %20.10g err= %.3e\n", j, k, i,
					num_diff, comp_diff, err);
			}
		}
	}
}

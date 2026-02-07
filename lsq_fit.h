/*
 * Somewhat generalized least square fit routines
 */

#ifndef	_LSQ_FIT_HEADER_
#define _LSQ_FIT_HEADER_

//
// Some functions are problematic to fit when there are large number of
// points and the fit is very good and there are large values involved.
// In this case the loss of precision from the accumulation causes
// problems. In particular, computing the residual (=error^2 sum) becomes
// very unreliable. By defining the extended precision flag, Kahan
// sumation is used in the accumulation phase. This incurres extra
// storage overhead and slows the data accumulation phase down. So
// If you don't care about the residual or the function to fit is
// well behaved or has a large fit error, then comment out the next line:
#define _LSQ_FIT_EXTENDED_PRECISION_

struct lsq_fit {
  int		n_var;					// #of variables
  int		n_func;					// #of functions

  double	(**F)(double *);		// Function pointer array
  double    (*FP)(double *, int);	// Pointer to parameterized function

  double	*A;						// accumulates data
  double	*D;						// dito
#ifdef _LSQ_FIT_EXTENDED_PRECISION_
  double	*Ac;					// Carries for Kahan summation
  double	*Dc;
#endif

  double    *C;						// Coefficients (solution)

  double    sum_f2;					// = sum (f*f)
  double    c_f2;					// carry (for summation)

  enum	    	{ALLOCATED=1,
		 INITIALIZED=2,
		 DATA_ADDED=3,
		 SOLVED=4} state;			// state of this fit
};

//
// Note: this is really plain vanila C-code
//       it could/should be a class, but then why bother - it works!
//
struct		lsq_fit	*new_lsq_fit (int n_var, int n_func, double (**F)(double *));
									// generates a new fit

struct		lsq_fit  *new_lsq_fit_p (int n_var, int n_para, double (*F)(double *, int));
									// dito, but with a parameterized function:
									// Functions are called in acending parameter order

int		init_lsq_fit (struct lsq_fit *lsq);
									// Initialize the fit

int		add_lsq_fit  (struct lsq_fit *lsq, double *x, double f);
									// Add a data-point

int		solve_lsq_fit(struct lsq_fit *lsq);
									// solve the Fit

double	eval_lsq_fit (struct lsq_fit *lsq, double *x);
									// Generate a fitted value

double	eval_function (double *x, double (**F)(double *), const double *C, int n_func);
									// Dito, but use a function array with coeficients

double  coeff_lsq_fit(struct lsq_fit *lsq, int i);
									// Get a coeficient from the solution

double  rsquare_lsq_fit(struct lsq_fit *lsq);
									// return R^2 (= sum (f(xi) - yi)^2) ) which measures the fit quality

void	free_lsq_fit(struct lsq_fit *lsq);
									// free the lsq_fit structure (destructor)


//////////////////////////////////////////////////////////////////////
//
// This part is used to invert one/multiple fits via Newton-Raphson
//
// Given a lsq-fit  G(x) = Sum(i) ci*gi(x)
//
// where: x     is a vector [x(0), ... ,x(n_var-1)]
//        gi(x) is a set of <n_func> functions of x
//        ci    are <n_func> constants that were computed via Least Square fit
//              such that Sum(i) (Si - G(xi))^2 is minimized for the set of
//              training tuples (Si, xi).
//
// Thus G(x) is a function that approximates S:   Si ~= G(xi)
//
// The machinery below tries to invert the function G, that is to find
// a vector <x> so that (G(x) - S)^2 is minimized for a given <S>
//
// Approach: start with a guess for xi and use NR iteration:
//
// F(x) is a vector [f(0)(x), ... ,f(n_var-1)(x)] with
//
//                                d
//                   f(i)(x) = ------- (S - G(x))^2
//                              d x(i)
//
// So the problem is to find <x> so that F(x) = [0, ... ,0]
//
// This is done by computing the vector <h> = [h(0), ... ,h(n_var-1)]
// By solving the linear equation system:
//
//                  JF(x)h = -F(x)
//
// Where JF is the Jacobian <n_var>x<n_var> matrix:
//
//                                   d
//                  jF(i)(j)(x) = ------- f(i)(x)
//                                 d x(j)
//
// Then one NR iteration does:
//
//                       x(n+1) = x(n) + h
//
// Subject to verify progress:  |F(x(n+1))| < |F(x(n))|
// If this test fails, reduce <h> linearily
//
// This code is set up that that F can be a sum of multiple lsq-fits,
// which arises if S is a vector that was fitted with one G(x) for
// each component where G(x) only differes in the constants <ci>, that
// is it uses the same generating functions gi(x).
//

typedef double (*func_ptr) (double *);

struct function_set {
	//
	// Function set: these functions are used to compute
	// the elements of the NR iteration
	//
	double		*val;					// Value vector of the functions
	double		(**gp)(double *);		// Function pointer array
	int			nf;						// #of defined functions
	int			nf_max;					// Allocation size
	int			n_var;					// #of variables (in argument)
};

struct function_1 {
	//
	// One function
	//
	int			n;						// #of components
	double		*coeff;					// Coefficient vector
	short		*index;					// Index vector
};

struct function_2 {
	//
	// Two functions
	//
	int			n;						// #of components
	double		*coeff;					// Coefficient vector
	short		*index1;				// Index vector 1
	short		*index2;				// Index vector 2
};

struct nr_invert {
	//
	// Data structure for inverting fitted functions
	//
	int			n_var;					// #of variables
	int			n_fits;					// #of fit systems = length of S (see above)
	void		(*update)(double *);	// optional update function to be called
										// before the derivatives are evaluated
	struct function_set	*fs;			// Pointer to the (shared) function set

	struct function_1	***Fs_comp;		// s-dependent component of F
	struct function_2	**Fd_comp;		// s-independent, partial derivative component of F

	struct function_1	***JFs_comp;	// s-dependent component of JF
	struct function_2	**JFd_comp;		// s-independent, partial derivative component of JF
};

//
// The function calls:
//

struct nr_invert *NRI_create	(		////// Constructor for the inversion function
		int n_sys,						// #of fit systems to be used
		struct lsq_fit **lsf_ptr,		// pointer to an array that points to the fit systems
		struct function_set **fsp,		// optional pointer to a function_set pointer:
										//    0		: a new function set will be allocated and used
										//   &(*fs) where *fs = 0	: fs will be set to a newly
										//							  allocated function set
										//   &(*fs) where *fs != 0  : the existing function set will be used
										//
										// Sharing function set can save space when a lot of fit-systems
										// need to be inverted
		double (***G_1d)(double *),		// Pointer to an array of pointer to arrays of the partial derivatives
		double (****G_2d)(double *),	// dito, but one more level of indirection to deal with the 2d
										// array of 2nd order partial derivatives wrt. to 2 variables
		void	(*update_f)(double *)	// either 0 or a pointer to an update function that is called
										// before the derivative functions (<G1_d> and <G2_d>) are called:
										// These functions (buf NOT the fit functions!) may share sub-expressions
										// that are time-consuming to compute. The purpose of the update function
										// is to evaluate these common sub-expressions and save them in some
										// shared memory for use by <G1_d> and <G2_d>. Thus the intended
										// side-effect of this function is to prepare values to be consumed by
										// the derivative functions. Proper encapsulation in a separate module
										// is advised for this hack to be used safely.
								);		// Returns a pointer to the system to invert

int NRI_iterate					(		////// Newton-Raphson iteration
		struct nr_invert *nri,			// pointer to the system to invert
		double *x,						// on call  : this vector must be initialized to the starting point
										// on return: this vector contains the solution
		const double *s,				// vector of the values the function ought to produce
		double eps,						// Termination condition: if the change to <x> in one iteration
										// falls below this value, the iteration will be stopped
										// Can be 0, in this case a fixed number of iterations will be done
		int n							// Max-number of iterations
								);		// 0 upon EPS-termination
										// 1 iteration failed to reduce residual after 3 tries with decreasing step size
										// 2 failure to solve the linear equation system that computes the next step
										// 3 Iteration limit exhausted without the step size falling below EPS

void verify_1st_order_diffs (int n_var, int n_func, double *x, double (**F)(double *), double (***G_1d)(double *), void	(*update_f)(double *));
void verify_2nd_order_diffs (int n_var, int n_func, double *x, double (**F)(double *), double (****G_2d)(double *), void (*update_f)(double *));
										// A debugging tool

//
// Export the linear equation solver interface
//
int lin_equ (int n, double *A, double *x, double *b);
double diff(double *x, double (*F)(double *), int n);

#endif

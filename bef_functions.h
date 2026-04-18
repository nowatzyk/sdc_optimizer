//
// non-linear fit functions for use in the binary ellipsoid fit sub-system
//

#ifndef _BEF_FUNCTIONS_DEF_
#define _BEF_FUNCTIONS_DEF_

// nl_lsq_fit expects array of function pointers:
extern double (*bef_function_ptr[1])(const double *x, const double *pa, int ip);
extern void   (*bef_diff_function_ptr[1])(double *d, const double *x, const double *pa, int ip);

int bef_param_ok (double *param, int ip);       // parameter acceptance function

//
// Misc.
//
void nl_check_nde(unsigned n_dims);             // Debug/diagnostics

int ellipsoid_intersect (
    const double *param,                        // Defines the ellipsoid (see above)
    unsigned n_dim,                             // #of dimensions
    const double *point,                        // Operating point *inside* the ellipsoid
    const double *dir,                          // Direction of interest: must be normalized to 1
    double &dist_p,                             // result: distance to the ellipsoid shell (> 0)
    double *dist_n = nullptr                    // Optional: distance to the shell in the negative direction
);

double vect_scalar_product(const double *va, const double *vb, unsigned n_dim);
void vect_normalize(double *vect, unsigned n_dim);

#endif

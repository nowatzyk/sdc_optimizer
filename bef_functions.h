//
// non-linear fit functions for use in the binary ellipsoid fit sub-system
//

#ifndef _BEF_FUNCTIONS_DEF_
#define _BEF_FUNCTIONS_DEF_

#include <eigen3/Eigen/Dense>

extern const double max_a;                      // sigmoid shape factor upper limit
extern const double min_a;                      //  .. and lower limit
                                                // purpose: dealing with fit parameter limits
extern const double min_foci_d;                 // minimal foci distance
extern const double gcs_min_scale;              // The gc scalefactor limits
extern const double gcs_max_scale;              //  ..

// bef_param_ok failure (not OK) explanation flags
const unsigned  bef_PnOK_min_a    = 0x0001;     // <a> too small
const unsigned  bef_PnOK_max_a    = 0x0002;     // <a> too large
const unsigned  bef_PnOK_f0_esc   = 0x0004;     // focal point 0 out of bound
const unsigned  bef_PnOK_f1_esc   = 0x0008;     // focal point 1 out of bound
const unsigned  bef_PnOK_cntr_esc = 0x0010;     // ellipsoid center out of bounds
const unsigned  bef_PnOK_f_merge  = 0x0020;     // foci too close
const unsigned  bef_PnOK_scale    = 0x0040;     // Scale factor out of bounds
const unsigned  bef_PnOK_fs2small = 0x0080;     // The focal sum is less than min_fs_excess larger
                                                // than the foci distance

// nl_lsq_fit expects array of function pointers:
extern double (*bef_function_ptr[1])(const double *x, const double *pa, int ip);
extern void   (*bef_diff_function_ptr[1])(double *d, const double *x, const double *pa, int ip);

int bef_param_ok (double *param, int ip);       // parameter acceptance function
int bef_param_ok_sa (double *param, int ip);    // Dito, but doesn't check <a>
unsigned get_bef_param_nOK_reason();            // Query why the last bef_param_ok() call failed

void set_private_a(double a);                   // Sets the sigmoid shape factor for the "*_sa" function set

// same as above, but without sigmoid shape factor
extern double (*bef_function_ptr_sa[1])(const double *x, const double *pa, int ip);
extern void   (*bef_diff_function_ptr_sa[1])(double *d, const double *x, const double *pa, int ip);

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

unsigned is_inside_ellipsoid (
    const double *param,                        // Defines the ellipsoid (see above)
    unsigned n_dim,                             // #of dimensions
    const double *point                         // n-dimensional point to be tested
);

double vect_scalar_product(const double *va, const double *vb, unsigned n_dim);
void vect_normalize(double *vect, unsigned n_dim);
void generate_random_dir(double *dir, unsigned n_dim);
double clip_to_unity(const double *pnt, const double *dir, unsigned n_dim);
Eigen::MatrixXd completeOrthogonalBasis(const Eigen::VectorXd &dir);

#endif

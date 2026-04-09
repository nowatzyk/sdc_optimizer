//
// The non-linear ellipsoid fit functions and their drivatives
//

#include "binary_ellipsoid_fit.h"


//
// Contraints that are ment to prevent the ellipsoid from degenerating:
//
const double min_focal_sum = 0.01;      // minimal focal_sum of the ellipsoid
// Note: focal_sum is defined as the sum of the distances
// to the two foci. This is twice of the radius of a sphere,
// if you think of a sphere as a degenerated ellipsoid
//  with 0 excentricity.

const double min_a = 0.001;             // The shape factor of the sigmoid function shall not 
// go negative or near 0: the gradient would vanish and LM is lost 

const double min_foci_d = 0.01;         // If the foci merge, their parameters loose independence and
// mayhem would ensure

const double GC_EPS = 1.0e-12;          // in gardiate correction this is a guard against a singularity
// that can happen when the point falls onto a focal point

const double gcs_min_scale = 0.1;       // The gc scale factor lower limit
const double gcs_max_scale = 10.0;      // The gc scale factor upper limit

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  A n-dimensional ellipsoid is defined by its two focal points and its focal-sum. However that makes
//  it rotational symmetric about its major axis. In the case of more than 2 dimensions, this
//  would mean that the scaling of the parameter is not independed across dimensions. However, each
//  dimension corresponds to one parameter of the design and there is no reason for all parameters
//  to have the same sensitivity or range. Thus for each dimension beyond 2, a independent scaling
//  parameter is addded. Thus the ellipsoid is defined as follows:
//
//  param[0] :          focal_sum the sum of the distances from a point on the 
//                      the ellipsoid shell to the two focal points. Must be > 0
//  param[1] : a        controls the slope of the sigmoid function. Must be > 0
//                      A large <a> makes the slope steeper (sharper transition between 0 and 1)
//
//  param[2, 1+n_dimensions]                 coordinates of first focal point
//  param[2+n_dimensions, 1+2*n_dimensions]  coordinates of second focal point
//  param[2+2*n_dimensions + i]              scale facor for each coordinate >2
//
// The ellipsoid mapps each point in the n-dimensional space to a probability for that point to be inside
// the ellipsoid (1) or outside (0). This function is the sigmoid function of an approximation to
// the distance of the point to the surface of the ellipsoid. The esact distance would findng the
// solution to a 4th order polynomial. While possible, that is way too cumbersome, so this distance
// is approximated by the difference between the focal-sum and the sum of the distances between
// the point and each focal point. Beacsue this approximarion is rather poor, it is corrected by
// a gradiant function that depends on the angle focal_point0 - point - focal_point1.
//

//
// Internal helper: apply the scaling transform to <pnt>.
// Returns the transformed point in <pnt_t> (caller-allocated, n_dimensions entries).
// C[] is the ellipsoid centre = (f0+f1)/2.
//
static void apply_scale(const double *pnt, const double *f0, const double *f1,
                        const double *s,   double *pnt_t, int n_dimensions)
{
    for (unsigned i = 0; i < n_dimensions; i++) {
        if (i < 2) {
            pnt_t[i] = pnt[i];
        } else {
            double Ci = (f0[i] + f1[i]) * 0.5;
            pnt_t[i] = (pnt[i] - Ci) * s[i-2] + Ci;
        }
    }
}

double n_dim_ellipsoid_gcs(const double *pnt, const double *param, int n_dimensions)
{
    double  fs = param[0];
    double   a = param[1];
    const double *f0 = param + 2;
    const double *f1 = f0 + n_dimensions;
    const double  *s = f1 + n_dimensions;   // n_dimensions-2 scale factors
    
    // Apply scaling transform
    double pnt_t[n_dimensions];
    apply_scale(pnt, f0, f1, s, pnt_t, n_dimensions);
    
    double d2f0 = 0.0, d2f1 = 0.0, t;
    for (unsigned i = 0; i < n_dimensions; i++) {
        t = pnt_t[i] - f0[i];  d2f0 += t * t;
        t = pnt_t[i] - f1[i];  d2f1 += t * t;
    }
    d2f0 = sqrt(d2f0);
    d2f1 = sqrt(d2f1);
    
    if ((d2f0 < GC_EPS) || (d2f1 < GC_EPS))
        return 1.0;
    
    double dr = fs - (d2f0 + d2f1);         // positive if point is inside ellipsoid
    
    double gc = 0.0;
    for (unsigned i = 0; i < n_dimensions; i++) {
        t = (pnt_t[i] - f0[i])/d2f0 + (pnt_t[i] - f1[i])/d2f1;
        gc += t * t;
    }
    gc = sqrt(gc);                          // gc = ||U0 + U1||
    
    return 1.0 / (1.0 + exp(-a * dr / gc));
}

//
// diff_n_dim_ellipsoid_gcs -- gradient-corrected derivatives, with scaling.
//
// Notation:
//   pt     = pnt'  (transformed point)
//   d0i    = pt[i] - f0[i]
//   d1i    = pt[i] - f1[i]
//   gc_v[] = U0 + U1  (sum of unit vectors)
//   gc     = ||gc_v||
//   dr     = fs - (d2f0 + d2f1)
//   sprime = phi*(1-phi)
//
// Chain rule for scaled dimensions (i >= 2):
//
//   d(pt[i])/d(f0[i]) = (1 - s[i-2]) / 2
//   d(pt[i])/d(f1[i]) = (1 - s[i-2]) / 2
//   d(pt[i])/d(s[i-2]) = pnt[i] - C[i]
//
// d(phi)/d(pt[j]) via quotient rule on dr/gc:
//   d(dr)/d(pt[j])   = -gc_v[j]
//   d(gc^2)/d(pt[j]) = 2*[ gc_v[j]*(1/d2f0+1/d2f1)
//                          - (gc_v.U0)*U0[j]/d2f0
//                          - (gc_v.U1)*U1[j]/d2f1 ]
//   d(gc)/d(pt[j])   = d(gc^2)/d(pt[j]) / (2*gc)
//   d(dr/gc)/d(pt[j])= d(dr)/d(pt[j])/gc - dr*d(gc)/d(pt[j])/gc^2
//
// d(phi)/d(f0[i]) = direct_f0[i]  +  dphi_dpt[i] * (1-s[i-2])/2   for i >= 2
// d(phi)/d(f1[i]) = direct_f1[i]  +  dphi_dpt[i] * (1-s[i-2])/2   for i >= 2
// d(phi)/d(s[k])  = dphi_dpt[k+2] * (pnt[k+2] - C[k+2])
//
void diff_n_dim_ellipsoid_gcs(double *diff, const double *pnt, const double *param, int n_dimensions)
{
    double  fs = param[0];
    double   a = param[1];
    const double *f0 = param + 2;
    const double *f1 = f0 + n_dimensions;
    const double  *s = f1 + n_dimensions;   // n-2 scale factors
    
    double *df0 = diff + 2;
    double *df1 = diff + 2 + n_dimensions;
    double *ds  = diff + 2 + 2*n_dimensions; // n-2 scale derivatives
    
    // --- Apply scaling transform ---
    double pnt_t[n_dimensions];
    apply_scale(pnt, f0, f1, s, pnt_t, n_dimensions);
    
    // --- Distances in transformed space ---
    double d2f0 = 0.0, d2f1 = 0.0, t;
    for (unsigned i = 0; i < n_dimensions; i++) {
        t = pnt_t[i] - f0[i];  d2f0 += t * t;
        t = pnt_t[i] - f1[i];  d2f1 += t * t;
    }
    d2f0 = sqrt(d2f0);
    d2f1 = sqrt(d2f1);
    
    // Handle singularity
    if ((d2f0 < GC_EPS) || (d2f1 < GC_EPS)) {
        unsigned n_params = 2 + 2*n_dimensions + (n_dimensions > 2 ? n_dimensions-2 : 0);
        for (unsigned i = 0; i < n_params; i++) diff[i] = 0.0;
        return;
    }
    
    // --- gc vector and gc, gcp3_2 via two-sqrt trick ---
    double gc_v[n_dimensions];
    double gc = 0.0;
    for (unsigned j = 0; j < n_dimensions; j++) {
        gc_v[j] = (pnt_t[j]-f0[j])/d2f0 + (pnt_t[j]-f1[j])/d2f1;
        gc += gc_v[j] * gc_v[j];               // accumulate gc^2
    }
    t             = sqrt(gc);                   // t      = gc
    double gcp3_2 = sqrt(gc * t);              // gcp3_2 = gc^(3/2)
    gc            = t;                          // gc     = gc
    
    double dr = fs - (d2f0 + d2f1);
    
    // --- Sigmoid and derivative factor ---
    double phi    = 1.0 / (1.0 + exp(-a * dr / gc));
    double sprime = phi * (1.0 - phi);
    double a_sp   = a * sprime;
    
    // --- d phi / d fs,  d phi / d a ---
    diff[0] = a_sp / gc;
    diff[1] = sprime * dr / gc;
    
    // --- Precompute dot products for direct f0/f1 derivatives ---
    double gc_dot_pnt_f0 = 0.0, gc_dot_pnt_f1 = 0.0;
    for (unsigned j = 0; j < n_dimensions; j++) {
        gc_dot_pnt_f0 += gc_v[j] * (pnt_t[j] - f0[j]);
        gc_dot_pnt_f1 += gc_v[j] * (pnt_t[j] - f1[j]);
    }
    double d2f0_cu    = d2f0 * d2f0 * d2f0;
    double d2f1_cu    = d2f1 * d2f1 * d2f1;
    double two_gc3    = 2.0 * gcp3_2 * gcp3_2;  // = 2*gc^3
    
    // --- d(phi)/d(pnt_t[j]) -- needed for scaled-dimension corrections ---
    // Only compute for i >= 2 (others need no correction).
    // d(dr/gc)/d(pt[j]) = -gc_v[j]/gc - dr*d(gc)/d(pt[j])/gc^2
    // d(gc^2)/d(pt[j])  = 2*[gc_v[j]*(1/d2f0+1/d2f1)
    //                        - (gc_v.U0)*U0[j]/d2f0
    //                        - (gc_v.U1)*U1[j]/d2f1]
    // Precompute gc.U0 and gc.U1
    double gc_dot_U0 = gc_dot_pnt_f0 / d2f0;   // = gc_v . U0  (since U0=(pt-f0)/d2f0)
    double gc_dot_U1 = gc_dot_pnt_f1 / d2f1;
    
    double dphi_dpt[n_dimensions];
    for (unsigned j = 0; j < n_dimensions; j++) {
        double U0j = (pnt_t[j]-f0[j])/d2f0;
        double U1j = (pnt_t[j]-f1[j])/d2f1;
        double dgc2 = 2.0*(gc_v[j]*(1.0/d2f0+1.0/d2f1)
        - gc_dot_U0*U0j/d2f0
        - gc_dot_U1*U1j/d2f1);
        // d(dr/gc)/d(pt[j]) = -gc_v[j]/gc - dr*(dgc2/(2*gc))/gc^2
        //                   = -gc_v[j]/gc - dr*dgc2/(2*gc^3)
        double d_ratio = -gc_v[j]/gc - dr*dgc2/two_gc3;
        dphi_dpt[j] = a_sp * d_ratio;
    }
    
    // --- Direct f0, f1 derivatives + scaling correction for i >= 2 ---
    for (unsigned i = 0; i < n_dimensions; i++) {
        double d0i = pnt_t[i] - f0[i];
        double d1i = pnt_t[i] - f1[i];
        
        // Direct d/df0[i] (pt treated as fixed)
        double dgc2_f0  = 2.0*d0i/d2f0_cu * gc_dot_pnt_f0 - 2.0*gc_v[i]/d2f0;
        double d_rat_f0 = d0i/(d2f0*gc) - dr*dgc2_f0/two_gc3;
        df0[i] = a_sp * d_rat_f0;
        
        // Direct d/df1[i]
        double dgc2_f1  = 2.0*d1i/d2f1_cu * gc_dot_pnt_f1 - 2.0*gc_v[i]/d2f1;
        double d_rat_f1 = d1i/(d2f1*gc) - dr*dgc2_f1/two_gc3;
        df1[i] = a_sp * d_rat_f1;
        
        // Scaling correction for i >= 2:
        // d(pt[i])/d(f0[i]) = d(pt[i])/d(f1[i]) = (1 - s[i-2]) / 2
        if (i >= 2) {
            double dpti_dfoi = (1.0 - s[i-2]) * 0.5;
            df0[i] += dphi_dpt[i] * dpti_dfoi;
            df1[i] += dphi_dpt[i] * dpti_dfoi;  // same correction for f1
        }
    }
    
    // --- d phi / d s[k]  for k = 0..n_dimensions-3 ---
    // d(pt[k+2])/d(s[k]) = pnt[k+2] - C[k+2]
    //                     = pnt[k+2] - (f0[k+2]+f1[k+2])/2
    if (n_dimensions > 2) {
        for (unsigned k = 0; k < n_dimensions-2; k++) {
            unsigned i   = k + 2;
            double Ci    = (f0[i] + f1[i]) * 0.5;
            double dpti_dsk = pnt[i] - Ci;
            ds[k] = dphi_dpt[i] * dpti_dsk;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Associated functions
//

int bef_param_ok (const double *param, int n_dimensions)
//
// Parameter check function to provide some guard rails against LM missbehaving
//
{
    if (param[0] < min_focal_sum) return 0; // Focal sum too small
    
    if (param[1] < min_a) return 0;         // sigmoid shape must remain positive (or mayhem ensures)
    
    double fd = 0.0;
    for (unsigned i = 0; i < n_dimensions; i++) {
        double t = param[2 + i] - param[2 + n_dimensions + i];
        fd += t * t;
    }
    
    if (fd < (min_foci_d * min_foci_d)) return 0; // foci too close
    
    if (n_dimensions > 2) {
        const double *s = param + (2 + 2*n_dimensions);
        for (unsigned i = 0; i < (n_dimensions - 2); i++)
            if ((s[i] < gcs_min_scale) || (s[i] > gcs_max_scale))
                return 0;                   // scale out of range
    }
    
    return 1;               // Parameters are OK
}

// nl_lsq_fit expects an array of function pointers:
double (*bef_function_ptr[1])(const double *x, const double *pa, int ip) =
    {n_dim_ellipsoid_gcs};
void   (*bef_diff_function_ptr[1])(double *d, const double *x, const double *pa, int ip) =
    {diff_n_dim_ellipsoid_gcs};

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Debugging aid
//
void nl_check_nde(unsigned n_dim)
//
// This compares the diff-function with a numerical differentiation. These values ought
// to agree to at least 5 significant digits. An error is likely otherwise.
// 
// Note: Make sure to re-run this test if you edit the above functions.
//
{
    srandom(123);
    assert(n_dim > 1);
    
    double *x = new double[n_dim];	// Some random point in the first quadrant unity n-cube
    // must not 0
    for (unsigned i = 0; i < n_dim; i++)
        x[i] = (double) random() / (double) 0x7fffffff; // the numerator is 2^31-1, max value of random
        
    double *param = new double[2 + 2*n_dim + ((n_dim > 2) ? n_dim - 2 : 0)];
    param[0] = 0.4123;                      // focal_sum > 0
    param[1] = 0.7777;                      // shape factor > 0
    
    for (unsigned i = 0; i < n_dim; i++) {
        // make up two foci within the unity n-cube
        param[i + 2]                = (double) random() / (double) 0x7fffffff;	
        param[i + 2 + n_dim] = (double) random() / (double) 0x7fffffff;	
    }
    
    printf("Checking Diff-function for NDE:\n");
    
    if (n_dim > 2) {
        unsigned n_scales = n_dim - 2;
        for (unsigned i = 0; i < n_scales; i++)
            param[2 + 2*n_dim + i] = 1.5 - ((double) random() / (double) 0x7fffffff); // in [0.5,1.5]
            
        check_diffs(2 + 2*n_dim + n_scales, 0, x, param, n_dim_ellipsoid_gcs, diff_n_dim_ellipsoid_gcs);
    } else
        check_diffs(2 + 2*n_dim, 0, x, param, n_dim_ellipsoid_gcs, diff_n_dim_ellipsoid_gcs);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Misc functions needed to deal with n-dimensional ellipsoids
//

double vect_scalar_product(const double *va, const double *vb, unsigned n_dim)
{
    double s = 0.0;
    for (unsigned i = 0; i < n_dim; i++)
        s += va[i] * vb[i];
    return s;
}

void vect_normalize(double *vect, unsigned n_dim)
// Normalizes vect so that ||vect|| = 1
{
    double sqs = 0.0;
    for (unsigned i = 0; i < n_dim; i++)
        sqs += vect[i]*vect[i];
    
    assert(sqs > 1e-20);
    sqs = 1.0 / sqrt(sqs);
    
    for (unsigned i = 0; i < n_dim; i++)
        vect[i] *= sqs;    
}

int ellipsoid_intersect (
    const double *param,                // Defines the ellipsoid (see above)
    unsigned n_dim,                     // #of dimensions
    const double *point,                // Operating point *inside* the ellipsoid
    const double *dir,                  // Direction of interest: must be normalized to 1
    double &dist_p,                     // result: distance to the ellipsoid shell in the positive direction
    double *dist_n)                     // Optional: distance to the shell in the negative direction
//
// This sfunction computes the distance from a point <point> within an ellipsoid to it's shell
// allong a given direction <dir>.
//
// It returns a 0 upon success and sets the desired distance.
// A return value of 1 signifies a problem (say start point is not within ellipsoid)
// inwhich case the <dist> variable will remain unchanged
//
{
    double fs  = param[0];
    const double *f0 = param + 2;
    const double *f1 = param + 2 + n_dim;
    
    double *a = new double[n_dim];      // = point - f0
    double *b = new double[n_dim];      // = point - f1
    for (unsigned i = 0; i < n_dim; i++) {
        a[i] = point[i] - f0[i];
        b[i] = point[i] - f1[i];
    }
    
    double nasq = vect_scalar_product(a, a, n_dim);     // = |a|^2
    nasq *= nasq;
    double nbsq = vect_scalar_product(b, b, n_dim);     // = |b|^2
    nbsq *= nbsq;
    
    for (unsigned i = 0; i < n_dim; i++)
        a[i] = b[i] - a[i];             // now a = b - a
        
    // intermediate variables
    double r4sq = fs*fs;                // = 4*focal_sum^2 (eventually)
    double L = 2.0 * vect_scalar_product(dir, a, n_dim);  // L = 2<dir>(b - a)
    double K = r4sq + nbsq - nasq;      // K = fs^2 + |b|^2 - |a|^2
    r4sq *= 4.0;
    
    // Coefficients for the quadratic equation  A*dist^2 + B*dist + C = 0
    double A = r4sq - L*L;                                          // A = 4fs^2 - L^2
    double B = 2.0*(r4sq*vect_scalar_product(b, dir, n_dim) - K*L); // B = 2(4fs^2*(b*dir) - K*L)
    double C = r4sq * nbsq - K*K;                                   // C = 4fs^2*|b|^2 - K^2
    
    // Now solve for dist
    double t = B*B - 4.0*A*C;
    if ((t < 0.0) || (fabs(A) < 1.0e-15)) return 1; // No good solution
    t = sqrt(t);
    double d = (-B + t) / (2.0 * A);
    if (d > 0.0) {                      // There is a solution
        dist_p = d;                     // Provide the positive solution
        if (dist_n != nullptr) {
            d = (-B - t) /  (2.0 * A);  // Second solution is for the other direction
            *dist_n = d;
        }
        return 0;
    }
    
    return 1;
}

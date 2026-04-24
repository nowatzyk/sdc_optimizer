//
// The non-linear ellipsoid fit functions and their drivatives
//

#include "binary_ellipsoid_fit.h"

//
// Contraints that are meant to prevent the ellipsoid from degenerating:
//
const double min_focal_sum = 0.01;      // minimal focal_sum of the ellipsoid
                                        // Note: focal_sum is defined as the sum of the distances
                                        // to the two foci. This is twice of the radius of a sphere,
                                        // if you think of a sphere as a degenerated ellipsoid
                                        //  with 0 excentricity.

const double min_a = 0.001;             // The shape factor of the sigmoid function shall not 
// go negative or near 0: the gradient would vanish and LM is lost
const double max_a = 100.0;             // The transition becomes too sharp, no more gradient
                                        // LM gets stuck on a solution can cannot converge anymore

const double min_foci_d = 0.01;         // If the foci merge, their parameters loose independence and
                                        // mayhem would ensure

const double GC_EPS = 1.0e-12;          // in gardiate correction this is a guard against a singularity
                                        // that can happen when the point falls onto a focal point

const double gcs_min_scale = 0.1;       // The gc scale factor lower limit
const double gcs_max_scale = 10.0;      // The gc scale factor upper limit

#define _POK_ELLIPSOID_CENTER_ONLY_     // When defined, the ellipsoid foci may wander out of the
                                        // parameter space as long as the ellipsoid center remains inside.
                                        // This does not seem to be a good idea
                                        // Revise after D2_latch in 5D: It is *much* better. Turns
                                        // out, after an inititial excursion of one focal point, in the
                                        // end, both foci moved back into p-space.

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
    //
    // Parameter layout (as offset into the <param> array:
    // 0                : focal-sum 
    // 1                : coordinates of focal point 0
    // 1+n_dimesnsions  : coordinates of focal point 1
    // 1+2*n_dimensions : scale factors for all dimensions > 1 (n_dimensions - 2) numbers
    // ...+1            : sigmoid shape factor <a>
    //
    // Note: This parameter layout differs from the initial version, which had <a> in second place
    //       This change was needed to allow fitting with *and* without <a> being subject to
    //       the LM LSQ fit procedure, by simply leaving out the last parameter.
    //
    double  fs = param[0];
    const double *f0 = param + 1;
    const double *f1 = f0 + n_dimensions;
    const double  *s = f1 + n_dimensions;   // n_dimensions-2 scale factors
    double   a = param[1 + 2*n_dimensions + ((n_dimensions > 2) ? (n_dimensions - 2) : 0)];
    
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

static double private_a;
void set_private_a(double a)
{
    private_a = a;
}

double n_dim_ellipsoid_gcs_sa(const double *pnt, const double *param, int n_dimensions)
{
    unsigned n = 2 + 2*n_dimensions + ((n_dimensions > 2) ? (n_dimensions - 2) : 0);
    double param_pa[n];
    for(unsigned i = 0; i < (n - 1); i++)
        param_pa[i] = param[i];
    param_pa[n - 1] = private_a;
    return n_dim_ellipsoid_gcs(pnt, param_pa, n_dimensions);
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
    const double *f0 = param + 1;
    const double *f1 = f0 + n_dimensions;
    const double  *s = f1 + n_dimensions;   // n_dimensions-2 scale factors
    double   a = param[1 + 2*n_dimensions + ((n_dimensions > 2) ? (n_dimensions - 2) : 0)];
    
    double *df0 = diff + 1;
    double *df1 = diff + 1 + n_dimensions;
    double *ds  = diff + 1 + 2*n_dimensions; // n-2 scale derivatives
    
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
    diff[1 + 2*n_dimensions + ((n_dimensions > 2) ? (n_dimensions - 2) : 0)] = sprime * dr / gc;
    
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

void diff_n_dim_ellipsoid_gcs_sa(double *diff, const double *pnt, const double *param, int n_dimensions)
{
    unsigned n = 2 + 2*n_dimensions + ((n_dimensions > 2) ? (n_dimensions - 2) : 0);
    double param_pa[n];
    for(unsigned i = 0; i < (n - 1); i++)
        param_pa[i] = param[i];
    param_pa[n - 1] = private_a;
    
    double diff_pa[n];
    diff_n_dim_ellipsoid_gcs(diff_pa, pnt, param_pa, n_dimensions);
    
    for (unsigned i = 0; i < (n - 1); i++)
        diff[i] = diff_pa[i];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Associated functions
//
static unsigned bef_pok_fail_reason = 0;
unsigned get_bef_param_nOK_reason()
    {return bef_pok_fail_reason;}

int bef_param_ok (double *param, int n_dimensions)
//
// Parameter check function to provide some guard rails against LM missbehaving
//
{
    unsigned nOK = 0;

    double a = param[1 + 2*n_dimensions + ((n_dimensions > 2) ? (n_dimensions - 2) : 0)];
    if (a < min_a)                          // sigmoid shape must remain positive (or mayhem ensures)
        nOK |= bef_PnOK_min_a;
    
    if (a > max_a)                          // shape is too sharp (LM can wander off in this direction
        nOK |= bef_PnOK_max_a;
    
    bef_param_ok_sa(param, n_dimensions);   // check the rest
    bef_pok_fail_reason |= nOK;             // combine failure codes
    
    return bef_pok_fail_reason == 0;        // No failure is good (returns 1)
}

int bef_param_ok_sa (double *param, int n_dimensions)
//
// Parameter check function to provide some guard rails against LM missbehaving
//
{
    bef_pok_fail_reason = 0;                // assume no failure

    if (param[0] < min_focal_sum) return 0; // Focal sum too small

    double fd = 0.0;
    for (unsigned i = 0; i < n_dimensions; i++) {
        double t0 = param[1 + i];
        double t1 = param[1 + n_dimensions + i];
#ifdef _POK_ELLIPSOID_CENTER_ONLY_
        double t = (t0 + t1) * 0.5;         // Center of ellipsoid coordinate
        if ((t < 0.0) || (1.0 < t))         // Must be within [0,1]
            bef_pok_fail_reason |= bef_PnOK_cntr_esc; 
#else
        if ((t0 < 0.0) || (1.0 < t0))       // focal point 0 must be in [0,1]
            bef_pok_fail_reason |= bef_PnOK_f0_esc;
        if ((t1 < 0.0) || (1.0 < t1))       // focal point 1 must be in [0,1]
            bef_pok_fail_reason |= bef_PnOK_f1_esc;      
#endif
        t0 -= t1;
        fd += t0 * t0;
    }
    
    if (fd < (min_foci_d * min_foci_d))     // foci too close ?
        bef_pok_fail_reason |= bef_PnOK_f_merge;
    
    if (n_dimensions > 2) {
        const double *s = param + (1 + 2*n_dimensions);
        for (unsigned i = 0; i < (n_dimensions - 2); i++)
            if ((s[i] < gcs_min_scale) || (s[i] > gcs_max_scale))
                bef_pok_fail_reason |= bef_PnOK_scale; // scale out of range
    }
    
    return bef_pok_fail_reason == 0;        // No failure is good (returns 1)
}

//
// nl_lsq_fit expects an array of function pointers:
//
double (*bef_function_ptr[1])(const double *x, const double *pa, int ip) =
    {n_dim_ellipsoid_gcs};
void   (*bef_diff_function_ptr[1])(double *d, const double *x, const double *pa, int ip) =
    {diff_n_dim_ellipsoid_gcs};
    
double (*bef_function_ptr_sa[1])(const double *x, const double *pa, int ip) =
    {n_dim_ellipsoid_gcs_sa};
void   (*bef_diff_function_ptr_sa[1])(double *d, const double *x, const double *pa, int ip) =
    {diff_n_dim_ellipsoid_gcs_sa};
    
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
    
    double *x = new double[n_dim];      // Some random point in the first quadrant unity n-cube
    // must not 0
    for (unsigned i = 0; i < n_dim; i++)
        x[i] = (double) random() / (double) 0x7fffffff; // the numerator is 2^31-1, max value of random
        
    double *param = new double[2 + 2*n_dim + ((n_dim > 2) ? n_dim - 2 : 0)];
    param[0] = 0.8123;                  // focal_sum > 0 (Note: should be larger than the foci distance)
    set_private_a(0.7777);
    param[1 + 2*n_dim + ((n_dim > 2) ? (n_dim - 2) : 0)] = 0.7777; // shape factor > 0

    
    for (unsigned i = 0; i < n_dim; i++) {
        // make up two foci within the unity n-cube
        param[i + 1]         = (double) random() / (double) 0x7fffffff;
        param[i + 1 + n_dim] = (double) random() / (double) 0x7fffffff;
    }
    
    printf("Checking Diff-function for NDE:\n");
    
    if (n_dim > 2) {
        unsigned n_scales = n_dim - 2;
        for (unsigned i = 0; i < n_scales; i++)
            param[1 + 2*n_dim + i] = 1.5 - ((double) random() / (double) 0x7fffffff); // in [0.5,1.5]
            
        check_diffs(2 + 2*n_dim + n_scales, n_dim, x, param, n_dim_ellipsoid_gcs_sa, diff_n_dim_ellipsoid_gcs_sa);
    } else
        check_diffs(2 + 2*n_dim, n_dim, x, param, n_dim_ellipsoid_gcs_sa, diff_n_dim_ellipsoid_gcs_sa);
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
    
    if (sqs < 1e-20) {
        //
        // vect is too short to normalize. Could throw an execption, fail an assertion
        // In order to be a bit more robust, a random direction vector is returned
        //
        fprintf(stderr, "vect_normalize: called on near 0 vector\n");
        generate_random_dir(vect, n_dim);
        return;
    }
    
    sqs = 1.0 / sqrt(sqs);
    
    for (unsigned i = 0; i < n_dim; i++)
        vect[i] *= sqs;    
}

double rnd_normal()
//
// Produces random numbers with gaussian distribution via Box-Mueller transform
//
{
    static unsigned valid = 0;
    static double rnd1;
    
    if (valid == 1) {
        valid = 0;
        return rnd1;
    }
    
    double s, c;
    sincos(2.0 * M_PI * rnd_01d(), &s, &c);
    
    double t;
    do {
        t = -2.0 * log(rnd_01d());
        t = sqrt(t);
    } while (isnan(t) || !isfinite(t));         // Just to be not annoyed by corner cases
    
    valid = 1;
    rnd1 = s * t;
    return c * t;
}

void generate_random_dir(double *dir, unsigned n_dim)
//
// Populated <dir> with a random direction, normalized to 1
//
{
    for (unsigned i = 0; i < n_dim; i++)
        dir[i] = rnd_normal();
    vect_normalize(dir, n_dim);
}

double clip_to_unity(const double *pnt, const double *dir, unsigned n_dim)
//
// Given a point <pnt> inside the unity <n_dim> dimensional hypercube (all coordinates in [0,1]) and
// a normalized direction vector <dir>, return the distance from <pnt> along the direction <dir> to a
// point on the surface of the unity cube.
//
// Returns -1, if the point is outside the unity cube, or <dir> has issues
//
{
    // Validate point is inside unit cube
    for (unsigned i = 0; i < n_dim; i++)
        if (pnt[i] < 0.0 || pnt[i] > 1.0)
            return -1.0;
        
    double d = __DBL_MAX__;
        
    for (unsigned i = 0; i < n_dim; i++) {
        if (dir[i] == 0.0)
            continue;                   // no movement in this dimension
                
        double d_i;
        if (dir[i] > 0.0) {
            if (pnt[i] >= 1.0) continue;// on this wall, ray runs along/away
            d_i = (1.0 - pnt[i]) / dir[i];
        } else {
            if (pnt[i] <= 0.0) continue;// on this wall, ray runs along/away
            d_i = -pnt[i] / dir[i];     // dir[i] < 0, result > 0
        }
            
        d = min(d, d_i);
    }
        
    if (d == __DBL_MAX__ || d < 0.0)
        return -1.0;    // degenerate: all dirs zero, or point on cube with ray pointing out
            
    return d;
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
// A return value of 1 signifies a problem (say the start point is not within the ellipsoid)
// in which case the <dist> variable will remain unchanged
//
{
    double fs  = param[0];
    const double *f0 = param + 1;
    const double *f1 = param + 1 + n_dim;
    
    double *a = new double[n_dim];      // = point - f0
    double *b = new double[n_dim];      // = point - f1
    for (unsigned i = 0; i < n_dim; i++) {
        a[i] = point[i] - f0[i];
        b[i] = point[i] - f1[i];
    }
    
    double nasq = vect_scalar_product(a, a, n_dim);     // = |a|^2
    double nbsq = vect_scalar_product(b, b, n_dim);     // = |b|^2
    
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
    if ((t < 0.0) || (fabs(A) < 1.0e-15))
        return 1;                       // No good solution

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

Eigen::MatrixXd completeOrthogonalBasis(const Eigen::VectorXd &e1)
//
// Given a unit vector e1 (the major axis direction),
// returns a matrix whose columns are n-1 orthonormal vectors
// perpendicular to e1, computed via Modified Gramm-Schmidt (MGS) with column pivoting
//
{
    const int n = e1.size();
    if (n < 2) throw std::invalid_argument("n must be >= 2");
    
    // --- Step 1: Build a pool of n candidate vectors.
    // Use the n standard basis vectors e_i. At least n-1 of them
    // are linearly independent of e1 (at most one can be parallel).
    // We will pick the n-1 best ones via pivoting.
    std::vector<Eigen::VectorXd> pool;
    pool.reserve(n);
    for (int i = 0; i < n; i++) {
        Eigen::VectorXd ei = Eigen::VectorXd::Zero(n);
        ei(i) = 1.0;
        // Project out e1 component immediately
        ei -= ei.dot(e1) * e1;
        pool.push_back(ei);
    }
    
    // --- Step 2: MGS with pivoting to extract n-1 orthonormal vectors.
    Eigen::MatrixXd basis(n, n - 1);
    
    for (int k = 0; k < n - 1; k++) {
        // Pivot: find the vector in the pool with the largest residual norm
        int pivot = -1;
        double bestNorm = -1.0;
        for (int j = k; j < n; j++) {
            double nm = pool[j].norm();
            if (nm > bestNorm) {
                bestNorm = nm;
                pivot = j;
            }
        }
        
        if (bestNorm < 1e-10)
            throw std::runtime_error("Numerical rank deficiency in basis completion");
        
        // Swap pivot into position k
        std::swap(pool[k], pool[pivot]);
        
        // Normalize to get the k-th basis vector
        basis.col(k) = pool[k] / bestNorm;
        
        // Project out this direction from all remaining vectors (MGS step)
        for (int j = k + 1; j < n; j++)
            pool[j] -= pool[j].dot(basis.col(k)) * basis.col(k);
    }
    
    return basis;  // n x (n-1), columns are orthonormal, all perp to e1
}

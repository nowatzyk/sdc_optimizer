//
// This is a sub-system to the BO functions. Specifically, it approxiamtes a multi-dimensional
// Shmoo plot with a n-dimensional hyper-ellipsoid. It tries to find the largest (most volume)
// such ellipsoid, reports iis center, and the marging for each parameter
//

#include "binary_ellipsoid_fit.h"

const unsigned n_lm_iteration = 150;            // #of of iteration for the LM solver
#define _BIN_EFIT_DEBUG_                        // When defined, produce extra debugging info on stdout

//#define _ADD_SYNTHETIC_FAILS_                 // When defined, it implements the idea that a fit is better
                                                // when it has about the same number of fail ponts as it does
                                                // pass points. So synthetic point are added when in the ray
                                                // search a fail point was encountered and the outward direction
                                                // has a search budget left. In this case, extra fail points
                                                // are distributed along the outward ray without actually running
                                                // a simulation. If the pass region is convex, these all must fail.
                                                // And if it isn't convex, that island of functionality should not 
                                                // be used anyway.
                                                // Well, it turns out that that is a bad idea: the results are
                                                // noticably worse. The pass region shrinks a litte, as expected
                                                // as the fit becomes more conservative. That is good. But it is also,
                                                // pushed in the wrong direction because of the too high influencs
                                                // of the clusters of fail point. I left the code here to document
                                                // a not so good idea. But also to allow more experimentation with
                                                // this idea: the nl-lsq-fit procedure allows weights being given
                                                // to data points. So by reducing the weight, it may be possible
                                                // to preserve the positive impact (more conservative extimate).
                                                
#define _BIN_EFIT_PASS_ONLY_OUTLIER_            // If defined, only pass points will be considered to be an outlier
                                                // The rationale is that for a fail to be considered to be
                                                // an outlier, the ellipsoid extends over regions where the
                                                // the simulation failed. This is is not a conservative thing to do

//#define _BIN_EFIT_EXCLUDE_SIG_SHAPE_            // If defined, the sigmoid shape factor is not part of the 
                                                // LM LSQ-fit procedure. It has a tendency to run away so that
                                                // the fitted function becoems a setp function, which has
                                                // undesirable consequences (lacks gradient to center the
                                                // ellipsoid)
                                                
///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  The points being expored are kept in an array to be replayed to the non-linear least square fit 
//  procedure. By managing this data collection explicitly (as opposed to let nl_lsq_fit handle it)
//  it is possible to filter out outliers, observe what the algorithm is doing, etc.
//

explored_points::explored_points(unsigned nd, unsigned sz) :
    n_dim(nd), n_alloc(sz), n_used(0), item_size(nd + 1)
{
    assert(sizeof(double) == sizeof(uint64_t)); // Just make sure the memory layout is sound

    assert(n_alloc > 1);
    data = (double*) std::malloc(item_size * n_alloc * sizeof(double));
    if (!data) throw std::bad_alloc();
}

void explored_points::resize(unsigned new_sz)
{
    assert(new_sz > n_alloc);
    n_alloc = new_sz;
    data = (double *) std::realloc(data, item_size * n_alloc * sizeof(double));
}

unsigned explored_points::add_pnt(double *pnt)
{
    if (n_used >= n_alloc)                      // need to add space
        resize(n_alloc + n_alloc/2);            // add 50%

    for (int i = 0; i < n_dim; i++)
        data[n_used * item_size + 1 + i] = pnt[i]; // copy coordinates
    
    meta(n_used) = 0ull;                        // Zero all attributes
    
    n_used += 1;
    return n_used - 1;                          // Returns the index to the new point
}

///////////////////////////////////////////////////////////////////////////////////////////////////


bin_ellipsoid_fit::bin_ellipsoid_fit(FILE* result_fp, vector<const_parameter *>& op, parameter *of_p,
                                     FILE *s_fp, unsigned n_itr, unsigned n_ray_mul, unsigned n_can,
                                     unsigned n_p_p_ray, double of) :
    result_fp(result_fp), opt_params(op), of_ptr(of_p), sum_fp(s_fp)
{
    n_dim = opt_params.size();
    
    // save start point:
    x_start = new double[n_dim];
    for (unsigned i = 0; i < n_dim; i++)
        x_start[i] = opt_params[i]->get_mapped_value();
    
    // Allocate data structures:
    ev_cache = new EvalCache(n_dim);
    exp_pnts = new explored_points(n_dim);

    n_rays = n_ray_mul * n_dim;
    n_candidates = n_can;
    n_probes_p_ray = n_p_p_ray;
    n_iterations = n_itr;
    outlier_frac = of;

    ray_dirs = new double[n_rays * n_dim];
    
    if (eval_pnt(x_start) == 0) {
        fprintf(stderr, "Error: The initial parameters must yield a passing configuration\n");
        exit(1);
    }
    
    explore(x_start);                           // Explore the starting configuration

#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_
    n_e_params = 2 + 2*n_dim + ((n_dim > 2) ? (n_dim - 2) : 0);
    e_params = new double[n_e_params];
    estimate_initial_e_params();
    
    e_fit = new nl_lsq_fit(n_e_params, n_dim, 1,
                           bef_function_ptr, bef_diff_function_ptr, bef_param_ok, n_dim, 0);
#else
    n_e_params = 2 + 2*n_dim + ((n_dim > 2) ? (n_dim - 2) : 0) - 1;
    e_params = new double[n_e_params];
    estimate_initial_e_params();
    
    e_fit = new nl_lsq_fit(n_e_params, n_dim, 1,
                           bef_function_ptr_1, bef_diff_function_ptr_1, bef_param_ok_1, n_dim, 0);
#endif
    
    //
    // Now the actual work is done here. It is modular, so it could beacome a separate module...
    //
    for (unsigned i = 0; i < n_iterations; i++) {
        if (i > 0) {
            explore(x_start);
            reject_outliers();
        }
        solve();
    }
    print_results();
#ifdef _BIN_EFIT_DEBUG_
    {
        FILE *of = fopen("bef_points.dat", "w");    // Dump all points
        assert(of != nullptr);
        fprintf(of, "#");                           // Mark header line as gnuplot comment
        for (unsigned i = 0; i < n_dim; i++)
            fprintf (of, " (%u):%s", i + 1, opt_params[i]->get_name());
        fprintf(of, " (%u):value (%u):fit_val (%u):outlier (%u):synt_pnt (%u):squashed\n",
                      n_dim+1,   n_dim+2,     n_dim+3,     n_dim+4,      n_dim+5);
        for(unsigned i = 0; i < exp_pnts->size(); i++) {
            for (unsigned j = 0; j < n_dim; j++)
                fprintf(of, "%.5lg ", opt_params[j]->map_01_to_parm(exp_pnts->value(i)[j]));

            fprintf(of, "%u %.5lg %u %u %u\n", (exp_pnts->meta(i) & bef_pnt_value) != 0ull,
                    e_fit->eval(exp_pnts->value(i)), (exp_pnts->meta(i) & bef_pnt_outlier) != 0ull,
                    (exp_pnts->meta(i) & bef_pnt_synth) != 0ull, (exp_pnts->meta(i) & bef_pnt_squashed) != 0ull);
        }
        fclose(of);
    }
#endif
}

unsigned bin_ellipsoid_fit::eval_pnt(double *pnt)
{
    double score = ev_cache->lookup(pnt);
    if (isfinite(score))
        return score > 0.0;                     // Cache hit: we are done!
        
    for (unsigned i = 0; i < n_dim; i++)
        opt_params[i]->set_mapped_value(pnt[i]);
    
    loop_complex.run_once(sum_fp);
    score = of_ptr->get_cur_value();
    if (!isfinite(score)) score = -DBL_MAX;     // NaN -> deep fail
    
    ev_cache->store(pnt, score);                // add datum to the cache
    
    unsigned i = exp_pnts->add_pnt(pnt);
    if (score > 0.0) {
        exp_pnts->meta(i) |= bef_pnt_value;
        return 1;
    }
    return 0;
}

void bin_ellipsoid_fit::explore(double* pnt)
    //
    // The region about the point <pnt> is explored by first casting rays in the direction of
    // the cartesian coordinates. Then <n_r_dirs> are cast in random directions.
    //
    // Note: There is a budget for how much computation is pratical, and probing in random
    //       directions is not really a good way to go about this. But finding a decoration of
    //       point on the surface of an n-dimensional sphere so that the minimal distance between
    //       any two points is maximized is hard.
    //
{
    // Set 1: let's explore the directions parallel to the Cartesian coordinate axis
    double dir[n_dim];
    for (unsigned i = 0; i < n_dim; i++) {
        for (unsigned j = 0; j < n_dim; j++)
            dir[j] = 0.0;
        dir[i] = 1.0;
        ray_search(pnt, dir);                   // Search in the positive direction
        dir[i] = -1.0;
        ray_search(pnt, dir);                   // dito in the negative direction
    }
    
    // Set 2: let's explore additional directions
    //        It is desirable to look in additional directions. Ideally it would be nice to
    //        use an algorithm to distribute points on the n-dimensional unity sphere so that
    //        the minimal distance between any two points is maximized (most uniform distribution
    //        of points on that shpere). But that is very hard to compute. This is a cheap,
    //        greedy approximation. For each direction to be added, <n_candidates> are randomly
    //        generated, and the cos-similarity to a prior directions is computed. The max-of these
    //        similarites is use as a measure on how similar these directions are to the existing ones.
    //        A CS of 1, means that the directions are identical, so that is not a good candidate.
    //
    for (unsigned i = 0; i < n_rays; i++) {
        double best_dir[n_dim];
        double best_score = 2.0;
        
        for (unsigned j = 0; j < n_candidates; j++) {
            for (unsigned k = 0; k < n_dim; k++)
                dir[k] = rnd_01d() - 0.5;
            vect_normalize(dir, n_dim);
            // <dir> is the candidate 
            
            double max_cs = 0.0;
            // first check the cos-similarity wrt. the coordinate axes
            for (unsigned k = 0; k < n_dim; k++)
                max_cs = fmax(max_cs, fabs(dir[k])); // Note: both axis directions are considered:
                                                     //       that is why the fabs is used
                
            // Now let's check wrt. to the prior extra directions:
            // Yes, this is an n_candidates^2 operation, so makeing n_candidates large is costly.
            for (unsigned k = 0; k < i; k++)
                max_cs = fmax(max_cs, vect_scalar_product(dir, ray_dirs + (k * n_dim), n_dim));
            
            if (max_cs < best_score) {
                // Got a better direction to explore:
                best_score = max_cs;
                for (unsigned k = 0; k < n_dim; k++)
                    best_dir[k] = dir[k];
            }
        }
        
        assert(best_score < 1.0);       // There has to be a new direction
        for(unsigned j = 0; j < n_dim; j++)
            dir[j] = best_dir[j];       // Need a copy because ray_search() will scale the direction
        ray_search(pnt, dir);           // explore it
        
        for(unsigned k = 0; k < n_dim; k++)
            ray_dirs[(i * n_dim) + k] = best_dir[k];    // save ray direction (wasted on the last iteration)
    }
}

void bin_ellipsoid_fit::ray_search(double *pnt, double *dir)
// Explore points along the direction <dir> from the origin point <pnt>
{
    //
    // scale <dir> so that <pnt>+<dir> is a point on the unity parameter bounding box
    //
    double s = 2.0;     // because ||dir|| = 1, <pnt>+2*<dir> is guaranteed to be outside. <s> can only shrink
    for (unsigned i = 0; i < n_dim; i++) {
        double d_min = -pnt[i];
        double d_max = 1.0 - pnt[i];
        // s*dir[i] must be in [d_min,d_max]
        if ((s * dir[i]) > d_max) s = d_max / dir[i];
        if ((s * dir[i]) < d_min) s = d_min / dir[i];
    }
    for (unsigned i = 0; i < n_dim; i++)
        dir[i] *= s;
    
    unsigned first = exp_pnts->size();      // This will be the frst point to add (will need this later)
    
    double probe_pnt[n_dim];                // Probe the end-point
    for (unsigned i = 0; i < n_dim; i++)
        probe_pnt[i] = fmin(1.0, fmax(0.0, pnt[i] + dir[i]));
        // Note: the clipping is needed because rounding errors can conspire to move the point
        //       slightly outside the parameter range of [0,1] that will cause trouble downstream

    if (eval_pnt(probe_pnt) > 0)            // The end point is a pass (un-expected!)
        ray_search_mode1(0.0, 1.0, pnt, dir, n_probes_p_ray);
    else                                    // The end point is a fail
        ray_search_mode0(0.0, 1.0, pnt, dir, n_probes_p_ray);
    
    //
    // At this point, the ray search is complete. But the added points may have encountered a non-convexity
    // issue. If the pass region was convex, the pass point and the fail points (if there are any) will
    // have just one transition from pass to fail when the points are odered by their distance to the origin.
    // Note that the added points do not need to be ordered by their distance. However, if the pass region
    // is not convex, then there may be failed points followed by pass points when ordered by distance from
    // the origin. Such points would degrade the ellipsoid fit. Therefore, the next part will look for this
    // potential issue and squash all the pass points following the closest fail point.
    //
    double min_d_fail = __DBL_MAX__;        // minimal distance^2 to a failed point
    for (unsigned i = first; i < exp_pnts->size(); i++) {
        if ((exp_pnts->meta(i) & bef_pnt_value) != 0ull)
            continue;                       // skip the passing points

        double s = 0.0;                     // accumulate the distance from the origin squared
        for (unsigned j = 0; j < n_dim; j++) {
            double t = pnt[j] - (exp_pnts->value(i))[j];
            s += t * t;
        }
        min_d_fail = fmin(min_d_fail, s);   // No need to apply sqrt(), only the order matters
    }
    
    if (min_d_fail != __DBL_MAX__) {
        //
        // There was at least one failed point, so let's squash all passing points beyond this fail
        //
        for (unsigned i = first; i < exp_pnts->size(); i++) {
            if ((exp_pnts->meta(i) & bef_pnt_value) == 0ull)
                continue;                   // skip the failed points, no need to squash those
                
            double s = 0.0;                 // accumulate the distance from the origin squared
            for (unsigned j = 0; j < n_dim; j++) {
                double t = pnt[j] - (exp_pnts->value(i))[j];
                s += t * t;
            }
            if (s >= min_d_fail) {          // Point to be squashed: pass beyond a fail as seen from origin
                exp_pnts->meta(i) &= ~bef_pnt_value;
                exp_pnts->meta(i) |=  bef_pnt_squashed;
            }
        }
    }
}

void bin_ellipsoid_fit::ray_search_mode0(double ds, double de, double *pnt, double *dir, unsigned n_probes)
    //
    // This function probes <n_probes> points along the ray from the origin point <pnt> to the end point
    // (<pnt> + 1.0 * <Dir>). <dir> is the direction vector that was scaled so that when it is added to
    // the point <pnt> it produces a point on the surface of the parameter bounding box. Both the 
    // start and end points have been evaluated: the start point is a pass and the end point is a fail.
    // <ds> and <de> are the distances from the origin to the start and end point respectively. Distances
    // are expressed as a fraction of the length of the direction vector <dir>: 0 <= {ds,de} <= 1.0
    // Obviously ds < de. <n_probes> are the number of probe points to be used.
    //
    // Given that the start point is a pass and the end point is a fail, the mode0 strategy is a binary search
    // pattern: find the half way point, rinse and repeat
    //
{
    assert((n_probes > 0) && (ds < de) && (ds >= 0) && (de <= 1.0));
    
    double d = (ds + de) * 0.5;
    double probe_pnt[n_dim];
    for (unsigned i = 0; i < n_dim; i++)
        probe_pnt[i] = pnt[i] + d * dir[i];
    
    n_probes--;                                 // About to use up one probe point
    if (eval_pnt(probe_pnt) > 0) {
        //
        // The half-way point is a pass:
        //
        if (n_probes <= 0)
            return;                             // Probe budget is up, we are done.
        unsigned n2 = (n_probes + 1) / 2;
        ray_search_mode0(d, de, pnt, dir, n2);
        n_probes -= n2;
        if (n_probes > 0)
            ray_search_mode1(ds, d, pnt, dir, n_probes);
    } else {
        //
        // The half-way point is a fail
        //
        if (n_probes <= 0)
            return;                             // Probe budget is up, we are done.
            
        ray_search_mode0(ds, d, pnt, dir, (n_probes + 1) / 2); // Spend 1/2 of probe budget going forward
#ifdef _ADD_SYNTHETIC_FAILS_
        n_probes /= 2;                          // spend the rest of it filling up with synthetic fail points
        if (n_probes <= 0)
            return;                             // Probe budget is up, we are done.

        //
        // Now we add n_probes synthetic fail points along the outboud portion of the ray.
        // This does not cost any simulation time, but it helps for the ellipsoid fit by
        // providing more known fail points
        //
        double dd = (de - d) / (double) (n_probes + 1);
        for (unsigned i = 0; i < n_probes; i++) {
            d += dd;
            for (unsigned j = 0; j < n_dim; j++)
                probe_pnt[j] = pnt[j] + d * dir[j];
            unsigned j = exp_pnts->add_pnt(probe_pnt);
            exp_pnts->meta(j) |= bef_pnt_synth; // Mark the extra points as being synthetic fails
        }
#endif
    }
}

void bin_ellipsoid_fit::ray_search_mode1(double ds, double de, double *pnt, double *dir, unsigned n_probes)
    //
    // This funtion is similar to the above function, but both end-points are known passes. In this case
    // a different strategy is used: <n_probes> point will be distributed about the ray. Each point is
    // evaluated, starting with the innermost. If the pass region is convex, all of these evaluations
    // will result in passes. However, if the region is not convex, a fail may happen: in that case,
    // all remaining points will be set to synthetic fails.
    //
{
    assert((n_probes > 0) && (ds < de) && (ds >= 0) && (de <= 1.0));
    
    double dd = (de - ds) / (double) (n_probes + 1);
    double probe_pnt[n_dim];
    unsigned fail_found = 0;
    
    for (unsigned i = 0; i < n_probes; i++) {
        ds += dd;
        for (unsigned j = 0; j < n_dim; j++)
            probe_pnt[j] = pnt[j] + ds * dir[j];
        
        if (!fail_found) {
            if (eval_pnt(probe_pnt) == 0) {     // Evidence of non-convexity: a interior fail
                fail_found = 1;
                if ((n_probes - i) > 0)         // there is probe-point buget left to locate the boundary better
                    ray_search_mode0(ds - dd, ds, pnt, dir, n_probes - i);
            }
        } else {
#ifdef _ADD_SYNTHETIC_FAILS_
            unsigned j = exp_pnts->add_pnt(probe_pnt);
            exp_pnts->meta(j) |= bef_pnt_synth; // Mark the extra points as being synthetic fails
#endif
        }
    }
}

void bin_ellipsoid_fit::estimate_initial_e_params()
    //
    // For the non-linear LSQ fit, a reasonable starting parameter configuration for the ellipsoid
    // is required. This funtion makes an initial guess for such a parameter set based on the probe points
    // accumulated so far.
    // The center of the ellipsoid is set to the average coordinate of all passing points.
    // The direction of the major axis is determined by the two passing points that have the largest distance.
    // The foci are set to be +/- 25% of this distance away from the center and the focal sum is set to be this
    // distance.
    //
{
    assert(exp_pnts->size() > n_e_params);      // Need to have enough points for a fit
    
    //
    // 1. Determine center of interior points
    //
    double cntr[n_dim];
    for (unsigned i = 0; i < n_dim; i++)
        cntr[i] = 0.0;
    
    unsigned n_cntr_pnts = 0;
    for (unsigned i = 0; i < exp_pnts->size(); i++) {
        if (!(exp_pnts->meta(i) & bef_pnt_value))
            continue;                           // skip all failed points
        n_cntr_pnts++;
        
        for (unsigned j = 0; j < n_dim; j++)
            cntr[j] += (exp_pnts->value(i))[j];
    }
    if (n_cntr_pnts < (n_dim + 2)) {
        fprintf(stderr, "Initial exploration yielded too few interior points: increase search budget\n");
        exit(1);
    }
    double t = 1.0 / (double) n_cntr_pnts;
    for (unsigned i = 0; i < n_dim; i++)
        cntr[i] *= t;                           // make it the average (center of the point cloud)

    //
    // 2. Find the two interior points that are most far apart
    //
    unsigned p0 = 0, p1 = 0;
    t = 0.0;
    for (unsigned i = 0; i < (exp_pnts->size() - 1); i++) {
        if (!(exp_pnts->meta(i) & bef_pnt_value))
            continue;                           // skip all failed points
        
        for (unsigned j = i + 1; j < exp_pnts->size(); j++) {
            if (!(exp_pnts->meta(j) & bef_pnt_value))
                continue;                       // skip all failed points

            double d = 0.0;
            for (unsigned k = 0; k < n_dim; k++) {
                double q = (exp_pnts->value(i))[k] - (exp_pnts->value(j))[k];
                d += q*q;
            }
            if (t < d) {
                t = d;
                p0 = i;
                p1 = j;
            }
        }
    }
    assert((t > 0.0) && (p0 != p1));            // we should have found a most distant point pair
    
    //
    // 3. Determine length and direction of the line between p0 and p1
    //
    double axis[n_dim];
    double len = 0.0;
    for (unsigned i = 0; i < n_dim; i++) {
        axis[i] = t = (exp_pnts->value(p0))[i] - (exp_pnts->value(p1))[i];
        len += t * t;
    }
    len = sqrt(len);                            // = length of the vector
    t = 1.0 / len;
    for (unsigned i = 0; i < n_dim; i++)
        axis[i] *= t;                           // Normalize the direction of the major axis to 1
        
    //
    // 3. Create guess for the initial parameter set
    //
    double *ep_ptr =  e_params;
    *ep_ptr++ = len;                            // Radius set to reach p0/p1
    
#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_
    *ep_ptr++ = 1.0;                            // Initial fuzziess
#else
    set_sigmoid_shape(1.0);
#endif
    
    for (unsigned i = 0; i < n_dim; i++) {
        // set initial foci to the center point +/- 1/4 of the axis
        ep_ptr[n_dim] = cntr[i] - 0.25 * len * axis[i];
        *ep_ptr++     = cntr[i] + 0.25 * len * axis[i];
    }
    if (n_dim > 2) {
        ep_ptr += n_dim;
        for (unsigned i = 0; i < (n_dim - 2); i++)
            *ep_ptr++ = 1.0;    // scaling factor for the higher dimensions is set to 1
    }
}

int bin_ellipsoid_fit::solve()
    //
    // This performs one LSQ fit of the n-dimensional ellipsoid using the data collected so far.
    // It is required that the nl_lsq_fit subsystem had been set up and was provided with a suitabe
    // initial guess for the ellipsoid parameters. The Levenberg-Marquardt algorithm used is not very
    // robust and can fail. Initial experimenst so far have not encountered any failures, but those test
    // were not particulrily exhausting.
    //
    // returns 0 on success, 1 on failure
    //
{
#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_
    e_params[1] = 1.0;              // LM tends to make the edge too sharp and paints itself into a corner
#else
    set_sigmoid_shape(1.0);
#endif
    
    int ec = e_fit->init(e_params);
    if (ec != 0) {
        //
        // Note: this means that the parameter check function failed the current e_params set.
        //       This set is either the result of a previous solve, which means the P_ok() should have
        //       checked them, whitch would expose a program bug, OR that the initial geuss for the 
        //       did not meet the P_ok() requirements, which is also a bug.
        //
        fprintf(stderr, "The initial ellipse parameter estimation was not viable: bad P_ok()\n");
        exit(1);
    }
    
    for (unsigned i = 0; i < n_lm_iteration; i++) { // One LM iteration
        
        for (unsigned j = 0; j < exp_pnts->size(); j++) {
            if ((exp_pnts->meta(j) & bef_pnt_outlier) != 0ull)
                continue;                       // skip outliers

            double f = 0.0;
            if ((exp_pnts->meta(j) & bef_pnt_value) != 0ull)
                f = 1.0;
            ec = e_fit->add_datum(exp_pnts->value(j), f);
            assert(ec == 0);                    // adding data failing is a sign of a bug
        }
        
        double cur_res, new_res;                // Current and new residuals
        ec = e_fit->solve_1s(cur_res, new_res);
        e_fit->get_params(e_params);            // retrieve the new ellipsoid parameters
        
#ifdef _BIN_EFIT_DEBUG_
        printf(" %2u : cr= %.6lg  nr=%.6lg ec=%d - ", i,  cur_res, new_res, ec);
        for (unsigned q = 0; q < n_e_params; q++)
            printf(" %.4lg", e_params[q]);
        printf("\n");
#endif
        
        if (ec == 1)
            break;
        
        if (ec == -3) {
            fprintf(stderr, "LM iteration %d max a issue\n", i);
            break;
        }
        
        if (ec < 0) {
            fprintf(stderr, "LM iteration %d failed: ec=%d\n", i, ec);
            return 1;
        }
#ifdef _BIN_EFIT_EXCLUDE_SIG_SHAPE_
        set_sigmoid_shape(1.0 + fmin(30.0, 0.5 * (double) i));
#endif
    }

#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_
    for (unsigned i = 0; i < n_dim; i++)        // relocate x_start to the new ellipsoid center
        x_start[i] = (e_params[2 + i] + e_params[2 + n_dim + i]) * 0.5;
#else
    for (unsigned i = 0; i < n_dim; i++)        // relocate x_start to the new ellipsoid center
        x_start[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
#endif
    
#ifdef _BIN_EFIT_DEBUG_
    //
    // Print all hyper eillipsoids
    //
    static unsigned n_solve = 0;
    for (unsigned i = 0; i < (n_dim - 1); i++)
        for (unsigned j = i + 1; j < n_dim; j++) {
            char buf[128];
            sprintf(buf, "bef_ellipsoid_s%u_%u_%u.dat", n_solve, i, j);
            bin_ellipsoid_fit::print_elliosoid(buf, i, j);
        }
    n_solve++;
#endif

    return 0;
}

struct outlier_rec {
    unsigned        ep_index;                   // Index to point in <exp_pnts>
    double          discrepancy;                // absolute difference between actual value and fitted value
};

int oulier_key(const void *a, const void *b)
{
    if (((outlier_rec *) a)->discrepancy < ((outlier_rec *) b)->discrepancy) return  1;
    if (((outlier_rec *) a)->discrepancy > ((outlier_rec *) b)->discrepancy) return -1;
    return 0;
}

void bin_ellipsoid_fit::reject_outliers()
    //
    // The pass region in the Shmoo plot can be non-convex and may have wierd shapes that not fit well
    // with the idea of approximating the pass region with an n-dimensional ellipsoid. The fit can be improved
    // by rejecting points that do not fit well. However this is a dicy preposition: it can easily succomb to
    // confirmation bias. Hence the fraction of data points to be rejected as outliners should be small.
    //
    // The outlier rejection ranks data point by the absolute difference between the fitted value and the 
    // actual value. It considers all data points, including those that previouslu were rejected as outliers.
    // The <otlier_fac> fraction of points will then be marked as outlier 
    //
    // Pre-requisite: there must have been a successfull lsq-fit
    //
{
    assert((0.0 <= outlier_frac) && (outlier_frac <= 1.0));
    unsigned n_out = (unsigned) floor((double)(exp_pnts->size()) * outlier_frac);
    if (n_out < 1)
        return;                                 // Nothing to do: no outliers
                                                // Note: the number of data points is increasing in each
                                                // iteration, so n_out will be monotonically increasing.
                                                // Thus if n_out is 0, then there cannot be any outliers so far,
                                                // hence it is unnecessary to clear the outlier flags.
    assert(n_out < exp_pnts->size());

    outlier_rec *o_rec = new outlier_rec[exp_pnts->size()];
    
    unsigned noc = 0;
    for (unsigned i = 0; i < exp_pnts->size(); i++) {
        exp_pnts->meta(i) &= ~bef_pnt_outlier;  // clear previous outlier designation
        
        if ((exp_pnts->meta(i) & (bef_pnt_synth | bef_pnt_squashed)) != 0ull)
            continue;                           // synthetic or squashed points are never outliers
            
#ifdef _BIN_EFIT_PASS_ONLY_OUTLIER_
        if ((exp_pnts->meta(i) & bef_pnt_value) == 0ull)
            continue;                           // Fails are not considered outliers
#endif

        //
        // Note: it may be a good idea to consider only points with a passing value for being an outlier
        //       on the grounds that that is the conservative thing to do and because the most likely
        //       source for ouliers are passing regions that are eithernot convex or are disconnected
        //       neither are good design points.
        //
        double f_act = (exp_pnts->meta(i) & bef_pnt_value) ? 1.0 : 0.0;  // actual value
        double f_fit = e_fit->eval(exp_pnts->value(i)); // fitted value
        
        o_rec[noc].discrepancy = fabs(f_act - f_fit);
        o_rec[noc].ep_index = i;
        noc++;                                  // one more outlier candidate
    }
    
    assert(noc > n_out);                        // It doesn't make sense to have more outliers to remove
                                                // than there are outlier candidates!
    
    qsort(o_rec, noc, sizeof(outlier_rec), oulier_key);
    
    for (unsigned i = 0; i < n_out; i++)        // Mark the new set of outliers
        exp_pnts->meta(o_rec[i].ep_index) |= bef_pnt_outlier;
    
    delete[] o_rec;
}

void bin_ellipsoid_fit::print_results()
{
    double dir[n_dim];

    for (unsigned i = 0; i < n_dim; i++) {
        
        for (unsigned j = 0; j < n_dim; j++)
            dir[j] = (i == j) ? 1.0 : 0.0;      // set up direction vector: point along the positive i-th axis
        
        double plus, minus;
#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_        
        int ec = ellipsoid_intersect(e_params, n_dim, x_start, dir, plus, &minus);
#else
        double ep_par[n_e_params + 1];
        sigmoid_shape_include(ep_par, e_params, n_e_params + 1);
        int ec = ellipsoid_intersect(ep_par, n_dim, x_start, dir, plus, &minus);
#endif
        assert(ec == 0);                        // e_intersect is broken if it cannot do this!
        
        if (i > 1) {                            // for dimensions 2, 3, ... undo scaling
#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_  
            plus  /= e_params[2 + 2*n_dim + (i - 2)];
            minus /= e_params[2 + 2*n_dim + (i - 2)];
#else
            plus  /= e_params[2 + 2*n_dim + (i - 2) - 1];
            minus /= e_params[2 + 2*n_dim + (i - 2) - 1];
#endif
        }
        
        //
        // Use the parameter mapping function to translate from fit space [0,1] to the actual parameter values
        //
        double min_value  = opt_params[i]->map_01_to_parm(fmax(0.0, x_start[i] + minus));
        double cntr_value = opt_params[i]->map_01_to_parm(          x_start[i]         );
        double max_value  = opt_params[i]->map_01_to_parm(fmin(1.0, x_start[i] + plus ));
        
        fprintf(result_fp, "%20s : %10.4lg [%10.4lg,%10.4lg] -%6.2lf%% +%6.2lf%%\n",
                opt_params[i]->get_name(), cntr_value, min_value, max_value,
                -(min_value/cntr_value - 1.0) * 100.0, (max_value/cntr_value - 1.0) * 100.0);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Misc. debugging and visualization functions
//

static const unsigned n_ellipsoid_pnts = 100;   // #of points to use when plotting an ellipsoid

void  bin_ellipsoid_fit::print_elliosoid(char *fn, unsigned x, unsigned y)
{
    assert ((x != y) && (x < n_dim) && (y < n_dim));
    FILE *of = fopen(fn, "w");
    assert(of != nullptr);
    
    double cntr[n_dim];
    double dir[n_dim];
    for (unsigned i = 0; i < n_dim; i++)  {             // Find center of ellipsoid
#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_
        cntr[i] = (e_params[2 + i] + e_params[2 + n_dim + i]) * 0.5;
#else
        cntr[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
#endif
        dir[i] = 0.0;
    }
    
    double t = 0.0;
    double dt = (2.0 * M_PI) / (double) (n_ellipsoid_pnts - 1);
    // It is intended that the first and last point overlap so that the line is closed
    
    for (unsigned i = 0; i < n_ellipsoid_pnts; i++) {
        sincos(t, dir + y, dir + x);
        double dist;
#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_
        int ec = ellipsoid_intersect(e_params, n_dim, cntr, dir, dist);
#else
        double params_p1[n_e_params + 1];
        sigmoid_shape_include(params_p1, e_params, n_e_params + 1);
        int ec = ellipsoid_intersect(params_p1, n_dim, cntr, dir, dist);
#endif
        assert(ec == 0);
        for (unsigned j = 0; j < n_dim; j++) {
            if (j) fprintf(of, " ");
            if (j < 2)
                fprintf(of, "%.6lg", opt_params[j]->map_01_to_parm(cntr[j] + dist * dir[j]));
            else
                fprintf(of, "%.6lg", opt_params[j]->map_01_to_parm(cntr[j] + dist * dir[j]
#ifndef _BIN_EFIT_EXCLUDE_SIG_SHAPE_                
                                     / e_params[2 + 2*n_dim + (j - 2)]));
#else
                                     / e_params[1 + 2*n_dim + (j - 2)]));
#endif
        }
        fprintf(of, "\n");
        t += dt;
    }
    fclose(of);
}


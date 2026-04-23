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
                                                
//#define _HP_FILTER_ELLIPSOID_NORMAL_            // If defined, the de-convexification process computes
                                                // the angle bisector to the two lines connecting the foci
                                                // to the fail point. If no defined, the normal is the line
                                                // connecting the fail point to the center of the ellipsoid.
                                                // I would have guessed that using the ellipsoid normal 
                                                // would work better becaue it is more like a tangential
                                                // hyperplane. However, in the D2 latch test case, it is
                                                // clearly worse.
                                                
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
    
    warn_high = 0;                              // No warnings issued so far
    warn_low = 0;
    
    explore(x_start);                           // Explore the starting configuration

    n_e_params = 2 + 2*n_dim + ((n_dim > 2) ? (n_dim - 2) : 0);
    e_params = new double[n_e_params];
    e_params_bak = new double[n_e_params];
    estimate_initial_e_params();
    for (unsigned i = 0; i < n_e_params; i++)
        e_params_bak[i] = e_params[i];          // Keep a copy
    
    e_fit = new nl_lsq_fit(n_e_params, n_dim, 1,
                           bef_function_ptr, bef_diff_function_ptr, bef_param_ok, n_dim, 0);
    e_fit_sa = new nl_lsq_fit(n_e_params - 1, n_dim, 1,
                           bef_function_ptr_sa, bef_diff_function_ptr_sa, bef_param_ok_sa, n_dim, 0);
    
    //
    // Now the actual work is done here. It is modular, so it could beacome a separate module...
    //
    for (unsigned i = 0; i < n_iterations; i++) {
        if (i > 0) {
            //explore(x_start);
            e_shell_search(300);  // <a> should be faily large!
            reject_outliers();
            if (i > 1)
                hp_filter();
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
        fprintf(of, " (%u):value (%u):fit_val (%u):outlier (%u):synt_pnt (%u):squashed (%u):convexified\n",
                      n_dim+1,   n_dim+2,     n_dim+3,     n_dim+4,      n_dim+5,       n_dim+6            );
        for(unsigned i = 0; i < exp_pnts->size(); i++) {
            for (unsigned j = 0; j < n_dim; j++)
                fprintf(of, "%.5lg ", opt_params[j]->map_01_to_parm(exp_pnts->value(i)[j]));

            fprintf(of, "%u %.5lg %u %u %u %u\n", (exp_pnts->meta(i) & bef_pnt_value) != 0ull,
                    e_fit->eval(exp_pnts->value(i)), (exp_pnts->meta(i) & bef_pnt_outlier) != 0ull,
                    (exp_pnts->meta(i) & bef_pnt_synth) != 0ull, (exp_pnts->meta(i) & bef_pnt_squashed) != 0ull,
                    (exp_pnts->meta(i) & bef_convexified) != 0ull);
        }
        fclose(of);
    }
#endif
}

unsigned bin_ellipsoid_fit::eval_pnt(double *pnt, int *ind)
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
    if (ind != nullptr)
        *ind = (int) i;                         // An int is used so that the caller can initialize it to -1
                                                // to see that a point was actually added by testing the sign
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
            generate_random_dir(dir, n_dim);
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
    double s = clip_to_unity(pnt, dir, n_dim);
    assert(s >= 0.0);                       // Could be 0 if the point is on the unity cube surface
    for (unsigned i = 0; i < n_dim; i++)
        dir[i] *= s;
    
    unsigned first = exp_pnts->size();      // This will be the frst point to add (will need this later)
    
    double probe_pnt[n_dim];                // Probe the end-point
    for (unsigned i = 0; i < n_dim; i++)
        probe_pnt[i] = fmin(1.0, fmax(0.0, pnt[i] + dir[i]));
        // Note: the clipping is needed because rounding errors can conspire to move the point
        //       slightly outside the parameter range of [0,1] that will cause trouble downstream

    int ind = -1;                           // may need to update meta-data
    if (eval_pnt(probe_pnt, &ind) > 0) {    // The end point is a pass (un-expected!)
        //
        // This can cause problems: if this happens, it means that there is at least one circuit parameter
        // that is constrained by the designer and not by the circuit failing. In other words, the Shmoo
        // pass region touches the boarder of the parameter space. Consequently, it is possible that the
        // fitting process of the pass ellipsoid can be pushed out of the parameter space. If the center of
        // the ellipsoid exits the parameter space, then the fit will fail because the center will be rejected
        // as it is no longer within the design space. To make this outcome less likeley, the pass point at
        // the parameter space boundary will be turned into a squashed fail point. Also, a warning will
        // be issued to alert the user that this happened.
        //
        if (ind >= 0) {                     // This was not a cache hit
            exp_pnts->meta(ind) |=  bef_pnt_squashed; // squash this point
            exp_pnts->meta(ind) &= ~bef_pnt_value;
        }
        for (unsigned i = 0; i < n_dim; i++) {
            double t = probe_pnt[i];
            if ((fabs(t) < 1.0e-4) && !(warn_low & (1u << i))) {
                warn_low |= 1u << i;        // Will never have more than 32 dimensions!
                fprintf(stderr, "Warning: Parameter %s has a passing region extending below its lower bound\n",
                        opt_params[i]->get_name());
            }
            if ((fabs(t - 1.0) < 1.0e-4) && !(warn_high & (1u << i))) {
                warn_high |= 1u << i;        // Will never have more than 32 dimensions!
                fprintf(stderr, "Warning: Parameter %s has a passing region extending above its upper bound\n",
                        opt_params[i]->get_name());
            }
        }
        ray_search_mode1(0.0, 1.0, pnt, dir, n_probes_p_ray);
    } else                                  // The end point is a fail
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

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// ellipsoid shell search: the ray search above is rather inefficient wrt information geathering about
// the shape of the passing region because most points are either inside or outside the ellipsoid. Hence
// they convey little information. The shell search only consideres points near the surface of the
// ellipsoid. Ideally, if the ellipsod matches the shape of the passing region, this means that there
// is a 50/50 chance that a point on the shell is inside or outside of the passing region, thus determining
// that design point via simulation yields the most information.
//
// The e-shell search requires that there is an approximate ellipsoid established. It will then select random
// points near the surface of the ellipsoid. The radial distribution of these points is governed by the 
// shape factor <a> to the sigmoid function. If <a> is large, then the points will be on the eillipsoid 
// surface. If <a> is small, the radial distribution is wider. The algorithm used here is drawing a
// random number from a uniform distribution over [0,1], and finding the intersect point with the sigmoid
// function of the ellipsoid function.
//
// There are a few notable approximations made here: the direction of the probe point from the ellipsoid center
// are drawn from a uniform distribution on a n-dimensional hyper sphere: it isn't uniform on the 
// ellipsoid, rather has lower density near the poles wrt the major axis. The secon approximation is that
// the distance to the ellipsoid surface is along the ray from the center, not perpendicular to the surface
// and also different from the gradient corrected distance that is used in the ellipsoid definition used
// by the fit. That is A-Ok, because all of this is a heuristic: there is no perfect way to do this, just
// inifinitelt many heuristics: the milage varies with the actual shape of the Shmoo plot.
//

void bin_ellipsoid_fit::e_shell_search(unsigned n)
{
    for (unsigned i = 0; i < n; i++) {
        double dir[n_dim];
        generate_random_dir(dir, n_dim);
        
        double d_e;                             // Distance to ellisoid surface
        int ec = ellipsoid_intersect(e_params, n_dim, x_start, dir, d_e);
        assert(ec == 0);                        // x_start ought to be the ellipsoid center!
        
        double t;
        do {
            do {
                t = rnd_01d();
            } while ((t <= 0.0) || (t >= 1.0)); // Make sure t is in (0,1) excluding 0 and 1
            t = log(1.0 / t - 1.0) / e_params[n_e_params - 1];
                                                // apply inverse sigmoid function to fuzzy shell
        } while ((d_e + t ) <= 0.0);            // There is a exponential tail, truncate it!
        d_e += t;
        
        if (n_dim > 2) {
            //
            // Now we need to take the scale of the dimensions >2 into account:
            //
            for (unsigned i = 2; i < n_dim; i++) // un-scale the direction vector, no longer normalized
                dir[i] /= e_params[1 + 2*n_dim + (i - 2)];
                
            d_e *= sqrt(vect_scalar_product(dir, dir, n_dim));  // correct distance
            vect_normalize(dir, n_dim);         // re-normalize dir (the direction may have changed a bit)
        }
        
        double d_c = clip_to_unity(x_start, dir, n_dim);
        assert(d_c > 0.0);                      // x_start ought to be within the legal parameter space
        
        double probe_pnt[n_dim];                // The point to be probed

        if (d_c < d_e) {
            //
            // The ellipsoid extends beyond the unity cube that limits the available design parameter
            // space. There are now 3 choices for what to do in this case:
            // 1. Ignore this direction and proceed with another direction
            // 2. Probe the intersection point with the unity cube via simulation and use the result.
            // 3. Declare the intersection point a fail (skip simulation) on the ground that it is undesirable
            //    to encourage the ellipsoid to extend beyond the space allowed by the desiger. This
            //    argument was used in the ray search to squash pass points on rays when they hit 
            //    the design space limits.
            //
            // Option *3* is implemented here:
            //
            for (unsigned j = 0; j < n_dim; j++)
                probe_pnt[j] = x_start[j] + d_c * dir[j];
            unsigned ind = exp_pnts->add_pnt(probe_pnt);
            exp_pnts->meta(ind) |= bef_pnt_synth;

        } else {
            //
            // Evaluate the point on the ellipsoid shell:
            //
            for (unsigned j = 0; j < n_dim; j++)
                probe_pnt[j] = x_start[j] + d_e * dir[j];
            eval_pnt(probe_pnt);
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
    
    for (unsigned i = 0; i < n_dim; i++) {
        // set initial foci to the center point +/- 1/4 of the axis
        ep_ptr[n_dim] = cntr[i] - 0.25 * len * axis[i];
        *ep_ptr++     = cntr[i] + 0.25 * len * axis[i];
    }
    ep_ptr += n_dim;
    if (n_dim > 2) {
        for (unsigned i = 0; i < (n_dim - 2); i++)
            *ep_ptr++ = 1.0;                    // scaling factor for the higher dimensions is set to 1
    }
    *ep_ptr = 1.0;                              // Initial fuzziess
}

void bin_ellipsoid_fit::hp_filter()
    //
    // This filter looks for failed interior points that and excludes all points that are outside
    // of a hyper plane that is defined by the failed point that is perpendicular to the angular
    // bisector of the two lines connecting the point to the foci of the current ellipsoid.
    //
    // Pre-requisite: there needs to be a valid ellipsoid fit
    //
{
    double *f0 = e_params + 1;                  // focal point 0
    double *f1 = f0 + n_dim;                    // focal point 1
    
    for (unsigned i = 0; i < exp_pnts->size(); i++)
        exp_pnts->meta(i) &= ~bef_convexified;  // clear flags from previous round
        
    for (unsigned i = 0; i < exp_pnts->size(); i++) {
        if (exp_pnts->meta(i) & (bef_pnt_value | bef_convexified | bef_pnt_synth | bef_pnt_squashed))
            continue;                           // skip passing, convexified, synthetic or squashed points
        
        double *pnt = exp_pnts->value(i);
        double normal[n_dim];
#ifdef _HP_FILTER_ELLIPSOID_NORMAL_
        double dir0[n_dim], dir1[n_dim];        // direction vectors from foci to point
        for (unsigned j = 0; j < n_dim; j++)  {
            dir0[j] = pnt[j] - f0[j];
            dir1[j] = pnt[j] - f1[j];
        }
        vect_normalize(dir0, n_dim);
        vect_normalize(dir1, n_dim);

        for (unsigned j = 0; j < n_dim; j++) 
            normal[j] = dir0[j] + dir1[j];      // normal is the angular bisector, albite not normalized to 1
#else
        for (unsigned j = 0; j < n_dim; j++)
            normal[j] = pnt[j] - (f0[j] + f1[j]) * 0.5;
#endif
        double nl_sq = vect_scalar_product(normal, normal, n_dim);
        if (nl_sq < 0.01)
            continue;                           // point is too close to the major axis: skip to avoid trouble

        
        for (unsigned j = 0; j < exp_pnts->size(); j++) {
            if ((exp_pnts->meta(j) & bef_pnt_value) == 0ull)
                continue;                       // skip failed points
            if (exp_pnts->meta(j) & (bef_convexified | bef_pnt_synth | bef_pnt_squashed))
                continue;                       // skip ...
            if (i == j)
                continue;                       // don't bother with self
                
            //
            // This is an n^2 process, but this is still small potatos compared to running a josim simulation.
            //
            double dir[n_dim];
            for (unsigned k = 0; k < n_dim; k++)
                dir[k] = exp_pnts->value(j)[k] - pnt[k];   // dir0 = direction from pnt to candidate
                
            if (vect_scalar_product(dir, normal, n_dim) > 0.0)
                exp_pnts->meta(j) |= bef_convexified;       // is on the far side of the hyper plane
        }
    }
}

int bin_ellipsoid_fit::solve()
    //
    // This performs one LSQ fit of the n-dimensional ellipsoid using the data collected so far.
    // It is required that the nl_lsq_fit subsystem had been set up and was provided with a suitabe
    // initial guess for the ellipsoid parameters. The Levenberg-Marquardt algorithm used is not very
    // robust and can fail. The primary failure modes observed are:
    //
    // 1. The sigmoid shape factor <a> increases beyond utility. This is actually sensible from the
    //    objective of a LSQ fit: the function being fitted is binay and the a large <a> will make the
    //    approxmination more like a step function. However that is undesirable because it eliminates
    //    any gradient towards a better fit and can lock in a sub-optimal fit. The P_ok() function
    //    can jeject such solution, but LM will not change direction and after the lne search retry
    //    count is exhausted, the LM fit will fail with a -3 return code. Often the solution at that
    //    point is pretty good, and the fit may resume as if convergence was achieved.
    //
    // 2. The ellipsoid is squeezed out of the parameter space: This can happen if the parameter space
    //    is constrained by external factors, for example JJ sizes smaller than x are not realizable,
    //    but the circuit would actually work OK (or better) if small JJ's were used. In this case the
    //    passing region extends beyond the limits of the parameter space. The LM solver then finds
    //    solutions where one focal point moves out of the p-space. Eventually, the center of the
    //    ellipsoid can move out of the p-space (= not all mapped parameters are within [0,1]) and
    //    P_ok() will fail the invalid solution.
    //
    //  Recovery strategies:
    //
    //  a) for a <a> run-away: if the solution was close to convergence: accept the solution as valid
    //                         if the solution has not converged, restart the LM fit but exclude <a>
    //                         from being fitted. instead, gradually increase <a> from 1 to 30.
    //  b) for ellipsoid escape: restart fit without <a>. However, this doesn't really address the underlying
    //                         problem. A better mitigation strategy is to add synthetic fail points
    //                         along the limiting hyperplane in the direction where the opbjective function
    //                         shows no failures. 
    //
    // returns 0 on success, 1 on failure
    //
{
    static unsigned n_solve = 0;    // Just count the solver invocations. Used in warning and error messages
    n_solve++;

    e_params[n_e_params - 1] = 1.0;  // LM tends to make the edge too sharp and paints itself into a corner
  
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
    
    double cur_res, new_res;                // Current and new residuals
    for (unsigned i = 0; i < n_lm_iteration; i++) { // One LM iteration
        
        //
        // Add the data points to the fit system:
        //
        for (unsigned j = 0; j < exp_pnts->size(); j++) {
            if ((exp_pnts->meta(j) & bef_pnt_outlier) != 0ull)
                continue;                       // skip outliers

            double f = 0.0;
            if (((exp_pnts->meta(j) & bef_pnt_value) != 0ull) &&    // If the value is 1
                ((exp_pnts->meta(j) & bef_convexified) == 0ull))    // AND this point is not convexified
                f = 1.0;
            ec = e_fit->add_datum(exp_pnts->value(j), f);
            assert(ec == 0);                    // adding data failing is a sign of a bug
        }
        
        ec = e_fit->solve_1s(cur_res, new_res);
        e_fit->get_params(e_params);            // retrieve the new ellipsoid parameters
        
#ifdef _BIN_EFIT_DEBUG_
        printf(" %2u : cr= %.6lg  nr=%.6lg ec=%d - ", i,  cur_res, new_res, ec);
        for (unsigned q = 0; q < n_e_params; q++)
            printf(" %.4lg", e_params[q]);
        printf("\n");
#endif
        
        if (ec != 0)
            break;
    }
    
    if (ec >= 0) {
        //
        // The LM has converged (ec == 1) or the number of LM iterations is exhausted (ec == 0), in
        // which case the solution is likely OK.
        //
        if (ec == 0)
            fprintf(stderr, "Solve %u: #of LM iterations exhausted. current/new residuals= %.6lg/%.6lg\n",
                    n_solve, cur_res, new_res);

        for (unsigned i = 0; i < n_dim; i++)    // relocate x_start to the new ellipsoid center
            x_start[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
        
        for (unsigned i = 0; i < n_e_params; i++)
            e_params_bak[i] = e_params[i];
        
#ifdef _BIN_EFIT_DEBUG_
        //
        // Print all hyper eillipsoids
        //
        for (unsigned i = 0; i < (n_dim - 1); i++)
            for (unsigned j = i + 1; j < n_dim; j++) {
                char buf[128];
                sprintf(buf, "bef_ellipsoid_s%u_%u_%u.dat", n_solve, i, j);
                bin_ellipsoid_fit::print_elliosoid(buf, i, j);
            }
#endif
        return 0;
    }
    
    if ((ec == -3) &&                                   // Failed due to P_ok() inervention
        (e_params[n_e_params - 1] >= (max_a - 1.0)) &&  // AND caused by the sigmoid shape factor ran away
        (new_res <= cur_res) &&                         // AND progress was being made (aka not recovery)
        ((cur_res - new_res) < 0.1) ) {                 // AND the expected residual change is pretty low
        //
        // Then this solution is good enough!
        //
        fprintf(stderr, "Solve %u: OK solution (<a> run-away recovery). current/new residuals= %.6lg/%.6lg\n",
                n_solve, cur_res, new_res);
        
        for (unsigned i = 0; i < n_dim; i++)    // relocate x_start to the new ellipsoid center
            x_start[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
        
        for (unsigned i = 0; i < n_e_params; i++)
            e_params_bak[i] = e_params[i];
        
#ifdef _BIN_EFIT_DEBUG_
        //
        // Print all hyper eillipsoids
        //
        for (unsigned i = 0; i < (n_dim - 1); i++)
            for (unsigned j = i + 1; j < n_dim; j++) {
                char buf[128];
                sprintf(buf, "bef_ellipsoid_s%u_%u_%u.dat", n_solve, i, j);
                bin_ellipsoid_fit::print_elliosoid(buf, i, j);
            }
#endif
        return 0;
    }
    
    //
    // Test if a foci was moved out of parameter space
    //
    ec = 0;
    for (unsigned i = 0; i < n_dim; i++)
        ec |= (e_params[ 1         + i] < 0.0) || (e_params[ 1         + i] > 1.0) ||
              (e_params[ 1 + n_dim + i] < 0.0) || (e_params[ 1 + n_dim + i] > 1.0) ;
    if (ec)
        build_wall(100);     // If so, try to fix this by adding synthetic fail points

    
    //
    // Plan A has failed, now trying recovery.
    //
    // Note: this recovery did not yet analyze which failure has happened. Ellipsoid escape need
    //       code to add a wall of synthetic fail points to address the problem. That is TBD for now.
    //
    for (unsigned i = 0; i < n_e_params; i++)
        e_params[i] = e_params_bak[i];
    set_private_a(1.0);
    ec = e_fit_sa->init(e_params);
    assert(ec == 0);                            // Worked above, must work now again!
    
    for (unsigned i = 0; i < n_lm_iteration; i++) { // One LM iteration
        
        //
        // Add the data points to the fit system:
        //
        for (unsigned j = 0; j < exp_pnts->size(); j++) {
            if ((exp_pnts->meta(j) & bef_pnt_outlier) != 0ull)
                continue;                       // skip outliers
                
            double f = 0.0;
            if (((exp_pnts->meta(j) & bef_pnt_value) != 0ull) &&    // If the value is 1
                ((exp_pnts->meta(j) & bef_convexified) == 0ull))    // AND this point is not convexified
            f = 1.0;
            ec = e_fit_sa->add_datum(exp_pnts->value(j), f);
            assert(ec == 0);                    // adding data failing is a sign of a bug
        }
        
        ec = e_fit_sa->solve_1s(cur_res, new_res);
        e_fit_sa->get_params(e_params);         // retrieve the new ellipsoid parameters
        
#ifdef _BIN_EFIT_DEBUG_
        printf("SA %2u : cr= %.6lg  nr=%.6lg ec=%d - ", i,  cur_res, new_res, ec);
        for (unsigned q = 0; q < (n_e_params - 1); q++)
            printf(" %.4lg", e_params[q]);
        printf("\n");
#endif
        if (ec != 0)
            break;
        
        set_private_a(1.0 + 30.0 * (double) i / (double) n_lm_iteration);
    }
    
    if (ec >= 0) {
        //
        // The LM has converged (ec == 1) or the number of LM iterations is exhausted (ec == 0), in
        // which case the solution is likely OK.
        //
        if (ec == 0)
            fprintf(stderr, "Solve %u: #of LM iterations exhausted. current/new residuals= %.6lg/%.6lg\n",
                    n_solve, cur_res, new_res);
            
        for (unsigned i = 0; i < n_dim; i++)    // relocate x_start to the new ellipsoid center
            x_start[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
        
        e_params[n_e_params - 1] = 1.0;
        for (unsigned i = 0; i < n_e_params; i++)
            e_params_bak[i] = e_params[i];
        
        #ifdef _BIN_EFIT_DEBUG_
        //
        // Print all hyper eillipsoids
        //
        for (unsigned i = 0; i < (n_dim - 1); i++)
            for (unsigned j = i + 1; j < n_dim; j++) {
                char buf[128];
                sprintf(buf, "bef_ellipsoid_s%u_%u_%u.dat", n_solve, i, j);
                bin_ellipsoid_fit::print_elliosoid(buf, i, j);
            }
            #endif
            return 0;
    }
    
    //
    // No luck today, will try anyway
    //
    fprintf(stderr, "Solve %u: recovery failed - trying anyway. current/new residuals= %.6lg/%.6lg\n",
            n_solve, cur_res, new_res);
    for (unsigned i = 0; i < n_e_params; i++)
        e_params[i] = e_params_bak[i];      // Go back to the back-up starting point
    for (unsigned i = 0; i < n_dim; i++)    // relocate x_start to the new ellipsoid center
        x_start[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
    
    return 1;
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
        int ec = ellipsoid_intersect(e_params, n_dim, x_start, dir, plus, &minus);
        assert(ec == 0);                        // e_intersect is broken if it cannot do this!
        
        if (i > 1) {                            // for dimensions 2, 3, ... undo scaling
            plus  /= e_params[1 + 2*n_dim + (i - 2)];
            minus /= e_params[1 + 2*n_dim + (i - 2)];
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

int bin_ellipsoid_fit::build_wall(unsigned n)
//
// This function addresses the problem caused by passing region that extends beyond the valid parameter space.
// It is meant to be called once the LM lsq fit process has pushed one focal point of the ellipsoid
// outside of the allowable parameter space. This function will then add <n> synthetic fail points along
// the parameter space boundary to build a wall of fail-points that should push pack the ellipsoid.
//
// Prerequisite: <e_params> has a solution where one focal point of the ellipsoid has moved outside of
//               the parameter space.
//
// What is done: The ellipsoid pole is cut off perpendicular to the major axis on the side that escaped
// the parameter space. The cut produces a "surface" (really a sub-space with one less dimension than that
// of the ellipsoid) that is populated with <n> points. These points are connected with lines to that other
// (interior) focal point that are used as directions to shoot rays at the boundary of the parameter space.
// The ray/parameter-space wall intersection points will then be added as synthetic fail points. Points
// are distributed to sample the the solid angle of the cut surface as seen from the interior focal point
// uniformly.
//
// Return 0 upon success
// error returns:
//    1: No focal point ouside the unity box
//    2: Both foci outside the unity box
//    3: Foci too close (distance < 1e-10 to prevent divide by 0 errors)
//    4: Bad focal sum
//
{
    assert((n > 0) && (n_dim > 1));         // This function doesn't make sense for 1 D
    
    //
    // Step 1: locate which focal point is the interior and which one is the exterior one
    //
    double *f0 = e_params + 1;
    double *f1 = f0 + n_dim;
    unsigned f0_is_ext = 0, f1_is_ext = 0;
    for (unsigned i = 0; i < n_dim; i++) {
        if ((f0[i] < 0.0) || (f0[i] > 1.0)) f0_is_ext = 1;
        if ((f1[i] < 0.0) || (f1[i] > 1.0)) f1_is_ext = 1;
    }

    if ((f0_is_ext ^ f1_is_ext) == 0) {
        if (f0_is_ext) {
            fprintf(stderr, "build_wall: double escape is not (yet) supported\n");
            return 2;
        }
        return 1;                           // No focal point escaped: nothing to do
    } else if (f0_is_ext) {
        f1 = e_params + 1;
        f0 = f1 + n_dim;
    }
    Eigen::Map<const Eigen::VectorXd> f_ext(f1, static_cast<Eigen::Index>(n_dim));
    Eigen::Map<const Eigen::VectorXd> f_int(f0, static_cast<Eigen::Index>(n_dim));        
    
    //
    // Step 2: compute the direction of the escape and the escape surface subspace
    //
    Eigen::VectorXd dir = f_ext - f_int;
    double d_foci = dir.norm();
    if (d_foci < 1.0e-10)
        return 3;
    dir *= 1.0 / d_foci;                    // Normalize: |dir| = 1
    Eigen::MatrixXd e_set = completeOrthogonalBasis(dir);  // e_set: escape surface sub-space basis vectors
    
    //
    // Step 3: decide on the radius for the target area
    //
    double d_w = clip_to_unity(f_int.data(), dir.data(), n_dim);  // distance from interior focal point to boundary
    assert((d_w > 0.0) && (d_w < d_foci));  // Must be true because f0 is inside and f1 is ouside the box.
    double fs = e_params[0];                // = focal sum
    double fs_sq = fs * fs;
    double t = fs_sq - d_foci*d_foci  + 2.0*d_foci*d_w;
    t = t*t - 4.0 * fs_sq * d_w*d_w;
    if (t < 0.0)
        return 4;
    double radius = sqrt(t) / (2.0 * fs);   // Radius of the target circle (in un-scaled coordinates)
    Eigen::VectorXd cap_center = f_int + d_w * dir;
    
    //
    // Step 4: add synthetic points
    //
    for (unsigned i = 0; i < n; i++) {
        double rd[n_dim - 1];                // A random sample in the target point subspace
        generate_random_dir(rd, n_dim - 1);
        Eigen::Map<Eigen::VectorXd> r_dir(rd, static_cast<Eigen::Index>(n_dim - 1));
        
        double r = radius * pow(rnd_01d(), 1.0 /((double) (n_dim - 1))); // corrected radial distribution
        Eigen::VectorXd aim = cap_center + r * (e_set * r_dir); // aim point
        if (n_dim > 2) {
            // take care of the scaling for the higher dimensions
            for (unsigned j = 2; j < n_dim; j++) {
                double ci = (f0[j] + f1[j]) * 0.5;  // Center point
                aim[j] = (aim[j] - ci) / e_params[2 + 2*n_dim + (j - 2)] + ci;
            }
        }
        Eigen::VectorXd ray_dir = aim - f_int;  // direction from f-int to aim
        ray_dir.normalize();
        d_w = clip_to_unity(f_int.data(), ray_dir.data(), n_dim); // Distance to the boundary point
        
        Eigen::VectorXd b_pnt = f_int + d_w*ray_dir;    // This ought to be the point we were looking for
        
        unsigned p_ind = exp_pnts->add_pnt(b_pnt.data());
        exp_pnts->meta(p_ind) |= bef_pnt_synth; // Add a synthetic fail point to buid the wall
    }
    
    return 0;
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
        cntr[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
        dir[i] = 0.0;
    }
    
    double t = 0.0;
    double dt = (2.0 * M_PI) / (double) (n_ellipsoid_pnts - 1);
    // It is intended that the first and last point overlap so that the line is closed
    
    for (unsigned i = 0; i < n_ellipsoid_pnts; i++) {
        sincos(t, dir + y, dir + x);
        double dist;
        int ec = ellipsoid_intersect(e_params, n_dim, cntr, dir, dist);
        assert(ec == 0);
        for (unsigned j = 0; j < n_dim; j++) {
            if (j) fprintf(of, " ");
            if (j < 2)
                fprintf(of, "%.6lg", opt_params[j]->map_01_to_parm(cntr[j] + dist * dir[j]));
            else
                fprintf(of, "%.6lg", opt_params[j]->map_01_to_parm(cntr[j] + dist * dir[j]
                                     / e_params[1 + 2*n_dim + (j - 2)]));
        }
        fprintf(of, "\n");
        t += dt;
    }
    fclose(of);
}


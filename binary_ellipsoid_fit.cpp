//
// This is a sub-system to the BO functions. Specifically, it approxiamtes a multi-dimensional
// Shmoo plot with a n-dimensional hyper-ellipsoid. It tries to find the largest (most volume)
// such ellipsoid, reports iis center, and the marging for each parameter
//

#include "binary_ellipsoid_fit.h"

const unsigned n_lm_iteration = 200;            // Max #of of iteration for the LM solver

const unsigned n_best_solutions = 16;           // Beam search limit: only this number of solutions will
                                                // be persued. If there are more, the worse ones will be
                                                // discarded.

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

void explored_points::print_stat(FILE *fp)
//
// Just a summary of the stored points
// (Note: needs handcrafting if meta-data is changed)
//
{
    const unsigned n_flags = 5;
    const char *flag_nemo[] =
        {"Value", "Synthetic", "Squashed", "Outlier", "Convexified"};
        
    unsigned *cnt = new unsigned[1 << n_flags];
    memset(cnt, 0, sizeof(unsigned) * (1 << n_flags));
    
    for (unsigned i = 0; i < n_used; i++)       // Count all flag combinations
        cnt[meta(i) & ((1 << n_flags) - 1)] += 1;
    
    fprintf(fp, " %u explored Points\n", n_used);
    for(unsigned i = 0; i < n_flags; i++)
        fprintf(fp, "%*sV%*s%s\n", 2 + 2*i, "", 4 + 2*(n_flags-i), "", flag_nemo[i]);
    for(unsigned i = 0; i < (1 << n_flags); i++) {
        if (cnt[i] == 0)
            continue;
        fprintf(fp, "  ");
        for (unsigned j = 0; j < n_flags; j++)
            fprintf(fp, "%c ", (i & (1 << j)) ? '1' : '.');
        fprintf(fp, "%u\n", cnt[i]);
    }
    
    delete[] cnt;
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

    estimate_initial_e_params();                // make up an initial fit parameter estimate

    
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
            e_shell_search(600);  // <a> should be faily large!
            reject_outliers();
            if (i > 1)
                hp_filter();
        }
        int ec = solve();
        if (ec < 0)
            printf("iteration %u completed with ec=%d\n", i + 1, ec);
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
// ellipsoid, rather has lower density near the poles wrt the major axis. The second approximation is that
// the distance to the ellipsoid surface is along the ray from the center, not perpendicular to the surface
// and also different from the gradient corrected distance that is used in the ellipsoid definition used
// by the fit. That is A-Ok, because all of this is a heuristic: there is no perfect way to do this, just
// inifinitely many heuristics: the milage varies with the actual shape of the Shmoo plot.
//
// Note: the sigmoid shape factor <a> is explicit and not taken from the <e_params> vector. <a> should be 
//       >= 50 initially and can be increased later when there there is a good fit.
//

void bin_ellipsoid_fit::e_shell_search(unsigned n, double a)
{
    assert(a > 0.1);
    
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
            t = log(1.0 / t - 1.0) / a;         // apply inverse sigmoid function to fuzzy shell
        } while ((d_e + t ) <= 0.0);            // There is a exponential tail, truncate it!
        d_e += t;
        
        double probe_pnt[n_dim];                // The point to be probed
        for (unsigned j = 0; j < n_dim; j++) {
            if (j <= 1)
                probe_pnt[j] = x_start[j] + d_e * dir[j];
            else
                probe_pnt[j] = x_start[j] + d_e * dir[j] /  e_params[1 + 2*n_dim + (j - 2)];
            dir[j] =  probe_pnt[j] - x_start[j];
        }
        d_e = sqrt(vect_scalar_product(dir, dir, n_dim));  // correct distance in scaled coordinates
        vect_normalize(dir, n_dim);             // re-normalize dir (including effect of scaling)
        
        double d_c = clip_to_unity(x_start, dir, n_dim);
        assert(d_c > 0.0);                      // x_start ought to be within the legal parameter space
        

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
    *ep_ptr++ = len * 1.5;                      // Radius set to reach p0/p1 + 50% slack to make an ellipse
    
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
                dir[k] = exp_pnts->value(j)[k] - pnt[k];    // dir0 = direction from pnt to candidate
                
            if (vect_scalar_product(dir, normal, n_dim) > 0.0)
                exp_pnts->meta(j) |= bef_convexified;   // is on the far side of the hyper plane
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The LM fit machinery
//

int solution_key(const void *a, const void *b)
    //
    // Solution sort key:
    //
    // 1. Solutions that completed the LM iteration cleanly come first
    // 2. Slolution with low residuals are better
    //
{
    // Primary key:
    switch ((((lm_solution *) a)->ec >= 0) * 2 + (((lm_solution *) b)->ec >= 0)) {
/*
        case 1:         // B completed without problem while A failed to do so:
            return 1;   // B should be placed first
        case 2:         // A completed without problem while B failed to do so:
            return -1;  // A should come first
*/
        default: {      // Undecided (both failed or both succeded)
            double t = (((lm_solution *) a)->end_res) - (((lm_solution *) b)->end_res);
            if (t < 0.0)
                return -1;  // A's residual is smaller: should come first
            else if (t > 0.0)
                return 1;
            else
                return 0;   // Unlikely! (if control gets here it is most likely due to a bug)
        }
    }
    return 0;           // not reachable
}

void bin_ellipsoid_fit::set_up_proto_solution(double ai, double da)
{
    lm_solution ps;
    ps.a_init = ai;
    ps.a_incr = da;
    ps.param_init = new double[n_e_params];
    ps.param_final = new double[n_e_params];
    for (unsigned i = 0; i < n_e_params; i++)
        ps.param_init[i] = e_params[i];
    solutions.push_back(ps);
}

void bin_ellipsoid_fit::derive_solution(lm_solution &sol)
    //
    // <sol> is a result from a LM LSQ fit attempt. This function takes the
    // end-point of <sol> and makes it the starting point of a new solution.
    //
    // Note: LM lsq-fit attempts are not guaranteed to be successful, thus the
    //       endpoint may *not* pass the bef_param_ok[_sa]() functions. However,
    //       starting points must pass this function. So fixing the parameter set
    //       might be required.
    //
{
    lm_solution ps;
    
    for (unsigned i = 0; i < 3; i++) {
        ps.a_init = sol.a_init;
        ps.a_incr = sol.a_incr;
        
        if (i > 0) {                                // No modify it
            if (ps.a_init > 0) {                    // sans_a case:
                ps.a_init *= (i == 1) ? 1.2 : 0.8;  // change a by +/- 20%
                if ((ps.a_init >= max_a) || (ps.a_init <= min_a))
                    continue;                       // just in case
            } else
                break;                              // fit all: may add turn it into a sans_a case,
                                                    //         but have those already, so why bother.
        }
        
        ps.param_init = new double[n_e_params];
        ps.param_final = new double[n_e_params];
        for (unsigned i = 0; i < n_e_params; i++)
            ps.param_init[i] = sol.param_final[i];
        
        if (sol.pok_fail != 0) {                    // Did the final parameter set pass bef_param_ok[_sa]()?
            //
            // NO: initiate recovery
            //
            // Note: there are other reasons why the LM LSQ-fit may have failed, but those don't really
            //       prevent trying again with more data points which may lead to a different outcome.
            //       The objective here is to create a viable starting point for a fit attempt.
            //
            ps.pok_fail = sol.pok_fail;
            ps.ec = sol.ec;
            if (param_recovery(ps) == 0) {
                // Parameter recovery failed
                fprintf(stderr, "derive_solution: parameter recovery failed - skipping this set\n");
                
                // give up:
                delete[] ps.param_init;
                delete[] ps.param_final;
                return;
            }
        }
    
        solutions.push_back(ps);
    }
}

unsigned  bin_ellipsoid_fit::param_recovery(lm_solution &sol)
    //
    // <sol> is a solution from a LM LSQ-fit attempt that ultimately failed due to an abort
    // because of bef_param_ok[_sa]() detecting a bad parameter configuration.
    // This function implements a number of mitigation strategies to get a better fit.
    // They may not lead to a perfect fit, but the parameters will pass bef_param_ok[_sa]().
    //
    // Returns 0 on failure and 1 on success
{
    assert(sol.ec == -3);                       // verify why we are here

    for (unsigned i = 0; i < n_e_params; i++)
        e_params[i] = sol.param_init[i];
    for (unsigned i = 0; i < n_dim; i++)        // relocate x_start to the new ellipsoid center
        x_start[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
    
    switch (sol.pok_fail) {
        
        case bef_PnOK_min_a:
        case bef_PnOK_max_a:
        case (bef_PnOK_min_a | bef_PnOK_max_a):
            assert(sol.a_init == 0.0);          // shoule happen only is fit all mode
            sol.param_init[n_e_params-1] = 50.0;// Back off <a>
            break;                              // That's all 
            
        case bef_PnOK_f0_esc:
        case bef_PnOK_f1_esc: {
            
            build_wall(150);                    // add synthetic fail points to discourage ellipsoid
                                                // drifting towards this p-space boundary
                                                // Note: this is why x_start/e_params were setup above
            
            double *f = nullptr;
            if ( sol.pok_fail == bef_PnOK_f0_esc)
                f = sol.param_init + 1;         // aka focal point 0
            else
                f = sol.param_init + (1 + n_dim);
            
            // back off focal point towards the center
            double dir[n_dim];
            for (unsigned i = 0; i < n_dim; i++)
                dir[i] = 0.5 - f[i];
            vect_normalize(dir, n_dim);
            // Let's move the offending focal point 0.1 unit towards the parameter space center
            for (unsigned i = 0; i < n_dim; i++)
                f[i] += 0.1 * dir[i];

            break;
        }
        
        case bef_PnOK_cntr_esc: {
            // This is the case of the ellipsoid center escape. Nudge center towards p-space center
            double dir[n_dim];
            for (unsigned i = 0; i < n_dim; i++)
                dir[i] = 0.5 - 0.5 * (sol.param_init[1 + i] + sol.param_init[1 + n_dim + i]);
            vect_normalize(dir, n_dim);
            // Let's move both focal points 0.1 unit towards the parameter space center
            for (unsigned i = 0; i < n_dim; i++) {
                sol.param_init[1 +         i] += 0.1 * dir[i];
                sol.param_init[1 + n_dim + i] += 0.1 * dir[i];
            }
            break;
        }
        
        case bef_PnOK_f_merge: {
            // Foci too close: nudge them apart
            double dir[n_dim];
            for (unsigned i = 0; i < n_dim; i++)
                dir[i] = sol.param_init[1 + i] - sol.param_init[1 + n_dim + i];
            vect_normalize(dir, n_dim);         // = direction from F1 to F0
            
            // Let's move both focal points 0.05 unit away from each other, along the major axis
            for (unsigned i = 0; i < n_dim; i++) {
                sol.param_init[1 +         i] += 0.05 * dir[i];
                sol.param_init[1 + n_dim + i] -= 0.05 * dir[i];
            }
            break;
        }
        
        case bef_PnOK_scale:                // scale factors out of range
            for (unsigned i = 0; i < (n_dim - 2); i++)
                sol.param_init[1 + 2*n_dim + i] = 1.0;  // reset them to 1
            break;
            
        case bef_PnOK_fs2small: {
            // focal sum too small (not an ellipsoid anymore)
            double d[n_dim];
            for (unsigned i = 0; i < n_dim; i++)
                d[i] = sol.param_init[1 + i] - sol.param_init[1 + n_dim + i];
            double f_dist = sqrt(vect_scalar_product(d, d, n_dim));
            sol.param_init[0] = 1.5 * f_dist; // set focal sum to 1.5 * foci distance
            break;
        }
        
        default:
            // can't handle multiple, failures (which should be unlikely)
            fprintf(stderr, "param_recovery failed due to multiple, concurrent problems.\n");
            return 0;
    }
    
    if (sol.a_init == 0.0)
        return bef_param_ok (sol.param_init, n_dim);
    else
        return bef_param_ok_sa (sol.param_init, n_dim);
}
unsigned  bin_ellipsoid_fit::min_param_recovery(double *p_ptr, unsigned sans_a)
//
// This is essentially the same as param_recovery(), but this version does only minimal
// changes while param_recovery() adds some margin that is meant to facilitate a LM LSQ-fit,
// while this version is meant to prepare a parmeter set that is used for the next exploration
// or to be presented as the final result.
//
// Returns 0 on failure and 1 on success
{
    if ((sans_a) ? bef_param_ok_sa (p_ptr, n_dim) : bef_param_ok (p_ptr, n_dim))
        return 1;                               // Nothing wrong with the parameters

    switch (get_bef_param_nOK_reason()) {
        
        case bef_PnOK_min_a:
        case bef_PnOK_max_a:
        case (bef_PnOK_min_a | bef_PnOK_max_a):
            p_ptr[n_e_params-1] = 50.0;         // Back off <a>
            break;                              // That's all 
            
        case bef_PnOK_f0_esc:
        case bef_PnOK_f1_esc: {
            for (unsigned i = 0; i < n_dim; i++) {
                p_ptr[1         + i] = fmin(1.0, fmax(0.0, p_ptr[1         + i]));
                p_ptr[1 + n_dim + i] = fmin(1.0, fmax(0.0, p_ptr[1 + n_dim + i]));
            }
            break;
        }
        
        case bef_PnOK_cntr_esc: {
            // This is the case of the ellipsoid center escape. Nudge away from limit
            for (unsigned i = 0; i < n_dim; i++) {
                double t = 0.5 * (p_ptr[1 + i] + p_ptr[1 + n_dim + i]);
                
                if (t <= 0.0)
                    t = 1.0e-10 - t;            // add a little margin to deal with fp-rounding
                else if (t >= 1.0)
                    t = 1.0 - (t + 1.0e-10);
                else
                    t = 0.0;
                    
                p_ptr[1         + i] += t;
                p_ptr[1 + n_dim + i] += t;
            }
            break;
        }
        
        case bef_PnOK_f_merge: {
            // Foci too close: nudge them apart
            double dir[n_dim];
            for (unsigned i = 0; i < n_dim; i++)
                dir[i] = p_ptr[1 + i] - p_ptr[1 + n_dim + i];
            double actual_fd = sqrt(vect_scalar_product(dir, dir, n_dim));
            vect_normalize(dir, n_dim);         // = direction from F1 to F0
            
            double t = min_foci_d - actual_fd;
            assert(t > 0.0);
            t = 0.5 * (t + 1e-10);

            for (unsigned i = 0; i < n_dim; i++) {
                p_ptr[1         + i] += t * dir[i];
                p_ptr[1 + n_dim + i] -= t * dir[i];
            }
            break;
        }
        
        case bef_PnOK_scale:                // scale factors out of range
            for (unsigned i = 0; i < (n_dim - 2); i++)
                p_ptr[1 + 2*n_dim + i] = fmin(gcs_max_scale, fmax(gcs_min_scale, p_ptr[1 + 2*n_dim + i]));
            break;
            
        case bef_PnOK_fs2small: {
            // focal sum too small (not an ellipsoid anymore)
            double d[n_dim];
            for (unsigned i = 0; i < n_dim; i++)
                d[i] = p_ptr[1 + i] - p_ptr[1 + n_dim + i];
            double f_dist = sqrt(vect_scalar_product(d, d, n_dim));
            p_ptr[0] = 1.01 * f_dist; // set focal 1% above min (Note: this case should not really happen in this context)
            break;
        }
        
        default:
            // can't handle multiple, failures (which should be unlikely)
            return 0;
    }
    
    // Verify that the problem was fixed
    if (sans_a)
        return bef_param_ok_sa (p_ptr, n_dim);
    else
        return bef_param_ok    (p_ptr, n_dim);
}

void bin_ellipsoid_fit::de_duplicate()
    //
    // Removes duplicated solutions
    //
    // Pre-requisit: the solutions must be solved: the final parameter set is compared
    //
    // Note: This is an O(n^2) operation, but n is rather small, so no worries.
    //
{
    for (unsigned i = 0; i < (solutions.size() - 1); i++) {
        unsigned i_sa = solutions[i].a_init > 0.0;
        for (unsigned j = i + 1; j < solutions.size(); j++) {
            if (i_sa != (solutions[j].a_init > 0.0))
                continue;                       // Don't compare apples to oranges
                
            unsigned match = 1;                 // Assume we have a match
            for (unsigned k = 0; k < (n_e_params - i_sa); k++) {
                double pi = solutions[i].param_final[k];
                double pj = solutions[j].param_final[k];
                double t = 0.5 * (fabs(pi) + fabs(pj));
                if (t < 1.0e-8)
                    continue;                   // average is ~zero: parameter matches
                if ((fabs(pi - pj) / t) < 0.001)
                    continue;                   // relative difference is less than 0.1% : match
                match = 0;                      // Not a match
                break;
            }
            
            if (match) {                        // Delete solution j
                delete[] solutions[j].param_init;
                delete[] solutions[j].param_final;
                solutions.erase(solutions.begin() + j);
            }
        }
    }
}

int bin_ellipsoid_fit::solve()
    //
    // This performs one LSQ fit of the n-dimensional ellipsoid using the data collected so far.
    // Because LM iterations are not very stable, a beam-search strategy is used: multiple starting
    // points and fit types are used. Each solve tries several starting points. The best ones are 
    // kept. The LM fitting process is far less expensive than doing data point collections, so
    // multipl attemps at the fit are made and the best ones are kept.
    //
{
    if (solutions.size() <= 0) {
        //
        // There are no prior solution to draw from, so 3 starting points are created from
        // the current, initial parameter set:
        //  1. A regular fit one (wth sigmod shape factor <a> subject to the LSQ fit
        //  2. A fit without <a> being fitted, rather kept it constant at 50
        //  3. A fit with <a> being changed from 10 to 100 (if max iterations)
        //
        set_up_proto_solution(0.0, 0.0);
        set_up_proto_solution(50.0, 0.0);
        set_up_proto_solution(50.0, 0.5);
    }
    
#ifdef _BIN_EFIT_DEBUG_
    exp_pnts->print_stat(stdout);
#endif
    
    for (unsigned i = 0; i < solutions.size(); i++)
        solve_one(solutions[i]);                // Solve all current starting points
        
    qsort(solutions.data(), solutions.size(), sizeof(lm_solution), solution_key);
    
    while (solutions.size() > (n_best_solutions - 4)) {
        //
        // Delete worst solutions (up to 3) to make room for new ones
        //
        delete[] solutions.back().param_init;
        delete[] solutions.back().param_final;
        solutions.pop_back();
    }
    
    unsigned n_data = exp_pnts->size();
    unsigned n_sol = solutions.size();
    for (unsigned i = 0; i < 4; i++)
        derive_solution(solutions[i]);          // add derived solutions: Note it may add more than 1 per call
    if (exp_pnts->size() > n_data)              // data points were added
        n_sol = 0;                              // need to re-do previous solves

    for (unsigned i = n_sol; i < solutions.size(); i++)
        solve_one(solutions[i]);
    de_duplicate();
    qsort(solutions.data(), solutions.size(), sizeof(lm_solution), solution_key);
    
#ifdef _BIN_EFIT_DEBUG_
    for (unsigned i = 0; i < solutions.size(); i++)
        printf(">>> %s%2u: cr=%10.6lg ni=%3u ec=%2d ai=%5.2lg da=%5.2lg\n",
               (solutions[i].a_init > 0.0) ? "SA" : "S ", i, solutions[i].end_res, solutions[i].n_iter,
               solutions[i].ec, solutions[i].a_init, solutions[i].a_incr);
#endif

    //
    // Use best solution to go forward:
    //
    double best_res = solutions[0].end_res;
    int s_sel = -1;
    for (unsigned i = 0; i < solutions.size(); i++) {
        if (solutions[i].end_res > (1.1 * s_sel))
            break;                              // Only consider solution within 10% of best
        if (solutions[i].ec >= 0) {             // Found a good solution that does not need fixing
            s_sel = i;
            break;
        }
    }
    if (s_sel < 0) {
        for (unsigned i = 0; i < solutions.size(); i++)
            if (min_param_recovery(solutions[i].param_final, solutions[i].a_incr > 0.0)) {
                s_sel = i;
                break;
            }
    }
    if (s_sel < 0) {
        fprintf(stderr, "solve: no viable solution found\n");
        exit(1);                                // There is only so much that this approach can do,
                                                // but I think that this case is very unlikely to be encountered
    }
    
    for (unsigned i = 0; i < n_e_params; i++)
        e_params[i] = solutions[s_sel].param_final[i];
    for (unsigned i = 0; i < n_dim; i++)        // relocate x_start to the new ellipsoid center
        x_start[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;

    return solutions[0].ec;
}

void bin_ellipsoid_fit::solve_one(lm_solution &sol)
//
// Performs one nl-lsq fit using LM iteration.
//
// The pre-requisit is that sol->param_init is initialized to a viable starting point. If
// <a_init> is > 0, then a fit without the sigmoid shape factor is performed.
// All other fields of <sol> will be updated. The parameter arrays in <sol> must have been
// allocated by the caller.
//
{
    for (unsigned i = 0; i < n_e_params; i++)   // Copy the parameter set
        e_params[i] = sol.param_init[i];
    
    double a = sol.a_init;
    unsigned sans_a = a > 0.0;

    int ec; 
    if (sans_a)                                  // Initialize the fit system
        ec = e_fit_sa->init(e_params);
    else
        ec = e_fit->init(e_params);

    if (ec != 0) {
        //
        // Note: this means that the parameter check function failed the current e_params set.
        //       This set is either the result of a previous solve, which means the P_ok() should have
        //       checked them, whitch would expose a program bug, OR that the initial guess for the 
        //       did not meet the P_ok() requirements, which is also a bug.
        //
        fprintf(stderr, "The initial ellipse parameter estimation was not viable: bad P_ok()\n");
        exit(1);
    }
    
    double cur_res, new_res;                    // Current and new (estimated) residuals
    for (unsigned i = 0; i < n_lm_iteration; i++) {
        
        set_private_a(a);                       // Note: this does not matter when !sans_a
        a += sol.a_incr;
        
        //
        // Add the data points to the fit system:
        //
        for (unsigned j = 0; j < exp_pnts->size(); j++) {
            if ((exp_pnts->meta(j) & bef_pnt_outlier) != 0ull)
                continue;                       // skip outliers
                
            double f = 0.0;                     // assume the datum is a fail
            if (((exp_pnts->meta(j) & bef_pnt_value) != 0ull) &&    // IF the value is 1
                ((exp_pnts->meta(j) & bef_convexified) == 0ull))    // AND this point is not convexified
                f = 1.0;                        // THEN the datum is a pass
            if (sans_a)
                ec = e_fit_sa->add_datum(exp_pnts->value(j), f);
            else
                ec = e_fit->add_datum(exp_pnts->value(j), f);
            assert(ec == 0);                    // adding data failing is a sign of a bug
        }
        
        if (sans_a) {                           // Perform one LM step and retrieve the new parameters
            ec = e_fit_sa->solve_1s(cur_res, new_res);
            e_fit_sa->get_params(sol.param_final); 
        } else {
            ec = e_fit->solve_1s(cur_res, new_res);
            e_fit->get_params(sol.param_final);
        }

        if (i == 0)                             // record info
            sol.begin_res = cur_res;
        sol.end_res = cur_res;
        sol.n_iter = i + 1;                    // +1 because iteration <i> just completed
        
#ifdef _BIN_EFIT_DEBUG_
        printf(" %2u : cr= %9.6lg  nr=%9.6lg ec=%2d p=", i,  cur_res, new_res, ec);
        for (unsigned q = 0; q < (n_e_params - sans_a); q++)
            printf(" %8.4lg", sol.param_final[q]);
        printf("\n");
#endif
        
        if (ec != 0)
            break;
    }
    
    //
    // LM fit completed (analysis and actions are performed by the caller)
    //
    sol.ec = ec;
    sol.pok_fail = get_bef_param_nOK_reason();
    
    for (unsigned i = 0; i < n_e_params; i++)   // Copy the new parameter set
        e_params[i] = sol.param_final[i];
    for (unsigned i = 0; i < n_dim; i++)        // Debug and other code need x_start to be the e-center
        x_start[i] = 0.5 * (e_params[1 + i] + e_params[1 + n_dim + i]);
    
#ifdef _BIN_EFIT_DEBUG1_
    static unsigned n_solve = 1;
    //
    // Print all hyper eillipsoids
    //
    {
        printf(">>>> ec=%d  PoK=0x%04x\n", sol.ec, sol.pok_fail);
        char buf[128];
        for (unsigned i = 0; i < (n_dim - 1); i++)
            for (unsigned j = i + 1; j < n_dim; j++) {
                sprintf(buf, "bef_ellipsoid_%s%u_%u_%u.dat",
                        (sans_a) ? "sa" : "s", n_solve, i, j);
                bin_ellipsoid_fit::print_elliosoid(buf, i, j);
            }
        sprintf(buf, "bef_ellipsoid_%s%u_axis.dat", (sans_a) ? "sa" : "s", n_solve);
        print_e_major_axis(buf);
    }
    n_solve++;
#endif
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Outlier filter:
//

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
    
    if (0) {    // Debug code to be removed
        printf("f_int: "); print_a_point(stdout, f_int.data());
        printf("f_ext: "); print_a_point(stdout, f_ext.data());
        printf("radius = %.6lg\n", radius);
        printf("cap-center: "); print_a_point(stdout, cap_center.data());
        
        FILE *of = fopen("q.dat", "w");
        double t = 0.0;
        double dt = (2.0 * M_PI) / 99.0;
        Eigen::VectorXd t_dir(2);
        Eigen::VectorXd t_pnt(3);
        
        for (unsigned i = 0; i < 100; i++) {
            sincos(t, &(t_dir[0]), &(t_dir[1]));
            t_pnt = -radius * (e_set * t_dir) + cap_center;
            print_a_point(of, t_pnt.data());
            t += dt;
        }
        
        fclose(of);
    }
    
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
                aim[j] = (aim[j] - ci) / e_params[1 + 2*n_dim + (j - 2)] + ci;
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

void  bin_ellipsoid_fit::print_a_point_ec(FILE *of, const double *pnt, const double *off, double d)
//
// Print a point with optional offset vector and distance (just debugging infrastructure)
//
{
    double point[n_dim];
    
    for (unsigned k = 0; k < n_dim; k++) {
        double t = pnt[k];
        if (off != nullptr)
            t += off[k] * d;
        if (k >= 2) {
            double cntr = 0.5 * (e_params[1 + k] + e_params[1 + n_dim + k]);
            t = (t - cntr) / e_params[1 + 2*n_dim + (k - 2)] + cntr;
        }
        point[k] = t;
    }
    print_a_point(of, point);
}

void  bin_ellipsoid_fit::print_a_point(FILE *of, const double *pnt, const double *off, double d)
//
// Print a point with optional offset vector and distance (just debugging infrastructure)
//
{
    for (unsigned k = 0; k < n_dim; k++) {
        if (k) fprintf(of, " ");
        double t = pnt[k];
        if (off != nullptr)
            t += off[k] * d;
        fprintf(of, "%.6lg", opt_params[k]->map_01_to_parm(t));
    }
    fprintf(of, "\n");
}

void  bin_ellipsoid_fit::print_elliosoid(char *fn, unsigned x, unsigned y)
{
    const unsigned n_slices = 5;                // should be an odd number

    assert ((x != y) && (x < n_dim) && (y < n_dim) && (n_dim >= 2));
    FILE *of = fopen(fn, "w");
    assert(of != nullptr);
    
    double cntr[n_dim];
    double dir[n_dim];
    for (unsigned i = 0; i < n_dim; i++)  {     // Find center of ellipsoid
        cntr[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
        dir[i] = 0.0;
    }
    
    // Find the lowest dimension that is not x or y
    int z = -1;                                 // assume that there is no such thing (n_dim == 2)
    for (int i = 0; i < n_dim; i++)
        if ((i != x) && (i != y)) {
            z = i;
            break;
        }

    double dz = 0.0;
    if (z >= 0) {                               // Determine perpendicular range
        dir[z] = 1.0;
        double z_min, z_max;
        int ec = ellipsoid_intersect(e_params, n_dim, cntr, dir, z_max, &z_min);
        assert(ec == 0);
        dz = (z_max - z_min) / (double) (n_slices + 1);
        dir[z] = 0.0;
        cntr[z] -= dz * (double) (n_slices / 2);
    }
    
    for (unsigned j = 0; j < ((z >= 0) ? n_slices : 1); j++) {
        if (j > 0)
            fprintf(of, "\n");
        
        double t = 0.0;
        double dt = (2.0 * M_PI) / (double) (n_ellipsoid_pnts - 1);
        // It is intended that the first and last point overlap so that the line is closed
        
        for (unsigned i = 0; i < n_ellipsoid_pnts; i++) {
            sincos(t, dir + y, dir + x);
            double dist;
            int ec = ellipsoid_intersect(e_params, n_dim, cntr, dir, dist);
            assert(ec == 0);
            print_a_point_ec(of, cntr, dir, dist);
            t += dt;
        }
        
        if (z >= 0)
            cntr[z] += dz;
    }
    
    fclose(of);
}

void  bin_ellipsoid_fit::print_e_major_axis(char *fn)
{
    FILE *of = fopen(fn, "w");
    assert(of != nullptr);
    
    double cntr[n_dim];
    double dir[n_dim];
    for (unsigned i = 0; i < n_dim; i++)  {         // Find center of ellipsoid
        cntr[i] = (e_params[1 + i] + e_params[1 + n_dim + i]) * 0.5;
        dir[i] = e_params[1 + n_dim + i] - e_params[1 + i];     // Dir: from f0 to f1
    }
    double fc_d = sqrt(vect_scalar_product(dir, dir, n_dim));   // Distance between foci

    vect_normalize(dir, n_dim);
    double d0, d1;
    int ec = ellipsoid_intersect(e_params, n_dim, cntr, dir, d0, &d1);
    assert(ec == 0);
    
    print_a_point(of, cntr);                        // Point 1: ellisoid center
    print_a_point_ec(of, e_params + 1);                // Point 2: f0
    print_a_point_ec(of, cntr, dir, -0.5 * fc_d);      // Point 3: should be the same!
    print_a_point_ec(of, e_params + 1 + n_dim);        // Point 4: f1
    print_a_point_ec(of, cntr, dir,  0.5 * fc_d);      // Point 5: should be the same!    
    print_a_point_ec(of, cntr, dir, d0);               // Point 6: Ellipsoid +pole
    print_a_point_ec(of, cntr, dir, 0.5* e_params[0]); // Point 7: should be the same
    print_a_point_ec(of, cntr, dir, d1);               // Point 8: Ellipsoid -pole
    print_a_point_ec(of, cntr, dir, -0.5*e_params[0]); // Point 9: should be the same

    for (unsigned i = 2; i < n_dim; i++)
        dir[i] /= e_params[1 + 2*n_dim + (i - 2)];
    vect_normalize(dir, n_dim);
    
    double d_c = clip_to_unity(cntr, dir, n_dim);
    print_a_point(of, cntr, dir, d_c);              // Point 10: Point on p-space boundary
    for (unsigned i = 0; i < n_dim; i++)
        dir[i] *= -1.0;
    d_c = clip_to_unity(cntr, dir, n_dim);
    print_a_point(of, cntr, dir, d_c);              // Point 11: Point on p-space boundary, -dir direction
    
    fclose(of);
}

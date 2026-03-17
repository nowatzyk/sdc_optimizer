#include "eval_cache.h"

#include <cmath>
#include <math.h>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <stdio.h>

#include "xrand.h"

constexpr unsigned MAX_PROBE = 8;       // max number of rehash/match attempts          

// ---------------------------------------------------------------------------
// Construction / destruction
// ---------------------------------------------------------------------------

EvalCache::EvalCache(size_t n_params, size_t size, double eps)
    : n_params(n_params)
    , entry_sz(sizeof(Entry) + n_params * sizeof(double))
    , n_entries(0)
    , n_hits(0)
    , n_overflows(0)
{
    assert((n_params > 0) && (eps > 0.0));
    eps_squared = eps * eps;
    
    // Find the first power of 2 that is >= size:
    for (log2_cap = 0; size > (1 << log2_cap); log2_cap++);
    cap = 1 << log2_cap;

    table = new char[cap * entry_sz]();     // value-initialised: all bytes zero
    for (unsigned i = 0; i < cap; i++)
        get_entry(i)->score = NAN;          // initialize cache to be all empty
    
    hyper = new double[log2_cap * n_params];// Allocate storage for the hyper-planes
    double *hp = hyper;
    for (unsigned i = 0; i < log2_cap; i++, hp += n_params)
        for (unsigned j = 0; j < n_params; j++)
            hp[j] = rnd_01d() - 0.5;        // Get a random number in [-0.5,0.5]
            
    vln_scale = (double) cap * pow(1.0/sqrt((double) n_params), (double) n_params);
}

EvalCache::~EvalCache()
{
    delete[] table;
}

// ---------------------------------------------------------------------------
// Hash
// ---------------------------------------------------------------------------
//
// FNV-1a over the raw bytes of the parameter vector.
// Any hash that computes to 0 is bumped to 1 to preserve the empty-slot sentinel.
//

//#define _GRID_LSH_HASH_FUNCTION_      // Don't use: it is bad - just an interesting experiment
//
// Note: Grid-hashing is clearly worse. I tested it with both random and clustered
//       artificial point distributions. The scale-factor (see below, left at 500)
//       relates to EPS and would need to be automated. There is an optimum that
//       depends on EPS, but I just tried a few and the optimum for EPS = 1e-3 is
//       somewhere near 200 to 500. No point in further expolration on what that
//       should be.
//

unsigned EvalCache::compute_hash(const double *p) const
    //
    // This is a Locality Sensitive Hash (LSH) function that uses the projection
    // method. For each bit of the hash function, a random hyper-plane is constructed.
    // The bit is set if the parameter vector is above the plane and cleared otherwise.
    //
{
    unsigned h = 0;
    
#ifdef _GRID_LSH_HASH_FUNCTION_
    for (unsigned i = 0; i < n_params; i++) {
        double t = 500.0 * p[i];
        h ^= (unsigned) nearbyint(t);
        h = h * 66049u + 3907864577u;
    }
    
    h &= cap - 1;
#else
    double *hp = hyper;
    for (unsigned i = 0; i < log2_cap; i++, hp += n_params) {
        double s = 0.0;
        for (unsigned j = 0; j < n_params; j++)
            s += (p[j] - 0.5) * hp[j];
//          s += p[j] * hp[j];          // This would perform worse!
        h = (h << 1) | (s >= 0.0);
    }
    
    //
    // The following computes the vector lengths sqrt(s) according to the Euclidian vector norm.
    // It then determins in which hypershere shell the vector ends. The hypersphere in n_param
    // dimensions with radius equal sqrt(n_params) is subdivided into (cap) shells of equal volume.
    // The shell # corresponding to the shel where the vector p[] ends is eclusive or-ed into
    // the hash: this means the cos-similarity is only applied to vectors of roughly equal length.
    // The equal volume shell constraint means that the hash space is utilized approximatley uniformly.
    //
    // Note: the number of possible hashes should correspond to 1/EPS. For example, if the hash uses
    //       11 bits (2048 possible hashes) good matching with an EPS of 1e-3 results. If the hash
    //       were to use 14 bits, the expected hit rate for matches within EPS drops to about 60%
    //       of optimal because the hash quantizes the vector space into smaller parcles than the
    //       EPS sphere. Points in different hash parcels will not be considered to match, even if their
    //       distance is less than EPS. So to increase cache capacity, associativity needs to be
    //       increased rather than using more bits for the hash.
    //
    double s = 0.0;
    for (unsigned i = 0; i < n_params; i++)
        s += p[i] * p[i];                   // Compute the length of the p-vector
    s = vln_scale * pow(s, 0.5 * (double) n_params);
    unsigned hs = (unsigned) floor(s);
    if (hs >= cap) hs = cap - 1;    // This can happen only if p[i] = 1.0 for all i
    
    h ^= hs;                        // EX-OR the quantified vecor length
//    printf(" %04x %4d", h, hs);
#endif
    return h;
}

// ---------------------------------------------------------------------------
// lookup
// ---------------------------------------------------------------------------

double EvalCache::lookup(const double *p) const
{
    unsigned h = compute_hash(p);

    for (int probe = 0; probe < MAX_PROBE; probe++) {
        unsigned idx = (h + ((probe + 1) * probe)/2) & (cap - 1);
        const Entry  *e = get_entry(idx);

        if (!isfinite(e->score))    // IF the score is a NAN then:
            return NAN;             //   return NAN to indicate a cache miss

        const double *ep = e->params();
        double dist_sq = 0.0;
        for (size_t k = 0; k < n_params; k++) {
            double d = ep[k] - p[k];
            dist_sq += d * d;
        }
        if (dist_sq < eps_squared) {
            n_hits++;
            return e->score;
        }
    }

    return NAN;             // Miss
}

// ---------------------------------------------------------------------------
// store
// ---------------------------------------------------------------------------

void EvalCache::store(const double *p, double score)
{
    if (!isfinite(score))
        return;                     // never cache invalid results

    unsigned h = compute_hash(p);
    for (int probe = 0; probe < MAX_PROBE; probe++) {
        unsigned idx = (h + ((probe + 1) * probe)/2) & (cap - 1);
        Entry  *e = get_entry(idx);
        double *ep = e->params();
        
        if (!isfinite(e->score)) {  // IF the score is a NAN then:
                                    // Add a new cache entry
            e->score = score;
            for (size_t k = 0; k < n_params; k++)
                ep[k] = p[k];
            
            n_entries++;            //an entry was added
            return;
        }
        
        double dist_sq = 0.0;
        for (size_t k = 0; k < n_params; k++) {
            double d = ep[k] - p[k];
            dist_sq += d * d;
        }
        if (dist_sq < eps_squared) {
            n_hits++;               // Hitting an existing entry:
            e->score = (e->score + score) / 2.0; // average scores
            //
            // Note: It might be adventageous to average the p-vector too
            //       in order to increase the chance of further hits.
            //
            return;
        }
    }
    
    // Did not find a place to store this datum
     n_overflows++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Misc
//

void EvalCache::print_stats(FILE *fp)
{
    fprintf(fp, "eval-cache: %lu/%lu used, %lu hits, %lu over-flow events\n", 
            n_entries, cap, n_hits, n_overflows);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Cache test/debugging function
//
//#define _RANDOM_POINT_TEST_           // When defined, just use points normally distributed
                                        // within the unity cube. If not defined, then points will
                                        // be clustered near a certain point in the unity cube.

void test_eval_cache(unsigned n_dim, unsigned n_tests, unsigned n_matches, double eps, double eps_fraction)
    //
    // n_dim        : number of dimensions
    // n_tests      : number of tests to be performed, each being one unique vector
    // n_matches    : number of vectors that ought to match (>= 2)
    // eps          : EPS for a match
    // eps_fraction : fraction of EPS to make up a similar vector
    //
{
    EvalCache *evc = new EvalCache(n_dim, 1000, eps);
    
    double *p = new double[n_dim];
    double *s = new double[n_dim];
    double *p1 = new double[n_dim];

    unsigned n_rejects = 0;
    
    for (unsigned i = 0; i < n_tests; i++) {
        
        // Make up one test vector:
#ifdef _RANDOM_POINT_TEST_
        for (unsigned j = 0; j < n_dim; j++) {
            p[j] = rnd_01d();
        }
#else
        // Clustered about (0.75, 0.75, ...., 0.75)
        
        // 1. Pick a random direction
        double sum = 0.0;
        for (unsigned j = 0; j < n_dim; j++) {
            double t = rnd_01d();
            sum += t;
            s[j] = t;
        }

        // 2. Normalize direction vector to length 1
        double scale_f = 1.0 / sqrt(sum);
        for (unsigned j = 0; j < n_dim; j++)
            s[j] *= scale_f;

        do {
            // 3. Pick a random distance from the cluster point using a negative exponential distribution
            double d = rnd_ned(1.0);    // Adjust lam (here 1.0) to adjust the density of points near the center
                                        // A larger lam causes more points near the cluster center. Smaller lam
                                        // means less clustering (and more rejects).
            unsigned j;
            for (j = 0; j < n_dim; j++) {
                double t = 0.75 + s[j] * d;
                if ((t < 0.0) || (t > 1.0))
                    break;      // outside of unity volume, try again
                p[j] = t;
            }
            if (j >= n_dim)
                break;          // Got a valid point
            n_rejects++;
        } while (1);
#endif

        evc->store(p, 1.0);
        
        for (unsigned k = 0; k < n_matches; k++) {
            double sum = 0.0;
            for (unsigned j = 0; j < n_dim; j++) {
                double t = rnd_01d();
                sum += t;
                s[j] = t;
            }
            double scale_f = eps * eps_fraction / sqrt(sum);
            for (unsigned j = 0; j < n_dim; j++)
                p1[j] = p[j] + scale_f * s[j];
            
            evc->store(p1, 2.0 + (double) k);
        }
//        printf("\n");
    }
    
    printf("n_rejects = %u\n", n_rejects);
    evc->print_stats(stdout);
}

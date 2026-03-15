#include "eval_cache.h"

#include <cmath>
#include <math.h>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <stdio.h>

#include "xrand.h"

constexpr unsigned MAX_PROBE = 4;       // max number of rehash/match attempts          

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

unsigned EvalCache::compute_hash(const double *p) const
    //
    // This is a Locality Sensitive Hash (LSH) function that uses the projection
    // method. For each bit of the hash function, a random hyper-plane is constructed.
    // The bit is set if the parameter vector is above the plane and cleared otherwise.
    //
{
    unsigned h = 0;
    
    double *hp = hyper;
    for (unsigned i = 0; i < log2_cap; i++, hp += n_params) {
        double s = 0.0;
        for (unsigned j = 0; j < n_params; j++)
            s += p[j] * hp[j];
        h = (h << 1) | (s >= 0.0);
    }

    return h ? h : 1;   // disallow 0: it is the empty-slot sentinel
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

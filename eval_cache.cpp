#include "eval_cache.h"

#include <cmath>
#include <cassert>
#include <cstring>
#include <cstdlib>

// ---------------------------------------------------------------------------
// Construction / destruction
// ---------------------------------------------------------------------------

EvalCache::EvalCache(size_t n_params, size_t capacity)
    : n_params(n_params)
    , cap(capacity)
    , entry_sz(sizeof(Entry) + n_params * sizeof(double))
    , n_entries(0)
    , n_hits(0)
{
    // capacity must be a power of 2 so we can use masking instead of modulo
    assert(cap > 0 && (cap & (cap - 1)) == 0);
    assert(n_params > 0);

    table = new char[cap * entry_sz]();  // value-initialised: all bytes zero
                                         // => all hashes == 0 == empty
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

uint64_t EvalCache::compute_hash(const double *p) const
{
    static constexpr uint64_t FNV_OFFSET = 14695981039346656037ULL;
    static constexpr uint64_t FNV_PRIME  = 1099511628211ULL;

    uint64_t h = FNV_OFFSET;
    const unsigned char *bp = reinterpret_cast<const unsigned char *>(p);
    const unsigned char *end = bp + n_params * sizeof(double);

    while (bp < end) {
        h ^= static_cast<uint64_t>(*bp++);
        h *= FNV_PRIME;
    }

    return h ? h : 1;   // disallow 0: it is the empty-slot sentinel
}

// ---------------------------------------------------------------------------
// lookup
// ---------------------------------------------------------------------------

std::optional<double> EvalCache::lookup(const double *p) const
{
    const uint64_t h   = compute_hash(p);
    const size_t   idx = h & (cap - 1);

    for (int probe = 0; probe < MAX_PROBE; probe++) {
        const size_t  i = (idx + static_cast<size_t>(probe * probe)) & (cap - 1);
        const Entry  *e = get_entry(i);

        if (e->hash == 0)
            return std::nullopt;    // empty slot: key is definitely not present

        if (e->hash == h) {
            // hash matches: verify with Euclidean distance (no sqrt needed)
            const double *ep = e->params();
            double dist_sq = 0.0;
            for (size_t k = 0; k < n_params; k++) {
                double d = ep[k] - p[k];
                dist_sq += d * d;
            }
            if (dist_sq < EPS_SQ) {
                n_hits++;
                return e->score;
            }
        }
        // hash collision or verify failed: continue probing
    }

    return std::nullopt;
}

// ---------------------------------------------------------------------------
// store
// ---------------------------------------------------------------------------

void EvalCache::store(const double *p, double score)
{
    if (std::isnan(score))
        return;     // never cache invalid results

    const uint64_t h   = compute_hash(p);
    const size_t   idx = h & (cap - 1);

    for (int probe = 0; probe < MAX_PROBE; probe++) {
        const size_t i = (idx + static_cast<size_t>(probe * probe)) & (cap - 1);
        Entry       *e = get_entry(i);

        if (e->hash == 0) {
            // empty slot: claim it
            e->hash  = h;
            e->score = score;
            std::memcpy(e->params(), p, n_params * sizeof(double));
            n_entries++;
            return;
        }

        if (e->hash == h) {
            // check if this is the same point (update in place)
            const double *ep = e->params();
            double dist_sq = 0.0;
            for (size_t k = 0; k < n_params; k++) {
                double d = ep[k] - p[k];
                dist_sq += d * d;
            }
            if (dist_sq < EPS_SQ) {
                e->score = score;   // update existing entry
                return;
            }
        }
        // slot occupied by a different key: keep probing
    }

    // All probe slots occupied: silently drop.
    // The cache is advisory so this is always correct.
}

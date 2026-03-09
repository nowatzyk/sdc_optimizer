#pragma once

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <optional>

//
// EvalCache -- advisory within-run deduplication cache for optimizer oracle calls.
//
// Stores (parameter_vector -> score) mappings to avoid re-evaluating the oracle
// at points already visited.  The cache is owned by the Optimizer instance,
// constructed once with a fixed capacity, and never persisted across runs.
//
// Layout of each table entry (all fields 8 bytes, naturally aligned, no padding):
//
//   offset  0 : uint64_t hash    -- 0 = empty slot, nonzero = occupied
//   offset  8 : double   score   -- fixed location, always accessible after hash match
//   offset 16 : double   params[n_params]  -- variable length, set at construction
//
// The hash serves as both a fast pre-check and an occupied/empty flag (sentinel trick).
// Full Euclidean verification is only performed on a hash match.
//
class EvalCache {
public:
    //
    // n_params  : number of optimizable parameters (known at construction time)
    // capacity  : number of slots; must be a power of 2.  Default 4096.
    //
    EvalCache(size_t n_params, size_t capacity = 4096);
    ~EvalCache();

    // Non-copyable: owns a raw allocation
    EvalCache(const EvalCache &)            = delete;
    EvalCache &operator=(const EvalCache &) = delete;

    //
    // lookup() -- query the cache.
    //
    // p        : pointer to n_params normalised [0,1] parameter values
    // returns  : the cached score if found, std::nullopt otherwise
    //
    std::optional<double> lookup(const double *p) const;

    //
    // store() -- insert a result.
    //
    // NaN scores are silently ignored (no point caching an invalid result).
    // If all probe slots are occupied (rare at sane load factors) the entry
    // is silently dropped -- the cache is advisory so this is always correct.
    //
    void store(const double *p, double score);

    size_t size()     const { return n_entries; }
    size_t hits()     const { return n_hits;    }
    size_t capacity() const { return cap;       }

private:
    // Entry layout descriptor
    struct Entry {
        uint64_t hash;   // 0 = empty
        double   score;  // fixed at offset 8
        // double params[n_params] follow immediately at offset 16
        double *params() {
            return reinterpret_cast<double *>(this + 1);
        }
        const double *params() const {
            return reinterpret_cast<const double *>(this + 1);
        }
    };

    static constexpr double   EPS      = 1e-6;
    static constexpr double   EPS_SQ   = EPS * EPS;
    static constexpr int      MAX_PROBE = 3;

    size_t   n_params;   // number of doubles per parameter vector
    size_t   cap;        // capacity in slots (power of 2)
    size_t   entry_sz;   // sizeof(Entry) + n_params * sizeof(double)
    size_t   n_entries;
    mutable size_t n_hits;

    char    *table;      // flat allocation: cap * entry_sz bytes, zero-initialised

    Entry *get_entry(size_t i) const {
        return reinterpret_cast<Entry *>(table + i * entry_sz);
    }

    uint64_t compute_hash(const double *p) const;
};

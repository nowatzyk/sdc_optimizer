#pragma once

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <optional>
#include <stdio.h>

//
// EvalCache -- advisory within-run deduplication cache for optimizer oracle calls.
//
// Stores (parameter_vector -> score) mappings to avoid re-evaluating the oracle
// at points already visited.  The cache is owned by the Optimizer instance,
// constructed once with a fixed capacity, and never persisted across runs.
//
// Layout of each table entry (all fields 8 bytes, naturally aligned, no padding):
//
//   offset  0 : double   score   -- fixed location, always accessible after hash match, NAN means empty
//   offset  8 : double   params[n_params]  -- variable length, set at construction
//
// The hash serves as both a fast pre-check and an occupied/empty flag (sentinel trick).
//
// The cache currently uses a locality sensitive hash (LSH) that uses the projection methods to
// clusted entities high cos-similarity together. The hash is then used in a direct mapped cache
// with quadartic probing. It is quite likely that this is not optimimal because LSH custers similar 
// entities together and an optimizer tends to proble near the optimum more frequently, thus it can be
// expected that the hash collision rate is high. Probing is limited to MAX_PROBE cycles and an
// overflow event is counted when this occurs. If the cache statistic show low cache utilization and 
// many overflow events, then it is time to use a better cache structure.  2-choice hashing will
// be next: have two different LSH's, two direct mapped caches as it is now, but insert in the cache
// where quadraic probing shows the most vacancies. If both have the same number of vacancies, always
// insert in the first. Deterministic choice is important here.
//

#define DEFAULT_EPS 1.0e-3          // This is the default Euclidian distance between two
                                    // parameter vectors to be consider the same
                                    // This is very conservative and can be overriden in the
                                    // EvalCache constructor.

class EvalCache {
public:
    //
    // n_params  : number of optimizable parameters (known at construction time)
    // capacity  : number of slots; will always be a power of 2
    //
    EvalCache(size_t n_params, size_t capacity = 4096, double eps = DEFAULT_EPS);
    ~EvalCache();

    // Non-copyable: owns a raw allocation
    EvalCache(const EvalCache &)            = delete;
    EvalCache &operator=(const EvalCache &) = delete;

    //
    // lookup() -- query the cache.
    //
    // p        : pointer to n_params normalised [0,1] parameter values
    // returns  : the cached score if found, NAN otherwise
    //
    double lookup(const double *p) const;   // Returns NAN on a miss

    //
    // store() -- insert a result.
    //
    // NaN scores are silently ignored (no point caching an invalid result).
    // If all probe slots are occupied an overflow event is counted and the entry
    // is silently dropped -- the cache is advisory so this is always correct.
    //
    void store(const double *p, double score);

    size_t size()     const { return n_entries; };
    size_t hits()     const { return n_hits;    };
    size_t capacity() const { return cap;       };
    
    void print_stats(FILE *fp); // Prints a the statistics

private:
    // Entry layout descriptor
    struct Entry {
        double   score;         // fixed at offset 0
        // double params[n_params] follow immediately at offset 8
        double *params() {
            return reinterpret_cast<double *>(this + 1);
        }
        const double *params() const {
            return reinterpret_cast<const double *>(this + 1);
        }
    };

    size_t   n_params;          // number of doubles per parameter vector
    size_t   cap;               // capacity in slots (= 2^log2_cap)
    size_t   log2_cap;          // = log2 of capacity (= number of bits in hashed address)
    size_t   entry_sz;          // sizeof(Entry) + n_params * sizeof(double)
    
    // Cache statistics
    size_t   n_entries;         // = how many entries are in use
    mutable size_t   n_hits;    // = how many hits occured (mutable because it is changed in const function)
    size_t   n_overflows;       // = how many store requestes were ignored due to bin-overflows
    
    double   eps_squared;       // cache hit match tolerance ^2

    double  *hyper;             // Hyperplanes for the LSH hash function
    double   vln_scale;         // Vector length normalization scale factor:
                                // See explanation in compute_hash()

    char    *table;             // flat allocation: cap * entry_sz bytes, NAN score initialized

    Entry *get_entry(size_t i) const {
        return reinterpret_cast<Entry *>(table + i * entry_sz);
    }

    unsigned compute_hash(const double *p) const;
};

void test_eval_cache(unsigned n_dim, unsigned n_tests, unsigned n_matches, double eps, double eps_fraction);

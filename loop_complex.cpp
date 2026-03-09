#include "loop_complex.h"
#include "csv_analyzer.h"       // for lsq_fit_function
#include "parameter.h"          // for update_assignments()
#include "lsq_fit_function.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Global instance
//
LoopComplex loop_complex;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// LoopComplex -- construction
//

LoopComplex::LoopComplex()
    : eval_expr(nullptr)
    , run_josim_cb(nullptr)
    , post_run_cb(nullptr)
    , current_level(0)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Registration
//

void LoopComplex::register_iterator(LoopIterator *it)
{
    assert(it);
    it->level = current_level++;    // post-increment: iterator owns this level
                                    // subsequently registered stateful objects get higher level
    iterators.push_back(it);
}

void LoopComplex::register_stateful(StatefulExpression *s)
{
    assert(s);
    s->level = current_level;       // level at time of registration = innermost active loop
    stateful.push_back(s);
}

void LoopComplex::register_lsq_fit(lsq_fit_function *f)
{
    assert(f);
    f->level = current_level;
    lsq_fits.push_back(f);
}

unsigned LoopComplex::n_optimizable_params() const
{
    // Delegate to parameter class for now; will move here during parameter strip-down
    unsigned n = 0;
    for (auto *it : iterators)
        (void)it;   // placeholder -- iterators don't expose optimizable flag yet
    // For SA/BO, the count comes from parameter::p_ptr.size()
    // This will be revisited when parameter is refactored
    return (unsigned) iterators.size();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Lifecycle helpers
//

void LoopComplex::initialize_level(unsigned lv)
{
    for (auto *s : stateful)
        if (s->level == lv) s->initialize();
    for (auto *f : lsq_fits)
        if (f->level == lv) f->initialize();
}

void LoopComplex::finalize_level(unsigned lv)
{
    // StatefulExpressions always before LsqFits: fits may depend on finalized expression values
    for (auto *s : stateful)
        if (s->level == lv) s->finalize();
    for (auto *f : lsq_fits)
        if (f->level == lv) f->finalize();
}

void LoopComplex::update_all()
{
    for (auto *s : stateful)
        s->update();
    for (auto *f : lsq_fits)
        f->add_datum();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// run_once() -- the evaluation oracle
//
// Executes the full nested loop, running JoSIM at each innermost point.
// Returns eval_expr->get_value() after the outermost loop completes,
// or NaN if no eval_expr is defined.
//
// Level convention: level 0 = innermost, level N-1 = outermost.
// iterators[] is in registration order, outermost first (highest index = level 0).
//
// Finalization order at each level:
//   StatefulExpressions finalize before LsqFits (separate vectors ensure this).
//

double LoopComplex::run_once()
{
    assert(run_josim_cb != nullptr);

    const unsigned n_levels = (unsigned) iterators.size();

    // --- Single pass if no iterators (plain expression evaluation) ---
    if (n_levels == 0) {
        for (auto *s : stateful)  s->initialize();
        for (auto *f : lsq_fits)  f->initialize();

// TBD        parameter::update_assignments();
        run_josim_cb();
        update_all();

        for (auto *s : stateful)  s->finalize();
        for (auto *f : lsq_fits)  f->finalize();

        if (post_run_cb) post_run_cb();

        return eval_expr ? eval_expr->get_value() : NAN;
    }

    // --- Initialize all levels outermost -> innermost ---
    for (unsigned lv = 0; lv < n_levels; lv++)
        iterators[lv]->initialize();
    for (unsigned lv = 0; lv < n_levels; lv++)
        initialize_level(lv);

    // --- Inner simulation loop ---
    for (;;) {
// TBD        parameter::update_assignments();
        run_josim_cb();
        update_all();

        if (post_run_cb) post_run_cb();

        // --- Advance levels innermost (index 0) -> outermost ---
        bool advanced = false;
        for (unsigned lv = 0; lv < n_levels; lv++) {
            if (iterators[lv]->next()) {
                // This level has more steps: restart all inner levels
                for (unsigned inner = 0; inner < lv; inner++) {
                    finalize_level(inner);
                    iterators[inner]->initialize();
                    initialize_level(inner);
                }
                advanced = true;
                break;
            }

            // Level lv exhausted: finalize it
            finalize_level(lv);

            if (lv == n_levels - 1) {
                // Outermost level done: entire loop complete
                return eval_expr ? eval_expr->get_value() : NAN;
            }
            // Continue to try advancing the next outer level
        }

        if (!advanced) {
            // No iterators advanced -- shouldn't happen with n_levels > 0
            // but guard against it
            break;
        }
    }

    return eval_expr ? eval_expr->get_value() : NAN;
}

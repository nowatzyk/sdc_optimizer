#include <string.h>
#include "lsq_fit_function.h"
#include "loop_complex.h"

extern "C" {
#include "fit_functions.h"
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Static registry (replaces root/next linked list)
//
vector<lsq_fit_function *> lsq_fit_function::registry;

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Construction
//

lsq_fit_function::lsq_fit_function(const char *nm, int o,
                                    expression *fv, expression *dv,
                                    LoopComplex &lc)
    : level(0)
    , name(nm)
    , free_var(fv)
    , dep_var(dv)
    , order((unsigned) o)
{
    assert(o >= 1 && o <= 4);

    lsq_fit_ptr = new_lsq_fit(1, order + 1, polynomial_o4);
    clear_results();
    init_lsq_fit(lsq_fit_ptr);

    // Self-register with LoopComplex -- level assigned here
    lc.register_lsq_fit(this);

    // Add to find() registry
    registry.push_back(this);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Lifecycle
//

void lsq_fit_function::initialize()
{
    clear_results();
    init_lsq_fit(lsq_fit_ptr);
}

void lsq_fit_function::add_datum()
{
    double x = free_var->get_value();
    if (!isfinite(x))
        return;     // NaN is a controlled way to exclude a data point (like 1/0 in gnuplot)

    double y = dep_var->get_value();
    if (!isfinite(y))
        return;

    add_lsq_fit(lsq_fit_ptr, &x, y);
}

void lsq_fit_function::finalize()
{
    int ie = solve_lsq_fit(lsq_fit_ptr);

    if (ie == 0) {
        residual = rsquare_lsq_fit(lsq_fit_ptr);
        for (unsigned i = 0; i <= order; i++)
            coefficients[i] = coeff_lsq_fit(lsq_fit_ptr, i);
    } else {
        clear_results();    // no solution: leave everything as NAN
    }

    init_lsq_fit(lsq_fit_ptr);     // ready for next cycle
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Helpers
//

void lsq_fit_function::clear_results()
{
    residual = NAN;
    for (unsigned i = 0; i < 5; i++)
        coefficients[i] = NAN;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// find() -- parser lookup by name
//

lsq_fit_function *lsq_fit_function::find(const char *nm)
{
    for (auto *f : registry)
        if (!strcmp(nm, f->name))
            return f;
    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Functions to be used in expressions
//


double test_lsq_fit(void *obj_ptr, double x)
{
    return ((lsq_fit_function *) obj_ptr)->fit_ok();
}

double lsq_fit_get_cx(void *obj_ptr, double x)
{
    int ci = (int) nearbyint(x);
    assert((0 <= ci) && (ci < 5));
    return ((lsq_fit_function *) obj_ptr)->get_ci(ci);
}

double lsq_fit_residual(void *obj_ptr, double x)
{
    return ((lsq_fit_function *) obj_ptr)->get_residual();
}

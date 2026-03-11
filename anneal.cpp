//
// Optimization via simulated annealing
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "loop_complex.h"
#include "anneal.h"


#define REJ_MAX     20                      // #of rejects before reducing change scale
                                            //   in simulated annealing

#define MAX_RESET   1000                    // patience factor: if opt wasn't improved in this

const double MIN_RANGE = 1.0e-6;            // Min. parameter change range
extern const char*  EVAL_PARAMETER;         // The evaluation function parameter


#define MAX_CHANGE_TRIES 50                 // #of times a parameter change can be rejected

#define EPS         1.0e-10                 // Greedy search cut-off

sim_anneal *sim_anneal_ptr = nullptr;       // The sim annealing instance (if used)

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The simulated annealer
//

sim_anneal::sim_anneal(char* log_file_name, parameter *ev_ptr, parameter *rj_ptr) :
    eval_ptr(ev_ptr), reject_ptr(rj_ptr), changed_param(nullptr), changed_value(0.0),
    best_eval(NAN), have_best(0), changed_flag(0), sum_fp(nullptr)
{
    lfp = fopen(log_file_name, "w");
    assert(lfp);
}

void sim_anneal::add_sa_sched(double t, unsigned n)
{
    assert((t >= 0.0) && (n > 0));
    an_sched ans;
    ans.temp = t;
    ans.n_steps = n;
    sa_schedule.push_back(ans);
}

void sim_anneal::add_parameter(const_parameter* p_ptr)
{
    sa_p_ptrs.push_back(new sa_parameter(p_ptr));
}

void sim_anneal::save_best(double best_ev)
{
    best_eval = best_ev;
    have_best = 1;
    for (auto *pp : sa_p_ptrs ) pp->save_best();
}

double sim_anneal::restore_best()
{
    assert(have_best == 1);
    for (auto *pp : sa_p_ptrs ) pp->restore_best();
    return best_eval;
}

unsigned long sim_anneal::n_steps()
{
    unsigned long tot = 0;

    for (auto sa : sa_schedule)
        tot += sa.n_steps;

    return tot;
}

void sim_anneal::save_result(FILE* fp)
{
    for (auto *pp : sa_p_ptrs ) pp->p_ptr->print_self(fp);
}

double sim_anneal::comp_eval()
//
// The actual simulation happens here
//
{
    loop_complex.run_once(sum_fp);
    return eval_ptr->get_cur_value();
}


void sim_anneal::change(double c_range)
      //
      // This function changes a parameter
      //
      // The extend of a change is controlled by 0 < <c_range> <= 0.5
      //
      // A range of 0.5 means that a variable may be changed
      // by half of the possible range for that variable.
      //
{
#ifdef _ADAPTIVE_RANGEING_
    if (changed_flag == 1) {           // change followed by change => last change was accepted
        changed_param->n_accept += 1;
        changed_param->acc_pos = changed_value < (changed_param->p_ptr->get_cur_value());
        changed_param->n_reject = 0;
        if (changed_param->n_accept > REJ_MAX) {
            changed_param->range *= 1.25;
            if (changed_param->range > 0.5)
                changed_param->range = 0.5;
            else if (lfp)
                fprintf(lfp, "# range of %s increased to %.4g\n", changed_param->p_ptr->get_name(),
                        changed_param->range);
            changed_param->n_accept = 0;
        }
    }
#endif

    for (int n_try = 0; n_try < 20; n_try++) {
        int i = rnd_ri(n_tune());           // pick a parameter to change

        changed_param = sa_p_ptrs[i];       // Remember
        changed_value = changed_param->p_ptr->get_mapped_value();   // Save old value
        changed_flag = 1;

        // Determine valid range for change
#ifdef _ADAPTIVE_RANGEING_
        double d = changed_param->range;
#else
        double d = c_range;
#endif
        double d_min = fmax(0.0, changed_value - d);
        double d_max = fmin(1.0, changed_value + d);

        for (int j = 0; j < MAX_CHANGE_TRIES; j++) {
            d = rnd_01d();                  // d=random number in [0,1)

            double v_new = d_min + d * (d_max - d_min);
            assert(v_new >= 0.0 && v_new <= 1.0);
            changed_param->p_ptr->set_mapped_value(v_new); // Perform the change

            if ((reject_ptr == nullptr) || (reject_ptr->get_cur_value() == 0.0)) {
                                            // No reject function or change accepetd
                return;
            }
        }

        changed_param->p_ptr->set_mapped_value(changed_value); // Unchange this one, before trying another!
        changed_param->n_reject += 1;       // This was essentially a reject!
        fprintf(stderr, "parameter::change_param - Change rejected too often, trying something else\n");
    }

    fprintf(stderr, "Excessive rejections - Aborting\n");
    assert(0);
    exit(1);
}


void sim_anneal::unchange()
//
// Un-does the effect of the last change
//
// Note: only that last change can be un-changed
//
{
    assert(changed_flag);                   // There had to be a not-unchanged change
    changed_flag = 0;
    
    #ifdef _ADAPTIVE_RANGEING_
    if (changed_param->acc_pos == (changed_value < (changed_param->p_ptr->get_mapped_value())))
        changed_param->n_accept = 0;        // Over-shoot: evidence that we should not increase the range!
    changed_param->n_reject += 1;
    if (changed_param->n_reject > REJ_MAX) {
        changed_param->range *= 0.8;
        if (changed_param->range < MIN_RANGE)
            changed_param->range = MIN_RANGE;
        else if (lfp)
            fprintf(lfp, "# range of %s reduced to %.4g\n",
                    changed_param->p_ptr->get_name(), changed_param->range);
        changed_param->n_reject = 0;
    }
    #endif
    
    changed_param->p_ptr->set_mapped_value(changed_value);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The core of the simulated anealing function:
//
// This follows the usual algorithm, but this the following changes:
// - It will never forget the best value encountereed during the search and will
//   return this solution.
// - The starting point may be the best value: this means that this code will never
//   resturn a solution worse than the initial starting point.
// - There is a rat-hole cut-off: if no progress has been made in a while, the solution
//   is rest to the preveiously best known solution.
// - The range of parameter changes is decreased as the search progresses. If the
//   parameters have widely different ranges, it is better to used adaptive parameter changes.
//

void sim_anneal::optimize()
     //
     // Use simulated annealing to optimize a parameter set
     //
{
    double t0, t1, dt;                      // annealing temperature
    double old_eval, new_eval, min;
    int i, j, k, nc, n_rej, t_max;
    double d_scale;                         // scales the range of changes
    double d, t, sigm;

    int out_cnt = 0;                        // Used to control the volume of the log-output

    assert(n_tune() > 0);                   // Make sure we are not wasting effort

    min = old_eval = comp_eval();           // eval initial parameter set

    save_best(min);                         // initial value set is best so far.

    t_max = MAX_RESET;
    n_rej = 0;
    d_scale = 0.25;
    nc = 0;

    for (i = 0; i < sa_schedule.size(); i++) {  // for each part of the annealing schedule

        //
        // t0 = temp at start of schedule interval, t1 = temp at end of interval
        //
        t0 = sa_schedule[i].temp;
        if ((i + 1) < sa_schedule.size()) t1 = sa_schedule[i + 1].temp;
        else                              t1 = 0.0;

        // Within one annealing step decrease temperature linearily
        dt = (t1 - t0) / (double) (sa_schedule[i].n_steps - 1);
        for (j = 0; j < sa_schedule[i].n_steps; j++, t0 += dt) {

            //
            // Try out a new parameter set
            //
            change(d_scale);

            //
            // Evaluate the change
            //
            new_eval = comp_eval();

            //
            // Check score
            //
            if (new_eval < min) {           // new best score
                min = new_eval;             // keep track of max-evaluation
                save_best(new_eval);        // .. and how it was achived
                t_max = MAX_RESET;          // reset the reset counter

            } else if (0 >= t_max--) {      // reset if we went on a wrong path for too long
                old_eval = new_eval = min;  // restore previous best
                restore_best();
                t_max = MAX_RESET;          // reset the reset counter
                continue;
            }

            t = (old_eval - new_eval) / old_eval;  // t: relative change (t>0 is better)
            nc++;

            //
            // See if we accept this change
            //
            if (t0 > EPS) {
                sigm = 1.0 / (1.0 + exp(-t/t0));// sigmoid function of t*temp
                d = rnd_01d();              // d= (-0, 1]
                k = (sigm >= d);
            } else
                k = (t > 0.0);              // once temp is too low, just be greedy

            if (k) {                        // Accept change
                old_eval = new_eval;
                n_rej = 0;
                if (lfp) {
                    out_cnt++;
                    if ( (out_cnt < 10000   )                       ||
                        ((out_cnt < 100000  ) && !(out_cnt % 10  )) ||
                        ((out_cnt < 10000000) && !(out_cnt % 100 )) ||
                                                 !(out_cnt % 1000)    )
                        fprintf(lfp, "%d %.5e\n", nc, new_eval);
                }
            } else {                        // Reject: undo change
                unchange();
                n_rej++;
            }

            //
            // too many rejects: scale back change range
            //
#ifndef _ADAPTIVE_RANGEING_
            if ((n_rej > REJ_MAX) && (d_scale > 1.0e-5)) {
                n_rej = 0;
                d_scale *= 0.8;
                if (d_scale < MIN_RANGE)
                    d_scale = MIN_RANGE;
                else if (of)
                    fprintf(of, "# range reduced to %.4g\n", d_scale);
            }
#endif
        }
    }

    restore_best();                         // Return new configuration
}    

//////////////////////////////////////////////////////////////////////////////////////
//
// Parameter implementation
//

sa_parameter::sa_parameter(const_parameter *cp_ptr) :
    p_ptr(cp_ptr), best_val(NAN)
    //
    // Parameter constructor:
    //
    // Changes from the original parameter code of the anneal.* package:
    //
    // 1. The following items are part of the (contant-)parameter class, so no longer here:
    // <name> : ascii string to identify this parameter
    // <ptr>  : pointer to the actual value
    // <min/max> : allowed tuning range
    //
    // 2. The reject function is now a one per SA fit instance and no longer a one 
    //    per parameter entity. This is more course, but unlikely to be ver useful.
    //    also, it no longer can have side-effects, which is probably a good thing here.
    // <reject> : reject function (optional), may have side-effects
    //
    // Note on the use of reject functions: The primary intent of this
    // facility is to provide a mechanism to express constraints that include
    // multiple parameters. For example if <x> and <y> where two independent
    // caratesian coordinates, it may be desirable to limit the value of
    // r = sqrt(x^2 + y^2).
    //
{
#ifdef _ADAPTIVE_RANGEING_
    n_accept = 0;                           // Reset counters
    acc_pos = -1;
    n_reject = 0;
    range = 0.25;                           // .. and initialized default range
#endif
}

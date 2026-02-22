//
//  This file contails the code for the parameter class.
//
//  This class is rather complex because it does many things that probably should not be mixed together:
//
//  1. It acts as the symbol table for entiies that the use in the modification of a JoSIM (spice) circuit
//     These entities have several functions:
//      a) just a simple variable that can be used in a substitution
//      b) as a variable that is used to collect some infor from the simulation and that is then used
//         in the report and plotted via gnuplot or some other post-processin utility
//      c) as accumulator for data extracted from multiple simulation runs
//      d) as the result of an expression using other variable. This can be used to parameterize the JoSIM
//         circuit is a more general way, because the papameter mechaism in JoSIM (.param) is rather
//         limited: for example it is not possible to construct subcircuits where several component
//         values depend on one papa,eter in different ways. A JJ + dumping resistor is one such example.
//  2. It controls looping, where the parameter value it iterated over multiple values. For example
//     in the scan-type parameter the value is iterated from a start to en end value in a certain number of steps.
//     Other example of looing constructs is the binary search  function.
//  3. It does optimization via simulated annealing
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <vector>
using namespace std;

#include "xrand.h"
#include "csv_analyzer.h"
#include "parameter.h"

extern "C" {
#include "parser.h"                         // Bison generated headers (in build dir)
#include "lex.yy.h"                         // Flex-generated header (in build dir)
#include "lsq_fit.h"
#include "fit_functions.h"
#include "parser_interf.h"                  // needed to integrate the parser
}

///////////////////////////////////////////////////////////////////////////////////////////////////

double parameter::get_cur_value()
{
    double rtn = NAN;
    
    switch (type) {
        case scan:
            rtn = (min_value + (double) step_cntr * d_value);
            break;
            
        case assignment:
            rtn =  expr->get_value();
            break;
            
        case sim_anneal:
            rtn = cur_value;
            break;
            
        default:
            assert (0);
    }
    
    return rtn;
}


int parameter::i_next()
{
    switch (type) {
        case scan:
            step_cntr++;
    
            if (step_cntr >= n_steps) {
                step_cntr = 0;
                return 1;
            } 
    
            return 0;
            break;
            
        case sim_anneal:
            return sim_anneal_cycle();
            
        default: assert(0);
    }
}

void parameter::print_name(FILE* fp)
{
    fprintf(fp, "%s", name);
}

void parameter::print_value(FILE* fp)
{
    fprintf(fp, " %.15lg", get_cur_value());
}

void parameter::reset()
{
    parameter *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->i_reset();
        p_ptr = p_ptr->next;
    }
}

unsigned parameter::list_names(FILE* fp)
{
    parameter *p_ptr = root;
    unsigned cnt = 1;
    while(p_ptr != nullptr) {
        fprintf(fp, " (%u):", cnt++);
        p_ptr->print_name(fp);
        p_ptr = p_ptr->next;
    }
    
    return cnt;
}

void parameter::list_c_val(FILE* fp)
{
    parameter *p_ptr = root;
    while(p_ptr != nullptr) {
        p_ptr->print_value(fp);
        p_ptr = p_ptr->next;
    }
}

void parameter::default_init(char* nm)
// Default initializer
{
    type = undefined;
    name = nm;
    next = root;
    root = this;
    
    min_value = - __DBL_MAX__;
    max_value = __DBL_MAX__;
    n_steps = 0;
    d_value = 0.0;
    step_cntr = 0;
    
    expr = nullptr;
    
    best_val = 0.0;
    cur_value = 0.0;
#ifdef _ADAPTIVE_RANGEING_
    n_accept = 0;                       // Reset counters
    acc_pos = -1;
    n_reject = 0;
    range = 0.25;                       // .. and initialized default range
#endif
}

int parameter::advance(unsigned &lvl)
//
// Advance the parameter values. Returns 1 when all combinations are done, 0 otherwise
//
// <lvl> is set to highest nesting level that was advanced. This is used in the summary report
// so that a plot function can use this info to only use data points that are a function of
// the innermost loop(s) having been completed.
//
// Parameters like scan, binary_search, and simulated annealing cause the JoSIM simulation to
// to be re-run, only with the paramters changing values in some way. The output of the simulation
// is then digested and analyzed is some way, which may change the way paramters are changed.
// This means that there are multiple levels of loops arround the simulation/analysis core.
// A san type loop will simply cycle one parameter through a prefefined series of values. Having a second
// paramter scan, just means that 2 paramters will cyles through their values, essentially an two-
// dimensional nested loop.
//
// Some analysis can depend on a particular loop being completed. For example, one scan type parameter
// has exhausetd all possible values, and some property of the circuit is computed over all these runs.
// The summ-operator is one such example: it sums up expressions that are functions of the outcome
// of each simulation, say to find an average, or a minmum, etc. If this is not the outermost loop,
// then the summation must be reset, once the loop is restarted. So there is some work that needs to
// to be done whenever one loop completes. The <iterator> structure keeps track of that.
//
// Where things get messy, is the order in which this housekeeping work is done. This is critical for
// the way this facility is used, but it is also obscure and not obvious from the annotated/instrumented
// circuit file. Thus this long-winded comment. In a nutshell once a loop is complete, this is done
// in this specific order:
//
// 1. Compute any least-squate fits: this means that expressions refering to the results of a LSQ fit
//    will now produce the correct values.
// 2. save the result of the eval-expression that is controlling the simulated annealing, which is
//    always the outermost, top-level loop. <eval> expressions, may depend on results from the LSQ fit,
//    which means that the above must be done first.
// 3. Zero (reset) any summation expression. This will destroy their current value. The eval-function
//    may depend on summation expressions, therefore, this maust be done last.
//
// This logic is hardwired into this code, but it is somthing of a mess, because this way of doing business
// is not obvios, but it is important for how to structure the pragmas that explore a circuit.
// Sorry, but I don't have a better idea right now.
//
{
    for (unsigned level = 0; level < nesting_level; level++) {

        if (iterator[level].param_ptr->i_next() == 0) {
            lvl = level;
            return 0;
        }
        
        // Level <level> is done, new it is time to compute any LSQ fits for this level
        for (fit_element *f_ptr = iterator[level].fit_ptr; f_ptr != nullptr; f_ptr = f_ptr->next)
            f_ptr->fit_ptr->perform_fit();
        
        if (enable_sim_anneal) {
            // Save the evaluation, which might depend on sums that are about to be zero-ed
            assert(eval_ptr != nullptr);
            new_eval = eval_ptr->get_cur_value();
        }
        
        // Level <level> completed. Need to zero any sums over this level:
        for (zel_element *z_ptr = iterator[level].zel_ptr; z_ptr != nullptr; z_ptr = z_ptr->next)
            z_ptr->e_ptr->zero_sum();
    }

    lvl = nesting_level - 1;
    return 1;
}

parameter::parameter(char* nm, double v_min, double v_max, unsigned n)
// Constructor for a scan-type parameter
{
    default_init(nm);
    type = scan;                        //This is a scan type parameter
    
     min_value = v_min;
     max_value = v_max;
     n_steps = n;
    
    if (enable_sim_anneal == 1) {
        fprintf(stderr, "Line %d: Can't iterate over simulated annealing\n", yylineno);
        //
        // Note: the sim annealing code currently is statically initialized, meaning
        //       it can run only once. So it must be the outermost loop/iterator.
        //       Adding a proper reset/init is simple, but pointless: iteratng over simulated
        //       annealing steps would vastly increase run time. Better to that manually...
        //
        yy_n_parse_err++;
        return;
    }

    if (n <= 1) d_value = 1.0;          // Does not matter, never used (only 1 step)
    else        d_value = (v_max - v_min) / (double) (n - 1);
    
    iterator[nesting_level++].param_ptr = this;  // This parameter is a looping one
}

parameter::parameter(char* nm, double v_min, double v_init, double v_max, parameter* ev_ptr, parameter *rj_ptr)
//
// Enables simulated annealing
//
{
    default_init(nm);
    type = sim_anneal;
    
    if (as.size() <= 0) {
        fprintf(stderr, "Line %d: Attempt to use simulated annealing without a schedule\n", yylineno);
        yy_n_parse_err++;
        return;
    }
    
     min_value = v_min;
     cur_value = v_init;
     max_value = v_max;
     
     assert ((eval_ptr == nullptr) || (eval_ptr == ev_ptr));
     if (eval_ptr == nullptr) {
         //
         // This is the first parameter that sets up the simulated annealing machinery
         //
         eval_ptr = ev_ptr;
         reject_ptr = rj_ptr;
         iterator[nesting_level++].param_ptr = this;  // This parameter is a looping one
         enable_sim_anneal = 1;
     }
     
     p_ptr.push_back(this);             // This parameter is part of the set that is being optimized
}


parameter::parameter(char* nm, class expression* expr_ptr)
// assignment type parameter
{
    default_init(nm);
    type = assignment;
    
    expr = expr_ptr;        // This is the only thing that matters!
}

parameter *parameter::find_parameter(const char* nm)
{
    parameter *pptr = root;
    while (pptr != nullptr) {
        if (!strcmp(nm, pptr->name))
            return pptr;
        pptr = pptr->next;
    }
    return nullptr;
}

void parameter::enqueue_zero_op(expression* expr_ptr)
{
    zel_element *z_ptr = new zel_element;
    z_ptr->e_ptr = expr_ptr;
    z_ptr-> next = iterator[nesting_level].zel_ptr;
    iterator[nesting_level].zel_ptr = z_ptr;
}

void parameter::enqueue_fit_function(lsq_fit_function *lsqf_ptr)
{
    fit_element *f_ptr = new fit_element;
    f_ptr->fit_ptr = lsqf_ptr;
    f_ptr-> next = iterator[nesting_level].fit_ptr;
    iterator[nesting_level].fit_ptr = f_ptr;
}

void parameter::save_best()
{
    for (unsigned i = 0; i < p_ptr.size(); i++) {
        p_ptr[i]->best_val = p_ptr[i]->cur_value;
    }
    have_best = 1;
}

void parameter::restore_best()
{
    assert(have_best == 1);
    for (unsigned i = 0; i < p_ptr.size(); i++) {
        p_ptr[i]->cur_value = p_ptr[i]->best_val;
    }
}

unsigned parameter::sim_anneal_cycle()
{
    if (n_sim_anneal == 0) {
        //
        // First time
        //
        if (!isfinite(new_eval)) {
            fprintf(stderr, "Aborting: sim_anneal needs a valid evaluation on the starting configuration\n");
            exit(1);
        }
        min_eval = new_eval;
        save_best();
        old_eval = new_eval;
    }
    
    if (sched_stp_cntr >= as[cur_sched].n_steps)  {
        // Done with this schedule step, move on to the next
        cur_sched++;
        if (cur_sched >= as.size()) {
            // Done with all schedule steps
            restore_best();
            return 1;
        }
        sched_stp_cntr = 0;
    }
    
    if (sched_stp_cntr == 0) {
        // First time for this annealing schedule step:
        //
        // t0 = temp at start of schedule interval, t1 = temp at end of interval
        //
        t0 = as[cur_sched].temp;
        if ((cur_sched + 1) < (int) as.size()) t1 = as[cur_sched + 1].temp;
        else		                           t1 = 0.0;

        // Within one annealing step decrease temperature linearily
        dt = (t1 - t0) / (double) (as[cur_sched].n_steps - 1);
    }
    
    if (n_sim_anneal > 0) {
        unsigned skip_anneal = 1;               // ... to avoid spagetti code
        
        // A change had been made, now decide if it  kept:

        //
        // Check score
        //
        if (new_eval < min_eval) {              // new best score
            min_eval = new_eval;                // keep track of max-evaluation
            save_best();                        // .. and how it was achived
            t_max = MAX_RESET;                  // reset the reset counter
        } else if (0 >= t_max--) {              // reset if we went on a wrong path for too long
            old_eval = new_eval = min_eval;		// restore previous best
            restore_best();
            t_max = MAX_RESET;                  // reset the reset counter
            skip_anneal = 0;
        }

        if (skip_anneal == 1) {
            unsigned accept_change = isfinite(new_eval); // A change can only be accepted if it resulted into a valid configuration
            double t = (old_eval - new_eval) / old_eval;  // t: relative change (t>0 is better)
            //
            // See if we accept this change
            //
            if (t0 > EPS) {
                double sigm = 1.0 / (1.0 + exp(-t/t0));// sigmoid function of t*temp
                double d = rnd_01d();           // d is in (0, 1]
                accept_change &= (sigm >= d);
            } else
                accept_change &= (t > 0.0);     // once temperature is too low, just be greedy

            if (accept_change) {                // Accept change
                old_eval = new_eval;
                n_rej = 0;
                if (sa_log != nullptr) {
                    out_cnt++;
                    if ( (out_cnt < 10000   )                       ||
                        ((out_cnt < 100000  ) && !(out_cnt % 10  )) ||
                        ((out_cnt < 10000000) && !(out_cnt % 100 )) ||
                                                 !(out_cnt % 1000)    )
                        fprintf(sa_log, "%u %.5e\n", n_sim_anneal, new_eval);
                }
            } else {							// Reject: undo change
                unchange_param();
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
                else if (sa_log)
                    fprintf(sa_log, "# range reduced to %.4g\n", d_scale);
            }
#endif
        }
    }

    change_param(d_scale);  // Make a change to a parameter

    n_sim_anneal++;
    
    sched_stp_cntr++;       // Done with this step
    t0 += dt;
    return 0;
}

unsigned parameter::open_sa_log_file(char *nm)
{
    sa_log = fopen(nm, "w");
    return (sa_log == nullptr);
}

unsigned parameter::read_sa_schedule(char *file_name)
{
	FILE *asf;
    char buf[128];
    an_sched t;     

    asf = fopen(file_name, "r");
    if (!asf) {
		fprintf(stderr, "read_an_schedule: failed to open schedule\n");
		return 1;
	}

    if (!fgets(buf, 120, asf) ||
		strncmp("Anealing Schedule", buf, 17)) {
	    fprintf(stderr, "read_an_schedule: file is not an anealing schedule\n");
	    return 1;
	}

	while(fgets(buf, 120, asf)) {

		if (buf[0] == '#' || buf[0] == '\n' || buf[0] == 0 || buf[0] == '\r')
			continue;							// skip comments or white lines

	    if (2 != sscanf(buf, "%lf%d", &(t.temp), &(t.n_steps)) ||
			t.temp < 0.0 || t.n_steps < 1) {
			fprintf(stderr, "read_an_schedule: bad line '%s' ignored\n", buf);
			continue;
	    }

	    as.push_back(t);
	}

    fclose(asf);
    return 0;
}

void parameter::change_param(double c_range)
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
    if (changed_flag == 1) {                        // change followed by change => last change was accepted
        changed_param->n_accept += 1;
        changed_param->acc_pos = changed_value < (changed_param->cur_value);
        changed_param->n_reject = 0;
        if (changed_param->n_accept > REJ_MAX) {
            changed_param->range *= 1.25;
            if (changed_param->range > 0.5)
                changed_param->range = 0.5;
            else if (sa_log)
                fprintf(sa_log, "# range of %s increased to %.4g\n", changed_param->name, changed_param->range);
            changed_param->n_accept = 0;
        }
    }
#endif

    for (int n_try = 0; n_try < 20; n_try++) {
        int i = rnd_ri(p_ptr.size());               // pick a parameter to change

        changed_param = p_ptr[i];                   // Remember
        changed_value = (changed_param->cur_value); // Save old value
        changed_flag = 1;

        // Determine valid range for change
#ifdef _ADAPTIVE_RANGEING_
        double d = (changed_param->max_value - changed_param->min_value) * changed_param->range;
#else
        double d = (changed_param->max_value - changed_param->min_value) * c_range;
#endif
        double d_min = fmax(changed_param->min_value, changed_value - d);
        double d_max = fmin(changed_param->max_value, changed_value + d);

        for (int j = 0; j < MAX_CHANGE_TRIES; j++) {
            d = rnd_01d();                          // d=random number in [0,1)

            double v_new = d_min + d * (d_max - d_min);
            assert(v_new >= changed_param->min_value && v_new <= changed_param->max_value);
            changed_param->cur_value = v_new;       // Perform the change
            
            if ((reject_ptr == nullptr) || (reject_ptr->get_cur_value() == 0.0))
                // We are done IF there is no reject function OR the reject function evaluates to 0 (signaling a good value)
                return;
        }

        changed_param->cur_value = changed_value;   // Unchange this one, before trying another!
        changed_param->n_reject += 1;               // This was essentially a reject!
        fprintf(stderr, "parameter::change_param - Change rejected too often, trying something else\n");
    }

    fprintf(stderr, "Excessive rejections - Aborting\n");
    assert(0);
    exit(1);
}

void parameter::unchange_param()
      //
      // Un-does the effect of the last change
      //
      // Note: only that last change can be un-changed
      //
{
    assert(changed_flag);                           // There had to be a not-unchanged change
    changed_flag = 0;

#ifdef _ADAPTIVE_RANGEING_
    if (changed_param->acc_pos == (changed_value < (changed_param->cur_value)))
        changed_param->n_accept = 0;                // Over-shoot: evidence that we should not increase the range!
    changed_param->n_reject += 1;
    if (changed_param->n_reject > REJ_MAX) {
        changed_param->range *= 0.8;
        if (changed_param->range < MIN_RANGE)
            changed_param->range = MIN_RANGE;
        else if (sa_log)
            fprintf(sa_log, "# range of %s reduced to %.4g\n", changed_param->name, changed_param->range);
        changed_param->n_reject = 0;
    }
#endif

    changed_param->cur_value = changed_value;
}

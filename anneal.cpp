//
// Optimization via simulated annealing
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "anneal.h"

#define REJ_MAX		20							// #of rejects before reducing change scale
												//   in simulated annealing

#define MAX_RESET	50000						// patience factor: if opt wasn't improved in this
												// many steps, reset to optimum

#define MAX_CHANGE_TRIES 50						// #of times a parameter change can be rejected

#define EPS 1e-10								// Greedy search cut-off

static vector<an_sched> as;						// The anealing schedule

///////////////////////////////////////////////////////////////////////////////////////////////////

///// Sigh, this should be in math...
inline double fmax(double a, double b)
	{return (a >= b) ? a : b; };
inline double fmin(double a, double b)
	{return (a <= b) ? a : b; };

int read_an_schedule(char *file_name)
      //
      // Read an anealing schedule from a file
      //
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

		if (buf[0] == '#' || buf[0] == '\n' || buf[0] == 0)
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

void optimize(FILE *of)
     //
     // Use simulated annealing to optimize a parameter set
     //
{
	double t0, t1, dt;							// annealing temperature
	double old_eval, new_eval, min;
	int i, j, k, nc, n_rej, t_max;
	double d_scale;								// scales the range of changes
	double d, t, sigm;

	int out_cnt = 0;							// Used to control the volume of the log-output

	assert(parameter::n_tune() > 0);			// Make sure we are not wasting effort

	min = old_eval = comp_eval();				// eval initial parameter set

    parameter::save_best();						// initial value set is best so far.

    t_max = MAX_RESET;
	n_rej = 0;
	d_scale = 0.25;
	nc = 0;

	for (i = 0; i < (int) as.size(); i++) {		// for each part of the annealing schedule

	    //
	    // t0 = temp at start of schedule interval, t1 = temp at end of interval
	    //
	    t0 = as[i].temp;
	    if ((i + 1) < (int) as.size()) t1 = as[i + 1].temp;
	    else		                   t1 = 0.0;

	    // Within one annealing step decrease temperature linearily
	    dt = (t1 - t0) / (double) (as[i].n_steps - 1);
	    for (j = 0; j < as[i].n_steps; j++, t0 += dt) {

			//
			// Try out a new parameter set
			//
			parameter::change_param(d_scale, of);

			//
			// Evaluate the change
			//
			new_eval = comp_eval();

			//
			// Check score
			//
			if (new_eval < min) {				// new best score
				min = new_eval;					// keep track of max-evaluation
				parameter::save_best();			// .. and how it was achived
				t_max = MAX_RESET;				// reset the reset counter

			} else if (0 >= t_max--) {			// reset if we went on a wrong path for too long
				old_eval = new_eval = min;		// restore previous best
				parameter::restore_best();
				t_max = MAX_RESET;				// reset the reset counter
				continue;
			}

			t = (old_eval - new_eval) / old_eval;  // t: relative change (t>0 is better)
			nc++;

			//
			// See if we accept this change
			//
			if (t0 > EPS) {
				sigm = 1.0 / (1.0 + exp(-t/t0));// sigmoid function of t*temp
				d = rnd_01d();					// d= (-0, 1]
				k = (sigm >= d);
			} else
				k = (t > 0.0);					// once temp is too low, just be greedy

			if (k) {							// Accept change
				old_eval = new_eval;
				n_rej = 0;
				if (of) {
					out_cnt++;
					if ( (out_cnt < 10000   )                       ||
						((out_cnt < 100000  ) && !(out_cnt % 10  )) ||
						((out_cnt < 10000000) && !(out_cnt % 100 )) ||
						                         !(out_cnt % 1000)    )
						fprintf(of, "%d %.5e\n", nc, new_eval);
				}
			} else {							// Reject: undo change
				parameter::unchange_param(of);
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

    parameter::restore_best();					// Return new configuration
}	

//////////////////////////////////////////////////////////////////////////////////////
//
// Parameter implementation
//

parameter*			parameter::root = 0;		// Root of parameter list
int					parameter::have_best = 0;	// Flag defines the validity of best value
vector<parameter*>	parameter::p_ptr;			// Vector of active parameters
int					parameter::disable_checks = 0; // Check bypass flag

int					parameter::changed_flag = 0; // Private stuff between change and unchange
parameter*			parameter::changed_param = 0;
double				parameter::changed_value;

void				(*parameter::restore)() = 0;// optional restore function

parameter::parameter(char *_name, double *ptr, double _min, double _max, int (*_reject)())
	//
	// Parameter constructor:
	//
	// <name> : ascii string to identify this parameter
	// <ptr>  : pointer to the actual value
	// <min/max> : allowed tuning range
	// <reject> : reject function (optional), may have side-effects
	//
	// Note on the use of reject functions: The primary intent of this
	// facility is to provide a mechanism to express constraints that include
	// multiple parameters. For example if <x> and <y> where two independent
	// caratesian coordinates, it may be desirable to limit the value of
	// r = sqrt(x^2 + y^2).
	// Another, more dubious use of this facility is to compute dependent
	// parameters: these are entities that cannot be changed directly, rather
	// they depend on one or more paramters that are tunable. So the reject
	// function may update the derived values whenever a parameter is changed.
	//
	// Warning: if the reject function has side-effects (like the updating
	// on derived variables) it is important that:
	//  1.) *all* side-effects are included in *all* update functions. For example
	//      if some dependent variable <a> depends on parameter <x> and <y>
	//      then <a> must be recomputed in the update functions for both
	//      <x> and <y>. (otherwise <a> will end up being inconsistent!)
	//  2.) Side-effects must not depend on the order in which parameters are
	//      changed.
	//  3.) Illegal configuration during parameter restore: When a set of
	//      parameters is restored from a file, this is done in an arbitrary
	//      order during which the partially updated parameter set may violate
	//      some constraints. To aviod this kind of problems, the reject function
	//      is called only after all parameters have been restored. This is
	//      still problematic if there are dependency chains: say variable
	//      <a> depends on <b> which depends parameters <x> and <y>, but also on
	//      <z>. So the reject function for <z> may only recompute <a>, but
	//      the update functions for <x> and <y> must recompute <b> and <a>
	//      in that order. While this is kind of obvious, this facility has
	//      plenty of rope to hang yourself.
	//
{
	assert((_min + EPS) < (_max - EPS));		// Must have some range

	if (find(_name)) {
		assert(0);
		fprintf(stderr, "Parameters names must be unique\n");
		exit(1);
	}

	next = root;								// Add to list of parameters
	root = this;

	val = ptr;									// Copy info
	min = _min;
	max = _max;
	name = strdup(_name);
	reject = _reject;

	is_active = 0;								// Not being tuned by default

#ifdef _ADAPTIVE_RANGEING_
	n_accept = 0;								// Reset counters
	acc_pos = -1;
	n_reject = 0;
	range = 0.25;								// .. and initialized default range
#endif
}

parameter* parameter::find(char *s)
	//
	// Find a named parameter
	//
{
	for (parameter *p = root; p; p = p->next)
		if (!strcmp(s, p->name))
			return p;							// Match: found

	return 0;									// Not found
}

void parameter::add()
	//
	// Add this to be tuned
	//
{
	if (is_active)
		return;									// Allready being tuned

	p_ptr.push_back(this);						// Add me to the set
	is_active = 1;
}

void parameter::add_all()
	//
	// Add them all...
	//
{
	for (parameter *p = root; p; p = p->next)
		p->add();
}

void parameter::set_tuneable(int is_t)
	//
	// Change the tunability of a parameter
	//
{
	assert(0 <= is_t && is_t <= 1);

	if (is_active == is_t)
		return;									// Nothing to do

	if (is_t) {									// Add to the list of tunable parameters
		add();
	}

	//
	// Remove from the list of tunable parameters
	//
	is_active = 0;
	for (int i = 0; i < p_ptr.size(); i++)
		if (this == p_ptr[i]) {
			p_ptr.erase(p_ptr.begin() + i);
			return;
		}

	assert(0);									// Inconsistent!
}

void parameter::change_range(double n_min, double n_max)
	//
	// Update range
	//
{
	assert(n_min < n_max);
	min = n_min;
	max = n_max;
}

int add_t_param(char *s, int is_t)
	//
	// Tries to add a tunable parameter
	//
{
	char name[128];
	double value, min, max;
	int nc;

	assert (0 <= is_t && is_t <= 1);			// is tunable flag

	if (2 != sscanf(s, "%s%lf%n", name, &value, &nc))
		return 1;

	parameter *p = parameter::find(name);
	if (!p) {
		fprintf(stderr, "add_t_param: '%s' is not a known parameter\n", name);
		return 1;
	}

	if (nc < strlen(s) && 2 == sscanf(s + nc, "%lf%lf", &min, &max)) {
		// New range supplied:
		if (min + EPS >= max - EPS || value < min || value > max) {
			fprintf(stderr, "add_t_param: '%s' min/max problem\n", name);
			return 1;
		}
		p->change_range(min, max);
	}

	if (p->set(value)) {
		fprintf(stderr, "Rejected %.5g for parameter '%s'\n", value, name);
		return 1;
	}

	p->set_tuneable(is_t);

	return 0;									// Changed successfully
}

int parameter::read(const char *s)
      //
      // Read the set of parameters (and their ranges)
      // that should be tuned
      //
{
	FILE *pf;
    char buf[128];
    int err = 0;

    pf = fopen(s, "r");
    if (!pf)
	    return 1;								// Failed to open file


	if (!fgets(buf, 120, pf) ||
		strncmp("Parameters:", buf, 11)) {
	    fprintf(stderr, "read_parameter_set: file is not a parameter set\n");
	    fclose(pf);
	    return 1;
	}

	disable_checks = 1;							// bypass parameter checking to
												// avoid problems when a new parameter set is
												// loaded incrementally so that
												// inconsistencies may happen during the update.

    while(fgets(buf, 120, pf)) {

	    switch(buf[0]) {
	    case '#':
	    case '\n':
	    case 0:
			break;								// skip comments

	    case 'p':
	    case 'P':								// Add a non-tunable parameter
			err |= add_t_param(buf + 1, 0);
			break;

	    case 't':
	    case 'T':								// Add a tunable parameter
			err |= add_t_param(buf + 1, 1);
			break;

	    default:
			fprintf(stderr, "Bad parameter directive: '%s'\n", buf);
			break;
	    }
	}    

	disable_checks = 0;

	if (restore != 0)
		restore();								// Call restore function if present

	err |= check_all();
	if (!err)
		err |= check_all();						// Do this twice: there may be some order dependencies!

    fclose(pf);

    return err;
}

int parameter::write(const char *s)
      //
      // Write the parameters to a file (text)
      //
{
	FILE *of = fopen(s, "w");

    if (!of)
	    return 1;								// No luck opening file

    fprintf(of, "Parameters:\n");
      
    fprintf(of, "#P name Value Min Max\n");
    for (parameter *p = root; p; p = p->next) {
	    fprintf(of, "%c %s %18.10e %.5g %.5g\n", (p->is_active) ? 'T' : 'P', 
		    p->name, *(p->val), p->min, p->max);
	}

    fclose(of);

    return 0;
}

void parameter::change_param(double c_range, FILE *of)
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
	if (changed_flag == 1) {					// change followed by change => last change was accepted
		changed_param->n_accept += 1;
		changed_param->acc_pos = changed_value < *(changed_param->val);
		changed_param->n_reject = 0;
		if (changed_param->n_accept > REJ_MAX) {
			changed_param->range *= 1.25;
			if (changed_param->range > 0.5)
				changed_param->range = 0.5;
			else if (of)
				fprintf(of, "# range of %s increased to %.4g\n", changed_param->name, changed_param->range);
			changed_param->n_accept = 0;
		}
	}
#endif

	for (int n_try = 0; n_try < 20; n_try++) {
		int i = rnd_ri(n_tune());					// pick a parameter to change

		changed_param = p_ptr[i];					// Remember
		changed_value = *(changed_param->val);		// Save old value
		changed_flag = 1;

		/* Debugging: useful to track down impropper side-effects of reject functions
		if (changed_param->reject) {				// Has reject function:
			int ec = (changed_param->reject)();		// Execute it (before change!)
			if (ec != 0) {
				fprintf(stderr, "Trying to change '%s', but prio state is rejected!\n", changed_param->name);
			}
		}
		//*/

		// Determine valid range for change
#ifdef _ADAPTIVE_RANGEING_
		double d = (changed_param->max - changed_param->min) * changed_param->range;
#else
		double d = (changed_param->max - changed_param->min) * c_range;
#endif
		double d_min = fmax(changed_param->min, changed_value - d);
		double d_max = fmin(changed_param->max, changed_value + d);

		for (int j = 0; j < MAX_CHANGE_TRIES; j++) {
			d = rnd_01d();							// d=random number in [0,1)

			double v_new = d_min + d * (d_max - d_min);
			assert(v_new >= changed_param->min && v_new <= changed_param->max);
			*(changed_param->val) = v_new;			// Perform the change

			if (!(changed_param->reject) || !(changed_param->reject)()) {
													// No reject function or change accepetd

				// added to debug un-change assertion failure
				//for (parameter *p = root; p; p = p->next)
				//	p->last_good = *(p->val);		// Save the last good parameter set

				return;
			}
		}

		*(changed_param->val) = changed_value;		// Unchange this one, before trying another!
		if (changed_param->reject)
			(changed_param->reject)();				// Call the reject function because it may have to undo some set-up side-effects
		changed_param->n_reject += 1;				// This was essentially a reject!
		fprintf(stderr, "parameter::change_param - Change rejected too often, trying something else\n");
	}

	fprintf(stderr, "Excessive rejections - Aborting\n");
	parameter::write("abort_params_dump.txt");
	assert(0);
	exit(1);
}

void parameter::unchange_param(FILE *of)
      //
      // Un-does the effect of the last change
      //
      // Note: only that last change can be un-changed
      //
{
    assert(changed_flag);							// There had to be a not-unchanged change
    changed_flag = 0;

#ifdef _ADAPTIVE_RANGEING_
	if (changed_param->acc_pos == (changed_value < *(changed_param->val)))
		changed_param->n_accept = 0;				// Over-shoot: evidence that we should not increase the range!
	changed_param->n_reject += 1;
	if (changed_param->n_reject > REJ_MAX) {
		changed_param->range *= 0.8;
		if (changed_param->range < MIN_RANGE)
			changed_param->range = MIN_RANGE;
		else if (of)
			fprintf(of, "# range of %s reduced to %.4g\n", changed_param->name, changed_param->range);
		changed_param->n_reject = 0;
	}
#endif

    *(changed_param->val) = changed_value;
	if (changed_param->reject) {
		int i = (changed_param->reject)();			// Undo the side-effects
		if (i) {
			/*
			for (parameter *p = root; p; p = p->next)
				if (*(p->val) != p->last_good)
					printf(">>> unchange miss-match: '%s' val=%.10g last_good=%.10g  diff=%.10g\n",
							p->name, *(p->val), p->last_good, *(p->val) - p->last_good);

			(changed_param->reject)();				// Do it again, so that it can be observed
			//*/
			parameter::write("q_param.txt");
		}
		assert(i == 0);								// That should have worked!
	}
}

void parameter::save_best()
      //
      // Makes a copy of the current parameter set
      //
{
    for (parameter *p = root; p; p = p->next)
	    p->best_val = *(p->val);
	have_best = 1;
}

void parameter::restore_best()
      //
      // Resets to the best value set
      //
{
	assert(have_best > 0);

    for (parameter *p = root; p; p = p->next) {
		if (!(p->is_active))
			continue;								// skip in-active parameters
	    *(p->val) = p->best_val;
	}

	//
	// Note: restoring the known best parameter set may lead to bad intermediate
	//       parameter configurations. It is therefor necessary to restore *all*
	//       parameters fisrt, and then verify consitency. The later is done because
	//       of paranoia.
	//
	if (restore != 0)
		(restore)();								// Call the restore function, it it exists

    for (parameter *p = root; p; p = p->next) {
		if (!(p->is_active))
			continue;								// skip in-active parameters

		if (p->reject) {
			int i = (p->reject)();					// Mind the side-effects
			assert(i == 0);
		}
	}
}

int parameter::set(double x)
	//
	// Tries to set a parameter to a new value
	//
{
	if (x < min || x > max)
		return 1;									// Rejected because it is out of range

	double t = *val;								// save current value

	if (t == x)
		return 0;									// No change: do nothing

	*val = x;
	if (!disable_checks && reject && (reject)()) {
		*val = t;									// Undo change
		(reject)();
		return 1;
	}

	return 0;
}

int parameter::check_all()
	//
	// Verify that all parameter reject functions are OK
	//
{
	int ec = 0;

	assert(disable_checks == 0);

	for (parameter *pp = root; pp; pp = pp->next) {
		if (!(pp->reject))
			continue;
		int e = (pp->reject)();
		if (e)
			fprintf(stderr, "parameter::check_all - '%s' reject-function = %d\n",
				pp->name, e);

		ec |= e;
	}

	return ec;
}

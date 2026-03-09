//
// A few functions to use simulated annealing to optimmize a function
//

#ifndef ANNEAL_DEFINED_
#define ANNEAL_DEFINED_

#include <vector>
using namespace std;

#include "xrand.h"

//
// The next function must be defined in the calleing program:
// Parameters will only changed one at a time, thus it is sufficient
// that Change_param() is prepared to keep just enough state so that
// unchange_param() can undo the last change.
//

double comp_eval();					// Objective function to tune!
									// Lower values are better, 0 is the best possible

#define _ADAPTIVE_RANGEING_  		// If defined, the change scale is adjusted
									// individually for each parameter

#define MIN_RANGE	1e-6			// Min. parameter change range

void optimize(FILE *of);			// The actual optimization function
int read_an_schedule(char *fn);		// Reads an annealing schedule

//////////////////////////////////////////////////
//
// parameter data structure
//
class parameter {					// Something that the anealing process may tweek

	char		*name;				// Name of this paramter (needed to access it at run-time)
    double		*val;				// Pointer to actual value (see above)
	int			(*reject)();		// Pointer to optional reject function
									// Note: see code for info about using this facility!!

    double		best_val;			// Best value encountered so far
	static int	have_best;			// is !0 if best value is defined

    double		min, max;			// Bound on the tuning range

#ifdef _ADAPTIVE_RANGEING_
	int			n_reject;			// #of of rejects in a row
	int			n_accept;			// #of accepts in a row
	int			acc_pos;			// last accept was a positive change
	double		range;				// Change-range
#endif

	static vector<parameter*> p_ptr;// Subject to anealing

	static parameter *root;			// List of all parameters
    parameter	*next;				// Next-pointer

	static int	changed_flag;		// Change-machinery: indicate validity
	static parameter *changed_param;// Which was changed last
	static double	changed_value;	// What was its pre-cahnge value

	unsigned	is_active:1;		// =1 if this parameter is being optimized

	void update_tune_list();		// Update the list of tunable parameters

	static int	disable_checks;		// If asserted (!= 0) kill the parameter check functions

	static void (*restore)();		// Optional restore function to be used after
									// the best parameter set is restored

	//double		last_good;			// Used for debugging //**

public:
	parameter (	char *name,			// Name of this parameter
				double *ptr,		// Pointer to actual value
				double min,			// Min value allowed
				double max,			// Max value allowed
				int (*reject)() = 0	// Optional reject function:
									// Will be called when parameter is changed and may be
									// used to adjust derived values. If the retrurned
									// value is != 0, the change will be undone and another
									// change will be attempted. After MAX_CHANGE_TRIES
									// attempts, another parameter will be tried and
									// and error message is produced on <stderr>
									// Note: see code for info about using this facility!!
			   );

	static void change_param(double range, FILE *of = 0);	// This function changes a parameter
	static void unchange_param(FILE *of = 0);				// Undo the last change

	static void save_best();		// Shall save all parameters
	static void restore_best();		// Shall reset all parameters to the most
									// recent set saved via save_best()
	static void set_restore_funct(void (*rf)())
		{restore = rf;};			// Defines a restore function

	static int write(const char *fn);// Write all parameters to a file
	static int read(const char *fn);// Read all parameters from a file

	static void add_all();			// Allow all parameters to be tuned
	static parameter *find(char *s);// Find a particular parameter

	void add();						// Add this parameter to the list of tuned parameters
	void set_tuneable(int is_t);	// Change tunable status

	int set(double);				// Set a parameter, returns 1 if rejected

	void change_range(double min, double max); // Updated range

	static int n_tune()				// Returns number of tunable parameters
		{ return (int) p_ptr.size();};

	static int check_all();			// Verifies that all parameters satisfy their constraints
};

struct an_sched {					// anneal schedule
	double		temp;				// Starting temp for this phase
	int			n_steps;			// #of steps
};

#endif

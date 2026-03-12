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

extern "C" {
#include "parser.h"                         // Bison generated headers (in build dir)
#include "lex.yy.h"                         // Flex-generated header (in build dir)
#include "parser_interf.h"                  // needed to integrate the parser
#include "lsq_fit.h"
#include "fit_functions.h"
}

#include "xrand.h"
#include "csv_analyzer.h"
#include "parameter.h"
#include "loop_complex.h"
#include "lsq_fit_function.h"
#include "expression.h"
#include "nodes_of_interest.h"
#include "spice_deck.h"
#include "parser_interf.h"
#include "anneal.h"
#include "bo_optimizer.h"

///////////////////////////////////////////////////////////////////////////////////////////////////

const char* csv_out_path = "/tmp/JoSim2csv_analyzer.csv";      // Output fifo from Josim to this program
const char* cir_inp_path = "/tmp/JoSim2csv_analyzer.cir";      // Input fifo to Josim
                                            
const char* EVAL_PARAMETER = "eval";        // Name of the evaluation expression
const char* REJECT_PARAMETER = "reject";    // Keyword fo a reject functon (name of a parameter)
const char* BAY_OPT_OBJECTIVE = "bo_objective"; // Keyword for baysian optimization objective function

const char spice_escape = '~';              // Escape character in spice decks for parameter substitutions
const unsigned spice_src_max_char = 1024;   // Max line length in spice source file

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Static/gobal things:
//

char *snapshot_file_name = nullptr;         // File name for snapshots (if enabled)
unsigned snapshot_first = 0;
unsigned snapshot_frequency = 0;

char josim_output_buf[max_ln_length];       // Where to put the output from Josim (large, thus not put on heap)

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The core simulation function
//

void run_josim()
{
    //
    // Do the actual work:
    //   1. Run JoSIM with the edited circuit file
    //   2. Digest the JoSim output
    //
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //  1   Run JoSim                                                                            //
    ///////////////////////////////////////////////////////////////////////////////////////////
    pid_t child_pid = fork();
    if (child_pid == 0) {
        //
        // This is the child process
        //
        freopen("run.log", "w", stdout);    // Keep this file for debugging purposes
        
        int ie = execl("/usr/local/bin/josim-cli", "-a", "1", "-o", csv_out_path, cir_inp_path, nullptr);
        
        fprintf(stderr, "execl failed with return of %d\n", ie);
        perror("Failed to stat Josim");
        exit(1);
    }
    
    if (circuit.write_cir_file(cir_inp_path)) {
        fprintf(stderr, "failed to write the spice deck to JoSim\n");
        exit(1);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //  2   Parse the output from JoSim                                                          //
    ///////////////////////////////////////////////////////////////////////////////////////////
    FILE *inpf = fopen(csv_out_path, "r");
    assert(inpf != nullptr);
        
    FILE *snap_of = nullptr;
    if ((snapshot_file_name != nullptr) && (n_josim_runs >= snapshot_first) &&
        ( ((n_josim_runs - snapshot_first) % snapshot_frequency) == 0)         ) {
        //
        // Time to take a snapshot
        //
        sprintf(josim_output_buf, "%s_%u.csv", snapshot_file_name, n_josim_runs);
        snap_of = fopen(josim_output_buf, "w");
        assert(snap_of);
    }

    //
    // Process first line of CVS file
    //
    if (!fgets(josim_output_buf, max_ln_length - 1, inpf) ||
        strlen(josim_output_buf) > (max_ln_length -4)) {            // Meant to detect buffer overflow
        fprintf(stderr, "Failed to read first line\n");
        exit(1);
    }

    if (snap_of != nullptr) fputs(josim_output_buf, snap_of);       // Copy first line to snapshot file
        
    char *buf_p;
    {
        //
        // Discard white space and double quotes
        //
        char *tmp = josim_output_buf;
        for (buf_p = josim_output_buf; *buf_p; buf_p++)
            if ((*buf_p != '"') && (*buf_p != ' ') && (*buf_p != '\t') && (*buf_p != '\n'))
                *tmp++ = *buf_p;
         *tmp = 0;
    }
    
    buf_p = josim_output_buf;
    for (unsigned ic = 0; 1; ic++) {
        char *tmp = buf_p;
        while (*tmp && (*tmp != ',')) tmp++;    // Find end of string
        int last = !*tmp;       // Set if this is the last column
        *tmp = 0;               // Replance comma with 0
        //
        // buf_p is now pointing to the name string for this column
        //
        nodes_of_interest *noi_ptr = nodes_of_interest::find(buf_p);
        if (n_josim_runs == 0) {
            // First Josim run: allocate the ts for the NOI's
            if (noi_ptr == nullptr) {
                // Don't care about this column
                josim_out_columns.push_back(nullptr);
            } else {
                josim_out_columns.push_back(new time_series());
                if (noi_ptr->set_col_index(ic)) {
                    fprintf(stderr, "Josim output has multiple columns named '%s'\n", tmp);
                    exit(1);
                }
            }
        } else {
            // Verify that the column order is the same as the first running
            if (noi_ptr != nullptr) {
                if (ic != noi_ptr->get_col_index()) {
                    fprintf(stderr, "Josim run %u: column order has changed\n", n_josim_runs + 1);
                    exit(1);
                }
                josim_out_columns[ic]->reset();         // Reset the TS to accept new data
            }
        }
            
        if (last) {
            if (n_josim_runs == 0) {
                printf("Detected %lu columns\n", josim_out_columns.size());
            }
            break;
        }
            
        buf_p = tmp + 1;
    }
        
    if (n_josim_runs == 0) {
        // Verify that all NOI's are defined
        nodes_of_interest *n_ptr = nodes_of_interest::find_undef();
        if (n_ptr != nullptr) {
            fprintf(stderr, "Node of interest '%s' not found in Josim output\n", n_ptr->get_name());
            exit(1);
        }
    }

    //
    // Now read the data (in CSV format)
    //
    int line_no = 1;
    do {
        line_no++;
        if (!fgets(josim_output_buf, max_ln_length - 1, inpf))
            break;                             // Input file exhausetd
        if (snap_of != nullptr) fputs(josim_output_buf, snap_of);       // Copy a row to the snapshot file
        buf_p = josim_output_buf;
        double val;
        int nc;
        for (int i = 0; i < josim_out_columns.size(); i++) {
            if (1 > sscanf(buf_p,"%lf%n", &val, &nc)) {
                fprintf(stderr, "Josim run %u: Failed to digest line %d\n", n_josim_runs + 1, line_no);
                exit(1);
            }
            if (i == 0) {
                // Verify time column
                if (fabs(val - sim_time((double)(line_no - 2))) > 0.1e-12) {
                    fprintf(stderr, "Josim output on line %d: time expected %.10lg time found %.10lg\n",
                            line_no, sim_time((double)(line_no - 2)), val);
                    exit(1);
                }
            }
            if (josim_out_columns[i] != nullptr)
                josim_out_columns[i]->add_datum(val);   // Selectively add data
            buf_p += nc +1;
        }
    } while (1);
        
    fclose(inpf);
    if (snap_of != nullptr) {               // Close the snapshot file (if there was one)
        fclose(snap_of);
        snap_of = nullptr;
    }
        
    int status;
    waitpid(child_pid, &status, 0);         // reap the child process
    
    if (n_josim_runs == 0) printf("Read %d rows of data, exit status=%d\n", line_no, status);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// The main program
//

int main(int argc, char *argv[])
{
    rnd_init(0);

    // Create the named pipes (FIFO) (unless already present)
    if ((access(csv_out_path, F_OK) != 0) && (mkfifo(csv_out_path, 0666) == -1)) {
        perror("mkfifo failed");
        exit(1);
    }
    
    if ((access(cir_inp_path, F_OK) != 0) && (mkfifo(cir_inp_path, 0666) == -1)) {
        perror("mkfifo failed");
        exit(1);
    }
    
    if (argc != 2) {
        fprintf(stderr, "csv_analyzer <base-file-name>\n");
        exit(1);
    }
    
    sprintf(josim_output_buf, "%s.cir", argv[1]);
    if (1 > circuit.read_cir_file(josim_output_buf)) {
        fprintf(stderr, "Failed to read spice deck from '%s' \n", josim_output_buf);
        exit(1);
    }
    sprintf(josim_output_buf, "%s_sum.dat", argv[1]);
    FILE *sum_fp = fopen(josim_output_buf, "w");
    assert(sum_fp != nullptr);
    fprintf(sum_fp, "#");
    fprintf(sum_fp, " (%u):level\n", parameter::list_names(sum_fp));

    if (sim_anneal_ptr != nullptr) {
        parameter::sa_p_export(sim_anneal_ptr);
        printf("Simulated Annealing to optimize %u parameters in %lu steps\n",
               sim_anneal_ptr->n_tune(), sim_anneal_ptr->n_steps());
        sim_anneal_ptr->specify_summary_file(sum_fp);
    }
    

    //
    // Initialization and setup done.
    //
    if (sim_anneal_ptr != nullptr) {
        sim_anneal_ptr->optimize();
        sprintf(josim_output_buf, "%s_opt_params.txt", argv[1]);
        FILE *ofp = fopen(josim_output_buf, "w");
        parameter::save_result(ofp);
        fclose(ofp);
    } else if (baysian_opt != nullptr) {
        baysian_opt->specify_summary_file(sum_fp);
        sprintf(josim_output_buf, "%s_bo_params.txt", argv[1]);
        FILE *ofp = fopen(josim_output_buf, "w");
        baysian_opt->run(ofp);
        fclose(ofp);
    } else
        loop_complex.run_once(sum_fp);
    
    fclose(sum_fp);

    exit(0);
}

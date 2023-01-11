/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/
/**
 * An example of benchmarking Nomad on a COCO suite.
 * This example is inspired by the grid search example provided with coco.
 *
 * Set the global parameter BUDGET_MULTIPLIER to suit your needs.
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "Nomad/nomad.hpp"

#include "coco.h"

#define max(a,b) ((a) > (b) ? (a) : (b))


/**
 * The maximal budget for evaluations done by an optimization algorithm equals dimension * BUDGET_MULTIPLIER.
 * Increase the budget multiplier value gradually to see how it affects the runtime.
 */
static const unsigned int BUDGET_MULTIPLIER = 400;

/**
 * The maximal number of independent restarts allowed for an algorithm that restarts itself.
 */
//static const long INDEPENDENT_RESTARTS = 1e5;
static const long INDEPENDENT_RESTARTS = 0;

/**
 * The random seed. Change if needed.
 */
static const uint32_t RANDOM_SEED = 0xdeadbeef;

/**
 * A function type for evaluation functions, where the first argument is the vector to be evaluated and the
 * second argument the vector to which the evaluation result is stored.
 */
typedef void (*evaluate_function_t)(const double *x, double *y);

/**
 * A pointer to the problem to be optimized (needed in order to simplify the interface between the optimization
 * algorithm and the COCO platform).
 */
static coco_problem_t *PROBLEM;

/**
 * Calls coco_evaluate_function() to evaluate the objective function
 * of the problem at the point x and stores the result in the vector y
 */
static void evaluate_function(const double *x, double *y) {
  coco_evaluate_function(PROBLEM, x, y);
}

/**
 * Calls coco_evaluate_constraint() to evaluate the constraints
 * of the problem at the point x and stores the result in the vector y
 */
static void evaluate_constraint(const double *x, double *y) {
  coco_evaluate_constraint(PROBLEM, x, y);
}

/* Declarations of all functions implemented in this file (so that their order is not important): */
void coco_Nomad_experiment(const char *suite_name,
                           const char * suite_instances,
                           const char *suite_options,
                           const char *observer_name,
                           const char *observer_options);


void nomadOpt(evaluate_function_t evaluate_func,
              evaluate_function_t evaluate_cons,
              const size_t dimension,
              const size_t number_of_objectives,
              const size_t number_of_constraints,
              const double *lower_bounds,
              const double *upper_bounds,
              const size_t number_of_integer_variables,
              const size_t max_budget);

/* Structure and functions needed for timing the experiment */
typedef struct {
	size_t number_of_dimensions;
	size_t current_idx;
	char **output;
	size_t previous_dimension;
	size_t cumulative_evaluations;
	time_t start_time;
	time_t overall_start_time;
} timing_data_t;
static timing_data_t *timing_data_initialize(coco_suite_t *suite);
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem);
static void timing_data_finalize(timing_data_t *timing_data);

// Definition of a specific evaluator for Coco C interface
class CocoCInterfaceEval : public NOMAD::Evaluator
{
private:
    evaluate_function_t _bb_objective;
    evaluate_function_t _bb_constraints;

public:
    // Constructor
    CocoCInterfaceEval(std::shared_ptr<NOMAD::EvalParameters> evalParams,
                       evaluate_function_t bb_objective,
                       evaluate_function_t bb_constraints)
        : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB),
          _bb_objective(bb_objective),
          _bb_constraints(bb_constraints)
    {
    }

    //Destructor
    ~CocoCInterfaceEval() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override
    {
        bool eval_ok = true;
        countEval = true;

        size_t n = x.size();
        
        double * bb_inputs = new double[n];
        // collect the inputs parameters
        for (size_t i = 0; i < n; ++i)
        {
            bb_inputs[i] = x[i].todouble();
        }

        auto bbOutputType = _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
        size_t nbCons = getNbConstraints(bbOutputType);
        double * bb_constraint_outputs = new double[nbCons];
        char buffer [50];
        try
        {
            // call function for objective (single objective only)
            double f;
            _bb_objective(bb_inputs, &f);

            
            // collect outputs parameters
            // THE ORDER MUST BE OBJ PB .... PB
            std::string bbo("");
            sprintf (buffer, "%.16f ", f);
            //bbo = std::to_string(f) + " ";
            bbo = buffer;
            
            // call function for constraints
            if (nbCons > 0 )
            {
                _bb_constraints(bb_inputs, bb_constraint_outputs);
                
                for (size_t i = 0; i < nbCons; ++i)
                {
                    sprintf(buffer,"%.16f ",bb_constraint_outputs[i]);
                    bbo += buffer;
                    //bbo += std::to_string(bb_constraint_outputs[i]) + " ";
                }
            }
            x.setBBO(bbo);
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }

        delete[] bb_inputs;
        delete[] bb_constraint_outputs;
        

        return eval_ok;
    }
};


/**
 * The main method initializes the random number generator and calls the example experiment on the
 * bbob suite.
 * The dimension to keep in the suite of test is given as first argument. If empty, all dimensions are
 * considered
 */
int main(int argc, char **argv)
{
    
    coco_random_state_t *random_generator = coco_random_new(RANDOM_SEED);
    
    /* Change the log level to "warning" to get less output */
    coco_set_log_level("info");
    
    printf("Running the example experiment... (might take time, be patient)\n");
    fflush(stdout);
    
    /**
     * Start the actual experiments on a test suite and use a matching logger, for
     * example one of the following:
     *   bbob                 24 unconstrained noiseless single-objective functions
     *   bbob-biobj           55 unconstrained noiseless bi-objective functions
     *   [bbob-biobj-ext       92 unconstrained noiseless bi-objective functions]
     *   [bbob-constrained*   54 constrained noiseless single-objective functions]
     *   bbob-largescale      24 unconstrained noiseless single-objective functions in large dimension
     *   bbob-mixint          24 unconstrained noiseless single-objective functions with mixed-integer variables
     *   bbob-biobj-mixint    92 unconstrained noiseless bi-objective functions with mixed-integer variables
     *
     * Suites with a star are partly implemented but not yet fully supported.
     *
     * Adapt to your need. Note that the experiment is run according
     * to the settings, defined in example_experiment(...) below.
     */
    coco_set_log_level("info");
    
    /**
     * For more details on how to change the default suite and observer options, see
     * http://numbbo.github.io/coco-doc/C/#suite-parameters and
     * http://numbbo.github.io/coco-doc/C/#observer-parameters. */
    
    /**
      Strings to select instances, dimensions and functions. These can be passed as arguments to select the runs of  a Coco experiment.
        By selecting small batch of runs we can split the work to be done in parallel.
     */
    std::string instancesStr = "1-"; // Default: 15 instances are considerd
    std::string dimensionStr = "2,3,5,10,20,40"; // Default: all dimensions are considered
    std::string functionStr = "1-" ;
    
   // Default: all functions are considered
    if (argc >= 2)
        instancesStr = argv[1];
    if (argc >= 3)
        dimensionStr = argv[2];
    if (argc >= 4)
        functionStr = argv[3];
    const std::string suiteInstanceStr = "instances: "+instancesStr;
    const std::string suiteOptionsStr = "dimensions: "+ dimensionStr + " function_indices: " +functionStr;
    const std::string observerOptionsStr = "result_folder: Nomad_on_bbob_constrained_dim"+dimensionStr+" algorithm_name: Nomad4_inst";
    coco_Nomad_experiment("bbob-constrained", suiteInstanceStr.c_str(), suiteOptionsStr.c_str(), "bbob-constrained", observerOptionsStr.c_str());
    
    printf("Done!\n");
    fflush(stdout);
    
    coco_random_free(random_generator);
    
    return 0;
}

/**
 * Nomad solver for single-objective optimization. The
 * problem's initial solution is evaluated first.
 *
 * @param evaluate_func The function used to evaluate the objective function.
 * @param evaluate_cons The function used to evaluate the constraints.
 * @param dimension The number of variables.
 * @param number_of_objectives The number of objectives.
 * @param number_of_constraints The number of constraints.
 * @param lower_bounds The lower bounds of the region of interested (a vector containing dimension values).
 * @param upper_bounds The upper bounds of the region of interested (a vector containing dimension values).
 * @param number_of_integer_variables The number of integer variables (if > 0, all integer variables come
 * before any continuous ones).
 * @param max_budget The maximal number of evaluations.

 */
void nomadOpt(evaluate_function_t evaluate_func,
              evaluate_function_t evaluate_cons,
              const size_t dimension,
              const size_t number_of_objectives,
              const size_t number_of_constraints,
              const double *lower_bounds,
              const double *upper_bounds,
              const size_t number_of_integer_variables,
              const size_t max_budget)
{
    
    // create Nomad problem parameters
    auto p = std::make_shared<NOMAD::AllParameters>();

    //
    // Fix pb parameters using Nomad syntax
    //
    p->setAttributeValue("DISPLAY_DEGREE", 0);
    
    // Option sans VNS (toutes dimensions)
    std::string stats_file_name = "stats_" + std::string(coco_problem_get_id(PROBLEM)) + ".txt";
    
    //Option avec VNS (toutes dimensions)
    // std::string stats_file_name = "stats_"+ std::string(coco_problem_get_id(PROBLEM)) + ".txt";
    
    
    NOMAD::ArrayOfString stats_file_format (stats_file_name);
    stats_file_format.add("BBE");
    stats_file_format.add("BBO");
    p->setAttributeValue("STATS_FILE",stats_file_format);
    
    p->setAttributeValue("DIMENSION",dimension);
    
    // Option VNS pour dim 2,3,5,10
    // p->setAttributeValue("VNS_MADS_SEARCH", true);
        
    //Option VNS pour dim 20,40
    // p->setAttributeValue("VNS_MADS_SEARCH", true);
    // p->setAttributeValue("EVAL_USE_CACHE", false);

        
    
    std::string outType = "BB_OUTPUT_TYPE";  // Handle single objective. If number_of_objectives > 2  -> error
    for (size_t i =0 ; i < number_of_objectives ; i++)
        outType += " OBJ";
    for (size_t i =0 ; i < number_of_constraints ; i++)
        outType += " PB";   // Use Progressive barrier for all constraints
    p->readParamLine(outType);
    
    std::string inType;
    inType = "BB_INPUT_TYPE (";
    for (size_t i =0 ; i < dimension ; i++)
    {
        if (i < number_of_integer_variables )
            inType += " I";
        else
            inType += " R";
    }
    inType += ")";
    p->readParamLine(inType);
    
    std::string x0S, lbS, ubS;
    double *x = coco_allocate_vector(dimension);
    coco_problem_get_initial_solution(PROBLEM, x);
    x0S = "X0 (";
    lbS = "lower_bound (";
    ubS = "upper_bound (";
    for (size_t i =0 ; i < dimension ; i++)
    {
        x0S += " " + std::to_string(x[i]);
        lbS += " " + std::to_string(lower_bounds[i]);
        ubS += " " + std::to_string(upper_bounds[i]);
    }
    x0S += ")";
    lbS += ")";
    ubS += ")";
    p->readParamLine(x0S);
    p->readParamLine(lbS);
    p->readParamLine(ubS);
    coco_free_memory(x);
    
    // for reproducibility (like disabling openMP)
    p->readParamLine("NB_THREADS_OPENMP 1");
    
    // and the number of blackbox allowed
    p->setAttributeValue("MAX_BB_EVAL",max_budget);
    // p->setAttributeValue("VNS_MADS_SEARCH", true);

    p->checkAndComply()
    
    // The main step algorithm
    NOMAD::MainStep TheMainStep;
    TheMainStep.setAllParameters(p);
    auto ev = std::make_unique<CocoCInterfaceEval>(p->getEvalParams(), evaluate_func, evaluate_cons );
    TheMainStep.addEvaluator(std::move(ev));

    try
    {
        TheMainStep.start();
        TheMainStep.run();
        TheMainStep.end();
    }
    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }
    
    NOMAD::OutputQueue::Flush();
    NOMAD::MainStep::resetComponentsBetweenOptimization();

}



/**
 * Benchmarking of Nomad on a given suite with default instances
 * that can serve also as a timing experiment.
 *
 * @param suite_name Name of the suite (e.g. "bbob" or "bbob-biobj").
 * @param suite_options Options of the suite (e.g. "dimensions: 2,3,5,10,20 instance_indices: 1-5").
 * @param observer_name Name of the observer matching with the chosen suite (e.g. "bbob-biobj"
 * when using the "bbob-biobj-ext" suite).
 * @param observer_options Options of the observer (e.g. "result_folder: folder_name")
 */
void coco_Nomad_experiment(const char *suite_name,
                        const char * suite_instances,
                        const char *suite_options,
                        const char *observer_name,
                        const char *observer_options)
{

  size_t run;
  coco_suite_t *suite;
  coco_observer_t *observer;
  timing_data_t *timing_data;
    
  /* Initialize the suite and observer. */
  suite = coco_suite(suite_name,suite_instances, suite_options);
  observer = coco_observer(observer_name, observer_options);

  /* Initialize timing */
  timing_data = timing_data_initialize(suite);


  /* Iterate over all problems in the suite */
  while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {

    size_t dimension = coco_problem_get_dimension(PROBLEM);
      
    printf("\nStarting Problem %s, id=%s, type=%s, dim=%d \n",
            coco_problem_get_name(PROBLEM),
            coco_problem_get_id(PROBLEM),
            coco_problem_get_type(PROBLEM),
           static_cast<int>(coco_problem_get_dimension(PROBLEM)));

    /* Run the algorithm at least once */
    for (run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {

      long evaluations_done = (long) (coco_problem_get_evaluations(PROBLEM) +
            coco_problem_get_evaluations_constraints(PROBLEM));
      long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluations_done;

      /* Break the loop if the target was hit or there are no more remaining evaluations */
      if ((coco_problem_final_target_hit(PROBLEM) &&
           coco_problem_get_number_of_constraints(PROBLEM) == 0)
           || (evaluations_remaining <= 0))
        break;

      /* Call the optimization algorithm for the remaining number of evaluations */
      nomadOpt(evaluate_function,
               evaluate_constraint,
               dimension,
               coco_problem_get_number_of_objectives(PROBLEM),
               coco_problem_get_number_of_constraints(PROBLEM),
               coco_problem_get_smallest_values_of_interest(PROBLEM),
               coco_problem_get_largest_values_of_interest(PROBLEM),
               coco_problem_get_number_of_integer_variables(PROBLEM),
               (size_t) evaluations_remaining);

      /* Break the loop if the algorithm performed no evaluations or an unexpected thing happened */
      if (coco_problem_get_evaluations(PROBLEM) == evaluations_done) {
        printf("WARNING: Budget has not been exhausted (%lu/%lu evaluations done)!\n",
        		(unsigned long) evaluations_done, (unsigned long) dimension * BUDGET_MULTIPLIER);
        break;
      }
      else if (coco_problem_get_evaluations(PROBLEM) < evaluations_done)
        coco_error("Something unexpected happened - function evaluations were decreased!");
    
        
    }
      
      long evaluations_done = (long) (coco_problem_get_evaluations(PROBLEM) +
            coco_problem_get_evaluations_constraints(PROBLEM));
      long evaluations_remaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluations_done;

    /* Keep track of time */
    timing_data_time_problem(timing_data, PROBLEM);
  }

  /* Output and finalize the timing data */
  timing_data_finalize(timing_data);

  coco_observer_free(observer);
  coco_suite_free(suite);

}







/**
 * Allocates memory for the timing_data_t object and initializes it.
 */
static timing_data_t *timing_data_initialize(coco_suite_t *suite) {

	timing_data_t *timing_data = (timing_data_t *) coco_allocate_memory(sizeof(*timing_data));
	size_t function_idx, dimension_idx, instance_idx, i;

	/* Find out the number of all dimensions */
	coco_suite_decode_problem_index(suite, coco_suite_get_number_of_problems(suite) - 1, &function_idx,
			&dimension_idx, &instance_idx);
	timing_data->number_of_dimensions = dimension_idx + 1;
	timing_data->current_idx = 0;
	timing_data->output = (char **) coco_allocate_memory(timing_data->number_of_dimensions * sizeof(char *));
	for (i = 0; i < timing_data->number_of_dimensions; i++) {
		timing_data->output[i] = NULL;
	}
	timing_data->previous_dimension = 0;
	timing_data->cumulative_evaluations = 0;
	time(&timing_data->start_time);
	time(&timing_data->overall_start_time);

	return timing_data;
}

/**
 * Keeps track of the total number of evaluations and elapsed time. Produces an output string when the
 * current problem is of a different dimension than the previous one or when NULL.
 */
static void timing_data_time_problem(timing_data_t *timing_data, coco_problem_t *problem) {

	double elapsed_seconds = 0;
    
    if (problem != NULL) {
        time_t now;
        time(&now);
        elapsed_seconds = difftime(now, timing_data->start_time) ;
        printf("Problem solved in: %04lds\n", std::lround(elapsed_seconds));
    }
    
	if ((problem == NULL) || (timing_data->previous_dimension != coco_problem_get_dimension(problem))) {

		/* Output existing timing information */
		if (timing_data->cumulative_evaluations > 0) {
			time_t now;
			time(&now);
			elapsed_seconds = difftime(now, timing_data->start_time) / (double) timing_data->cumulative_evaluations;
			timing_data->output[timing_data->current_idx++] = coco_strdupf("d=%lu done in %.2e seconds/evaluation\n",
					timing_data->previous_dimension, elapsed_seconds);
            printf("Problem solved in: %04lds\n", std::lround(difftime(now, timing_data->start_time)));
            
		}

		if (problem != NULL) {
			/* Re-initialize the timing_data */
			timing_data->previous_dimension = coco_problem_get_dimension(problem);
			timing_data->cumulative_evaluations = coco_problem_get_evaluations(problem);
			time(&timing_data->start_time);
		}

	} else {
		timing_data->cumulative_evaluations += coco_problem_get_evaluations(problem);
	}
}

/**
 * Outputs and finalizes the given timing data.
 */
static void timing_data_finalize(timing_data_t *timing_data) {

	/* Record the last problem */
	timing_data_time_problem(timing_data, NULL);

  if (timing_data) {
  	size_t i;
  	double elapsed_seconds;
		time_t now;
		int hours, minutes, seconds;

		time(&now);
		elapsed_seconds = difftime(now, timing_data->overall_start_time);

  	printf("\n");
  	for (i = 0; i < timing_data->number_of_dimensions; i++) {
    	if (timing_data->output[i]) {
				printf("%s", timing_data->output[i]);
				coco_free_memory(timing_data->output[i]);
    	}
    }
  	hours = (int) elapsed_seconds / 3600;
  	minutes = ((int) elapsed_seconds % 3600) / 60;
  	seconds = (int)elapsed_seconds - (hours * 3600) - (minutes * 60);
  	printf("Total elapsed time: %dh%02dm%02ds\n", hours, minutes, seconds);

    coco_free_memory(timing_data->output);
    coco_free_memory(timing_data);
  }
}

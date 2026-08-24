# NomadBBO: Integration of Nomad in Python
# Based on work done by Christophe Tribes in NOMAD 3

#cython: language_level=3

from libcpp cimport bool
from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.vector cimport vector

from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np

from .python_src.problem_definition import ProblemDefinition
from .python_src.optimizers import (
    MixedOptimizer,
    QuantitativeOptimizer,
)
from .python_src.data_managers import (
    MixedDataManager,
    QuantitativeDataManager,
)


# Module-level global variable for the optimizer to be accessible in all scripts
cdef object optimizer_as_global_variable = None


__version__ = "1.0b5"

_RUN_FLAGS = """
            NomadBBO termination run flags
            ------------------------------
            1  Objective target reached or MADS converged to a feasible point.
            0  Feasible point found; budget or max iteration reached.
            -1  Mesh converged but no feasible point found.
            -2  No feasible point found; budget or max iteration reached.
            -3  Initial point failed to evaluate.
            -4  Time limit reached.
            -5  CTRL-C or user stopped.
            -6  Stop on feasible point.
            """

def version():
    return f"NomadBBO {__version__}"

def usage():
    print(f"""
    NomadBBO usage
    --------------

    Example
    -------
    from NomadBBO import ProblemDefinition, MixedOptimizer, optimize

    problem = ProblemDefinition(...)
    optimizer = MixedOptimizer()

    result = optimize(
        problemDefinition=problem,
        optimizer=optimizer,
        budget=1000
    )

    {_RUN_FLAGS}
    """)

def info():
    print(f"""
    NomadBBO {__version__}

    Python interface to NOMAD for mixed-variable blackbox optimization.

    Main objects:
        - ProblemDefinition
        - BaseOptimizer, either QuantitativeOptimizer or MixedOptimizer
        - optimize()

    For parameter help:
        NomadBBO.help("MAX_BB_EVAL")

    {_RUN_FLAGS}
    """)


# Define the interface function to get nomad help
def help(about):
    about = about.encode(u"ascii")
    printNomadHelp(about)


# ------------------ NEW ------------------------------- #
cdef class OptimizationRunState:

    cdef NomadBBOBlock uFeas
    cdef NomadBBOBlock uInfeas

    cdef vector[string] eParams

    cdef int runFlag
    cdef size_t nbIters
    cdef size_t nbEvals
    cdef string stopReason

    cdef double fReturn
    cdef double hReturn

    cdef list xReturn
    cdef list xBestFeas
    cdef list xBestFeasBBO
    cdef list xBestInfeas
    cdef list xBestInfeasBBO

    cdef int seed
    cdef int verbose
    cdef object budget
    cdef str fileLogName

    def __cinit__(self, int seed=0, object budget=None, int verbose=3, str fileLogName=""):
        self.uFeas = NomadBBOBlock()
        self.uInfeas = NomadBBOBlock()

        self.eParams = vector[string]()

        self.runFlag = 0
        self.nbIters = 0
        self.nbEvals = 0
        self.stopReason = string()

        self.fReturn = float("inf")
        self.hReturn = float("inf")

        self.xReturn = []
        self.xBestFeas = []
        self.xBestFeasBBO = []
        self.xBestInfeas = []
        self.xBestInfeasBBO = []

        self.seed = seed
        self.budget = budget
        self.verbose = verbose
        self.fileLogName = fileLogName


cdef void _appendDisplayParameters(OptimizationRunState runState, list params) except *:

    if runState.verbose == 0 or runState.verbose == 1:
        # 0: completely silent
        # 1: NOMAD is silent; NomadBBO prints the final result after optimization
        params.extend([
            "DISPLAY_DEGREE 0",
            "DISPLAY_ALL_EVAL no",
        ])

    elif runState.verbose == 2:
        # Display improvements only
        params.extend([
            "DISPLAY_DEGREE 2",
            "DISPLAY_ALL_EVAL no",
            "DISPLAY_STATS BBE (SOL) OBJ CONS_H",
        ])

    elif runState.verbose == 3:
        # Detailed display: all evaluations
        params.extend([
            "DISPLAY_DEGREE 2",
            "DISPLAY_ALL_EVAL yes",
            "DISPLAY_STATS BBE (SOL) OBJ CONS_H",
        ])
    elif runState.verbose == 4:
        # Very detailed display. More info from Nomad algorithm steps.
        params.extend([
            "DISPLAY_DEGREE 3",
        ])

    else:
        raise ValueError(
            f"verbose must be one of {{0, 1, 2, 3, 4}}, got {runState.verbose}."
        )


cdef void _appendNomadParameters(
    OptimizationRunState runState,
    object problemDefinition,
    list params,
    bint includeCatGroup):

    cdef size_t i
    cdef size_t nbParams

    cat_group = ""
    bb_input_type = "( "
    lower_bounds = "( "
    upper_bounds = "( "
    index = 0

    variable_names = problemDefinition.getVariableNames()

    for var_name in variable_names:
        var = problemDefinition.getVariableByName(var_name)

        if runState.verbose >= 3:
            print(f"\nVariable: {var_name}")
            print(f"  Type: {var.type}")
            print(f"  Bounds: {var.getBounds()}")

        if var.type == "binary":
            bb_input_type += "B "
            lower_bounds += "- "
            upper_bounds += "- "

        elif var.type == "choice":

            if runState.verbose >= 3:
                print(f"  Options: {var.options}")

            num_options = len(var.options)

            lower_bounds += "0 "
            upper_bounds += f"{num_options - 1} "
            bb_input_type += "I "
            cat_group += f"{index} "

        elif var.type in ["integer", "real"]:

            if runState.verbose >= 3:
                print(f"  Range: [{var.bounds[0]}, {var.bounds[1]}]")

            lower_bounds += f"{var.bounds[0]} "
            upper_bounds += f"{var.bounds[1]} "

            if var.type == "integer":
                bb_input_type += "I "
            else:
                bb_input_type += "R "

        else:
            raise ValueError(
                f"Unknown variable type '{var.type}' for variable '{var_name}'."
            )

        index += 1

    lower_bounds += " )"
    upper_bounds += " )"
    bb_input_type += " )"

    nbParams = len(params)
    for i in range(nbParams):
        runState.eParams.push_back(params[i].encode(u"ascii"))

    runState.eParams.push_back(
        f"DIMENSION {index}".encode(u"ascii")
    )
    runState.eParams.push_back(
        f"BB_INPUT_TYPE {bb_input_type}".encode(u"ascii")
    )
    runState.eParams.push_back(
        f"LOWER_BOUND {lower_bounds}".encode(u"ascii")
    )
    runState.eParams.push_back(
        f"UPPER_BOUND {upper_bounds}".encode(u"ascii")
    )

    if includeCatGroup and len(cat_group.strip()) > 0:
        runState.eParams.push_back(
            f"CAT_GROUP {cat_group}".encode(u"ascii")
        )

    for i, x0 in enumerate(problemDefinition.X0s):

        if runState.verbose >= 3:
            print(f"Initial point {i}: {x0}")

        x0_str = str(i) + " ( "
        x0_str += " ".join(str(coord) for coord in x0)
        x0_str += " )"

        runState.eParams.push_back(
            f"X0 {x0_str}".encode(u"ascii")
        )

    bb_output_types = " ".join(problemDefinition.bbot)
    runState.eParams.push_back(
        f"BB_OUTPUT_TYPE {bb_output_types}".encode(u"ascii")
    )



cdef dict _buildOptimizationResult(OptimizationRunState runState, object problemDefinition):
    stopReasonU = runState.stopReason.decode("utf-8")

    if runState.uFeas.size() > 0:
        p = runState.uFeas.get_x(0)
        runState.fReturn = p.getF()
        runState.hReturn = p.getH()

        for i in range(p.size()):
            runState.xReturn.append(p.get_coord(i))

        for j in range(runState.uFeas.size()):
            p = runState.uFeas.get_x(j)
            x_tmp = [p.get_coord(k) for k in range(p.size())]
            runState.xBestFeas.append(problemDefinition.convertPointToMixedVariableInput(x_tmp))
            runState.xBestFeasBBO.append([p.getBBO().encode("utf-8").decode("utf-8")])

    if runState.uInfeas.size() > 0:
        p = runState.uInfeas.get_x(0)
        runState.fReturn = p.getF()
        runState.hReturn = p.getH()
        runState.xReturn = [p.get_coord(i) for i in range(p.size())]

        for j in range(runState.uInfeas.size()):
            p = runState.uInfeas.get_x(j)
            x_tmp = [p.get_coord(k) for k in range(p.size())]
            runState.xBestInfeas.append(problemDefinition.convertPointToMixedVariableInput(x_tmp))
            runState.xBestInfeasBBO.append([p.getBBO().encode("utf-8").decode("utf-8")])

    x_return = problemDefinition.convertPointToMixedVariableInput(runState.xReturn)

    return {
        "x_single_best": x_return,
        "f_single_best": runState.fReturn,
        "h_single_best": runState.hReturn,
        "x_best_feas": runState.xBestFeas,
        "x_best_feas_bbo": runState.xBestFeasBBO,
        "x_best_infeas": runState.xBestInfeas,
        "x_best_infeas_bbo": runState.xBestInfeasBBO,
        "nb_evals": runState.nbEvals,
        "nb_iters": runState.nbIters,
        "run_flag": runState.runFlag,
        "stop_reason": stopReasonU}


# Standard NOMAD optimization for problems with binary, integer and continuous variables
def _optimize_quantitative(problemDefinition, optimizer, budget=None, verbose=3, fileLogName=""):
    # Algorithmic options from the optimizer with default implementations
    seed = int(optimizer.seed)

    # Initialization of Nomad run-state optimization paramters
    cdef OptimizationRunState runState = OptimizationRunState(
        seed=seed,
        budget=budget,
        verbose=int(verbose),
        fileLogName=fileLogName)

    # Bind optimizer with problemDefinition & DataManager
    optimizer.problemDefinition = problemDefinition
    optimizer.dataManager = QuantitativeDataManager(
        problemDefinition,
        seed=seed,
        isModelRequired=optimizer.isModelRequired,
        isModelUpdated=optimizer.updateModel,
        modelUpdateHardCap=optimizer.modelUpdateHardCap)

    global optimizer_as_global_variable
    optimizer_as_global_variable = optimizer

    # ----------------------------------------- DOE ------------------------------------------ #
    # Use DOE generator to create/complement initial points when a surrogate model is required.
    #
    # Target DOE size:
    #   min(20% of the budget, 35 * (n_continuous + n_integer + n_binary))
    #
    # If X0s are provided by the user, they count toward the target DOE size.

    # Points provided by the user
    x0s = problemDefinition.getX0s()
    if x0s is None:
        x0s = []
    X0s_size = len(x0s)

    if optimizer.isModelRequired:

        # Retrieve problem information
        idx_by_type = problemDefinition.getVariableIndicesByType()
        n_continuous = len(idx_by_type.get("real", []) or [])
        n_integer = len(idx_by_type.get("integer", []) or [])
        n_binary = len(idx_by_type.get("binary", []) or [])

        n_model_variables = n_continuous + n_integer + n_binary

        # Problem-based DOE size
        problem_doe_size = int(35 * n_model_variables)

        # Target DOE size
        if budget is not None:
            budget_doe_size = int(0.2 * budget)
            target_doe_size = min(budget_doe_size, problem_doe_size)
        else:
            target_doe_size = problem_doe_size

        # Number of additional DOE points required after accounting for user-provided X0s
        doe_sample_size = max(0, target_doe_size - X0s_size)

        # Minimum number of points recommended for GP construction
        min_required_points = n_model_variables + 1

        # Generate the missing DOE points
        if doe_sample_size > 0:

            doe_generator = optimizer.initialDOEgenerator
            generated_points = doe_generator.generate_raw_points(
                problemDefinition,
                doe_sample_size,
            )

            if generated_points:
                x0s.extend(generated_points)

            # Persist back because getX0s() returns a copy.
            problemDefinition.setX0(x0s)

        # Warn if the resulting DOE is too small for the GP
        if len(x0s) < min_required_points:
            print(
                f"Warning: only {len(x0s)} initial points are available, "
                f"but at least {min_required_points} are recommended "
                f"to construct the surrogate model."
            )

    # Store the complete initial DoE in the DataManager
    optimizer.dataManager.setInitialDoEPoints(x0s)

    # ---------------------------------------------------------------------------------------- #
    # Parameters
    params = [
        "USER_SEARCH yes",  # User-defined search method is always passed to NOMAD.
        "QUAD_MODEL_SEARCH yes",
        f"SEED {runState.seed}",
        "RNG_ALT_SEEDING true"
    ]

    if optimizer.quantPollType == "ADS":
        params.append("ADS_OPTIMIZATION yes")

    _appendDisplayParameters(runState, params)

    if runState.budget is not None:
        params.append(f"MAX_BB_EVAL {int(runState.budget)}")

    if runState.fileLogName != "":
        params.append(f"STATS_FILE {runState.fileLogName} bbe sol obj cons_h")

    # Use helper function common for mixed and quantitative optimization cases
    _appendNomadParameters(runState=runState, problemDefinition=problemDefinition, params=params, includeCatGroup=False)

    fBB = problemDefinition.evaluate

    # ------------------------ Standard Nomad run and return helper --------------------------- #
    runState.runFlag = runNomad(
                        cb, 
                        cbL, 
                        <void*> fBB,
                        # --- Callbacks --- #
                        NULL,  # For now, let's use the quad model not the custom model.   
                        modelCacheUpdate,                   
                        modelBuild,                         
                        NULL,                            # catFreePoll is not used for quantitative/standard Nomad
                        customSearch,               
                        # ----------------- #
                        runState.eParams,
                        runState.uFeas.c_block_ptr,
                        runState.uInfeas.c_block_ptr,
                        runState.nbEvals,
                        runState.nbIters,
                        runState.stopReason)
                        
    result = _buildOptimizationResult(runState, problemDefinition)

    if runState.verbose == 1:
        print(result)

    return result


def _optimize_mixed(problemDefinition, optimizer, budget=None, verbose=3, fileLogName=""):

    # Algorithmic options from the optimizer with default implementations
    seed = int(optimizer.seed)
    catKernelType = optimizer.catKernelType
    extendedPollTrigger = optimizer.extendedPollTrigger
    if extendedPollTrigger < 0:
        raise ValueError(f"ExtendedPollTrigger must be > 0, got {extendedPollTrigger}.")

    # Initialization of Nomad run-state optimization paramters
    cdef OptimizationRunState runState = OptimizationRunState(
        seed=seed,
        budget=budget,
        verbose=int(verbose),
        fileLogName=fileLogName)

    # ------------- Bind optimizer with problemDefinition & DataManager  ------------------ # 
    # Bind first: nbOfNeighbors default depends on it
    optimizer.problemDefinition = problemDefinition

    # Automatic default implementation for nbOfNeighbors
    if optimizer.nbOfNeighbors is None:
        nb_categories_per_variable = (problemDefinition.getNbCategoriesPerCategoricalVariable() or [])

        l_cat = 1
        for nb_cat in nb_categories_per_variable:
            l_cat *= int(nb_cat)

        # Default rule, at least 2 to avoid cycling between the same categorical components
        optimizer.nbOfNeighbors = int(max(2, 2 * int(np.sqrt(l_cat))))
    else:
        optimizer.nbOfNeighbors = int(optimizer.nbOfNeighbors)

    # Check if nbOfNeighbors is a positive integer
    if optimizer.nbOfNeighbors <= 0:
        raise ValueError(f"optimizer.nbOfNeighbors must be > 0, got {optimizer.nbOfNeighbors}.")

    optimizer.dataManager = MixedDataManager(
        problemDefinition,
        catKernelType=catKernelType,
        seed=seed,
        isModelRequired=optimizer.isModelRequired,
        isModelUpdated=optimizer.updateModel,
        modelUpdateHardCap=optimizer.modelUpdateHardCap)

    # Retrieve the effective kernel after constructing MixedDataManager (it has default construction rules)
    catKernelType = optimizer.dataManager.getCatKernelType()
    # ---------------------------------------------------------------------------------------------#

    # ------ Optimizer global variable for NOMAD callbacks: must be after binding above ---------- #
    global optimizer_as_global_variable
    optimizer_as_global_variable = optimizer
    # -------------------------------------------------------------------------------------------- #

    # ----------------------------------------- DOE ------------------------------------------ #
    # Use DOE generator to create/complement initial points.
    #
    # Target DOE size:
    #   min(20% of the budget, 35 * (n_continuous + n_integer + n_binary + ceil(sqrt(l_cat)))),
    # where l_cat is the number of categorical combinations.
    #
    # If X0s are provided by the user, they count toward the target DOE size.

    # Points provided by the user
    x0s = problemDefinition.getX0s()
    if x0s is None:
        x0s = []
    X0s_size = len(x0s)

    # Retrieve problem information
    idx_by_type = problemDefinition.getVariableIndicesByType()
    n_continuous = len(idx_by_type.get("real", []) or [])
    n_integer = len(idx_by_type.get("integer", []) or [])
    n_binary = len(idx_by_type.get("binary", []) or [])

    nb_categories_per_variable = (problemDefinition.getNbCategoriesPerCategoricalVariable() or [])

    l_cat = 1
    for nb_cat in nb_categories_per_variable:
        l_cat *= int(nb_cat)

    # Problem-based DOE size
    problem_doe_size = int(35 * (n_continuous + n_integer + n_binary + np.ceil(np.sqrt(l_cat))))

    # Target DOE size
    if budget is not None:
        budget_doe_size = int(0.2 * budget)
        target_doe_size = min(budget_doe_size, problem_doe_size)
    else:
        target_doe_size = problem_doe_size

    # Number of additional DOE points required after accounting for user-provided X0s
    doe_sample_size = max(0, target_doe_size - X0s_size)


    # If a model is required (catPoll=SURROGATE or isBOSearchUsed=True), 
    # determine the minimum number of points recommended for GP construction
    if optimizer.isModelRequired:

        n_model_variables = n_continuous + n_integer + n_binary
        n_cat = len(nb_categories_per_variable)

        min_required_points = 0

        if catKernelType == "GOWER":
            min_required_points = n_model_variables + n_cat + 1

        elif catKernelType == "CONT_RELAX":
            sum_cat = sum(nb_categories_per_variable)
            min_required_points = n_model_variables + sum_cat + 1

        elif catKernelType == "HOMO_HSPHERE":
            n_hsphere_params = sum(
                L * (L - 1) // 2
                for L in nb_categories_per_variable
            )
            min_required_points = n_model_variables + n_hsphere_params + 1


    # Generate the missing DOE points
    if doe_sample_size > 0:

        doe_generator = optimizer.initialDOEgenerator
        generated_points = doe_generator.generate_raw_points(problemDefinition, doe_sample_size)

        if generated_points:
            x0s.extend(generated_points)

        # Persist back because getX0s() returns a copy.
        problemDefinition.setX0(x0s)


    # Warn if the resulting DOE is too small for the GP
    if optimizer.isModelRequired and len(x0s) < min_required_points:
        print(
            f"Warning: only {len(x0s)} initial points are available, "
            f"but at least {min_required_points} are recommended for "
            f"catKernelType='{catKernelType}'."
        )


    # Store the complete initial DoE in the DataManager
    optimizer.dataManager.setInitialDoEPoints(x0s)

    # FOR NOW LET'S KEEP IT SIMPLE:
    # The DOE is given to NOMAD as initial points (X0) without BBO information.
    # In the future we will have several options here:
    # - we could also use the DOE to initialize the CatMads model cache before launching the optimization
    # - we could pass the DOE points with their BBO information in NOMAD cache
    # - ...
    # ---------------------------------------------------------------------------------------- #


    # ----------------------------------- Parameters ----------------------------------------- #    
    params = [
        f"CAT_EXTENDED_POLL_TRIGGER {extendedPollTrigger}",
        "DIRECTION_TYPE ORTHO 2N",
        "QUAD_MODEL_SEARCH yes",
        "NM_SEARCH no",
        f"SEED {runState.seed}",
        "RNG_ALT_SEEDING true"
    ]

    if optimizer.quantPollType == "ADS":
        params.append("CATADS_OPTIMIZATION yes")
    elif optimizer.quantPollType == "MADS":
        params.append("CATMADS_OPTIMIZATION yes")

    _appendDisplayParameters(runState, params)

    if runState.budget is not None:
        params.append(f"MAX_BB_EVAL {int(runState.budget)}")

    if runState.fileLogName != "":
        params.append(f"STATS_FILE {runState.fileLogName} bbe sol obj cons_h")

    # Use helper function common for mixed and quantitative optimization cases
    _appendNomadParameters(runState=runState, problemDefinition=problemDefinition, params=params, includeCatGroup=True)

    fBB = problemDefinition.evaluate
    # ---------------------------------------------------------------------------------------- #


    # ------------------------ CatMADS Nomad run and return helper --------------------------- #
    runState.runFlag = runNomad(
                        cb, 
                        cbL, 
                        <void*> fBB,
                        # --- Callbacks --- #
                        completeTrialPointsInformation,     
                        modelCacheUpdate,                   
                        modelBuild,                         
                        catMadsFreePoll,
                        customSearch,                  
                        # ----------------- #
                        runState.eParams,
                        runState.uFeas.c_block_ptr,
                        runState.uInfeas.c_block_ptr,
                        runState.nbEvals,
                        runState.nbIters,
                        runState.stopReason)

    result = _buildOptimizationResult(runState, problemDefinition)

    if runState.verbose >= 1:
        print(result)

    return result


# Interface function to perform optimization
def optimize(problemDefinition, optimizer=None, budget=None, verbose = 3, fileLogName=""):
    # If no budget is provided, then the terminiation is "Min mesh index reached" (Algo)

    # Verbose can be set as 0,1,2,3. Default value is 3.
    # 0 : Silent. The result is returned without displaying optimization output.
    # 1 : Final result only. The returned result is also displayed.
    # 2 : Progress output. Successful/improving evaluations are displayed.
    # 3 : Detailed output. All evaluations and additional run information are displayed.
    # 4 : Very detailed output. More info from Nomad algorithm steps.

    # Sanity checks: run options validation 
    verbose = int(verbose)
    if verbose not in (0, 1, 2, 3, 4):
        raise ValueError(
            f"verbose must be one of {{0, 1, 2, 3, 4}}, got {verbose}."
        )

    if budget is not None:
        budget = int(budget)
        if budget <= 0:
            raise ValueError(f"budget must be > 0, got {budget}.")


    # Detect problem type
    idx_by_type = problemDefinition.getVariableIndicesByType()
    has_categorical = len(idx_by_type.get("choice", []) or []) > 0


    # Dispatch into the problem optimizer and data_manager type
    if has_categorical:
        if optimizer is None:
            optimizer = MixedOptimizer()
        elif not isinstance(optimizer, MixedOptimizer):
            raise TypeError(
                "Problems with categorical variables require a MixedOptimizer, "
                f"got {type(optimizer).__name__}.")

        return _optimize_mixed(problemDefinition=problemDefinition, optimizer=optimizer, budget=budget, verbose=verbose, fileLogName=fileLogName)

    else:
        if optimizer is None:
            optimizer = QuantitativeOptimizer()
        elif not isinstance(optimizer, QuantitativeOptimizer):
            raise TypeError(
                "Problems only integer, continuous and binary variables require a QuantitativeOptimizer, "
                f"got {type(optimizer).__name__}.")

        return _optimize_quantitative(problemDefinition=problemDefinition, optimizer=optimizer, budget=budget, verbose=verbose, fileLogName=fileLogName)



# Define callback function for a single EvalPoint and Mixed Variables (or regular variables)  ---> link with Python
cdef int cb(void *f, shared_ptr[EvalPoint] x) noexcept:

    cdef NomadBBOEvalPoint u = NomadBBOEvalPoint()

    u.c_ep_ptr = x

    # NOTE: For now, we assume that the problem definition is stored in a global variable cat_mads_problem_definition that can be accessed in this callback.
    # f is not use.

    # Use the problem_definition from the global variable
    # set in CatMadsModelProblemInitialization to evaluate the blackbox function.
    global optimizer_as_global_variable
    if optimizer_as_global_variable is None:
        print("Error: CatMads manager is not set. Cannot evaluate the blackbox function.")
        return 0  # Échec

    problem_definition = optimizer_as_global_variable.problemDefinition
    if problem_definition is None:
        print("Error: CatMads problem definition is not set in the manager. Cannot evaluate the blackbox function.")
        return 0  # Échec

    # Extract coordinates
    point = []
    for i in range(u.size()):
        point.append(u.get_coord(i))

    # Convert for evaluate()
    eval_input = problem_definition.convertPointToMixedVariableInput(point)

    # Evaluate
    try:
        result = problem_definition.evaluate(eval_input)
        f_value = result["OBJ"]
        pb_values = result.get("PB", [])

        # The values for OBJ and PB must be set in the same order as defined in the BB_OUTPUT_TYPE parameter.
        # Loop on BB_OUTPUT_TYPE to set the corresponding values in NOMAD
        rawBBO =""
        pb_index = 0
        for output_type in problem_definition.bbot:
            if output_type == "OBJ":
                rawBBO += f"{f_value} "
            elif output_type.startswith("PB"):
                # We can have PB OBJ PB PB! So we need to keep track of the index of the PB values to set the right value for each PB output.
                if pb_index < len(pb_values):
                    pb_value = pb_values[pb_index]
                    # Set the PB value in NOMAD. We can use the getBBO/setBBO mechanism to store additional values, but we need to encode the index in the string. For example
                    rawBBO += f"{pb_value} "
                else:
                    rawBBO += f"INF "  # If not enough PB values are returned, set to INF or some default value
                pb_index += 1
            else:
                print(f"Warning: Unrecognized output type '{output_type}' in BB_OUTPUT_TYPE.")
        u.setBBO(rawBBO.encode("UTF-8"))
        return 1  # Succes
    except Exception as e:
        print(f"Evaluation error: {e}")
        return 0  # Fail



# Define callback function for a block (vector) of EvalPoints and Mixed Variables (or regular variables)  ---> link with Python
cdef vector[int] cbL(void *f, shared_ptr[Block] block) noexcept:
    cdef NomadBBOBlock u = NomadBBOBlock()

    u.c_block_ptr = block

    # A vector of int to return the success/failure of each evaluation
    cdef vector[int] results
    results.resize(u.size())

    # Use the problem_definition from the global variable
    # set in CatMadsModelProblemInitialization to evaluate the blackbox function.
    global optimizer_as_global_variable
    if optimizer_as_global_variable is None:
        print("Error: CatMads manager is not set. Cannot evaluate the blackbox function.")
        # set all results to failure
        for i in range(u.size()):
            results[i] = 0
        return results  # Échec

    problem_definition = optimizer_as_global_variable.problemDefinition
    if problem_definition is None:
        print("Error: CatMads problem definition is not set in the manager. Cannot evaluate the blackbox function.")
        # set all results to failure
        for i in range(u.size()):
            results[i] = 0
        return results  # Échec

    # Loop over each EvalPoint in the block
    cdef NomadBBOEvalPoint u_i = NomadBBOEvalPoint()
    for i in range(u.size()):
        u_i.c_ep_ptr = deref(u.c_block_ptr)[i]

        # Extraire les coordonnées
        point = []
        for k in range(u_i.size()):
            point.append(u_i.get_coord(k))

        # Convert before evaluate()
        eval_input = problem_definition.convertPointToMixedVariableInput(point)

        # Evaluate
        try:
            result = problem_definition.evaluate(eval_input)
            f_value = result["OBJ"]
            pb_values = result.get("PB", [])

            # The values for OBJ and PB must be set in the same order as defined in the BB_OUTPUT_TYPE parameter.
            # Loop on BB_OUTPUT_TYPE to set the corresponding values in NOMAD
            rawBBO =""
            pb_index = 0
            for output_type in problem_definition.bbot:
                if output_type == "OBJ":
                    rawBBO += f"{f_value} "
                elif output_type.startswith("PB"):
                    # We can have PB OBJ PB PB! So we need to keep track of the index of the PB values to set the right value for each PB output.
                    if pb_index < len(pb_values):
                        pb_value = pb_values[pb_index]
                        # Set the PB value in NOMAD. We can use the getBBO/setBBO mechanism to store additional values, but we need to encode the index in the string. For example
                        rawBBO += f"{pb_value} "
                    else:
                        rawBBO += f"INF "  # If not enough PB values are returned, set to INF or some default value
                    pb_index += 1
                else:
                    print(f"Warning: Unrecognized output type '{output_type}' in BB_OUTPUT_TYPE.")
            u_i.setBBO(rawBBO.encode("UTF-8"))
            results[i] = 1  # Succès
        except Exception as e:
            print(f"Evaluation error: {e}")
            results[i] = 0  # Échec

    return results


# Cython function called by a nomadCySimpleInterface wrapper.
# Let's use ModelManager class to completeTrialPointsInformation using model
# Used only for problems with categorical variables (MixedOptimizer)
cdef vector[vector[double]] completeTrialPointsInformation(vector[vector[double]]& trialPoints):
    """Gets a vector<vector<double>> from C++, converts to 2D numpy array, passes to model manager, returns result as 2D numpy array."""

    # Default empty vector to return in case of error or if model is not ready.
    cdef vector[vector[double]] cpp_vec


    cdef Py_ssize_t nbPoints = trialPoints.size()
    if nbPoints == 0:
        return cpp_vec  # Return empty vector if no trial points

    cdef Py_ssize_t n= trialPoints[0].size() if nbPoints > 0 else 0
    cdef np.ndarray[np.double_t, ndim=2] arrTP = np.empty((nbPoints, n), dtype=np.float64)
    cdef int i, j
    for i in range(nbPoints):
        for j in range(n):
            arrTP[i, j] = trialPoints[i][j]

    # Let's use a global model manager so that we can update the model when required
    global optimizer_as_global_variable
    if optimizer_as_global_variable is None:
        raise TypeError("ModelManager must be initialized first.")

    problem_definition = optimizer_as_global_variable.problemDefinition
    if problem_definition is None:
        raise TypeError("ModelManager must have a problem definition.")

    if not optimizer_as_global_variable.dataManager.getModelIsReady():
        # Return empty vector if no model is defined
        return cpp_vec

    # Pass the numpy array to the cat mads model manager. Get the result as 1D numpy array.
    result = optimizer_as_global_variable.dataManager.completeTrialPointsInformation(arrTP)

    if not isinstance(result, np.ndarray):
        raise TypeError("CompleteTrialPointsInformation Callback did not return a numpy array.")

    if result.ndim != 2:
        raise ValueError("CompleteTrialPointsInformation Callback: result is not a 2D numpy array")

    # Check that the first dimension matches nbPoints
    if result.shape[0] != nbPoints:
        raise ValueError(f"Expected result with first dimension {nbPoints}, got {result.shape[0]}")

    # Optionally ensure contiguous and dtype float64
    if not result.flags['C_CONTIGUOUS']:
        result = np.ascontiguousarray(result, dtype=np.float64)

    # Convert np array result into a vector of vectors of double
    cdef double[:, :] result_mv = result  # Create a 2D memoryview
    cdef vector[double] row_vec
    cdef Py_ssize_t result_cols = result.shape[1]

    for i in range(nbPoints):
        row_vec.clear()  # Clear the vector for reuse
        for j in range(result_cols):
            row_vec.push_back(result_mv[i, j])
        cpp_vec.push_back(row_vec)

    return cpp_vec


# Cython function called by a nomadCySimpleInterface wrapper to update Model cache.
# Let's use ModelManager class to update model with an eval point
cdef bool modelCacheUpdate(vector[double]& evalPoint):
    """Gets a vector<double> from C++, converts to 1D numpy array, passes to Model Manager, returns bool for success."""

    cdef Py_ssize_t nAndM= evalPoint.size()
    cdef np.ndarray[np.double_t, ndim=1] arrEP = np.empty(nAndM, dtype=np.float64)
    cdef int i
    for i in range(nAndM):
        arrEP[i] = evalPoint[i]

    # Let's use a global cat_mads model manager so that we can update the model when required
    global optimizer_as_global_variable
    if optimizer_as_global_variable is None:
        raise TypeError("ModelManager must be initialized first.")

    problem_definition = optimizer_as_global_variable.problemDefinition
    if problem_definition is None:
        raise TypeError("ModelManager must have a problem definition.")

    # Pass the numpy array to the cat mads model manager. Get the result as 1D numpy array.
    success = optimizer_as_global_variable.dataManager.addToCache(arrEP)

    return success

# Cython function called by a nomadCySimpleInterface wrapper to build Model.
# Let's use ModelManager class to build model
cdef bool modelBuild():
    """Calls Model Manager to build the model, returns bool for success."""

    # Let's use a global model manager so that we can update the model when required
    global optimizer_as_global_variable
    if optimizer_as_global_variable is None:
        raise TypeError("ModelManager must be initialized first.")

    problem_definition = optimizer_as_global_variable.problemDefinition
    if problem_definition is None:
        raise TypeError("ModelManager must have a problem definition.")

    # No surrogate model is required for this optimization
    if not optimizer_as_global_variable.isModelRequired:
        return True

    # Call the build model function
    success = optimizer_as_global_variable.dataManager.constructModel()

    return success


# Cython function called by a nomadCySimpleInterface wrapper to generate CatMads poll on categorical variables.
# Let's use CatMads ModelManager class to generate CatMads poll direction with an eval point
cdef vector[vector[double]] catMadsFreePoll(vector[double]& evalPoint):
    """Gets a vector<double> from C++, converts to 1D numpy array, passes to CatMads Model Manager, returns a vector of vector of double (trial points directions)."""

    # print('Hello world of cb')
    cdef Py_ssize_t nAndM= evalPoint.size()
    cdef np.ndarray[np.double_t, ndim=1] arrEP = np.empty(nAndM, dtype=np.float64)
    cdef int i, j

    # Initialize an empty list to hold the evaluation point values
    cdef list ep_list = []

    for i in range(nAndM):
        # print(evalPoint[i])
        arrEP[i] = evalPoint[i]

    print("CatMadsFreePoll callback called with evalPoint: \n", arrEP  )

    # Let's use a global cat_mads model manager so that we can update the model when required
    global optimizer_as_global_variable
    if optimizer_as_global_variable is None:
        print("Error: CatMads manager is not set. Cannot generate cat poll directions.")
        raise TypeError("CatMads ModelManager must be initialized first.")

    problem_definition = optimizer_as_global_variable.problemDefinition
    if problem_definition is None:
        raise TypeError("CatMads ModelManager must have a problem definition.")

    # Pass the numpy array to the cat mads model manager.
    # Get the result as 2D numpy array.
    cdef np.ndarray[np.double_t, ndim=2] tpNP
    if hasattr(optimizer_as_global_variable, 'categoricalPoll'):
        tpNP = optimizer_as_global_variable.categoricalPoll(arrEP)
    else:
        raise ModuleNotFoundError("CatMads ModelManager must have a function categoricalPoll.")

    #print("CatMads poll directions generated by categoricalPoll:", tpNP)
    #print("Number of Cat poll directions :", tpNP.shape[0])
    #print("Dimension of a Cat poll direction :", tpNP.shape[1])

    # Convert np array result into a vector of vector of double
    cdef vector[vector[double]] directions
    cdef vector[double] direction
    # Pointer to the underlying data of the NumPy array
    cdef double* data_ptr
    for i in range(tpNP.shape[0]):
        direction.resize(tpNP.shape[1])
        data_ptr= &tpNP[i, 0]  # Pointer to the start of the i-th row
        for j in range(tpNP.shape[1]):
            direction[j] = data_ptr[j]
        directions.push_back(direction)

    return directions


# Cython function called by a nomadCySimpleInterface wrapper to generate search points.
# Let's use ModelManager class to generate search points with a feasible and an infeasible eval point
cdef vector[vector[double]] customSearch(vector[double]& feasEvalPoint, vector[double]& infeasEvalPoint):
    """Gets one or two vector<double> from C++,
    converts to 1D numpy array,
    passes to Model Manager,
    returns a vector of vector of double (trial points).
    NOTE: User-defined custom search method is always passed to NOMAD. 
          In BOSearch, a flag enables/disables the production of points
    """

    # print('Hello world of cb')
    cdef Py_ssize_t nFeas= feasEvalPoint.size()
    cdef Py_ssize_t nInfeas= infeasEvalPoint.size()

    if nFeas != 0 and nInfeas != 0 and nFeas != nInfeas:
        raise ValueError(f"Feasible and infeasible evaluation points must have the same size. Got {nFeas} and {nInfeas}.")

    if nFeas == 0 and nInfeas == 0:
        raise ValueError("At least one of the evaluation points must be non-empty.")

    cdef np.ndarray[np.double_t, ndim=1] feasEP = np.empty(nFeas, dtype=np.float64)
    cdef np.ndarray[np.double_t, ndim=1] infeasEP = np.empty(nInfeas, dtype=np.float64)
    cdef int i, j

    if nFeas > 0:
        for i in range(nFeas):
            feasEP[i] = feasEvalPoint[i]
    if nInfeas > 0:
        for i in range(nInfeas):
            infeasEP[i] = infeasEvalPoint[i]


    # Let's use a global model manager so that we can update the model when required
    global optimizer_as_global_variable
    if optimizer_as_global_variable is None:
        print("Error: ModelManager is not set. Cannot generate search points.")
        raise TypeError("ModelManager must be initialized first.")

    problem_definition = optimizer_as_global_variable.problemDefinition
    if problem_definition is None:
        raise TypeError("ModelManager must have a problem definition.")

    # Pass the numpy array to the model manager.
    # Get the result as 2D numpy array.
    cdef np.ndarray[np.double_t, ndim=2] tpNP = np.empty((0, 0), dtype=np.float64)
    if hasattr(optimizer_as_global_variable, 'BOSearch'):
        tpNP = optimizer_as_global_variable.BOSearch(feasEP, infeasEP)
    else:
        raise ModuleNotFoundError("ModelManager must have a user search function called BOSearch.")

    print("User search (BO) points:", tpNP)

    # Convert np array result into a vector of vector of double
    cdef vector[vector[double]] trialPoints
    cdef vector[double] tp
    # Pointer to the underlying data of the NumPy array
    cdef double* data_ptr
    for i in range(tpNP.shape[0]):
        tp.resize(tpNP.shape[1])
        data_ptr= &tpNP[i, 0]  # Pointer to the start of the i-th row
        for j in range(tpNP.shape[1]):
            tp[j] = data_ptr[j]
        trialPoints.push_back(tp)

    return trialPoints




cdef class NomadBBOPoint:

    def __cinit__(self):
        self.c_p = Point()

    def __cinit__(self, vector[double] & v):
        self.c_p = Point(v)

    def get_coord(self, size_t i):
      cdef NomadBBODouble coord = NomadBBODouble()
      # coord.c_d = deref(self.c_p_ptr)[i]
      coord.c_d = self.c_p[i]
      cdef double coord_d
      if (coord.isDefined()):
        coord_d = coord.todouble()
      else:
        coord_d = float("inf")
      return coord_d

    def size(self):
      cdef size_t n
      #n = deref(self.c_p_ptr).size()
      n = self.c_p.size()
      return n

cdef class NomadBBOArrayOfDouble:

    def __cinit__(self):
      self.c_aod = ArrayOfDouble()

    def __cinit__(self, vector[double] & v):
      self.c_aod = ArrayOfDouble(v)

    def get_coord(self, size_t i):
      cdef NomadBBODouble coord = NomadBBODouble()
      coord.c_d = self.c_aod[i]
      cdef double coord_d
      if (coord.isDefined()):
        coord_d = coord.todouble()
      else:
        coord_d = float("inf")
      return coord_d

    def size(self):
      cdef size_t n
      n = self.c_aod.size()
      return n

cdef class NomadBBODouble:

    def __cinit__(self, double v):
        self.c_d = Double(v)
    def __cinit__(self):
        pass
    def todouble(self):
        return self.c_d.todouble()
    def isDefined(self):
        return self.c_d.isDefined()


cdef extern from *:
    """
    inline NOMAD::FHComputeType get_bb_compute_type() {
        NOMAD::FHComputeTypeS cts { NOMAD::ComputeType::STANDARD, NOMAD::HNormType::L2 };
        NOMAD::FHComputeType ct { NOMAD::EvalType::BB, cts };
        return ct;
    }
    inline NOMAD::FHComputeType get_model_compute_type() {
        NOMAD::FHComputeTypeS cts { NOMAD::ComputeType::STANDARD, NOMAD::HNormType::L2 };
        NOMAD::FHComputeType ct { NOMAD::EvalType::MODEL, cts };
        return ct;
    }
    inline NOMAD::FHComputeType get_surrogate_compute_type() {
        NOMAD::FHComputeTypeS cts { NOMAD::ComputeType::STANDARD, NOMAD::HNormType::L2 };
        NOMAD::FHComputeType ct { NOMAD::EvalType::SURROGATE, cts };
        return ct;
    }
    """
    FHComputeType get_bb_compute_type()
    FHComputeType get_model_compute_type()
    FHComputeType get_surrogate_compute_type()


cdef class NomadBBOEvalPoint:
    cdef void _require_ptr(self):
        if self.c_ep_ptr.get() == NULL:
            raise RuntimeError("NomadBBOEvalPoint: null EvalPoint pointer")

    def get_coord(self, size_t i):
        self._require_ptr()
        cdef NomadBBODouble coord = NomadBBODouble()
        coord.c_d = deref(self.c_ep_ptr)[i]
        cdef double coord_d
        if (coord.isDefined()):
            coord_d = coord.todouble()
        else:
            coord_d = float("inf")
        return coord_d


    def setBBO(self, string bbo):
        self._require_ptr()
        deref(self.c_ep_ptr).setBBO(bbo)

    def getBBO(self):
        self._require_ptr()
        return deref(self.c_ep_ptr).getBBO(EvalType.BB).decode()

    def getF(self):
        self._require_ptr()
        cdef NomadBBODouble f = NomadBBODouble()
        cdef FHComputeTypeS defaultFHComputeTypeS = FHComputeTypeS(ComputeType.STANDARD,
                                                                   HNormType.L2)
        cdef FHComputeType defaultFHComputeType = FHComputeType(EvalType.BB,
                                                                FHComputeTypeS)
        f.c_d = deref(self.c_ep_ptr).getF(defaultFHComputeType)
        cdef double f_d
        if ( f.isDefined() ):
            f_d = f.todouble()
        else:
            f_d = float('inf')
        return f_d

    def getH(self):
        self._require_ptr()
        cdef NomadBBODouble h = NomadBBODouble()
        cdef FHComputeTypeS defaultFHComputeTypeS = FHComputeTypeS(ComputeType.STANDARD,
                                                                   HNormType.L2)
        cdef FHComputeType defaultFHComputeType = FHComputeType(EvalType.BB,
                                                                FHComputeTypeS)
        h.c_d = deref(self.c_ep_ptr).getH(defaultFHComputeType)
        cdef double h_d
        if ( h.isDefined() ):
            h_d = h.todouble()
        else:
            h_d = 0
        return h_d

    def getModelF(self):
        self._require_ptr()
        cdef NomadBBODouble f = NomadBBODouble()
        f.c_d = deref(self.c_ep_ptr).getF(get_model_compute_type())
        cdef double f_d
        if ( f.isDefined() ):
            f_d = f.todouble()
        else:
            f_d = float('inf')
        return f_d

    def getModelH(self):
        self._require_ptr()
        cdef NomadBBODouble h = NomadBBODouble()
        h.c_d = deref(self.c_ep_ptr).getH(get_model_compute_type())
        cdef double h_d
        if ( h.isDefined() ):
            h_d = h.todouble()
        else:
            h_d = 0.0
        return h_d
    def getSurrogateF(self):
        self._require_ptr()
        cdef NomadBBODouble f = NomadBBODouble()
        f.c_d = deref(self.c_ep_ptr).getF(get_surrogate_compute_type())
        cdef double f_d
        if ( f.isDefined() ):
            f_d = f.todouble()
        else:
            f_d = float('inf')
        return f_d

    def getSurrogateH(self):
        self._require_ptr()
        cdef NomadBBODouble h = NomadBBODouble()
        h.c_d = deref(self.c_ep_ptr).getH(get_surrogate_compute_type())
        cdef double h_d
        if ( h.isDefined() ):
            h_d = h.todouble()
        else:
            h_d = 0.0
        return h_d

    def size(self):
        self._require_ptr()
        cdef size_t n
        n = deref(self.c_ep_ptr).size()
        return n

    def displayFullNomad(self):
        self._require_ptr()
        return deref(self.c_ep_ptr).display()

    def displayX(self):
        cdef str ret =''
        for i in range(self.size()):
             ret += str(self.get_coord(i)) + ' '
        return ret


cdef class NomadBBOBlock:
    # Define what is needed to use blocks
    def size(self):
        if self.c_block_ptr.get() == NULL:
            return 0
        return deref(self.c_block_ptr).size()

    def get_x(self, size_t i):
        cdef size_t n
        cdef NomadBBOEvalPoint x_i = NomadBBOEvalPoint()
        if self.c_block_ptr.get() == NULL:
            raise RuntimeError("NomadBBOBlock.get_x called on null block")
        n = deref(self.c_block_ptr).size()
        if i >= n:
            raise IndexError("NomadBBOBlock.get_x index out of range")
        x_i.c_ep_ptr = deref(self.c_block_ptr)[i]
        if x_i.c_ep_ptr.get() == NULL:
            raise RuntimeError("NomadBBOBlock.get_x returned null EvalPoint pointer")
        return x_i


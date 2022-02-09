#ifndef _NOMADSTDCINTERFACE_H_
#define _NOMADSTDCINTERFACE_H_

// strongly inspired by Ipopt C interface

#ifdef __cplusplus
extern "C"
{
#endif

    // Contains the most important things to solve a Nomad problem
    struct NomadProblemInfo;

    // Pointer to a Nomad problem
    typedef struct NomadProblemInfo *NomadProblem;

    // To pass other informations required in the blackbox (can be useful for other interfaces)
    typedef void *NomadUserDataPtr;

    // Blackbox functions types (TODO add for blackbox block evaluation functions ?)
    typedef bool (*Callback_BB_single)(int, double *, int, double *, bool *, NomadUserDataPtr);

    NomadProblem createNomadProblem(
        Callback_BB_single bb_single, // black box function
        int nb_inputs,                // number of inputs
        int nb_outputs                // number of outputs
    );

    void freeNomadProblem(NomadProblem nomad_problem);



    // parameters settings
    bool addNomadParam(NomadProblem nomad_problem, char *keyword_value_pair);

    bool addNomadValParam(NomadProblem nomad_problem, char *keyword, int value);

    bool addNomadDoubleParam(NomadProblem nomad_problem, char *keyword, double value);

    bool addNomadBoolParam(NomadProblem nomad_problem, char *keyword, bool value);

    bool addNomadStringParam(NomadProblem nomad_problem, char *keyword, char *param_str);

    bool addNomadArrayOfDoubleParam(NomadProblem nomad_problem, char *keyword, double *array_param);

    // TODO precise a return status, i.e. a stop reason

    // For the moment, do not allow the warm start
    bool solveNomadProblem(NomadProblem nomad_problem,
                           int nb_starting_points,          // number of starting points
                           double *x0s,                     // starting points
                           bool *exists_feas_sol,           // indicates if the algorithm finds a feasible solution
                           double *bb_best_x_feas,          // At the end, contains feas solution if this last one exists
                           double *bb_best_feas_outputs,    // At the end, contains feas output solution if this last one exists
                           bool *exists_inf_sol,            // Indicates if the algorithm finds an infeasible solution
                           double *bb_best_x_inf,           // At the end, contains infeas input solution if this last one exists
                           double *bb_best_inf_outputs,     // At the end, contains infeas output solution if this last one exists
                           NomadUserDataPtr user_data_ptr); // Anything, responsability is on you

#ifdef __cplusplus
}
#endif
#endif /* ifndef SYMBOL */

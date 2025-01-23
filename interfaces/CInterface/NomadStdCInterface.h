#ifndef _NOMADSTDCINTERFACE_H_
#define _NOMADSTDCINTERFACE_H_

#ifdef _MSC_VER
  #ifdef NOMAD_INTERFACE_C_DLL
    #define DLL_EXPORT_API __declspec(dllexport)
  #else
    #define DLL_EXPORT_API __declspec(dllimport)
  #endif
#else
  #define DLL_EXPORT_API
#endif

#include <stdbool.h>
#include "nomad_platform.hpp"
// strongly inspired by Ipopt C interface

#ifdef __cplusplus
extern "C"
{
#endif

// Contains the most important things to solve a Nomad problem
struct NomadProblemInfo;

// Pointer to a Nomad problem
typedef struct NomadProblemInfo *NomadProblem;

// To pass other information required in the blackbox (can be useful for other
// interfaces)
typedef void *NomadUserDataPtr;

// Blackbox functions types
typedef bool (*Callback_BB_single)(int, double *, int, double *, bool *,
                                   NomadUserDataPtr);
typedef void (*Callback_BB_block)(int, int, double *, int, double *, bool *, bool *,
                                  NomadUserDataPtr);

// Contains information concerning the solution return by Nomad
struct NomadResultInfo;

// Pointer to a Nomad result
typedef struct NomadResultInfo *NomadResult;

// Nomad problem functions API

DLL_EXPORT_API NomadProblem createNomadProblem(
    Callback_BB_single bb_single, // blackbox function (single evaluation)
    Callback_BB_block bb_block,   // blackbox function (evalution per block)
    const int nb_inputs,          // number of inputs
    const int nb_outputs          // number of outputs
);

DLL_EXPORT_API void freeNomadProblem(NomadProblem nomad_problem);

// parameters settings
DLL_EXPORT_API bool addNomadParam(const NomadProblem nomad_problem,
                                  const char *keyword_value_pair);

DLL_EXPORT_API bool addNomadValParam(const NomadProblem nomad_problem,
                                     const char *keyword,
                                     const int value);

DLL_EXPORT_API bool addNomadDoubleParam(const NomadProblem nomad_problem,
                                        const char *keyword,
                                        const double value);

DLL_EXPORT_API bool addNomadBoolParam(const NomadProblem nomad_problem,
                                      const char *keyword,
                                      const bool value);

DLL_EXPORT_API bool addNomadStringParam(const NomadProblem nomad_problem,
                                        const char *keyword,
                                        const char *param_str);

DLL_EXPORT_API bool addNomadArrayOfDoubleParam(const NomadProblem nomad_problem,
                                               const char *keyword,
                                               const double *array_param);


// Solve the problem and return the following flag
// *  1 - Objective target reached OR Mads converged (mesh criterion) to a feasible point (true problem).
// *  0 - At least one feasible point obtained and evaluation budget (single bb or block of bb) spent
//        or max iteration (user option) reached.
// * -1 - Mads mesh converged but no feasible point obtained (only infeasible) for the true problem.
// * -2 - No feasible point obtained (only infeasible) and evaluation budget (single bb or block of bb)
//        spent or max iteration (user option) reached
// * -3 - Initial point failed to evaluate
// * -4 - Time limit reached (user option)
// * -5 - CTRL-C or user stopped (callback function)
// * -6 - Stop on feasible point (user option)
// * -7 - Wrong parameters
// * -8 - Something has gone wrong with the evaluation
// NB: this function does not allow a warm-start.
DLL_EXPORT_API int solveNomadProblem(
    const NomadResult result,
    const NomadProblem nomad_problem,
    const int nb_starting_points, // number of starting points
    const double *x0s,            // starting points
    NomadUserDataPtr user_data_ptr); // Anything, responsibility is on you

// Nomad result API

DLL_EXPORT_API NomadResult createNomadResult(void);
DLL_EXPORT_API void freeNomadResult(NomadResult result);

DLL_EXPORT_API int nbInputsNomadResult(const NomadResult result);
DLL_EXPORT_API int nbOutputsNomadResult(const NomadResult result);
DLL_EXPORT_API int nbSolutionsNomadResult(const NomadResult result);
DLL_EXPORT_API bool feasibleSolutionsFoundNomadResult(const NomadResult result);

// Return true if all solutions have been loaded, false otherwise.
// WARNING: it is the responsibility of the user to allocate the correct
// memory size (nb_solutions * nb_inputs) before calling this function.
// Note that nb_solutions must be inferior or equal to the number of
// solutions returned by Nomad.
DLL_EXPORT_API bool loadInputSolutionsNomadResult(double *input_solutions,
                                                  const int nb_solutions,
                                                  const NomadResult result);

// Return true if all solutions have been loaded, false otherwise.
// WARNING: it is the responsibility of the user to allocate the correct
// memory size (nb_solutions * nb_outputs) before calling this function.
// Note that nb_solutions must be inferior or equal to the number of
// solutions returned by Nomad.
DLL_EXPORT_API bool loadOutputSolutionsNomadResult(double *output_solutions,
                                                   const int nb_solutions,
                                                   const NomadResult result);

#ifdef __cplusplus
}
#endif
#endif /* _NOMADSTDCINTERFACE_H_ */

#include "../Interfaces/NomadStdCInterface.h"

#include "../Algos/EvcInterface.hpp"
#include "../Math/RNG.hpp"
#include "../Nomad/nomad.hpp"
#include "../Param/AllParameters.hpp"
#include "../Type/LHSearchType.hpp"

#include <string.h>
#include <iostream>

struct NomadProblemInfo
{

    // parameters of the Nomad problem
    std::shared_ptr<NOMAD::AllParameters> p;

    // pointer to the black box function, whose arguments are:
    // nb_inputs: number of blackbox inputs
    // bb_inputs: array of bb_inputs
    // nb_outputs: number of blackbox outputs
    // bb_outputs: array of bb_outputs
    // count_eval: indicates if the black box counter has to be incremented or not
    // must return true if works, otherwise false.
    // WARNING: all arrays must be allocated before and deallocated after.
    Callback_BB_single bb_single;

    // TODO function of blocks of inputs
    // bool* (*bb_multiple)(int, int, double**, int double**)

    int nb_inputs;  // number of inputs
    int nb_outputs; // number of outputs

    double *x_lb; // lower bounds (can be null)
    double *x_ub; // upper bounds (can be null)

    char *type_bb_outputs; // follow the conventions of Nomad.
    char *type_bb_inputs;  // If null, all variables will be considered as real;
    // otherwise, follow the conventions of Nomad.

    double *granularity_bb_inputs; // by default, set to 0

    // display attributes
    int display_degree;
    bool display_all_eval;
    bool display_infeasible;
    bool display_unsuccessful;

    // eval parameters
    int max_bb_eval;         // maximum number of evaluations allowed.
    bool opportunistic_eval; // boolean fixing if the search is opportunistic or not
    bool use_cache;          // uses the cache or not

    // run parameters
    int lh_search_init;      // number of points evaluated by the LH before starting the optimization
    int lh_search_iter;      // number of points evaluated by the LH at each iteration
    bool speculative_search; // boolean fixing if applies speculative search or not
    bool nm_search;          // uses the nelder mead search flag

    // TODO: add models
};

NomadProblem createNomadProblem(
    Callback_BB_single bb_single, // black box function
    int nb_inputs,                // number of inputs
    int nb_outputs,               // number of outputs
    double *x_lb,                 // lower bounds (can be null)
    double *x_ub,                 // upper bounds (can be null)
    char *type_bb_inputs,         // follow the convention of Nomad (can be null)
    char *type_bb_outputs,        // follow the conventions of Nomad.
    int max_bb_eval               // maximum number of evaluations allowed
)
{
    // check inputs; redundant but no choice
    if (nb_inputs < 1 || nb_outputs < 1 || bb_single == nullptr || type_bb_outputs == nullptr || max_bb_eval < 1)
    {
        return nullptr;
    }
    NomadProblem retval = new NomadProblemInfo;

    // problems parameters
    retval->bb_single = bb_single;
    retval->nb_inputs = nb_inputs;
    retval->nb_outputs = nb_outputs;

    if (x_lb != nullptr)
    {
        retval->x_lb = new double[nb_inputs];
        for (int i = 0; i < nb_inputs; ++i)
        {
            retval->x_lb[i] = x_lb[i];
        }
    }
    else
    {
        retval->x_lb = nullptr;
    }

    if (x_ub != nullptr)
    {
        retval->x_ub = new double[nb_inputs];
        for (int i = 0; i < nb_inputs; ++i)
        {
            retval->x_ub[i] = x_ub[i];
        }
    }
    else
    {
        retval->x_ub = nullptr;
    }

    if (type_bb_inputs != nullptr)
    {
        retval->type_bb_inputs = new char[strlen(type_bb_inputs) + 1];
        strcpy(retval->type_bb_inputs, type_bb_inputs);
    }
    else
    {
        retval->type_bb_inputs = nullptr;
    }

    retval->type_bb_outputs = new char[strlen(type_bb_outputs) + 1];
    strcpy(retval->type_bb_outputs, type_bb_outputs);

    retval->granularity_bb_inputs = new double[nb_inputs];
    for (int i = 0; i < nb_inputs; ++i)
    {
        retval->granularity_bb_inputs[i] = 0.0;
    }

    // display parameters: default options
    retval->display_degree = 2;
    retval->display_all_eval = false;
    retval->display_infeasible = false;
    retval->display_unsuccessful = true;

    // eval parameters
    retval->max_bb_eval = max_bb_eval;
    retval->opportunistic_eval = true; // boolean fixing if the search is opportunistic or not
    retval->use_cache = true;          // uses the cache or not

    // run parameters
    retval->lh_search_init = 0;        // number of points evaluated by the LH before starting the optimization
    retval->lh_search_iter = 0;        // number of points evaluated by the LH at each iteration
    retval->speculative_search = true; // by default, applies speculative search
    retval->nm_search = true;          // uses the nelder mead search flag

    retval->p = std::make_shared<NOMAD::AllParameters>();

    return retval;
}

void freeNomadProblem(NomadProblem nomad_problem)
{
    // problem parameters
    if (!nomad_problem->x_lb)
    {
        delete[] nomad_problem->x_lb;
    }
    if (!nomad_problem->x_ub)
    {
        delete[] nomad_problem->x_ub;
    }
    nomad_problem->bb_single = nullptr;
    delete[] nomad_problem->type_bb_outputs;

    if (!nomad_problem->type_bb_inputs)
    {
        delete[] nomad_problem->type_bb_inputs;
    }

    delete[] nomad_problem->granularity_bb_inputs;

    // problem parameters
    nomad_problem->p = nullptr;
}

// Problem parameters
bool setNomadGranularityBBInputs(NomadProblem nomad_problem, double *granularity_bb_inputs)
{
    bool is_valid = true;
    for (int i = 0; i < nomad_problem->nb_inputs; ++i)
    {
        if (granularity_bb_inputs[i] < 0)
        {
            std::cerr << "The granularity must be positive" << std::endl;
            is_valid = false;
            break;
        }
    }
    if (is_valid)
    {
        for (int i = 0; i < nomad_problem->nb_inputs; ++i)
        {
            nomad_problem->granularity_bb_inputs[i] = granularity_bb_inputs[i];
        }
    }
    return is_valid;
}

// Display parameters
bool setNomadDisplayDegree(NomadProblem nomad_problem, int display_degree)
{
    bool is_valid = false;
    if ((display_degree >= 0) && (display_degree <= 3)) 
    {
        is_valid = true;
        nomad_problem->display_degree = display_degree;
    }
    else
    {
        std::cerr << "DISPLAY_DEGREE can only take values between 0 and 3 included" << std::endl;
    }
    return is_valid;
}

bool setNomadDisplayAllEval(NomadProblem nomad_problem, bool display_all_eval)
{
    nomad_problem->display_all_eval = display_all_eval;
    return true;
}

bool setNomadDisplayInfeasible(NomadProblem nomad_problem, bool display_infeasible)
{
    nomad_problem->display_infeasible = display_infeasible;
    return true;
}

bool setNomadDisplayUnsuccessful(NomadProblem nomad_problem, bool display_unsuccessful)
{
    nomad_problem->display_unsuccessful = display_unsuccessful;
    return true;
}

// Eval parameters
bool setNomadOpportunisticEval(NomadProblem nomad_problem, bool opportunistic_eval)
{
    nomad_problem->opportunistic_eval= opportunistic_eval;
    return true;
}

bool setNomadUseCache(NomadProblem nomad_problem, bool use_cache)
{
    nomad_problem->use_cache = use_cache;
    return true;
}

// Run parameters
bool setNomadLHSearchParams(NomadProblem nomad_problem, int lh_search_init, int lh_search_iter)
{
    bool is_valid = false;
    if ((lh_search_init < 0) || (lh_search_iter < 0)) 
    {
        std::cerr << "The LH_SEARCH must be positive or null" << std::endl;
    }
    else
    {
        is_valid = true;
        nomad_problem->lh_search_init = lh_search_init;
        nomad_problem->lh_search_iter = lh_search_iter;
    }
    return is_valid;
}

bool setNomadSpeculativeSearch(NomadProblem nomad_problem, bool speculative_search)
{
    nomad_problem->speculative_search = speculative_search;
    return true;
}

bool setNomadNMSearch(NomadProblem nomad_problem, bool nm_search)
{
    nomad_problem->nm_search = nm_search;
    return true;
}

// Solve the problem

// Definition of a specific evaluator for C interface
class CInterfaceEval : public NOMAD::Evaluator
{
private:
    Callback_BB_single _bb_single;
    NomadUserDataPtr _data_user_ptr;
    int _nbInputs;
    int _nbOutputs;
    bool _hasSgte;

public:
    // Constructor
    CInterfaceEval(std::shared_ptr<NOMAD::EvalParameters> evalParams,
                   Callback_BB_single bb_single,
                   int nbInputs,
                   int nbOutputs,
                   bool hasSgte,
                   NomadUserDataPtr user_data_ptr)
        : NOMAD::Evaluator(evalParams),
          _bb_single(bb_single),
          _nbInputs(nbInputs),
          _nbOutputs(nbOutputs),
          _hasSgte(hasSgte),
          _data_user_ptr(user_data_ptr)
    {
        if (_hasSgte)
        {
            std::cerr << "Warning: Surrogate evaluations are not currently supported." << std::endl;
            _hasSgte = false;
        }
    }

    //Destructor
    ~CInterfaceEval() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override
    {
        bool eval_ok = false;

        size_t n = x.size();

        double *bb_inputs = new double[_nbInputs];
        double *bb_outputs = new double[_nbOutputs];

        // collect the inputs parameters
        for (size_t i = 0; i < _nbInputs; ++i)
        {
            bb_inputs[i] = x[i].todouble();
        }

        try
        {
            // call function
            eval_ok = _bb_single(_nbInputs, bb_inputs, _nbOutputs, bb_outputs, &countEval, _data_user_ptr);

            // collect outputs parameters
            auto bbOutputType = _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
            std::string bbo("");
            for (size_t i = 0; i < _nbOutputs; ++i)
            {
                bbo += std::to_string(bb_outputs[i]) + " ";
            }

            const NOMAD::EvalType &evalType = getEvalType();
            x.setBBO(bbo, bbOutputType, evalType);
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }

        delete[] bb_inputs;
        delete[] bb_outputs;

        return eval_ok;
    }
};

bool solveNomadProblem(NomadProblem nomad_problem,
                       double *x0,
                       bool *exists_feas_sol,
                       double *bb_best_x_feas,
                       double *bb_best_feas_outputs,
                       bool *exists_inf_sol,
                       double *bb_best_x_inf,
                       double *bb_best_inf_outputs,
                       NomadUserDataPtr data_user_ptr)
{
    if (!x0 || !exists_feas_sol || !bb_best_x_feas || !bb_best_feas_outputs || !exists_inf_sol || !bb_best_x_inf || !bb_best_inf_outputs)
    {
        std::cerr << "All parameters must not be null" << std::endl;
        return 1;
    }
    // Configure main parameters

    // 1- pb parameters
    nomad_problem->p->getPbParams()->setAttributeValue("DIMENSION", nomad_problem->nb_inputs);

    // TODO : check according to C string C++string
    std::string type_bb_outputs_wrap(nomad_problem->type_bb_outputs);
    nomad_problem->p->getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE",
                                                         NOMAD::stringToBBOutputTypeList(type_bb_outputs_wrap));

    if (nomad_problem->type_bb_inputs != nullptr)
    {
        std::string type_bb_inputs_wrap(nomad_problem->type_bb_inputs);
        nomad_problem->p->getEvalParams()->setAttributeValue("BB_INPUT_TYPE",
                                                             NOMAD::stringToBBInputType(type_bb_inputs_wrap));
    }

    NOMAD::ArrayOfDouble granular_params(nomad_problem->nb_inputs);
    for (size_t i = 0; i < nomad_problem->nb_inputs; ++i)
    {
        granular_params[i] = nomad_problem->granularity_bb_inputs[i];
    }
    nomad_problem->p->getPbParams()->setAttributeValue("GRANULARITY", granular_params);

    // Fix bounds
    if (nomad_problem->x_lb != nullptr)
    {
        NOMAD::ArrayOfDouble lb(nomad_problem->nb_inputs);
        for (size_t i = 0; i < nomad_problem->nb_inputs; ++i)
        {
            lb[i] = nomad_problem->x_lb[i];
        }
        nomad_problem->p->getPbParams()->setAttributeValue("LOWER_BOUND", lb);
    }

    if (nomad_problem->x_ub != nullptr)
    {
        NOMAD::ArrayOfDouble ub(nomad_problem->nb_inputs);
        for (size_t i = 0; i < nomad_problem->nb_inputs; ++i)
        {
            ub[i] = nomad_problem->x_ub[i];
        }
        nomad_problem->p->getPbParams()->setAttributeValue("UPPER_BOUND", ub);
    }

    // 2- starting points
    if (x0 != nullptr)
    {
        NOMAD::Point start_x0(nomad_problem->nb_inputs);
        for (size_t i = 0; i < nomad_problem->nb_inputs; ++i)
        {
            start_x0[i] = x0[i];
        }
        nomad_problem->p->getPbParams()->setAttributeValue("X0", start_x0);
    }

    // 3- display attributes
    nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_DEGREE", nomad_problem->display_degree);
    nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_ALL_EVAL", nomad_problem->display_all_eval);
    nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_INFEASIBLE", nomad_problem->display_infeasible);
    nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_UNSUCCESSFUL", nomad_problem->display_unsuccessful);
    nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("EVAL ( SOL ) OBJ CONS_H H_MAX"));

    // 4- eval parameters
    nomad_problem->p->getEvaluatorControlParams()->setAttributeValue("MAX_BB_EVAL", nomad_problem->max_bb_eval);
    nomad_problem->p->getEvaluatorControlParams()->setAttributeValue("OPPORTUNISTIC_EVAL", nomad_problem->opportunistic_eval);
    nomad_problem->p->getEvaluatorControlParams()->setAttributeValue("USE_CACHE", nomad_problem->opportunistic_eval);
    // TODO : for the moment allow only one blackbox call.
    nomad_problem->p->getEvaluatorControlParams()->setAttributeValue<size_t>("BB_MAX_BLOCK_SIZE", 1);

    // 5- run parameters
    std::string lh_search_str = std::to_string(nomad_problem->lh_search_init) + " " + std::to_string(nomad_problem->lh_search_iter);
    nomad_problem->p->getRunParams()->setAttributeValue("LH_SEARCH", NOMAD::LHSearchType(lh_search_str));
    nomad_problem->p->getRunParams()->setAttributeValue("SPECULATIVE_SEARCH", nomad_problem->speculative_search);
    nomad_problem->p->getRunParams()->setAttributeValue("NM_SEARCH", nomad_problem->nm_search);

    // do not allow restart
    nomad_problem->p->getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    nomad_problem->p->getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);

    // parameters validation
    nomad_problem->p->checkAndComply();

    // Initialization of main results
    *exists_feas_sol = false;
    *exists_inf_sol = false;

    auto bestFeasEvalPoint = std::make_shared<NOMAD::EvalPoint>();
    auto bestInfEvalPoint = std::make_shared<NOMAD::EvalPoint>();
    bestFeasEvalPoint = nullptr;
    bestInfEvalPoint = nullptr;

    // Resolution
    try
    {
        int stopflag = 0;

        NOMAD::MainStep TheMainStep;
        TheMainStep.setAllParameters(nomad_problem->p);

        std::unique_ptr<CInterfaceEval> ev(new CInterfaceEval(nomad_problem->p->getEvalParams(),
                                                              nomad_problem->bb_single,
                                                              nomad_problem->nb_inputs,
                                                              nomad_problem->nb_outputs,
                                                              false,
                                                              data_user_ptr));
        TheMainStep.setEvaluator(std::move(ev));

        TheMainStep.start();
        stopflag = TheMainStep.run();
        TheMainStep.end();

        // Set the best feasible and best infeasible solutions ; TODO maybe change
        // for the moment, we do not consider fixed variables
        std::vector<NOMAD::EvalPoint> evalPointFeasList, evalPointInfList;
        auto nbFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointFeasList, NOMAD::Point(), NOMAD::EvalType::BB, nullptr);
        auto nbInf = NOMAD::CacheBase::getInstance()->findBestInf(evalPointInfList, NOMAD::INF, NOMAD::Point(), NOMAD::EvalType::BB, nullptr);

        if (nbInf > 0)
        {
            // Use first infeasible point. This could be generalized to show
            // all infeasible points.
            // Smart pointers should also be used.
            // For now, keep the same logic as with PyNomad 1.
            // One of the pointer is set to null to identify the type of solution
            NOMAD::EvalPoint evalPointInf = evalPointInfList[0];
            bestInfEvalPoint = std::make_shared<NOMAD::EvalPoint>(evalPointInf);
            if (0 == nbFeas)
            {
                bestFeasEvalPoint = nullptr;
            }
        }
        if (nbFeas > 0)
        {
            NOMAD::EvalPoint evalPointFeas = evalPointFeasList[0];
            bestFeasEvalPoint = std::make_shared<NOMAD::EvalPoint>(evalPointFeas);
            bestInfEvalPoint = nullptr;
        }

        if (bestInfEvalPoint != nullptr)
        {
            *exists_inf_sol = true;
            for (size_t i = 0; i < nomad_problem->nb_inputs; ++i)
            {
                bb_best_x_inf[i] = (*bestInfEvalPoint->getX())[i].todouble();
            }
            for (size_t i = 0; i < nomad_problem->nb_outputs; ++i)
            {
                bb_best_inf_outputs[i] = bestInfEvalPoint->getEval()->getBBOutput().getBBOAsArrayOfDouble()[i].todouble();
            }
        }

        if (bestFeasEvalPoint != nullptr)
        {
            *exists_feas_sol = true;
            for (size_t i = 0; i < nomad_problem->nb_inputs; ++i)
            {
                bb_best_x_feas[i] = (*bestFeasEvalPoint->getX())[i].todouble();
            }
            for (size_t i = 0; i < nomad_problem->nb_outputs; ++i)
            {
                bb_best_feas_outputs[i] = bestFeasEvalPoint->getEval()->getBBOutput().getBBOAsArrayOfDouble()[i].todouble();
            }
        }

        // reset parameters in case someone wants to restart an optimization again
        nomad_problem->p->resetToDefaultValues();
        
        // clear cache
        NOMAD::OutputQueue::Flush();
        NOMAD::CacheBase::getInstance()->clear();
        // set counter to 0
        NOMAD::CacheBase::getInstance()->resetNbCacheHits();
        NOMAD::EvcInterface::getEvaluatorControl()->setNbEval(0);
        // set seed to default for deterministic option
        NOMAD::RNG::resetPrivateSeedToDefault();

        return 0;
    }

    catch (std::exception &e)
    {
        printf("NOMAD exception (report to developper):\n%s\n", e.what());
    }

    // clean cache at the end of the iteration
    NOMAD::OutputQueue::Flush();
    NOMAD::CacheBase::getInstance()->clear();
    
    // set seed to 0 for deterministic option
    NOMAD::RNG::resetPrivateSeedToDefault();

    return -1;
}
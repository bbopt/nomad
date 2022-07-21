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
#include "../CInterface/NomadStdCInterface.h"

#include "../Algos/EvcInterface.hpp"
#include "../Cache/CacheBase.hpp"
#include "../Math/RNG.hpp"
#include "../Nomad/nomad.hpp"
#include "../Param/AllParameters.hpp"
#include "../Type/EvalSortType.hpp"
#include "../Type/DirectionType.hpp"
#include "../Type/LHSearchType.hpp"
#include "../Type/SgtelibModelFeasibilityType.hpp"
#include "../Type/SgtelibModelFormulationType.hpp"

#include <string.h>
#include <iostream>


struct NomadProblemInfo
{

    // parameters of the problem
    std::shared_ptr<NOMAD::AllParameters> p;

    // pointer to the black box function, whose arguments are:
    // nb_inputs: number of blackbox inputs
    // bb_inputs: array of bb_inputs
    // nb_outputs: number of blackbox outputs
    // bb_outputs: array of bb_outputs
    // count_eval: set to true if the evaluation has to be counted.
    // must return true if works, otherwise false.
    // WARNING: all arrays must be allocated before and deallocated after.
    Callback_BB_single bb_single;

    // TODO function of blocks of inputs
    // bool* (*bb_multiple)(int, int, double**, int double**)

    int nb_inputs;  // number of inputs
    int nb_outputs; // number of outputs
};

NomadProblem createNomadProblem(
    Callback_BB_single bb_single, // black box function
    int nb_inputs,                // number of inputs
    int nb_outputs)               // number of outputs
{
    // check inputs; redundant but no choice
    if (nb_inputs < 1 || nb_outputs < 1 || bb_single == nullptr)
    {
        return nullptr;
    }
    NomadProblem retval = new NomadProblemInfo;

    retval->bb_single = bb_single;
    retval->nb_inputs = nb_inputs;
    retval->nb_outputs = nb_outputs;

    retval->p = std::make_shared<NOMAD::AllParameters>();

    return retval;
}

void freeNomadProblem(NomadProblem nomad_problem)
{
    nomad_problem->bb_single = nullptr;
    nomad_problem->p = nullptr;
}

bool addNomadParam(NomadProblem nomad_problem, char *keyword_value_pair)
{
    nomad_problem->p->readParamLine(std::string(keyword_value_pair));
    return true;
}

bool addNomadValParam(NomadProblem nomad_problem, char *keyword, int value)
{
    nomad_problem->p->setAttributeValue(std::string(keyword), value);
    return true;
}

bool addNomadDoubleParam(NomadProblem nomad_problem, char *keyword, double value)
{
    nomad_problem->p->setAttributeValue(std::string(keyword), NOMAD::Double(value));
    return true;
}

bool addNomadBoolParam(NomadProblem nomad_problem, char *keyword, bool value)
{
    nomad_problem->p->setAttributeValue(std::string(keyword), value);
    return true;
}

bool addNomadStringParam(NomadProblem nomad_problem, char *keyword, char *param_str)
{
    // particular values
    if (std::string(keyword) == "BB_INPUT_TYPE")
    {
        nomad_problem->p->getPbParams()->setAttributeValue("BB_INPUT_TYPE",
                                                           NOMAD::stringToBBInputTypeList(std::string(param_str)));
    }
    else if (std::string(keyword) == "BB_OUTPUT_TYPE")
    {
        nomad_problem->p->getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE",
                                                             NOMAD::stringToBBOutputTypeList(std::string(param_str)));
    }
    else if (std::string(keyword) == "EVAL_QUEUE_SORT")
    {
        nomad_problem->p->getEvalParams()->setAttributeValue("EVAL_QUEUE_SORT",
                                                             NOMAD::stringToEvalSortType(std::string(param_str)));
    }
    else if (std::string(keyword) == "DIRECTION_TYPE")
    {
        nomad_problem->p->setAttributeValue("DIRECTION_TYPE", NOMAD::stringToDirectionType(std::string(param_str)));
    }
    else if (std::string(keyword) == "DIRECTION_TYPE_SECONDARY_POLL")
    {
        nomad_problem->p->setAttributeValue("DIRECTION_TYPE_SECONDARY_POLL", NOMAD::stringToDirectionType(std::string(param_str)));
    }
    else if (std::string(keyword) == "LH_SEARCH")
    {
        nomad_problem->p->getRunParams()->setAttributeValue("LH_SEARCH", NOMAD::LHSearchType(std::string(param_str)));
    }
    else if (std::string(keyword) == "SGTELIB_MODEL_FEASIBILITY")
    {
        nomad_problem->p->setAttributeValue("SGTELIB_MODEL_FEASIBILITY",
                                                            NOMAD::stringToSgtelibModelFeasibilityType(std::string(param_str)));
    }
    else if (std::string(keyword) == "SGTELIB_MODEL_FORMULATION")
    {
        nomad_problem->p->setAttributeValue("SGTELIB_MODEL_FORMULATION",
                                                            NOMAD::stringToSgtelibModelFormulationType(std::string(param_str)));
    }
    else if (std::string(keyword) == "DISPLAY_STATS")
    {
        nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString(std::string(param_str)));
    }
    // other values
    else
    {
        nomad_problem->p->setAttributeValue(std::string(keyword), std::string(param_str));
    }
    return true;
}

bool addNomadArrayOfDoubleParam(NomadProblem nomad_problem, char *keyword, double *array_param) 
{
    NOMAD::ArrayOfDouble param(nomad_problem->nb_inputs);
    for (size_t i = 0; i < (size_t)nomad_problem->nb_inputs; ++i) 
    {
        param[i] = array_param[i];
    }
    nomad_problem->p->setAttributeValue(std::string(keyword), param);
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
        for (size_t i = 0; i < (size_t)_nbInputs; ++i)
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
            for (size_t i = 0; i < (size_t)_nbOutputs; ++i)
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
                       int nb_starting_points,
                       double *x0s,
                       bool *exists_feas_sol,
                       double *bb_best_x_feas,
                       double *bb_best_feas_outputs,
                       bool *exists_inf_sol,
                       double *bb_best_x_inf,
                       double *bb_best_inf_outputs,
                       NomadUserDataPtr data_user_ptr)
{
    if ((nb_starting_points < 1) || !x0s || !exists_feas_sol || !bb_best_x_feas || !bb_best_feas_outputs || !exists_inf_sol || !bb_best_x_inf || !bb_best_inf_outputs)
    {
        std::cerr << "All parameters must not be null" << std::endl;
        return 1;
    }

    // Configure main parameters
    nomad_problem->p->getPbParams()->setAttributeValue("DIMENSION", nomad_problem->nb_inputs);

    // starting points
    if (x0s != nullptr) {
        NOMAD::ArrayOfPoint start_x0s;
        for (size_t pt_index = 0; pt_index < (size_t)nb_starting_points; ++pt_index)
        {
            NOMAD::Point start_x0(nomad_problem->nb_inputs);
            for (size_t i = 0; i < (size_t)nomad_problem->nb_inputs; ++i)
            {
                start_x0[i] = x0s[i + pt_index * nomad_problem->nb_inputs];
            }
            start_x0s.push_back(start_x0);
        }
        nomad_problem->p->setAttributeValue<NOMAD::ArrayOfPoint>("X0", start_x0s);
    }

    // TODO : for the moment allow only one blackbox call.
    nomad_problem->p->getEvaluatorControlGlobalParams()->setAttributeValue<size_t>("BB_MAX_BLOCK_SIZE", 1);

    // TODO: For the moment, let these attributes
    nomad_problem->p->getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    nomad_problem->p->getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);
    nomad_problem->p->setAttributeValue("HOT_RESTART_ON_USER_INTERRUPT", false);

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

        // Must perform check and comply before creating the evaluator
        nomad_problem->p->checkAndComply();

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
        std::vector<NOMAD::EvalPoint> evalPointFeasList, evalPointInfList;
        const NOMAD::EvalType &evalType = NOMAD::EvalType::BB;
        auto nbFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointFeasList, NOMAD::Point(), evalType, NOMAD::ComputeType::STANDARD, nullptr);
        auto nbInf = NOMAD::CacheBase::getInstance()->findBestInf(evalPointInfList, NOMAD::INF, NOMAD::Point(), evalType, NOMAD::ComputeType::STANDARD, nullptr);

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
            for (size_t i = 0; i < (size_t)nomad_problem->nb_inputs; ++i)
            {
                bb_best_x_inf[i] = (*bestInfEvalPoint->getX())[i].todouble();
            }
            for (size_t i = 0; i < (size_t)nomad_problem->nb_outputs; ++i)
            {
                bb_best_inf_outputs[i] = bestInfEvalPoint->getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble()[i].todouble();
            }
        }

        if (bestFeasEvalPoint != nullptr)
        {
            *exists_feas_sol = true;
            for (size_t i = 0; i < (size_t)nomad_problem->nb_inputs; ++i)
            {
                bb_best_x_feas[i] = (*bestFeasEvalPoint->getX())[i].todouble();
            }
            for (size_t i = 0; i < (size_t)nomad_problem->nb_outputs; ++i)
            {
                bb_best_feas_outputs[i] = bestFeasEvalPoint->getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble()[i].todouble();
            }
        }

        // reset parameters in case someone wants to restart an optimization again
        // nomad_problem->p->resetToDefaultValues();

        NOMAD::OutputQueue::Flush();
        NOMAD::MainStep::resetComponentsBetweenOptimization();

        return 0;
    }

    catch (std::exception &e)
    {
        printf("NOMAD exception (report to developper):\n%s\n", e.what());
    }

    NOMAD::OutputQueue::Flush();
    NOMAD::CacheBase::getInstance()->clear();

    return true;
}

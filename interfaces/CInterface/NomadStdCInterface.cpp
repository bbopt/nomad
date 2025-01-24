/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
#include "NomadStdCInterface.h"

#include "Algos/EvcInterface.hpp"
#include "Cache/CacheBase.hpp"
#include "Math/RNG.hpp"
#include "Nomad/nomad.hpp"
#include "Param/AllParameters.hpp"
#include "Type/DirectionType.hpp"
#include "Type/DMultiMadsSearchStrategyType.hpp"
#include "Type/EvalSortType.hpp"
#include "Type/LHSearchType.hpp"
#include "Type/SgtelibModelFeasibilityType.hpp"
#include "Type/SgtelibModelFormulationType.hpp"

#include <string>
#include <iostream>

struct NomadProblemInfo
{
    // parameters of the problem
    std::shared_ptr<NOMAD::AllParameters> p;

    // pointer to the blackbox function, whose arguments are:
    // - nb_inputs: number of blackbox inputs
    // - bb_inputs: array of bb_inputs
    // - nb_outputs: number of blackbox outputs
    // - bb_outputs: array of bb_outputs
    // - count_eval: set to true if the evaluation has to be counted.
    // must return true if works, otherwise false.
    // WARNING: all arrays must be allocated before and deallocated after.
    Callback_BB_single bb_single;

    // pointer to the black box function (evaluation per block), whose arguments are:
    // - block_size: size of the block
    // - nb_inputs: number of blackbox inputs
    // - bb_inputs: array of bb_inputs of size block_size * nb_inputs
    // - nb_outputs: number of blackbox outputs
    // - bb_outputs: array of bb_outputs of size block_size * nb_outputs
    // - count_eval: array of bool of size block_size,
    //               set to true if the evaluation for an element
    //               of the block has to be counted
    // - eval_ok: array of bool of size block_size,
    //            set to true if the evaluation for an element of the block
    //            has succeeded, false otherwise.
    // WARNING: all arrays must be allocated before and deallocated after.
    Callback_BB_block bb_block;

    int nb_inputs;  // number of inputs
    int nb_outputs; // number of outputs
};

struct NomadResultInfo
{
    int nb_inputs;
    int nb_outputs;
    int nb_solutions;
    bool feasible_solutions_found;
    double *solutions_inputs;
    double *solutions_outputs;
};

// NomadProblem API functions

NomadProblem createNomadProblem(
    Callback_BB_single bb_single, // blackbox function
    Callback_BB_block bb_block,   // blackbox function (evaluation per block)
    const int nb_inputs,                // number of inputs
    const int nb_outputs)               // number of outputs
{
    // check inputs; redundant but no choice
    if (nb_inputs < 1 || nb_outputs < 1 || (bb_single == nullptr))
    {
        return nullptr;
    }
    auto retval = new(std::nothrow) NomadProblemInfo;
    if (retval == nullptr)
        return nullptr;

    retval->bb_single = bb_single;
    retval->bb_block = bb_block;
    retval->nb_inputs = nb_inputs;
    retval->nb_outputs = nb_outputs;

    try
    {
        retval->p = std::make_shared<NOMAD::AllParameters>();
    }
    catch (std::exception& e)
    {
        delete retval;
        return nullptr;
    }

    return retval;
}

void freeNomadProblem(const NomadProblem nomad_problem)
{
    if (nomad_problem == nullptr)
        return;

    delete nomad_problem;
}

bool addNomadParam(const NomadProblem nomad_problem,
                   const char *keyword_value_pair)
{
    if (nomad_problem == nullptr)
        return false;

    try
    {
        nomad_problem->p->readParamLine(std::string(keyword_value_pair));
    }
    catch (std::exception& e)
    {
        return false;
    }
    return true;
}

bool addNomadValParam(const NomadProblem nomad_problem,
                      const char *keyword,
                      const int value)
{
    if (nomad_problem == nullptr)
        return false;

    try
    {
        nomad_problem->p->setAttributeValue(std::string(keyword), value);
    }
    catch (std::exception& e)
    {
        return false;
    }
    return true;
}

bool addNomadDoubleParam(const NomadProblem nomad_problem,
                         const char *keyword,
                         const double value)
{
    if (nomad_problem == nullptr)
        return false;

    try
    {
        nomad_problem->p->setAttributeValue(std::string(keyword), NOMAD::Double(value));
    }
    catch (std::exception& e)
    {
        return false;
    }
    return true;
}

bool addNomadBoolParam(const NomadProblem nomad_problem,
                       const char *keyword,
                       const bool value)
{
    if (nomad_problem == nullptr)
        return false;

    try
    {
        nomad_problem->p->setAttributeValue(std::string(keyword), value);
    }
    catch (std::exception& e)
    {
        return false;
    }
    return true;
}

bool addNomadStringParam(const NomadProblem nomad_problem,
                         const char *keyword,
                         const char *param_str)
{
    if (nomad_problem == nullptr)
        return false;

    // particular values
    if (std::string(keyword) == "BB_INPUT_TYPE")
    {
        try
        {
            nomad_problem->p->getPbParams()->setAttributeValue("BB_INPUT_TYPE",
                                                               NOMAD::stringToBBInputTypeList(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "BB_OUTPUT_TYPE")
    {
        try
        {
            nomad_problem->p->getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE",
                                                                 NOMAD::stringToBBOutputTypeList(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "DIRECTION_TYPE")
    {
        try
        {
            nomad_problem->p->setAttributeValue("DIRECTION_TYPE", NOMAD::stringToDirectionType(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "DIRECTION_TYPE_SECONDARY_POLL")
    {
        try
        {
            nomad_problem->p->setAttributeValue("DIRECTION_TYPE_SECONDARY_POLL", NOMAD::stringToDirectionType(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "DISPLAY_STATS")
    {
        try
        {
            nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "DMULTIMADS_NM_STRATEGY")
    {
        try
        {
            nomad_problem->p->setAttributeValue("DMULTIMADS_NM_STRATEGY",
                                                NOMAD::stringToDMultiMadsNMSearchType(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "DMULTIMADS_QUAD_MODEL_STRATEGY")
    {
        try
        {
            nomad_problem->p->setAttributeValue("DMULTIMADS_QUAD_MODEL_STRATEGY",
                                                NOMAD::stringToDMultiMadsQuadSearchType(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "EVAL_QUEUE_SORT")
    {
        try
        {
            nomad_problem->p->getEvalParams()->setAttributeValue("EVAL_QUEUE_SORT",
                                                                 NOMAD::stringToEvalSortType(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "LH_SEARCH")
    {
        try
        {
            nomad_problem->p->getRunParams()->setAttributeValue("LH_SEARCH", NOMAD::LHSearchType(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    else if (std::string(keyword) == "SGTELIB_MODEL_FEASIBILITY")
    {
        try
        {
            nomad_problem->p->setAttributeValue("SGTELIB_MODEL_FEASIBILITY",
                                                NOMAD::stringToSgtelibModelFeasibilityType(std::string(param_str)));
        }
        catch (std::exception&e )
        {
            return false;
        }
    }
    else if (std::string(keyword) == "SGTELIB_MODEL_FORMULATION")
    {
        try
        {
            nomad_problem->p->setAttributeValue("SGTELIB_MODEL_FORMULATION",
                                                NOMAD::stringToSgtelibModelFormulationType(std::string(param_str)));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    // other values
    else
    {
        try
        {
            nomad_problem->p->setAttributeValue(std::string(keyword), std::string(param_str));
        }
        catch (std::exception& e)
        {
            return false;
        }
    }
    return true;
}

bool addNomadArrayOfDoubleParam(const NomadProblem nomad_problem,
                                const char *keyword,
                                const double *array_param)
{
    if (nomad_problem == nullptr)
        return false;

    NOMAD::ArrayOfDouble param(nomad_problem->nb_inputs);
    for (size_t i = 0; i < (size_t)nomad_problem->nb_inputs; ++i)
    {
        param[i] = array_param[i];
    }
    try
    {
        nomad_problem->p->setAttributeValue(std::string(keyword), param);
    }
    catch (std::exception& e)
    {
        return false;
    }
    return true;
}

// Solve the problem

// Definition of specific evaluators for C interface
class CInterfaceEvalSingle : public NOMAD::Evaluator
{
private:
    Callback_BB_single _bb_single;
    NomadUserDataPtr _data_user_ptr;
    int _nbInputs;
    int _nbOutputs;
    bool _hasSgte;

public:
    // Constructor
    CInterfaceEvalSingle(const std::shared_ptr<NOMAD::EvalParameters>& evalParams,
                         Callback_BB_single bb_single,
                         const int nbInputs,
                         const int nbOutputs,
                         const bool hasSgte,
                         NomadUserDataPtr user_data_ptr)
        : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB),
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
    ~CInterfaceEvalSingle() override = default;

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override
    {
        bool eval_ok = false;

        auto bb_inputs = new double[_nbInputs];
        auto bb_outputs = new double[_nbOutputs];

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
            std::string bbo;
            for (size_t i = 0; i < (size_t)_nbOutputs; ++i)
            {
                bbo += std::to_string(bb_outputs[i]) + " ";
            }
            x.setBBO(bbo,_bbOutputTypeList, _evalType);
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

class CInterfaceEvalBlock : public NOMAD::Evaluator
{
private:
    Callback_BB_single _bb_single;
    Callback_BB_block _bb_block;
    NomadUserDataPtr _data_user_ptr;
    int _nbInputs;
    int _nbOutputs;
    bool _hasSgte;

public:
    // Constructor
    CInterfaceEvalBlock(const std::shared_ptr<NOMAD::EvalParameters>& evalParams,
                        Callback_BB_single bb_single,
                        Callback_BB_block bb_block,
                        const int nbInputs,
                        const int nbOutputs,
                        const bool hasSgte,
                        NomadUserDataPtr user_data_ptr)
            : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB),
              _bb_single(bb_single),
              _bb_block(bb_block),
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
    ~CInterfaceEvalBlock() override = default;

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override
    {
        bool eval_ok = false;

        auto bb_inputs = new double[_nbInputs];
        auto bb_outputs = new double[_nbOutputs];

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
            std::string bbo;
            for (size_t i = 0; i < (size_t)_nbOutputs; ++i)
            {
                bbo += std::to_string(bb_outputs[i]) + " ";
            }
            x.setBBO(bbo,_bbOutputTypeList, _evalType);
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

    std::vector<bool> eval_block(std::vector<std::shared_ptr<NOMAD::EvalPoint>> &block,
                                 const NOMAD::Double& hMax,
                                 std::vector<bool> &countEval) const override
    {
        const int block_size = (int)block.size();
        std::vector<bool> evalOk(block_size, false);

        auto bb_inputs = new double[_nbInputs * block_size];
        auto bb_outputs = new double[_nbOutputs * block_size];

        // collect the inputs parameters
        for (size_t index = 0; index < block_size; ++index)
        {
            const auto& x = block[index];
            for (size_t i = 0; i < (size_t)_nbInputs; ++i)
            {
                bb_inputs[index * _nbInputs + i] = x->operator[](i).todouble();
            }
        }

        auto eval_ok = new bool[block_size];
        auto count_eval = new bool[block_size];

        try
        {
            // call function
            _bb_block(block_size, _nbInputs, bb_inputs,
                      _nbOutputs, bb_outputs,
                      count_eval, eval_ok,
                      _data_user_ptr);

            // collect outputs parameters
            for (size_t index = 0; index < block_size; ++index)
            {
                std::string bbo;
                for (size_t i = 0; i < (size_t)_nbOutputs; ++i)
                {
                    bbo += std::to_string(bb_outputs[index * _nbOutputs + i]) + " ";
                }
                block[index]->setBBO(bbo, _bbOutputTypeList, _evalType);
            }

            // Set evalOk and countEval
            for (size_t index = 0; index < block_size; ++index)
            {
                evalOk[index] = eval_ok[index];
                countEval[index] = count_eval[index];
            }
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }


        delete[] bb_inputs;
        delete[] bb_outputs;
        delete[] eval_ok;
        delete[] count_eval;

        return evalOk;
    }
};


int solveNomadProblem(const NomadResult result,
                      const NomadProblem nomad_problem,
                      const int nb_starting_points,
                      const double *x0s,
                      NomadUserDataPtr data_user_ptr)
{
    if ((nb_starting_points < 1) || !x0s || !result || !nomad_problem)
    {
        std::cerr << "All parameters must not be null" << std::endl;
        return -7;
    }

    // Configure main parameters
    nomad_problem->p->getPbParams()->setAttributeValue("DIMENSION", nomad_problem->nb_inputs);

    // starting points
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

    // For the moment, do not allow warm-start.
    nomad_problem->p->getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    nomad_problem->p->getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);
    nomad_problem->p->setAttributeValue("HOT_RESTART_ON_USER_INTERRUPT", false);

    result->nb_inputs = 0;
    result->nb_outputs = 0;
    result->nb_solutions = 0;
    result->feasible_solutions_found = false;
    if (!result->solutions_inputs)
        delete[] result->solutions_inputs;
    if (!result->solutions_outputs)
        delete[] result->solutions_outputs;

    // Must perform a checking and compliance of parameters before the resolution
    try
    {
        nomad_problem->p->checkAndComply();
    }
    catch (std::exception& e)
    {
        printf("NOMAD exception (report to developer):\n%s\n", e.what());
        return -7;
    }

    // Resolution
    int runFlag;
    try
    {
        NOMAD::MainStep TheMainStep;
        TheMainStep.setAllParameters(nomad_problem->p);

        if (nomad_problem->bb_block == nullptr)
        {
            auto ev = std::make_unique<CInterfaceEvalSingle>(
                    nomad_problem->p->getEvalParams(), nomad_problem->bb_single,
                    nomad_problem->nb_inputs, nomad_problem->nb_outputs, false,
                    data_user_ptr);
            TheMainStep.addEvaluator(std::move(ev));
        }
        else
        {
            auto ev = std::make_unique<CInterfaceEvalBlock>(
                    nomad_problem->p->getEvalParams(),
                    nomad_problem->bb_single, nomad_problem->bb_block,
                    nomad_problem->nb_inputs, nomad_problem->nb_outputs, false,
                    data_user_ptr);
            TheMainStep.addEvaluator(std::move(ev));
        }

        TheMainStep.start();
        TheMainStep.run();
        TheMainStep.end();

        runFlag = TheMainStep.getRunFlag();

        // Access to h norm type
        const auto hNormType = nomad_problem->p->getAttributeValue<NOMAD::HNormType>("H_NORM");

        // Use default BB eval type and STANDARD compute type
        const auto evalType = NOMAD::EvalType::BB;
        const NOMAD::FHComputeType computeType= {evalType, {NOMAD::ComputeType::STANDARD, hNormType}};

        // Set the best feasible and best infeasible solutions;
        std::vector<NOMAD::EvalPoint> evalPointFeasList, evalPointInfList;
        const size_t nbFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointFeasList, NOMAD::Point(), computeType);
        const size_t nbInf = NOMAD::CacheBase::getInstance()->findBestInf(evalPointInfList, NOMAD::INF, NOMAD::Point(), computeType);

        // If there are no feasible points, return a set of infeasible solutions with minimal h value.
        if (nbFeas > 0)
        {
            result->nb_inputs = nomad_problem->nb_inputs;
            result->nb_outputs = nomad_problem->nb_outputs;
            result->nb_solutions = (int)nbFeas;
            result->feasible_solutions_found = true;
            result->solutions_inputs = new double[nbFeas * nomad_problem->nb_inputs];
            result->solutions_outputs = new double[nbFeas * nomad_problem->nb_outputs];
            for (size_t index = 0; index < (size_t)result->nb_solutions; ++index)
            {
                const auto bestXFeas = evalPointFeasList[index].getX();
                for (size_t i = 0; i < (size_t)result->nb_inputs; ++i)
                {
                    result->solutions_inputs[index * result->nb_inputs + i] = bestXFeas->operator[](i).todouble();
                }
                const auto bestFeasOutputs = evalPointFeasList[index].getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble();
                for (size_t i = 0; i < (size_t)result->nb_outputs; ++i)
                {
                    result->solutions_outputs[index * result->nb_outputs + i] =  bestFeasOutputs[i].todouble();
                }
            }
        }
        else if (nbInf > 0)
        {
            result->nb_inputs = nomad_problem->nb_inputs;
            result->nb_outputs = nomad_problem->nb_outputs;
            result->nb_solutions = (int) nbInf;
            result->feasible_solutions_found = false;
            result->solutions_inputs = new double[nbInf * nomad_problem->nb_inputs];
            result->solutions_outputs = new double[nbInf * nomad_problem->nb_outputs];
            for (size_t index = 0; index < (size_t)result->nb_solutions; ++index)
            {
                const auto bestXInf = evalPointInfList[index].getX();
                for (size_t i = 0; i < (size_t)result->nb_inputs; ++i)
                {
                    result->solutions_inputs[index * result->nb_inputs + i] = bestXInf->operator[](i).todouble();
                }
                const auto bestInfOutputs = evalPointInfList[index].getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble();
                for (size_t i = 0; i < (size_t)result->nb_outputs; ++i)
                {
                    result->solutions_outputs[index * result->nb_outputs + i] =  bestInfOutputs[i].todouble();
                }
            }
        }

        // reset parameters in case someone wants to restart an optimization again
        // nomad_problem->p->resetToDefaultValues();

        NOMAD::OutputQueue::Flush();
        NOMAD::MainStep::resetComponentsBetweenOptimization();

        return runFlag;
    }
    catch (std::exception &e)
    {
        printf("NOMAD exception (report to developer):\n%s\n", e.what());
        runFlag = -8;
    }

    NOMAD::OutputQueue::Flush();
    NOMAD::CacheBase::getInstance()->clear();

    return runFlag;
}

NomadResult createNomadResult()
{
    auto result = new(std::nothrow) NomadResultInfo;
    if (result == nullptr)
        return nullptr;

    result->nb_inputs = 0;
    result->nb_outputs = 0;
    result->nb_solutions = 0;
    result->feasible_solutions_found = false;
    result->solutions_inputs = nullptr;
    result->solutions_outputs = nullptr;
    return result;
}

void freeNomadResult(NomadResult result)
{
    if (result == nullptr)
        return;

    delete [] result->solutions_inputs;
    delete [] result->solutions_outputs;

    delete result;
}

int nbInputsNomadResult(const NomadResult result)
{
    if (result == nullptr)
        return 0;
    return result->nb_inputs;
}

int nbOutputsNomadResult(const NomadResult result)
{
    if (result == nullptr)
        return 0;
    return result->nb_outputs;
}


int nbSolutionsNomadResult(const NomadResult result)
{
    if (result == nullptr)
        return 0;
    return result->nb_solutions;
}

bool feasibleSolutionsFoundNomadResult(const NomadResult result)
{
    if (result == nullptr)
        return false;
    return result->feasible_solutions_found;
}

bool loadInputSolutionsNomadResult(double *input_solutions,
                                   const int nb_solutions,
                                   const NomadResult result)
{
    if (nb_solutions <= 0)
        return false;
    if (result == nullptr)
        return false;
    if (result->nb_solutions == 0)
        return false;
    if (result->nb_solutions < nb_solutions)
        return false;
    if (result->nb_inputs == 0)
        return false;
    if (result->solutions_inputs == nullptr)
        return false;

    const int nb_inputs = result->nb_inputs;
    for (int i = 0; i < nb_solutions; ++i)
    {
        for (int j = 0; j < nb_inputs; ++j)
        {
            input_solutions[i * nb_inputs + j] = result->solutions_inputs[i * nb_inputs + j];
        }
    }
    return true;
}

bool loadOutputSolutionsNomadResult(double *output_solutions,
                                    const int nb_solutions,
                                    const NomadResult result)
{
    if (nb_solutions <= 0)
        return false;
    if (result == nullptr)
        return false;
    if (result->nb_solutions == 0)
        return false;
    if (result->nb_solutions < nb_solutions)
        return false;
    if (result->nb_outputs == 0)
        return false;
    if (result->solutions_outputs == nullptr)
        return false;

    const int nb_outputs = result->nb_outputs;
    for (int i = 0; i < nb_solutions; ++i)
    {
        for (int j = 0; j < nb_outputs; ++j)
        {
            output_solutions[i * nb_outputs + j] = result->solutions_outputs[i * nb_outputs + j];
        }
    }
    return true;
}

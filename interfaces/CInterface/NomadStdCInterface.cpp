/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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
#include "../Nomad/nomad.hpp"
#include "../Param/AllParameters.hpp"

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

    int max_bb_eval; // maximum number of evaluations allowed
};

NomadProblem createNomadProblem(
    Callback_BB_single bb_single, // black box function
    int nb_inputs,                // number of inputs
    int nb_outputs,               // number of outputs
    double *x_lb,                 // lower bounds (can be null)
    double *x_ub,                 // upper bounds (can be null)
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

    retval->bb_single = bb_single;
    retval->nb_inputs = nb_inputs;
    retval->nb_outputs = nb_outputs;
    retval->max_bb_eval = max_bb_eval;

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

    retval->type_bb_outputs = new char[strlen(type_bb_outputs) + 1];
    strcpy(retval->type_bb_outputs, type_bb_outputs);

    retval->p = std::make_shared<NOMAD::AllParameters>();

    return retval;
}

void freeNomadProblem(NomadProblem nomad_problem)
{
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
    nomad_problem->p = nullptr;
}

// TODO think about conversion string c c++

// Display parameters
bool addNomadBoolDispParam(NomadProblem nomad_problem,
                           char *keyword,
                           bool flag)
{
    nomad_problem->p->getDispParams()->setAttributeValue(keyword, flag);
    return true;
}

bool addNomadValDispParam(NomadProblem nomad_problem,
                          char *keyword,
                          int flag)
{
    nomad_problem->p->getDispParams()->setAttributeValue(keyword, flag);
    return true;
}

bool addNomadStringDispParam(NomadProblem nomad_problem,
                             char *keyword,
                             char *param_str)
{
    nomad_problem->p->getDispParams()->setAttributeValue(keyword, param_str);
    return true;
}

// Run parameters
bool addNomadBoolRunParam(NomadProblem nomad_problem,
                          char *keyword,
                          bool flag)
{
    nomad_problem->p->getRunParams()->setAttributeValue(keyword, flag);
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
            eval_ok = _bb_single(_nbInputs, bb_inputs, _nbOutputs, bb_outputs, _data_user_ptr);

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

        countEval = true;
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
    nomad_problem->p->getPbParams()->setAttributeValue("DIMENSION", nomad_problem->nb_inputs);

    // TODO : check according to C string C++string
    std::string type_bb_outputs_wrap(nomad_problem->type_bb_outputs);
    nomad_problem->p->getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE",
                                                         NOMAD::stringToBBOutputTypeList(type_bb_outputs_wrap));

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


    // starting points
    if (x0 != nullptr) {
        NOMAD::Point start_x0(nomad_problem->nb_inputs);
        for (size_t i = 0; i < nomad_problem->nb_inputs; ++i) {
            start_x0[i] = x0[i];
        }
        nomad_problem->p->getPbParams()->setAttributeValue("X0", start_x0);
    }

    // Fix max bb evaluations
    nomad_problem->p->getEvaluatorControlParams()->setAttributeValue("MAX_BB_EVAL", nomad_problem->max_bb_eval);

    // TODO : for the moment allow only one blackbox call.
    nomad_problem->p->getEvaluatorControlParams()->setAttributeValue<size_t>("BB_MAX_BLOCK_SIZE", 1);

    // TODO: For the moment, let these attributes
    nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_DEGREE", 2);
    nomad_problem->p->getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("EVAL ( SOL ) OBJ CONS_H H_MAX"));

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
        std::vector<NOMAD::EvalPoint> evalPointFeasList, evalPointInfList;
        const NOMAD::Eval* refevalFeas;
        const NOMAD::Eval* refevalInfeas;
        auto nbFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointFeasList, NOMAD::Point(), NOMAD::EvalType::BB,refevalFeas);
        auto nbInf = NOMAD::CacheBase::getInstance()->findBestInf(evalPointInfList, NOMAD::INF, NOMAD::Point(), NOMAD::EvalType::BB,refevalInfeas);

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

        NOMAD::OutputQueue::Flush();
        NOMAD::CacheBase::getInstance()->clear();

        return 0;
    }

    catch (std::exception &e)
    {
        printf("NOMAD exception (report to developper):\n%s\n", e.what());
    }

    NOMAD::OutputQueue::Flush();
    NOMAD::CacheBase::getInstance()->clear();

    return -1;
}

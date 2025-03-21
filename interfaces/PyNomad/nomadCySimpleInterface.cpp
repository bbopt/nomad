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
// June 2019
// Version 1.0 is with NOMAD 3.
#define NOMAD_PYTHON_VERSION "2.2"

#include "Algos/EvcInterface.hpp"
#include "Math/RNG.hpp"
#include "Nomad/nomad.hpp"
#include "Param/AllParameters.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/EvalType.hpp"

#include <Python.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>

typedef int (*Callback)(void * apply,
                        std::shared_ptr<NOMAD::EvalPoint> x);
typedef std::vector<int> (*CallbackL)(void * apply,
                                      std::shared_ptr<NOMAD::Block> block);


static void printPyNomadVersion()
{
    std::cout << "-----------------------------------------------------------"  << std::endl;
    std::cout << " Python Interface to NOMAD version " << NOMAD_PYTHON_VERSION  << std::endl;
    std::cout << " NOMAD version " << NOMAD_VERSION_NUMBER                      << std::endl;
    std::cout << "-----------------------------------------------------------"  << std::endl;
}


static void printPyNomadUsage()
{
    std::cout << "--------------------------------------------------------------"       << std::endl;
    std::cout << " PyNomad interface usage"                                             << std::endl;
    std::cout << "--------------------------------------------------------------"       << std::endl;
    std::cout << "  Run NOMAD : result = PyNomad.optimize(bb, x0, lb, ub, param)"       << std::endl;
    std::cout << "--------------------------------------------------------------"       << std::endl;
    std::cout << "    Info    : PyNomad.info()"                                         << std::endl;
    std::cout << "    Help    : PyNomad.help(\"keywords\") or PyNomad.help()"           << std::endl;
    std::cout << "    Version : PyNomad.version()"                                      << std::endl;
    std::cout << "    Usage   : PyNomad.usage()"                                        << std::endl;
    std::cout << "--------------------------------------------------------------"       << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << " PyNomad.optimize input parameters:"                                  << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << "  bb      ---> name of the blackbox function - mandatory"             << std::endl;
    std::cout << "  x0      ---> list of values for initial point - optional"           << std::endl;
    std::cout << "               ex: x0 = [0, 0, 0, 0, 0]"                              << std::endl;
    std::cout << "  lb      ---> list of values for initial point - optional"           << std::endl;
    std::cout << "               ex: lb = [-1, -2, -3, -4, -5])"                        << std::endl;
    std::cout << "  ub      ---> list of values for initial point - optional"           << std::endl;
    std::cout << "               ex: ub = [1, 2, 3, 4, 5]"                              << std::endl;
    std::cout << "  params  ---> list of Nomad parameters - some are mandatory"         << std::endl;
    std::cout << "               ex: param =[\"DIMENSION 5\", \"BB_OUTPUT_TYPE OBJ\"]"  << std::endl;
    std::cout << "               Consult Nomad documentation, or type PyNomad.help(),"  << std::endl;
    std::cout << "               to see all parameters"                                 << std::endl;
    std::cout                                                                           << std::endl;

    std::cout << " PyNomad.optimize output parameters:"                                 << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << "  x_best    ---> list of values for the best feasible or infeasible"  << std::endl;
    std::cout << "               points at the end of Nomad optimization"               << std::endl;
    std::cout << "  f_best    ---> Objective function value for the best point obtained"<< std::endl;
    std::cout << "  h_best    ---> Infeasibility measure value for best point obtained" << std::endl;
    std::cout << "               (0 if feasible)"                                       << std::endl;
    std::cout << "  nb_evals   --> Number of blackbox evaluations"                      << std::endl;
    std::cout << "  nb_iters   --> * Currently not supported *"                         << std::endl;
    std::cout << "               (would be: Number of iterations of the Mads algorithm)"<< std::endl;
    std::cout << "  run_flag   --> Run flag for Nomad termination (see details below) " << std::endl;
    std::cout << "  stop_reason -> Nomad termination stop reason "                      << std::endl;
    std::cout                                                                           << std::endl;

    std::cout << "-----------------------------------------------------------"          << std::endl;
    std::cout << " Blackbox definition examples"                                        << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << " Form 1: A single point is passed to blackbox function"               << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << " def bb(x):"                                                          << std::endl;
    std::cout << "     # size() gets the number of variables"                           << std::endl;
    std::cout << "     dim = x.size()"                                                  << std::endl;
    std::cout << "     # get_coord(i) accesses the variable i"                          << std::endl;
    std::cout << "     # Compute objective f"                                           << std::endl;
    std::cout << "     f = sum([x.get_coord(i)**2 for i in range(dim)])"                << std::endl;
    std::cout << "     # Set output on x"                                               << std::endl;
    std::cout << "     # Converting to bytes is necessary"                              << std::endl;
    std::cout << "     x.setBBO(str(f).encode(\"UTF-8\"))"                              << std::endl;
    std::cout << "     # return 1 if evaluation is successful, 0 if it failed"          << std::endl;
    std::cout << "     return 1"                                                        << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << " Form 2: A block (list) of points is passed to blackbox function"     << std::endl;
    std::cout << "         This form is used when parameter BB_MAX_BLOCK_SIZE is"       << std::endl;
    std::cout << "         greater than 1. If BB_MAX_BLOCK_SIZE is equal to 1, or"      << std::endl;
    std::cout << "         not specified, then form 1 (above) is used."                 << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << " def bb_block(block):"                                                << std::endl;
    std::cout << "     # size() gets the number of points to evaluate"                  << std::endl;
    std::cout << "     nbPoints = block.size()"                                         << std::endl;
    std::cout << "     # evalOk is a list of booleans"                                  << std::endl;
    std::cout << "     evalOk = [False for i in range(nbPoints)]"                       << std::endl;
    std::cout << "     for k in range(nbPoints):"                                       << std::endl;
    std::cout << "         # eval each point"                                           << std::endl;
    std::cout << "         x = block.get_x(k)"                                          << std::endl;
    std::cout << "         dim = x.size()"                                              << std::endl;
    std::cout << "         f = sum([x.get_coord(i)**2 for i in range(dim)])"            << std::endl;
    std::cout << "         x.setBBO(str(f).encode(\"UTF-8\"))"                          << std::endl;
    std::cout << "         evalOk[k] = True"                                            << std::endl;
    std::cout << "     # return a list where 1 is success, 0 is a failed evaluation"    << std::endl;
    std::cout << "     return evalOk"                                                   << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << "-----------------------------------------------------------"          << std::endl;
    
    std::cout << "-----------------------------------------------------------"          << std::endl;
    std::cout << " Nomad termination run flags"                                         << std::endl;
    std::cout                                                                           << std::endl;
    std::cout << "   1 - Objective target reached OR Mads converged (mesh criterion) "  << std::endl;
    std::cout << "       to a feasible point (true problem)."                           << std::endl;
    std::cout << "   0 - At least one feasible point obtained and evaluation budget "   << std::endl;
    std::cout << "       (single bb or block of bb) spent or max iteration (user "      << std::endl;
    std::cout << "       option) reached."                                              << std::endl;
    std::cout << "  -1 - Mads mesh converged but no feasible point obtained (only      "<< std::endl;
    std::cout << "       infeasible) for the true problem."                             << std::endl;
    std::cout << "  -2 - No feasible point obtained (only infeasible) and evaluation"   << std::endl;
    std::cout << "       budget (single bb or block of bb) spent or max iteration (user"<< std::endl;
    std::cout << "       option) reached"                                               << std::endl;
    std::cout << "  -3 - Initial point failed to evaluate"                              << std::endl;
    std::cout << "  -4 - Time limit reached (user option)"                              << std::endl;
    std::cout << "  -5 - CTRL-C or user stopped (callback function)"                    << std::endl;
    std::cout << "  -6 - Stop on feasible point (user option)"                          << std::endl;
    
}


static void printPyNomadInfo()
{
    std::cout << "-----------------------------------------------------------"          << std::endl;
    std::cout << " NOMAD: Nonlinear Optimization using the MADS Algorithm"              << std::endl;
    std::cout << " Released under the GNU Lesser General Public License:"               << std::endl;
    std::cout << " http://www.gnu.org/copyleft/lesser.html"                             << std::endl;
    //std::cout << " Source available from: https://www.gerad.ca/nomad/" << std::endl;

    printPyNomadVersion();

    std::cout << " NOMAD solves a Global MINLP/NLP in the form  " << std::endl;
    std::cout << "          min f(x)                            " << std::endl;
    std::cout << "          subject to:                         " << std::endl;
    std::cout << "                   nlcon(x) <= 0              " << std::endl;
    std::cout << "                   lb <= x <= ub              " << std::endl;
    std::cout << "                   x in R                     " << std::endl;

    printPyNomadUsage();

}


static void printNomadHelp(std::string about)
{
    NOMAD::AllParameters allParameters;
    NOMAD::toupper(about);
    allParameters.displayHelp(about, false, std::cout);
}


//Python Evaluator Class
class PyEval : public NOMAD::Evaluator
{
private:
    Callback  _cb;
    CallbackL _cbL;
    void*     _apply;

public:
    // Constructor
    PyEval(std::shared_ptr<NOMAD::EvalParameters> evalParams,
           Callback cb,
           CallbackL cbL,
           void * apply,
           NOMAD::EvalType evalType)
      : NOMAD::Evaluator(evalParams, evalType),
        _cb(cb),
        _cbL(cbL),
        _apply(apply)
    {
    }

    //Destructor
    ~PyEval() {}


    std::vector<bool> eval_block(NOMAD::Block &block,
                        const NOMAD::Double& NOMAD_UNUSED(hMax),
                        std::vector<bool> &countEval) const override
    {
        size_t nbPoints = block.size();
        std::vector<bool> evalOk(nbPoints, false);
        countEval.resize(nbPoints, false);
        
        // eval_block is always called.
        // if cbL is NULL, this means that the block must be of size 1, and that
        // cb should be used.
        if (nullptr == _cbL)
        {
            NOMAD::EvalPointPtr x_ptr = block[0];
            PyGILState_STATE state = PyGILState_Ensure();
            evalOk[0] = _cb(_apply, x_ptr);
            PyGILState_Release(state);
            countEval[0] = true;   // Always true in Python
        }
        else
        {
            // Call Python blackbox function on a vector of EvalPoints
            try
            {
                // Call callback
                std::shared_ptr<NOMAD::Block> block_ptr = std::make_shared<NOMAD::Block>(block);
                PyGILState_STATE state = PyGILState_Ensure();
                std::vector<int> cblret = _cbL(_apply, block_ptr);
                PyGILState_Release(state);
                block = *block_ptr;
                for (size_t i = 0; i < nbPoints; i++)
                {
                    evalOk[i]    = cblret[i];
                    countEval[i] = true;   // Always true in Python
                }
            }
            //If these errors occur, it is due to errors in python code
            catch(...)
            {
                printf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
                //Force exit
                raise(SIGINT);
                return evalOk;
            }
        }
        return evalOk;
    }
};


// Helper function for runNomad
static void initAllParams(std::shared_ptr<NOMAD::AllParameters> allParams,
                          std::vector<double> X0,
                          std::vector<double> LB,
                          std::vector<double> UB,
                          const std::vector<std::string> &params)
{
    size_t dimension;
    bool dimensionDefined = false;

    try
    {
        size_t ndec = X0.size();
        if ( ndec != 0 )
        {
            NOMAD::Point px0(ndec);

            dimension = ndec;
            dimensionDefined = true;

            for (size_t i = 0; i < ndec; i++)
            {
                px0[i] = X0[i];
            }

            allParams->setAttributeValue("X0", px0);
        }

        size_t nlb = LB.size();
        if ( nlb != 0 )
        {
            NOMAD::ArrayOfDouble plb(nlb);

            if ( ! dimensionDefined )
            {
                dimension = nlb;
                dimensionDefined = true;
            }

            if ( ndec != 0 && nlb !=ndec )
            {
                throw NOMAD::Exception("",0,"The lower bound size is inconsistent with X0 size");
            }


            for (size_t i = 0; i < nlb; i++)
            {
                plb[i] = LB[i];
            }

            allParams->setAttributeValue("LOWER_BOUND", plb);
        }

        size_t nub = UB.size();
        if ( nub != 0 )
        {
            NOMAD::ArrayOfDouble pub(nub);

            if ( ! dimensionDefined )
            {
                dimension = nub;
                dimensionDefined = true;
            }

            if ( nlb != 0 && nub != nlb )
            {
                throw NOMAD::Exception("",0,"The upper bound size is inconsistent with lower bound size");
            }

            if ( ndec != 0 && nub != ndec )
            {
                throw NOMAD::Exception("",0,"The upper bound size is inconsistent with X0 size");
            }

            for (size_t i = 0; i < nub; i++)
            {
                pub[i] = UB[i];
            }

            allParams->setAttributeValue("UPPER_BOUND", pub);

        }

        if ( dimensionDefined )
        {
            allParams->setAttributeValue("DIMENSION", dimension);
        }

        // The seed will always be to its default value
        NOMAD::RNG::resetPrivateSeedToDefault();

        for (size_t i = 0; i < params.size(); i++)
        {
            // Elements of the params array look like this:
            // "PARAM_NAME VALUE"
            // We can interpret them with readParamLine.
            allParams->readParamLine(params[i]);
        }

        allParams->checkAndComply();
    }

    catch(std::exception &e)
    {
        printf("NOMAD Parameter Error:\n%s\n",e.what());
    }

}


static int runNomad(Callback cb,
                    CallbackL cbL,
                    void * apply,
                    std::vector<double> X0,
                    std::vector<double> LB,
                    std::vector<double> UB,
                    const std::vector<std::string> &params,
                    std::shared_ptr<NOMAD::EvalPoint>& bestFeasSol,
                    std::shared_ptr<NOMAD::EvalPoint>& bestInfeasSol,
                    size_t &nbEvals,
                    size_t &nbIters,
                    std::string & stopReason)
{
    auto allParams = std::make_shared<NOMAD::AllParameters>();
    initAllParams(allParams, X0, LB, UB, params);
    bestFeasSol = nullptr;
    bestInfeasSol = nullptr;
    int runFlag = -3;
    stopReason ="No stop reason";

    try
    {
        
        
        Py_BEGIN_ALLOW_THREADS

        NOMAD::MainStep TheMainStep;
        TheMainStep.setAllParameters(allParams);

        // Set cbL to NULL if blocks are not used.
        // Set cb to NULL if blocks are used.
        if (allParams->getEvaluatorControlGlobalParams()->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE") > 1)
        {
            // Using blocks
            cb = nullptr;
        }
        else
        {
            cbL = nullptr;
        }
        
        auto ev = std::make_unique<PyEval>(allParams->getEvalParams(), cb, cbL, apply, NOMAD::EvalType::BB);
        TheMainStep.addEvaluator(std::move(ev));

        TheMainStep.start();
        TheMainStep.run();
        TheMainStep.end();

        nbEvals = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
        nbIters = 0; // Not supported in this version of NOMAD 4
        // Keeping the value for compatibility with PyNomad 1

        const auto hNormType = allParams->getAttributeValue<NOMAD::HNormType>("H_NORM");
        const auto evalType = NOMAD::EvalType::BB;
        const NOMAD::FHComputeType computeType= {evalType, {NOMAD::ComputeType::STANDARD, hNormType}};
        
        // Set the best feasible and best infeasible solutions
        std::vector<NOMAD::EvalPoint> evalPointFeasList, evalPointInfList;
        auto nbFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointFeasList, NOMAD::Point(), computeType);
        auto nbInf  = NOMAD::CacheBase::getInstance()->findBestInf(evalPointInfList, NOMAD::INF, NOMAD::Point(), computeType);
        
        // For now
        // If nbFeas > 0 we return a single best feasible point (no infeasible point)
        // Else (if nbFeas == 0) we return the least infeasible point with the smallest f (index 0, see findBestInf)
        // Maybe this could be generalized to show the best feasible point and all undominated infeasible points.
        // The same logic for Nomad C and Matlab interfaces and for PyNomad.
        if (nbFeas > 0)
        {
            NOMAD::EvalPoint evalPointFeas = evalPointFeasList[0];
            bestFeasSol = std::make_shared<NOMAD::EvalPoint>(evalPointFeas);
            bestInfeasSol = nullptr;
        }
        else if (0 == nbFeas)
        {
            bestFeasSol = nullptr;
            if (nbInf > 0)
            {
                // One of the pointer is set to null to identify the type of solution
                NOMAD::EvalPoint evalPointInf = evalPointInfList[0];
                bestInfeasSol = std::make_shared<NOMAD::EvalPoint>(evalPointInf);
            }
            else if(0 == nbInf)
            {
                bestInfeasSol = nullptr;
            }
        }
        
        runFlag = TheMainStep.getRunFlag();
        
        stopReason = TheMainStep.getAllStopReasons()->getStopReasonAsString();
        
        NOMAD::MainStep::resetComponentsBetweenOptimization();
        
        Py_END_ALLOW_THREADS
        return runFlag;
    }

    catch(std::exception &e)
    {
        printf("NOMAD exception (report to developer):\n%s\n",e.what());
    }

    return runFlag; // Default is for Nomad error

}

static int runNomad(Callback cb,
                    CallbackL cbL,
                    void * applyBB,
                    void * applySurrogate,
                    std::vector<double> X0,
                    std::vector<double> LB,
                    std::vector<double> UB,
                    const std::vector<std::string> &params,
                    std::shared_ptr<NOMAD::EvalPoint>& bestFeasSol,
                    std::shared_ptr<NOMAD::EvalPoint>& bestInfeasSol,
                    size_t &nbEvals,
                    size_t &nbIters,
                    std::string & stopReason)
{
    auto allParams = std::make_shared<NOMAD::AllParameters>();
    initAllParams(allParams, X0, LB, UB, params);
    bestFeasSol = nullptr;
    bestInfeasSol = nullptr;
    stopReason = "No stop reason";
    
    int runFlag = -3 ;

    std::cout<<"Run nomad with surrogate"<<std::endl;
    try
    {
        Py_BEGIN_ALLOW_THREADS

        NOMAD::MainStep TheMainStep;
        TheMainStep.setAllParameters(allParams);

        // Set cbL to NULL if blocks are not used.
        // Set cb to NULL if blocks are used.
        if (allParams->getEvaluatorControlGlobalParams()->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE") > 1)
        {
            // Using blocks
            cb = nullptr;
        }
        else
        {
            cbL = nullptr;
        }
        
        auto evBB = std::make_unique<PyEval>(allParams->getEvalParams(), cb, cbL, applyBB, NOMAD::EvalType::BB);
        TheMainStep.addEvaluator(std::move(evBB));

        auto evSurrogate = std::make_unique<PyEval>(allParams->getEvalParams(), cb, cbL, applySurrogate, NOMAD::EvalType::SURROGATE);
        TheMainStep.addEvaluator(std::move(evSurrogate));
        
        TheMainStep.start();
        TheMainStep.run();
        TheMainStep.end();

        nbEvals = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
        nbIters = 0; // Not supported in this version of NOMAD 4
                     // Keeping the value for compatibility with PyNomad 1

        const auto hNormType = allParams->getAttributeValue<NOMAD::HNormType>("H_NORM");
        const auto evalType = NOMAD::EvalType::BB;
        const NOMAD::FHComputeType computeType= {evalType, {NOMAD::ComputeType::STANDARD, hNormType}};

        // Set the best feasible and best infeasible solutions
        std::vector<NOMAD::EvalPoint> evalPointFeasList, evalPointInfList;
        auto nbFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointFeasList, NOMAD::Point(), computeType);
        auto nbInf  = NOMAD::CacheBase::getInstance()->findBestInf(evalPointInfList, NOMAD::INF, NOMAD::Point(), computeType);

        // For now
        // If nbFeas > 0 we return a single best feasible point (no infeasible point)
        // Else (if nbFeas == 0) we return the least infeasible point with the smallest f (index 0, see findBestInf)
        // Maybe this could be generalized to show the best feasible point and all undominated infeasible points.
        // The same logic for Nomad C and Matlab interfaces and for PyNomad.
        if (nbFeas > 0)
        {
            NOMAD::EvalPoint evalPointFeas = evalPointFeasList[0];
            bestFeasSol = std::make_shared<NOMAD::EvalPoint>(evalPointFeas);
            bestInfeasSol = nullptr;
        }
        else if (0 == nbFeas)
        {
            bestFeasSol = nullptr;
            if (nbInf > 0)
            {
                // One of the pointer is set to null to identify the type of solution
                NOMAD::EvalPoint evalPointInf = evalPointInfList[0];
                bestInfeasSol = std::make_shared<NOMAD::EvalPoint>(evalPointInf);
            }
            else if(0 == nbInf)
            {
                bestInfeasSol = nullptr;
            }
        }

        runFlag = TheMainStep.getRunFlag();
        
        stopReason = TheMainStep.getAllStopReasons()->getStopReasonAsString();
        
        NOMAD::MainStep::resetComponentsBetweenOptimization();

        Py_END_ALLOW_THREADS
        return runFlag;
    }

    catch(std::exception &e)
    {
        printf("NOMAD exception (report to developper):\n%s\n",e.what());
    }

    return runFlag;

}

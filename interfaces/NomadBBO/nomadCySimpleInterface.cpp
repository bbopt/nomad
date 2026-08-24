// June 2019
// Version 1.0 is with NOMAD 3.
#define NOMAD_PYTHON_VERSION "2.2"

#include "Algos/EvcInterface.hpp"
#include "Algos/Mads/MadsMegaIteration.hpp"
#include "Algos/CatMads/CatMads.hpp"
#include "Algos/CatAds/CatAds.hpp"
#include "Algos/CatMads/CatCustomOrder.hpp"

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

// The boolean return is for stopping (MegaIter or else)
typedef bool (*CallbackU)(void * apply,
                          std::shared_ptr<NOMAD::Block> block);

// These functions typedef must be the same as the one defined in PyNomad.pyx
typedef std::vector<std::vector<double>> (*CustomOrderCompPYXFunc)(std::vector<std::vector<double>>& trialPoints) ;

// For CustomOrder completeTrialPointsInformation function
CustomOrderCompPYXFunc __customOrderCompPYXFunc = [](std::vector<std::vector<double>>& trialPoints){return std::vector<std::vector<double>>();};
NOMAD::CustomOrderCompWrapper customOrderCompWrapper = [](std::vector<std::vector<double>> & trialPoints) -> std::vector<std::vector<double>>
{
    // This wrapper function calls a cython defined function (from PyNomad.pyx) (trial points vectors conversion into numpy array is done by cython)
    // The wrapper is important to ensure gil state/release while calling python

    PyGILState_STATE gstate = PyGILState_Ensure();
    std::vector<std::vector<double>> trialPointsModelEval = __customOrderCompPYXFunc(trialPoints);
 
   // Display the result
   std::cout << "Result from CustomOrderCompWrapper: " << std::endl;
   for (size_t i = 0; i < trialPointsModelEval.size(); i++)
   {
         std::cout << "Trial point " << i << ": ";
         for (size_t j = 0; j < trialPointsModelEval[i].size(); j++)
         {
                std::cout << trialPointsModelEval[i][j] << " ";
         }
         std::cout << std::endl;
   }

    PyGILState_Release(gstate);
    return trialPointsModelEval;
};


// For updating Model cache
typedef bool (*ModelCacheUpdatePYXFunc)(std::vector<double>& trialPoint) ;
ModelCacheUpdatePYXFunc __modelCacheUpdatePYXFunc = [](std::vector<double>& ep){return false;};
bool modelCacheUpdateWrapper(NOMAD::EvalPointPtr evalPoint)
{
    bool success = false;

    // This wrapper function calls a cython defined function (from PyNomad.pyx) (trial point vector conversion into numpy array is done by cython)
    // The wrapper is important to ensure gil state/release while calling python

    // Convert eval point into vector of double
    auto eval = evalPoint->getEval(NOMAD::EvalType::BB);
    if (nullptr != eval)
    {
        std::vector<double> ep;
        const auto x = evalPoint->getX();
        // std::cout << "Updating Model cache for point: " << *x << ", ";

        for (size_t i =0; i < x->size(); i++ )
        {
            ep.push_back((*x)[i].todouble());
        }
        auto bbo = eval->getBBOutput().getBBOAsArrayOfDouble();
        // std::cout << " BBO=(" << eval->getBBO() << ")" << std::endl;
        for (size_t i =0 ; i < bbo.size(); i++)
        {
            // Do not update Model cache if evaluation has no BBO values
            if ( bbo[i] == NOMAD::INF )
            {
                //std::cout << "Can not update model with this point " << std::endl;
                return false;
            }

            ep.push_back(bbo[i].todouble());
        }

        PyGILState_STATE gstate = PyGILState_Ensure();
        success = __modelCacheUpdatePYXFunc(ep);
        PyGILState_Release(gstate);

    }

    return success;
};

// For CatMads/CatAds poll direction generation
typedef std::vector<std::vector<double>> (*CatPollFreePYXFunc)(std::vector<double>& frameCenter) ;
CatPollFreePYXFunc __catPollFreePYXFunc = [](std::vector<double>& ep){return std::vector<std::vector<double>>() ;};
std::vector<NOMAD::Direction> catPollFreeWrapper(NOMAD::EvalPointPtr evalPoint)
{
    std::vector<NOMAD::Direction> directions;

    // This wrapper function calls a cython defined function (from PyNomad.pyx) (direction vector conversion into numpy array is done by cython)
    // The wrapper is important to ensure gil state/release while calling python

    // Convert eval point into vector of double
    auto eval = evalPoint->getEval(NOMAD::EvalType::BB);
    if (nullptr != eval)
    {
        std::vector<double> ep;
        const auto x = evalPoint->getX();
        if (! x->NOMAD::ArrayOfDouble::isDefined())
        {
            throw NOMAD::Exception("PyNomad",0,"Error in CatMads/CatAds poll callback: frame center point is not defined.");
        }

        for (size_t i =0; i < x->size(); i++ )
        {
            ep.push_back((*x)[i].todouble());
        }

        PyGILState_STATE gstate = PyGILState_Ensure();
        std::vector<std::vector<double>> ve = __catPollFreePYXFunc(ep);

        // Convert to Direction
        for (const auto & v : ve)
        {
            NOMAD::Point p(v);
            directions.push_back(NOMAD::Direction(p));
        }

        PyGILState_Release(gstate);

    }

    return directions;
};

// Custom search. Can be used for CatMads/CatAds poll direction generation or as a Mads user-defined search method. 
typedef std::vector<std::vector<double>> (*CustomSearchPYXFunc)(std::vector<double>& vFeas, std::vector<double>& vInfeas) ;
CustomSearchPYXFunc __customSearchPYXFunc = [](std::vector<double>& vFeas, std::vector<double>& vInfeas){return std::vector<std::vector<double>>() ;};
std::vector<NOMAD::EvalPoint> customSearchWrapper(NOMAD::EvalPointPtr feasEvalPoint, NOMAD::EvalPointPtr infeasEvalPoint )
{

    // This wrapper function calls a cython defined function (from PyNomad.pyx) (direction vector conversion into numpy array is done by cython)
    // The wrapper is important to ensure gil state/release while calling python

    std::vector<NOMAD::EvalPoint> trialPoints;

    // Convert eval point into vector of double
    if (feasEvalPoint == nullptr && infeasEvalPoint == nullptr)
    {
        return trialPoints; // No information to send to python, return empty trial points list
    }


    NOMAD::Eval* feasEval = nullptr;
    NOMAD::Eval* infeasEval = nullptr;
    if (feasEvalPoint != nullptr)
    {
        feasEval = feasEvalPoint->getEval(NOMAD::EvalType::BB);
    }

    std::vector<double> epFeas;
    if (nullptr != feasEval)
    {
        const auto x = feasEvalPoint->getX();
        if (! x->NOMAD::ArrayOfDouble::isDefined())
        {
            throw NOMAD::Exception("PyNomad",0,"Error in custom search callback: frame center point is not defined.");
        }
        for (size_t i =0; i < x->size(); i++ )
        {
            epFeas.push_back((*x)[i].todouble());
        }
    }

    if (infeasEvalPoint != nullptr)
    {
        infeasEval = infeasEvalPoint->getEval(NOMAD::EvalType::BB);
    }

    std::vector<double> epInfeas;
    if (nullptr != infeasEval)
    {
        const auto xInfeas = infeasEvalPoint->getX();
        if (! xInfeas->NOMAD::ArrayOfDouble::isDefined())
        {
            throw NOMAD::Exception("PyNomad",0,"Error in custom search callback: frame center point is not defined.");
        }

        for (size_t i =0; i < xInfeas->size(); i++ )
        {
            epInfeas.push_back((*xInfeas)[i].todouble());
        }
    }

    PyGILState_STATE gstate = PyGILState_Ensure();

    std::vector<std::vector<double>> tps = __customSearchPYXFunc(epFeas, epInfeas);

    PyGILState_Release(gstate);

    // Convert to trial points
    for (const auto & tp : tps)
    {
        NOMAD::Point p(tp);
        trialPoints.push_back(NOMAD::EvalPoint(p));
    }

    return trialPoints;
};

// For building python Model from cache after DOE
typedef bool (*ModelBuildPYXFunc)() ;
ModelBuildPYXFunc __modelBuildPYXFunc = [](){return false;};
bool modelBuildWrapper()
{
    bool success = false;

    PyGILState_STATE gstate = PyGILState_Ensure();
    success = __modelBuildPYXFunc();

    PyGILState_Release(gstate);

    return success;
};

bool madsBuildModelCustomSearchCallbackFn(const NOMAD::Step& step,
                          NOMAD::EvalPointSet & trialPoints)
{
    // Important: by default USER_CALLS are disabled when doing quad model optimization
    // -> NO call to this function when doing quad model search.

    // Several NOMAD::Algorithm are used by NOMAD.
    // For now, we do search only on Mads.

    auto mads = dynamic_cast<const NOMAD::Mads*>(step.getRootAlgorithm());
    if (nullptr == mads)
    {
        return false; // Not a Mads algorithm, do not call user-defined search method
    }

    auto barrier = mads->getMegaIterationBarrier();
    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"No barrier available.");
    }
    NOMAD::EvalPointPtr feasEvalPoint = barrier->getCurrentIncumbentFeas();
    NOMAD::EvalPointPtr infeasEvalPoint = barrier->getCurrentIncumbentInf();

    bool modelBuilt = modelBuildWrapper(); // Build model if available before calling user-defined search method

    if (!modelBuilt)
    {
        return false; // No model built, do not call user-defined search method
    }

    std::vector<NOMAD::EvalPoint> trialPointsVec = customSearchWrapper(feasEvalPoint, infeasEvalPoint);

    for (const auto & tp : trialPointsVec)
    {
        trialPoints.insert(tp);
    }
    if (trialPoints.empty())
    {
        return false; // No trial points generated by user-defined search method
    }

    return true;
};


static void printNomadHelp(std::string about)
{
    NOMAD::AllParameters allParameters;
    NOMAD::toupper(about);
    allParameters.displayHelp(about, false, std::cout);
}

// NOTE: This is not available in the NomadBBO interface
// Let's keep it until we have a callback passed to runNomad for that.
// This is for MEGA_ITER_END callback function.
// Maybe set or not.
CallbackU  _megaIterEndCb;
void*     _megaIterEndApply;
static void setCustomMegaIterEndCallbackFunction(CallbackU cbU,
                                                void * apply)
{
    _megaIterEndCb = cbU;
    _megaIterEndApply = apply;
}
/*---------------------------------------------------------------*/
/* After each mega iteration verify if stop (may not be used).   */
/*---------------------------------------------------------------*/
void customMegaIterEndCbFunc(const NOMAD::Step& step,
                  bool &stop)
{
        // Important: by default USER_CALLS are disabled when doing quad model optimization
        // -> NO call to this function when doing quad model search.

        // Several NOMAD::Algorithm are used by NOMAD.
        // We are interested only on the main Mads (Mega) Iteration.
        // Use a dynamic cast to make sure with have the Mads (Mega) Iteration.
        auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>(&step);

        if (nullptr != megaIter)
        {

            // Set the best feasible solutions
            // A single best feasible solution should be sufficient. Let's pass all of them
            // in case, the user wants access to all best feasible points.
            std::vector<NOMAD::EvalPoint> evalPointFeasList;
            auto nbFeas = NOMAD::CacheBase::getInstance()->findBestFeas(evalPointFeasList, NOMAD::Point(), NOMAD::defaultFHComputeType /* for BB and default compute type */);

            if (nbFeas == 0)
            {
                stop = false;
                return;
            }

            try
            {
                if (_megaIterEndCb)
                {

                    NOMAD::Block block;
                    for (const auto & ep: evalPointFeasList)
                    {
                        block.push_back(std::make_shared<NOMAD::EvalPoint>(ep));
                    }
                    std::shared_ptr<NOMAD::Block> block_ptr = std::make_shared<NOMAD::Block>(block);
                    PyGILState_STATE state = PyGILState_Ensure();
                    stop = _megaIterEndCb(_megaIterEndApply, block_ptr);
                    PyGILState_Release(state);
                }
            }
            //If these errors occur, it is due to errors in python code
            catch(...)
            {
                printf("Unrecoverable error in User Mega Iteration Callback, Exiting NOMAD...\n\n");
                //Force exit
                raise(SIGINT);
            }
        }
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
                        [[maybe_unused]] const NOMAD::Double& hMax,
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

            // Update CatMads model if available
            if (evalOk[0])
            {
                bool catModelUpdateSuccess = modelCacheUpdateWrapper(x_ptr);
            }
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

                    // Update CatMads model if available
                    if (evalOk[i])
                    {
                        bool catModelUpdateSuccess = modelCacheUpdateWrapper(block[i]);
                    }
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
            // "PARAM_NAME VALUE [VALUE2 ...]".
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
                    CustomOrderCompPYXFunc customOrderCompPYXFunc,
                    ModelCacheUpdatePYXFunc modelCacheUpdatePYXFunc,
                    ModelBuildPYXFunc modelBuildPYXFunc,
                    CatPollFreePYXFunc catPollFreePYXFunc,
                    CustomSearchPYXFunc customSearchPYXFunc,
                    const std::vector<std::string> &params,
                    std::shared_ptr<NOMAD::Block>& bestIncFeasSol,
                    std::shared_ptr<NOMAD::Block>& bestIncInfeasSol,
                    size_t &nbEvals,
                    size_t &nbIters,
                    std::string & stopReason)
{

    auto allParams = std::make_shared<NOMAD::AllParameters>();

    std::vector<double> X0, LB, UB; // Undefined vectors. All infos are in params for mixed variables
    initAllParams(allParams, X0, LB, UB, params);
    bestIncFeasSol = nullptr;
    bestIncInfeasSol = nullptr;
    int runFlag = -3;
    stopReason ="No stop reason";

    Py_Initialize();

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

        // Need to start the main step to have an evaluator control ready.
        TheMainStep.start();

        // Callbacks for both CatMads and regular Mads.
        __modelCacheUpdatePYXFunc = modelCacheUpdatePYXFunc;
        __customSearchPYXFunc = customSearchPYXFunc; 
        __modelBuildPYXFunc = modelBuildPYXFunc; 
        
        bool isCatOptim = (allParams->getAttributeValue<bool>("CATMADS_OPTIMIZATION") || allParams->getAttributeValue<bool>("CATADS_OPTIMIZATION") );
        if (isCatOptim)
        {

            try {

                // Callbacks for CatMads/CatAds
                __customOrderCompPYXFunc = customOrderCompPYXFunc;
                __catPollFreePYXFunc = catPollFreePYXFunc;



                // Set the custom comparison function for categorical variable order in CatMads/CatAds.
                auto bbot = allParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
                auto customOrder = std::make_shared<NOMAD::CatCustomOrder>(customOrderCompWrapper, bbot);

                NOMAD::EvcInterface::getEvaluatorControl()->setCatCompMethod(customOrder);

                // Pass the method to build the model to CatMads/CatAds.
                auto catAlgo =  std::dynamic_pointer_cast<NOMAD::CatAlgoUtils>(TheMainStep.getAlgo(NOMAD::StepType::ALGORITHM_CATALGO));
                if ( !catAlgo)
                {
                    throw NOMAD::Exception("NomadBBO",0,"Cannot access CatMads/CatAds algorithm.");
                }
                
                bool success = catAlgo->setBuildModelCallback(modelBuildWrapper);
                if (!success)
                    throw NOMAD::Exception("NomadBBO",0,"Cannot build the model.");

                success = catAlgo->setPollFreeCallback(catPollFreeWrapper);
                if (!success)
                    throw NOMAD::Exception("NomadBBO",0,"Cannot set the CatMads/CatAds poll callback.");

                success = catAlgo->setCategoricalSearchCallback(customSearchWrapper);
                if (!success)
                    throw NOMAD::Exception("NomadBBO",0,"Cannot set the custom search callback.");

            }
            catch (const std::exception& e)
            {
                std::cerr << "Error setting CatMads/CatAds: " << e.what() << std::endl;
                throw NOMAD::Exception("NomadBBO",0,"Cannot continue.");
            }
        }
        else
        {
            // Regular Mads optimization
            
            // Registering the callback function to generate Mads BO search trial points
             // USER_SEARCH is always enabled but a flag controlling the use of this callback is managed in the python code (optimizer).
            TheMainStep.addCallback<NOMAD::MadsCallbackType::USER_METHOD_SEARCH>(madsBuildModelCustomSearchCallbackFn);
        }

        // NOTE: This is not available in the NomadBBO interface
        // Let's keep it until we have a callback passed to runNomad for that.
        // // Link MEGA_ITER_END callback function with user python function defined locally
        // if (_megaIterEndCb)
        // {
        //     TheMainStep.addCallback<NOMAD::AlgoCallbackType::MEGA_ITERATION_END>(customMegaIterEndCbFunc);
        // }



        TheMainStep.run();

        // Get the points before ending the main step, as end() will reset some components of NOMAD.
        //
        // For single objective optimization, we have usually single feas/infeas incumbent points.
        // For multi-objective optimization, we have several feasible incumbent points but infeasible incumbent are not supported.
        const auto feasIncPoints = TheMainStep.getBarrierIncumbentPoints(true);
        const auto infeasIncPoints = TheMainStep.getBarrierIncumbentPoints(false);

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

        // Convert to list of coordinates for Python interface
        if (!feasIncPoints.empty())
        {
            bestIncFeasSol = std::make_shared<NOMAD::Block>();
            bestIncFeasSol->reserve(feasIncPoints.size());
            for (const auto &ep : feasIncPoints)
            {
                bestIncFeasSol->push_back(std::make_shared<NOMAD::EvalPoint>(ep));
            }
        }
        else
        {
            bestIncFeasSol = nullptr;
        }

        if (!infeasIncPoints.empty())
        {
            bestIncInfeasSol = std::make_shared<NOMAD::Block>();
            bestIncInfeasSol->reserve(infeasIncPoints.size());
            for (const auto &ep : infeasIncPoints)
            {
                bestIncInfeasSol->push_back(std::make_shared<NOMAD::EvalPoint>(ep));
            }
        }
        else
        {
            bestIncInfeasSol = nullptr;
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

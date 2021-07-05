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

// Generic
#include "../Algos/Algorithm.hpp"
#include "../Algos/EvcInterface.hpp"
#include "../Algos/Iteration.hpp"
#include "../Algos/MegaIteration.hpp"
#include "../Algos/Step.hpp"
#include "../Cache/CacheBase.hpp"
#include "../Output/OutputQueue.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
bool NOMAD::Step::_userInterrupt = false;
bool NOMAD::Step::_userTerminate = false;

NOMAD::StepEndCbFunc NOMAD::Step::_cbIterationEnd = defaultStepEnd;
NOMAD::StepEndCbFunc NOMAD::Step::_cbMegaIterationEnd = defaultStepEnd;
NOMAD::HotRestartCbFunc NOMAD::Step::_cbHotRestart = defaultHotRestart;

bool NOMAD::Step::_showWarnings = true;


/*---------------------------------------------------------*/
/*  user interruption (static, called by pressing ctrl-c)  */
/*---------------------------------------------------------*/
void NOMAD::Step::userInterrupt(int signalValue)
{
    std::cout << std::endl << "NOMAD caught User interruption." << std::endl;
    if (_userInterrupt)
    {
        // This is the 2nd CTRL-C. Terminate.
        std::cout << "Terminate NOMAD." << std::endl;
        setUserTerminate();
        throw NOMAD::UserTerminateException(__FILE__, __LINE__, "User termination");
    }
    else
    {
        std::cout << "Please wait..." << std::endl;
        // Some steps will check for _userInterrupt and then call
        // hotRestartOnUserInterrupt(). Here we are in a static method
        // so we cannot call it.
    }

    // Set this stop reason to be tested by EvaluatorControl
    NOMAD::AllStopReasons::set(NOMAD::BaseStopType::HOT_RESTART);

    NOMAD::Step::_userInterrupt = true;
}


void NOMAD::Step::debugSegFault(int signalValue)
{
    NOMAD::OutputQueue::Flush();
#ifdef _OPENMP
    #pragma omp critical
#endif
    std::cerr << "Caught seg fault in thread " << NOMAD::getThreadNum() << std::endl;
    throw NOMAD::Exception(__FILE__,__LINE__,"Caught seg fault");
}


void NOMAD::Step::init()
{
    if (nullptr != _parentStep)
    {
        // If the parent is ROOT the params will be null
        if (nullptr == _runParams)
        {
            _runParams = _parentStep->_runParams;
        }
        if (nullptr == _pbParams)
        {
            _pbParams = _parentStep->_pbParams;
        }
    }
}


NOMAD::Step::~Step()
{
    NOMAD::OutputQueue::Flush();
}


void NOMAD::Step::addCallback(const NOMAD::CallbackType& callbackType,
                              const NOMAD::StepEndCbFunc& stepEndCbFunc)
{
    switch (callbackType)
    {
        case NOMAD::CallbackType::ITERATION_END:
            _cbIterationEnd = stepEndCbFunc;
            break;
        case NOMAD::CallbackType::MEGA_ITERATION_END:
            _cbMegaIterationEnd = stepEndCbFunc;
            break;
        default:
            break;
    }
}


void NOMAD::Step::addCallback(const NOMAD::CallbackType& callbackType,
                              const NOMAD::HotRestartCbFunc& hotRestartCbFunc)
{
    if (NOMAD::CallbackType::HOT_RESTART == callbackType)
    {
        _cbHotRestart = hotRestartCbFunc;
    }
}


void NOMAD::Step::runCallback(NOMAD::CallbackType callbackType,
                              const NOMAD::Step& step,
                              bool &stop)
{
    switch(callbackType)
    {
        case NOMAD::CallbackType::ITERATION_END:
            _cbIterationEnd(step, stop);
            break;
        case NOMAD::CallbackType::MEGA_ITERATION_END:
            _cbMegaIterationEnd(step, stop);
            break;
        default:
            break;
    }
}


void NOMAD::Step::runCallback(NOMAD::CallbackType callbackType,
                              std::vector<std::string>& paramLines)
{
    if (NOMAD::CallbackType::HOT_RESTART == callbackType)
    {
        _cbHotRestart(paramLines);
    }
}


void NOMAD::Step::AddOutputInfo(const std::string& s, bool isBlockStart, bool isBlockEnd) const
{
    // NB. Set the output level as LEVEL_INFO by default.
    OUTPUT_INFO_START
    NOMAD::OutputInfo outputInfo(getName(), s, NOMAD::OutputLevel::LEVEL_INFO, isBlockStart, isBlockEnd);
    NOMAD::OutputQueue::Add(std::move(outputInfo));
    OUTPUT_INFO_END
}


void NOMAD::Step::AddOutputInfo(const std::string& s, NOMAD::OutputLevel outputLevel) const
{
    if (NOMAD::OutputQueue::GoodLevel(outputLevel))
    {
        NOMAD::OutputInfo outputInfo(getName(), s, outputLevel);
        NOMAD::OutputQueue::Add(std::move(outputInfo));
    }
}


void NOMAD::Step::AddOutputError(const std::string& s) const
{
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_ERROR);
}


void NOMAD::Step::AddOutputWarning(const std::string& s) const
{
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_WARNING);
}


void NOMAD::Step::AddOutputVeryHigh(const std::string& s) const
{
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_VERY_HIGH);
}


void NOMAD::Step::AddOutputHigh(const std::string& s) const
{
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_HIGH);
}


void NOMAD::Step::AddOutputDebug(const std::string& s) const
{
    OUTPUT_DEBUG_START
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END
}


void NOMAD::Step::AddOutputInfo(NOMAD::OutputInfo outputInfo) const
{
    NOMAD::OutputQueue::Add(std::move(outputInfo));
}


/// Implementation of virtual functions : default start
void NOMAD::Step::start()
{
    defaultStart();
    startImp();
}


void NOMAD::Step::end()
{
    endImp();
    defaultEnd();
}


bool NOMAD::Step::run()
{
    return runImp();
}


void NOMAD::Step::observe(const std::vector<NOMAD::EvalPoint>& evalPointList)
{
    // Should not be called if it is not reimplemented.
    throw NOMAD::StepException(__FILE__,__LINE__,"Observe is not implemented in step " + getName(), this);
}


void NOMAD::Step::defaultStart()
{
    // Test shared_ptr here because MainStep has no stopReason
    if ( _stopReasons && ! _stopReasons->checkTerminate() )
    {
        _stopReasons->setStarted();
    }

    AddOutputInfo("Start step " + getName() , true, false);
}


void NOMAD::Step::defaultEnd()
{
    AddOutputInfo("End step " + getName(), false, true);
    // Flush because the step is done.
    NOMAD::OutputQueue::Flush();
}


void NOMAD::Step::verifyParentNotNull()
{
    if (nullptr == _parentStep)
    {
        std::string err = "Parent step for \"" + getName() + "\" should not be NULL";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::Step::verifyGenerateAllPointsBeforeEval(const std::string& method, const bool expected) const
{
    bool actual = _runParams->getAttributeValue<bool>("MEGA_SEARCH_POLL");

    if (expected != actual)
    {
        std::string err = "Error: " + method + " should only be called if ";
        err += " parameter MEGA_SEARCH_POLL is ";
        err += (expected ? "true" : "false");
        throw NOMAD::StepException(__FILE__,__LINE__,err, this);
    }
}


bool NOMAD::Step::isAnAlgorithm() const
{
    NOMAD::Step* step = const_cast<Step*>(this);
    return (nullptr != dynamic_cast<NOMAD::Algorithm*>(step));
}


const NOMAD::Algorithm* NOMAD::Step::getRootAlgorithm() const
{
    auto algo = isAnAlgorithm() ? dynamic_cast<const NOMAD::Algorithm*>(this)
                                : getParentOfType<NOMAD::Algorithm*>();

    auto parentAlgo = algo->getParentOfType<NOMAD::Algorithm*>();
    while (nullptr != parentAlgo)
    {
        algo = parentAlgo;
        parentAlgo = algo->getParentOfType<NOMAD::Algorithm*>();
    }

    return algo;
}


const NOMAD::Algorithm* NOMAD::Step::getFirstAlgorithm() const
{
    auto algo = isAnAlgorithm() ? dynamic_cast<const NOMAD::Algorithm*>(this)
                                : getParentOfType<NOMAD::Algorithm*>();

    return algo;
}


std::string NOMAD::Step::getAlgoName() const
{
    std::string s = "";
    auto algo = getFirstAlgorithm();
    if (nullptr != algo)
    {
        s = algo->getName();
    }

    // Append a space for easiness of use
    if (!s.empty())
    {
        s += " ";
    }

    return s;
}


// Get MeshBase from the Iteration ancestor.
const std::shared_ptr<NOMAD::MeshBase> NOMAD::Step::getIterationMesh() const
{
    std::shared_ptr<NOMAD::MeshBase> mesh = nullptr;
    const NOMAD::Iteration* iteration = getParentOfType<NOMAD::Iteration*>();

    if (nullptr != iteration)
    {
        mesh = iteration->getMesh();
    }
    return mesh;
}


// Get Barrier from the MegaIteration ancestor.
const std::shared_ptr<NOMAD::Barrier> NOMAD::Step::getMegaIterationBarrier() const
{
    std::shared_ptr<NOMAD::Barrier> barrier = nullptr;
    NOMAD::MegaIteration* megaIter = nullptr;

    if (isAnAlgorithm())
    {
        // An Algorithm has its own MegaIteration member. Get it.
        auto algo = dynamic_cast<const NOMAD::Algorithm*>(this);
        megaIter = algo->getMegaIteration().get();
    }
    else
    {
        // Is current Step a MegaIteration?
        auto constMegaIter = dynamic_cast<const NOMAD::MegaIteration*>(this);
        if (nullptr == constMegaIter)
        {
            // Get first parent of type MegaIteration.
            constMegaIter = getParentOfType<NOMAD::MegaIteration*>();
        }
        megaIter = const_cast<NOMAD::MegaIteration*>(constMegaIter);
    }

    if (nullptr != megaIter)
    {
        barrier = megaIter->getBarrier();
    }

    return barrier;
}


bool NOMAD::Step::solHasFeas() const
{
    bool hasFeas = NOMAD::CacheBase::getInstance()->hasFeas(NOMAD::EvalType::BB);

    if (!hasFeas)
    {
        // No feasible point in cache, but possibly in MegaIteration ancestor's barrier.
        auto barrier = getMegaIterationBarrier();
        if (nullptr != barrier)
        {
            for (auto xFeas : barrier->getAllXFeas())
            {
                if (xFeas.isEvalOk(NOMAD::EvalType::BB) && xFeas.isFeasible(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD))
                {
                    hasFeas = true;
                    break;
                }
            }
        }
    }

    return hasFeas;
}


bool NOMAD::Step::hasPhaseOneSolution() const
{
    bool hasPhaseOneSol = false;

    // A phase one solution has a PHASE_ONE Eval with f = 0.
    std::vector<NOMAD::EvalPoint> evalPointList;
    NOMAD::CacheBase::getInstance()->find(NOMAD::EvalPoint::isPhaseOneSolution, evalPointList);

    // Points have to verify hMax.
    auto barrier = getMegaIterationBarrier();
    NOMAD::Double hMax = (nullptr != barrier) ? barrier->getHMax() : _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
    for (auto evalPoint : evalPointList)
    {
        NOMAD::Double h;
        if (NOMAD::EvalStatusType::EVAL_OK == evalPoint.getEvalStatus(NOMAD::EvalType::BB))
        {
            h = evalPoint.getH(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD);
        }
        if (h.isDefined() && h <= hMax)
        {
            hasPhaseOneSol = true;
            break;
        }
    }

    if (!hasPhaseOneSol)
    {
        // No feasible point in cache, but possibly in MegaIteration ancestor's barrier.
        if (nullptr != barrier)
        {
            auto xFeas = barrier->getFirstXFeas();
            if (nullptr != xFeas && NOMAD::EvalStatusType::EVAL_OK == xFeas->getEvalStatus(NOMAD::EvalType::BB))
            {
                NOMAD::Double h = xFeas->getH(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD);
                hasPhaseOneSol = NOMAD::EvalPoint::isPhaseOneSolution(*xFeas) && (h <= hMax);
            }
        }
    }

    return hasPhaseOneSol;
}


void NOMAD::Step::hotRestartOnUserInterrupt()
{
    if (   !_stopReasons->testIf(NOMAD::BaseStopType::HOT_RESTART)
        && _stopReasons->checkTerminate())
    {
        return;
    }
    hotRestartBeginHelper();

    hotRestartEndHelper();
}


void NOMAD::Step::hotRestartBeginHelper()
{
    if (nullptr != _runParams
        && !_runParams->getAttributeValue<bool>("HOT_RESTART_ON_USER_INTERRUPT"))
    {
        setUserTerminate();
        _stopReasons->set(NOMAD::BaseStopType::CTRL_C);
    }
}


void NOMAD::Step::hotRestartEndHelper()
{
    // Call function on parent step.
    if (nullptr != _parentStep)
    {
        (const_cast<Step*>(_parentStep))->hotRestartOnUserInterrupt();
    }

    // _userInterrupt must be set to false only when the hot restart process is done.
    // Reset base stop reason (remove ctrl-c)
    if (!_userTerminate && _userInterrupt)
    {
        _userInterrupt = false;
        _stopReasons->set(NOMAD::BaseStopType::STARTED);
    }
}


void NOMAD::Step::debugShowCallStack() const
{
    std::vector<std::string> stepNameStack;
    NOMAD::Step* step = const_cast<NOMAD::Step*>(this);
    while (nullptr != step)
    {
        stepNameStack.push_back(step->getName());
        step = const_cast<NOMAD::Step*>(step->getParentStep());
    }

    if (stepNameStack.empty())
    {
        return;
    }

    // Show the steps in order, this is why we created the stack.
    std::cout << "Call stack:" << std::endl;
    // NOTE: Using "i < stepNameStack.size()" as condition for loop,
    // since i is a size_t (it is always >= 0).
    for (size_t i = stepNameStack.size()-1; i < stepNameStack.size(); i--)
    {
        for (size_t j = 0; j < (stepNameStack.size()-i-1); j++)
        {
            // indentation
            std::cout << "  ";
        }
        std::cout << stepNameStack[i] << std::endl;
    }
    std::cout << std::endl;
}



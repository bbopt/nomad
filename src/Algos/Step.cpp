
// Generic
#include "../Algos/Algorithm.hpp"
#include "../Algos/EvcInterface.hpp"
#include "../Algos/Iteration.hpp"
#include "../Algos/MainStep.hpp"
#include "../Algos/MegaIteration.hpp"
#include "../Algos/Step.hpp"
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
    NOMAD::AllStopReasons::set( NOMAD::BaseStopType::CTRL_C );

    NOMAD::Step::_userInterrupt = true;
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
    NOMAD::OutputInfo outputInfo(_name, s, NOMAD::OutputLevel::LEVEL_INFO, isBlockStart, isBlockEnd);
    NOMAD::OutputQueue::Add(std::move(outputInfo));
    OUTPUT_INFO_END
}


void NOMAD::Step::AddOutputInfo(const std::string& s, NOMAD::OutputLevel outputLevel) const
{
    if (NOMAD::OutputQueue::GoodLevel(outputLevel))
    {
        NOMAD::OutputInfo outputInfo(_name, s, outputLevel);
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
        std::string err = "Parent step for \"" + _name + "\" should not be NULL";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::Step::verifyGenerateAllPointsBeforeEval(const std::string& method, const bool expected) const
{
    bool actual = _runParams->getAttributeValue<bool>("GENERATE_ALL_POINTS_BEFORE_EVAL");

    if (expected != actual)
    {
        std::string err = "Error: " + method + " should only be called if ";
        err += " parameter GENERATE_ALL_POINTS_BEFORE_EVAL is ";
        err += (expected ? "true" : "false");
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }
}


bool NOMAD::Step::isAnAlgorithm() const
{
    NOMAD::Step* step = const_cast<Step*>(this);
    return (nullptr != dynamic_cast<NOMAD::Algorithm*>(step));
}


std::string NOMAD::Step::getAlgoName() const
{
    std::string s = "";
    if (isAnAlgorithm())
    {
        s = getName();
    }
    else
    {
        auto algo = getParentOfType<NOMAD::Algorithm*>();
        if (nullptr != algo)
        {
            s = algo->getName();
        }
    }

    // Append a space for easiness of use
    if (!s.empty())
    {
        s += " ";
    }

    return s;
}


std::string NOMAD::Step::getAlgoComment() const
{
    std::string algoComment;
    const NOMAD::MainStep* mainstep = getParentOfType<NOMAD::MainStep*>(false);
    if (nullptr != mainstep)
    {
        algoComment = mainstep->getAlgoComment();
    }
    return algoComment;
}


void NOMAD::Step::setAlgoComment(const std::string& algoComment, const bool force)
{
    NOMAD::MainStep* mainstep = getParentOfType<NOMAD::MainStep*>(false);
    if (nullptr != mainstep)
    {
        mainstep->setAlgoComment(algoComment, force);
    }
}


void NOMAD::Step::resetPreviousAlgoComment(const bool force)
{
    NOMAD::MainStep* mainstep = getParentOfType<NOMAD::MainStep*>(false);
    if (nullptr != mainstep)
    {
        mainstep->resetPreviousAlgoComment(force);
    }
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


// Get the frame center from the iteration ancestor.
const std::shared_ptr<NOMAD::EvalPoint> NOMAD::Step::getIterationFrameCenter() const
{
    std::shared_ptr<NOMAD::EvalPoint> frameCenter = nullptr;
    const NOMAD::Iteration* iteration = getParentOfType<NOMAD::Iteration*>();

    if (nullptr != iteration)
    {
        frameCenter = iteration->getFrameCenter();
    }
    return frameCenter;
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
        // Get first parent of type MegaIteration.
        auto constMegaIter = getParentOfType<NOMAD::MegaIteration*>();
        megaIter = const_cast<NOMAD::MegaIteration*>(constMegaIter);
    }

    if (nullptr != megaIter)
    {
        barrier = megaIter->getBarrier();
    }
    return barrier;
}


void NOMAD::Step::hotRestartOnUserInterrupt()
{
    hotRestartBeginHelper();

    hotRestartEndHelper();
}


void NOMAD::Step::hotRestartBeginHelper()
{
    if (nullptr != _runParams
        && !_runParams->getAttributeValue<bool>("HOT_RESTART_ON_USER_INTERRUPT"))
    {
        setUserTerminate();
        _stopReasons->set( NOMAD::BaseStopType::CTRL_C);
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



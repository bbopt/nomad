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

#ifndef __NOMAD400_STEP__
#define __NOMAD400_STEP__

#include "../Algos/MeshBase.hpp"
#include "../Eval/Barrier.hpp"
#include "../Output/OutputInfo.hpp"
#include "../Param/PbParameters.hpp"
#include "../Param/RunParameters.hpp"
#include "../Type/CallbackType.hpp"
#include "../Type/EvalType.hpp"
#include "../Util/AllStopReasons.hpp"

#include "../nomad_nsbegin.hpp"

class Step;
class Algorithm;

typedef std::function<void(const Step& step, bool &stop)> StepEndCbFunc;  ///< Type definitions for callback functions at the end of a step.
typedef std::function<void(std::vector<std::string>& paramLines)> HotRestartCbFunc; ///< Type definitions for callback functions for hot restart.

/// Base class of all types of steps (Iteration, Termination, Initialization, Poll, Mads,...).
class Step
{

protected:
    static bool _userInterrupt; ///< Interrupt NOMAD if Ctrl-C is pressed.
    static bool _userTerminate; ///< Terminate NOMAD if Ctrl-C is pressed again.

    const Step*  _parentStep; ///< The parent of this step.
    std::string         _name;  ///< The name of this step.

    std::shared_ptr<AllStopReasons>      _stopReasons; ///< The stop reasons of an algorithm.


    std::shared_ptr<RunParameters>       _runParams; ///< The run parameters that control a step.
    std::shared_ptr<PbParameters>        _pbParams;  ///< The problem parameters that control a step.

    // Callbacks that may be re-implemented by the user
    static StepEndCbFunc    _cbIterationEnd;
    static StepEndCbFunc    _cbMegaIterationEnd;
    static HotRestartCbFunc _cbHotRestart;

    // By default, always show warnings.
    // Some warnings do not need to be shown in some cases, ex. unit tests.
    static bool _showWarnings;

public:

    /// Constructor #1 for MainStep (no parent)
    /**
     */
    explicit Step()
      : _parentStep(nullptr),
        _name("Main Step"),
        _stopReasons(nullptr),
        _runParams(nullptr),
        _pbParams(nullptr)
    {
        init();
    }


    /// Constructor #2 for child step of a parent sharing the same stopReason
    /**
     \param parentStep      The parent of this step (cannot be nullptr).
     \param runParams       The run parameters that control this step (null by default).
     \param pbParams        The problem parameters that control this step (null by default).
     */
    explicit Step(const Step* parentStep,
                  const std::shared_ptr<RunParameters>   &runParams = nullptr,
                  const std::shared_ptr<PbParameters> &pbParams = nullptr)
      : _parentStep(parentStep),
        _name("Step"),
        _runParams(runParams),
        _pbParams(pbParams)
    {
        if (_parentStep == nullptr)
        {
            throw Exception(__FILE__, __LINE__, "Parent step is NULL. This constructor is for child steps having a parent only.");
        }
        else
        {
            _name = "Child step";
            _stopReasons = parentStep->getAllStopReasons();
        }
        init();
    }

    /// Constructor #3: for a child Step with a provided stopReason (such as an algorithm)
    /**
     \param parentStep      The parent of this step (can be nullptr if child of a MainStep).
     \param stopReasons     The stop reasons for all the steps of an algo (cannot be nullptr)
     \param runParams       The run parameters that control this step (null by default).
     \param pbParams        The problem parameters that control this step (null by default).
     */
    explicit Step(const Step* parentStep,
                  std::shared_ptr<AllStopReasons> stopReasons,
                  const std::shared_ptr<RunParameters> &runParams = nullptr,
                  const std::shared_ptr<PbParameters> &pbParams = nullptr)
      : _parentStep(parentStep),
        _name("Step"),
        _stopReasons(stopReasons),
        _runParams(runParams),
        _pbParams(pbParams)
    {
        if (nullptr == _stopReasons)
        {
            throw Exception(__FILE__, __LINE__, "StopReason is NULL. Must be provided for this child step.");
        }

        init();
    }


    /// Destructor
    /**
     Upon destruction of a step the output queue is flushed. Time to print.
     */
    virtual ~Step();

    // Get / Set

    /// Interruption call by user.
    /**
     Called by pressing Ctrl-C.
     */
    static bool getUserTerminate() { return _userTerminate; }

    /// Interruption requested
    static void setUserTerminate() { _userTerminate = true; }

    /// Get the parent step.
    /**
     * There is no setParentStep(). We should not change parent step externally.

     \return The parent step of this step.
     */
    const Step* getParentStep() const { return _parentStep; }

    /// Get the name of this step
    /**
     \return A /c string containing the name of this step.
     */
    virtual const std::string& getName() const { return _name; }

    /// Set the name of this step
    /**
     \param name    The name is provided as a \c string -- \b IN.
     */
    void setName(const std::string& name) { _name = name; }

    const std::shared_ptr<AllStopReasons>& getAllStopReasons() const { return _stopReasons ; }

    const std::shared_ptr<RunParameters>& getRunParams() const { return _runParams; }
    const std::shared_ptr<PbParameters>& getPbParams() const { return _pbParams; }

    /// Interruption call by user.
    /**
     * Called when the user pressed Ctrl-C.
     \param signalValue Signal value -- \b IN.
     */
    static void userInterrupt(int signalValue);
    static void debugSegFault(int signalValue);

    static bool getUserInterrupt() { return _userInterrupt; }

    /// \brief Set user callback
    void addCallback(const CallbackType& callbackType,
                     const StepEndCbFunc& stepEndCbFunc);
    void addCallback(const CallbackType& callbackType,
                     const HotRestartCbFunc& hotRestartCbFunc);

    /// \brief Run user callback
    static void runCallback(CallbackType callbackType,
                            const Step& step,
                            bool &stop);
    static void runCallback(CallbackType callbackType,
                            std::vector<std::string>& paramLines);

    static void disableWarnings() { _showWarnings = false; }

    /// \brief display output
    void AddOutputInfo(const std::string& s, bool isBlockStart, bool isBlockEnd) const;
    void AddOutputInfo(const std::string& s, OutputLevel outputLevel = OutputLevel::LEVEL_INFO) const;
    void AddOutputError(const std::string& s) const;
    void AddOutputWarning(const std::string& s) const;
    void AddOutputVeryHigh(const std::string& s) const;
    void AddOutputHigh(const std::string& s) const;
    void AddOutputDebug(const std::string& s) const;
    void AddOutputInfo(OutputInfo outputInfo) const;

    /// Template function to get the parent of given type.
    /**
     * Starting with parent of current Step, and going through ancestors,
    get first Step that is of type T.
     * By default, stop if an Algorithm is found. Returned Step could be
     irrelevant otherwise. To go further up than an Algorithm, set optional
     parameter stopAtAlgo to false.
     */
    template<typename T>
    T getParentOfType(const bool stopAtAlgo = true) const
    {
        Step* retStep = nullptr;

        Step* step = const_cast<Step*>(_parentStep);
        while (nullptr != step)
        {
            if (nullptr != dynamic_cast<T>(step))
            {
                retStep = step;
                break;
            }
            else if (stopAtAlgo && step->isAnAlgorithm())
            {
                break;
            }
            step = const_cast<Step*>(step->getParentStep());
        }

        return dynamic_cast<T>(retStep);
    }


    bool isAnAlgorithm() const;

    /// Get Algorithm ancestor that has no Algorithm ancestor.
    const Algorithm* getRootAlgorithm() const;

    /**
     \return the name of the first Algorithm ancestor of this Step,
     or the Step itself, if it is an Algorithm.
     \note If the Algorithm ancestor exists, a blank space is added at
     the end of the string for easier use. This method is mostly used
     to compute Step names as sub-steps of algorithms.
     */
    std::string getAlgoName() const;

    /**
      Get comment that will be shown in normal display for additional information.
    */
    virtual std::string getAlgoComment() const;
    /**
      Set comment that will be shown in normal display for additional information.
    */
    virtual void setAlgoComment(const std::string& algoComment, const bool force = false);

    virtual void resetPreviousAlgoComment(const bool force = false);

    /**
     \return The MeshBase for the first Iteration ancestor of this Step.
     */
    const std::shared_ptr<MeshBase> getIterationMesh() const;

    /**
     \return The frameCenter for the first Iteration ancestor of this Step.
     */
    const std::shared_ptr<EvalPoint> getIterationFrameCenter() const;

    /**
     \return The Barrier for the main MegaIteration ancestor of this Step.
     */
    const std::shared_ptr<Barrier> getMegaIterationBarrier() const;

    /**
    Start of the Step. Initialize values for the run.
    */
    void start() ;

    /**
    Placeholder to be implemented in derived classes. Called by start.
    */
    virtual void startImp() = 0 ;

    /**
     * Perform main step task.
     * Main part of the Step
     \return \c true if the Step was positive, for instance, a success was found;
     \c false if there was no success running this step
    */
    bool run();

    /**
    Placeholder to be implemented in derived classes. Called by run.
    */
    virtual bool runImp() = 0 ;

    /**
     * End of the Step. Clean up structures, flush output.
    */
    void end();

    /**
    Placeholder to be implemented by derived classes. Called by end.
    */
    virtual void endImp() = 0 ;

    /// Helper for hot restart functionalities
    virtual void hotRestartOnUserInterrupt();

    /// For debugging purposes. Show the stack of Steps for this step.
    void debugShowCallStack() const;

protected:
    /// Helper for constructors.
    /**
     Throw Exception when not verified.
     */
    void verifyParentNotNull();

    /// Helper for validating steps depending on parameter GENERATE_ALL_POINTS_BEFORE_EVAL
    void verifyGenerateAllPointsBeforeEval(const std::string& method, const bool expected) const;

    /// Helpers for hot restart, to be called at the start and end of any override.
    void hotRestartBeginHelper();
    /// Helpers for hot restart, to be called at the start and end of any override.
    void hotRestartEndHelper();

private:

    /// Helper for constructor
    void init();

    // Default callbacks. They do nothing.
    static void defaultStepEnd(const Step& step  __attribute__((unused)), bool &stop) { stop = false; }
    static void defaultHotRestart(std::vector<std::string>& paramLines  __attribute__((unused))) {};

    /**
     Default task always executed when start() is called
     */
    void defaultStart();

    /**
     Default task always executed when end() is called
     */
    void defaultEnd();

};


class StepException : public Exception
{
public:
    /// Constructor
    StepException(const std::string& file, const size_t line, const std::string & msg, const Step* step)
      : Exception(file, line, msg)
    {
        if (nullptr != step)
        {
            step->debugShowCallStack();
        }
    }
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_STEP__

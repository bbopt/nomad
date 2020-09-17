/**
  \file   MainStep.hpp
  \brief  Main Step to hold MADS, or other Algorithms
  \author Viviane Rochon Montplaisir
  \date   June 2018
*/
#ifndef __NOMAD400_MAINSTEP__
#define __NOMAD400_MAINSTEP__

#include "../Algos/Algorithm.hpp"
#include "../Eval/Evaluator.hpp"
#include "../Param/AllParameters.hpp"

#include "../nomad_nsbegin.hpp"

/// Main step to manage an algorithm execution.
/**
* MainStep is the root step of an executable running an algorithm.
* Think of it as a wrapper around the algorithm.
* MainStep manages all algorithm needs: Parameters, Evaluator
* and has utility functions (ex. displayUsage).
* MainStep takes care of the OpenMP parallelism.
* An algorithm can call other algorithms during its execution.
*/
class MainStep: public Step
{
private:
    std::string                         _paramFileName;  ///< Name of the file containing the parameters.
    std::shared_ptr<AllParameters>      _allParams;
    std::shared_ptr<Evaluator>          _evaluator; ///< Used in library running mode (not batch mode)
    std::vector<std::shared_ptr<Algorithm>>  _algos;
    std::string                         _algoComment;   ///< Comment to appear in the stats, e.g. "Phase One"
    std::vector<std::string>            _prevAlgoComment; ///< Pile of previous comments, used when going back to the main algo after running a sub-algo.
    bool                                _forceAlgoComment; ///< When true, do not change comment until reset is called


public:
    /// Constructor
    explicit MainStep()
    : Step(),
        _paramFileName(""),
        _evaluator(nullptr),
        _algos(),
        _algoComment(""),
        _prevAlgoComment(),
        _forceAlgoComment(false)
    {
        init();
    }

    /// Destructor
    virtual ~MainStep();

    /*---------*/
    /* Get/Set */
    /*---------*/

    /**
     In batch mode: Set the parameter file name, which will be read in start().
     In library mode: Set the parameters directly using set_PARAM_NAME(...). In library mode, it is also possible to read a parameter file.
     */
    void setParamFileName(const std::string& paramFileName) { _paramFileName = paramFileName;}

    void setAllParameters(const std::shared_ptr<AllParameters> &allParams);

    /**
     The evaluator may be shared between main threads.
     */
    void setEvaluator(std::shared_ptr<Evaluator> ev) { _evaluator = ev;}

    /// Set comment to be added at the end of the display stats, e.g., "Phase One"
    void setAlgoComment(const std::string& algoComment, const bool force = false) override;
    /// Reset comment to the previous in the stack
    void resetPreviousAlgoComment(const bool force = false) override;
    /// Get current comment on the top of the stack
    std::string getAlgoComment() const override { return _algoComment; }


    /*---------*/
    /* Others  */
    /*---------*/

    /**
     Once all algos have been added are executed (call start, run and end for each algo) in the MainStep::run.
     */
    void addAlgo(const std::shared_ptr<Algorithm> algo) { _algos.push_back(algo); }

    void clearAlgos() { _algos.clear(); }

    /// Helper function called by the code main function if necessary.
    void displayUsage(const char* exeName);

    /// Helper function called by the code main function if necessary.
    void displayVersion();

    /// Helper function called by the code main function if necessary.
    void displayInfo();

    /// Helper function called by the code main function if necessary.
    void displayHelp(const std::string& helpSubject = "all", bool devHelp = false);

    /**
     The user has requested a hot restart. Update the parameters with the changes requested by the user (read file or set inline).
     */
    void hotRestartOnUserInterrupt() override;

    /// Helper to reset some components (used by the runner when running multiple optimization)
    static void resetComponentsBetweenOptimization();

protected:
    /// Specific implementation to start NOMAD
    /**
     Implementation called by Step::start.
     During the main step start, the parameters are read (if parameter file is available), the cache, the evaluator and the evaluator control are created. If an algorithm is set for run in the parameter file, it is added (MainStep::addAlgo).
     */
    virtual void startImp() override;

    /**
     Implementation called by Step::run.
     Once all algos have been added are executed (call start, run and end for each algo) in the MainStep::run. \n

     If a stop reason (not the default STARTED) is propagated to the MainStep, the sequence of algos is stopped.
     */
    virtual bool runImp() override;

    /**
     Implementation called by Step::end.
     */
    virtual void endImp() override;

    /// Helper function when creating the evaluator
    int getNumThreads() const;

    /// Set the number of threads for the evaluator
    void setNumThreads() const;

    void printNumThreads() const;

    /// Detect if a Phase One search is required
    /**
     A phase one search is required if an EB type constraint
     is not feasible for X0.
     */
    bool detectPhaseOne();

    /// Helper for start
    void createCache() const;

    /// Helper for start
    void updateX0sFromCache() const;


private:
    /// Helper for constructor
    void init();
};


#include "../nomad_nsend.hpp"


#endif // __NOMAD400_MAINSTEP__

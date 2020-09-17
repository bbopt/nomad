
#ifndef __NOMAD400_ALGORITHM__
#define __NOMAD400_ALGORITHM__

#include "../Algos/Initialization.hpp"
#include "../Algos/MegaIteration.hpp"
#include "../Algos/Step.hpp"
#include "../Algos/Termination.hpp"

#include "../nomad_nsbegin.hpp"


/// Generic class for any direct search optimizer algorithm
/**
  \note: AllParameters and EvaluatorControl are held by MainStep.
 \note: Cache is a singleton all by itself.
 \note MegaIteration holds the algorithm-related structures.
 */
class Algorithm: public Step
{
protected:

    std::unique_ptr<Initialization>  _initialization;   ///< To initialize the algorithm (X0)
    std::unique_ptr<Termination>     _termination;      ///< To verify termination conditions
    std::shared_ptr<MegaIteration>   _megaIteration;    ///< MegaIteration used to keep information between steps

    bool _endDisplay;

#ifdef TIME_STATS
    double _startTime;
    double _totalRealAlgoTime;
    double _totalCPUAlgoTime;
#endif // TIME_STATS

public:
    /// Constructor
    /**
     \param parentStep          The parent of this Step -- \b IN.
     \param stopReasons         The stop reasons of this algo -- \b IN.
     \param runParams           The run parameters that control the algorithm -- \b IN.
     \param pbParams            The problem parameters that control the algorithm -- \b IN.
     */
    explicit Algorithm(const Step* parentStep,
                       std::shared_ptr<AllStopReasons> stopReasons,
                       const std::shared_ptr<RunParameters>& runParams,
                       const std::shared_ptr<PbParameters>& pbParams )
      : Step(parentStep, stopReasons, runParams, pbParams),
        _initialization(nullptr),
        _termination(nullptr),
        _megaIteration(nullptr),
        _endDisplay(true)
#ifdef TIME_STATS
        ,_startTime(0.0),
        _totalRealAlgoTime(0.0),
        _totalCPUAlgoTime(0.0)
#endif // TIME_STATS
    {
        init();
    }

    /// Destructor
    virtual ~Algorithm();

    /*---------*/
    /* Get/Set */
    /*---------*/
    const std::shared_ptr<MegaIteration> getMegaIteration() const { return _megaIteration; }
    void setMegaIteration(const std::shared_ptr<MegaIteration> megaIteration) { _megaIteration = megaIteration; }

    void setEndDisplay( bool endDisplay ) {_endDisplay = endDisplay; }


protected:
    ///  Helper for Constructor.
    void init();

    /// Default implementation of the start tasks of an algorithm
    /**
     If doing a hot restart get the algorithm ready to continue. \n
     If starting a new algorithm, reset the stop reason, the lap evaluation counter, and perform initialization.
     */
    virtual void startImp() override;

    /// Default implementation of the end tasks of an algorithm
    /**
     Display some information, reset the lap counters and save information for a potential hot restart.
     */
    virtual void endImp() override;

    /// Each algorithm must implement its run tasks.
    /**
     Run algorithm execution for single-objective.
     \return \c true
     */
    virtual bool runImp() override = 0;

    /// Helper for start() when doing a hot restart.
    virtual void readInformationForHotRestart() = 0;

    /// Helper for end()
    void saveInformationForHotRestart() const;
    /// Helper for end()
    void displayBestSolutions() const;
    /// Helper for end()
    void displayEvalCounts() const;

    /// Helper for hot restart
    void hotRestartOnUserInterrupt() override;

public:
    /**
     Sub-algo: an algorithm can be part of an algorithm.
     */
    bool isSubAlgo() const;
    bool isMainAlgo() const { return !isSubAlgo(); }

    /*---------*/
    /* Others  */
    /*---------*/
    /// Verify if this Algorithm is ready to be terminated
    bool terminate(size_t iteration);

    virtual void read(std::istream& is);
    virtual void display(std::ostream& os) const;

};

/// Operator to write parameters used for hot restart.
std::ostream& operator<<(std::ostream& os, const Algorithm& algo);

/// Operator to read parameters used for hot restart.
std::istream& operator>>(std::istream& is, Algorithm& algo);


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_ALGORITHM__

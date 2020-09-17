#ifndef __NOMAD400_POLL__
#define __NOMAD400_POLL__

#include <set>

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Mads/PollMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the poll (step 3) of MADS algorithm.
/**
 Generate the trial points (Poll::startImp), launch evaluation (Poll::runImp) and postprocecssing (Poll::endImp).
 */
class Poll: public Step, public IterationUtils
{
private:
    std::shared_ptr<PollMethodBase> _pollMethod; ///< Unlike for search, a single Poll method is executed
#ifdef TIME_STATS
    static double  _pollTime;        ///< Total time spent running the poll
    static double  _pollEvalTime;    ///< Total time spent evaluating poll points
#endif // TIME_STATS


public:
    /// Constructor
    /**
     \param parentStep The parent of this poll step
     */
    explicit Poll(const Step* parentStep)
      : Step(parentStep),
        IterationUtils(parentStep)
    {
        init();
    }
    virtual ~Poll() {}

    /// Generate new points to evaluate
    /**
     The trial points are obtained by:
        - adding poll directions (Poll::setPollDirections) to the poll center (frame center).
        - snaping points (and directions) to bounds.
        - projecting points on mesh.
     */
    void generateTrialPoints() override ;

#ifdef TIME_STATS
    /// Time stats
    static std::vector<double> getPollTime()       { return _pollTime; }
    static std::vector<double> getPollEvalTime()   { return _pollEvalTime; }
#endif // TIME_STATS



private:
    /// Helper for constructor
    void init();

    /// Implementation for start tasks for MADS poll.
    /**
     Call to generate trial points and test for mesh precision
     */
    virtual void    startImp() override ;

    /// Implementation for run tasks for MADS poll.
    /**
     Start trial points evaluation.
     \return Flag \c true if found better solution \c false otherwise.
     */
    virtual bool    runImp() override;

    /// Implementation for end tasks for MADS poll.
    /**
     Call the IterationUtils::postProcessing of the points.
     */
    virtual void    endImp() override ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_POLL__

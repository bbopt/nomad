#ifndef __NOMAD400_MEGASEARCHPOLL__
#define __NOMAD400_MEGASEARCHPOLL__

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the mega search and poll of MADS
/**
 Calling the start function generates search and poll trial points at the same time before starting evaluation.
 Calling the run function starts the evaluaions.
 The postprocessing is performed when calling the end funcion.
 */
class MegaSearchPoll: public Step, public IterationUtils
{
private:
    /**
     Hash table to remember which iteration generated this point.
     I tried working around it, but in the end it is easier to just
    remember the iteration.
    Mutable because it is updated in generateTrialPoints().
     */
    mutable std::map<EvalPoint, std::shared_ptr<MadsIteration>, EvalPointCompare> _iterForPoint;

public:
    /// Constructor
    /**
     \param parentStep The parent of this step
     */
    explicit MegaSearchPoll(const Step* parentStep )
      : Step( parentStep ),
        IterationUtils( parentStep ),
        _iterForPoint()
    {
        init();
    }

    // Destructor
    virtual ~MegaSearchPoll()
    {
        _iterForPoint.clear();
    }

    /**
     Get which iteration generated a point. This is used by evaluator control interface. Having the iteration, gives access to the mesh and some of its attributes.
     */
    const std::shared_ptr<MadsIteration> getIterForPoint(const EvalPoint& point) const;

private:

    /// Generate the trial poins for the search and poll steps.
    /**
     Call MegaSearchPoll::generateTrialPoints.
     */
    virtual void    startImp() override;

    ///Start evaluations
    virtual bool    runImp() override;

    /**
     Call for postprocessing: computation of a new hMax and update of the barrier.
     */
    virtual void    endImp() override ;

    void init();

    /// Generate new points to evaluate
    /**
     The trial points are produced using poll and search. The duplicates are removed and they are merged all together.
     */
    void generateTrialPoints() override ;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_MEGASEARCHPOLL__

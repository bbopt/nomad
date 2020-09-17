#ifndef __NOMAD400_POLLMETHODBASE__
#define __NOMAD400_POLLMETHODBASE__

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Step.hpp"
#include "../../Math/Direction.hpp"

#include "../../nomad_nsbegin.hpp"


/// Class for generic poll method of MADS. Run by Poll.
/**
 */
class PollMethodBase: public Step , public IterationUtils
{
private:

    std::string _comment; ///<  Comment shown when a poll method is used

public:
    /// Constructor
    /**
     /param parentStep      The parent of this poll step -- \b IN.
     */
    explicit PollMethodBase( const Step* parentStep )
      : Step( parentStep ),
        IterationUtils ( parentStep )
    {
        init();
    }

    /// Implementation of endImp. Intermediate function called to generate the trial points
    /**
     - Call for the intermediate base PollMethodBase::generateTrialPoints (call generateTrialPointsImp, snap on bounds and mesh).
     - Sanity check on generated trial points
     - Update the points with frame center
     */
    void startImp() override;

    /// Implementation of endImp. Function called to evaluate the trial points
    /**
     - Evaluate the trial points and update the barrier.
     - The projection of trial points on bounds and on mesh is performed before this function is called and after the function PollMethodBase::generateTrialPointsImp is called.
     */
    bool runImp() override;

    /// Implementation of endImp
    /**
        Call to the postProcessing function to update the Barrier
    */
    void endImp() override ;

    /// Intermediate function (not yet the  implementation that generate the trial points)
    /**
     - Display before and after generation comments.
     - Launches the implementation of the poll method to generate the trial points (::generateTrialPointsImp).
     - Snap the points to bounds and mesh.
     */
    void generateTrialPoints() override;

    /// Generate poll directions on a unitary frame. See derived classes (Ortho2nPollMethod, Np1UniPollMethod,...) for implementations.
    virtual void generateUnitPollDirections(std::list<Direction> &directions, size_t dim) const = 0 ;

protected:
    void init();

private:

    /// Scale and project on mesh poll directions.
    /**
     /param dirs      The unit directions to be scaled and projected on mesh -- \b IN/OUT.
     */
    void scaleAndProjectOnMesh(std::list<Direction> & dirs);


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_POLLMETHODBASE__


#ifndef __NOMAD400_SEARCHMETHODSIMPLE__
#define __NOMAD400_SEARCHMETHODSIMPLE__

#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for generic search method of MADS. Run by Search.
/**
 - Pure virtual class. Final class derived of this must implement ::generateTrialPointsImp.
 - The evaluation of derived class is performed when ::runImp is called.
 - Projection on mesh and bounds is performed after ::generateTrialPointsImp is called by the base class SearchMethodBase.
 */
class SearchMethodSimple: public SearchMethodBase
{
public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SearchMethodSimple( const Step* parentStep )
      : SearchMethodBase( parentStep )
    {
    }

    /// Intermediate function called to generate the trial points
    /**
     - Call for the intermediate base SearchMethodBase::generateTrialPoints (call generateTrialPointsImp, snap on bounds and mesh).
     - Sanity check on generated trial points
     - Update the points with frame center
     */
    void startImp() override;

    /// Function called to evaluate the trial points
    /**
     - Evaluate the trial points and update the barrier.
     - The projection of trial points on bounds and on mesh is performed before this function is called and after the function SearchMethodBase::generateTrialPointsImp is called.
     */
    bool runImp() override;

    /**
     - Pure virtual function.
     - The derived class must provide the implementation that generate the trial point.
     */
    virtual void generateTrialPointsImp() override = 0 ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SEARCHMETHODSIMPLE__


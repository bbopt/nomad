#ifndef __NOMAD400_SSDMADSMEGAITERATION__
#define __NOMAD400_SSDMADSMEGAITERATION__

#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Math/RandomPickup.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the iterations of SSD MADS
/**
Manager for Mads iterations.

*/
class SSDMadsMegaIteration: public MadsMegaIteration
{
private:
    std::vector<std::shared_ptr<Mads>> _madsList; ///< A collection of mads on subproblems.

    RandomPickup _randomPickup;

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling -- \b IN.
     \param mesh            Mesh on which other Iteration meshes are based -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit SSDMadsMegaIteration(const Step* parentStep,
                                  size_t k,
                                  std::shared_ptr<Barrier> barrier,
                                  std::shared_ptr<MeshBase> mesh,
                                  SuccessType success)
      : MadsMegaIteration(parentStep, k, barrier, mesh, success),
        _randomPickup(_pbParams->getAttributeValue<size_t>("DIMENSION"))
    {
        _randomPickup.reset();
    }
    // No Destructor needed - keep defaults.

    /// Implementation of the start tasks for MADS mega iteration.
    /**
     Creates a MadsIteration for each frame center and each desired mesh size.
     Use all xFeas and xInf available.
     For now, not using other frame centers.
     */
    virtual void startImp() override ;

    /// Implementation of the run tasks for MADS mega iteration.
    /**
     Manages the generation of points: either all poll and search points are generated all together before starting evaluation using the MegaSearchPoll or they are generated using a MadsIteration with search and poll separately. A run parameter controls the behavior.
     */
    virtual bool runImp() override;

private:
    void setupSubproblemParams(std::shared_ptr<PbParameters> & subProblemPbParams, std::shared_ptr<RunParameters> & subProblemRunParams, const Point & bestPoint, bool isPollster );

};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SSDMADSMEGAITERATION__

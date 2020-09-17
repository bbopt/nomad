#ifndef __NOMAD400_QUAD_MODEL_SINGLE_PASS__
#define __NOMAD400_QUAD_MODEL_SINGLE_PASS__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelIterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/**
 Final class to generate points for single pass on a given frame center. Used as Search method.
 The start, run and end tasks are empty. No evaluations are performed.
 The QuadModelSinglePass::generateTrialPoints function manages the creation process. The sample set to build the quad model is created by calling QuadModelIteration::startImp(). The points are not projected on mesh (done in SearchMethodBase).
 */
class QuadModelSinglePass final: public QuadModelIteration, public QuadModelIterationUtils
{

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param frameCenter    The "center" around which to construct the training set  -- \b IN.
     \param madsMesh           The Mads mesh for constructing the training set (no snap on mesh is performed) -- \b IN.
     */
    explicit QuadModelSinglePass(const Step* parentStep,
                                 const std::shared_ptr<EvalPoint>& frameCenter,
                                 const std::shared_ptr<MeshBase>& madsMesh)
      : QuadModelIteration(parentStep, frameCenter, 0, madsMesh),
        QuadModelIterationUtils(parentStep)
    {
        _stopReasons = std::make_shared<AlgoStopReasons<ModelStopType>>();
    }
    // No Destructor needed - keep defaults.

    /// Implementation of start task. Nothing to do.
    void startImp() override {}

    /// Implementation of run task. Nothing to do.
    bool runImp() override { return  false;}

    /// Implementation of run task. Nothing to do.
    void endImp() override {}

    /// Generate trial points
    /**
     - Update the quadratic model.
     - Optimize the quadratic model problem.
     - Insert the best feasible and best infeasible (if available) as trial points.
     */
    void generateTrialPoints() override;

};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_NMALLREFLECTIVE__

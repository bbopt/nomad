#ifndef __NOMAD400_QUAD_MODEL_ITERATION__
#define __NOMAD400_QUAD_MODEL_ITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Algos/MeshBase.hpp"
#include "../../Eval/EvalPoint.hpp"
#include "../../../ext/sgtelib/src/Surrogate.hpp"
#include "../../../ext/sgtelib/src/TrainingSet.hpp"

#include "../../nomad_nsbegin.hpp"

///
class QuadModelIteration: public Iteration
{
private:

    void init();

    /**
     - The center point of the model.
     - Cache points used to build the model are taken around this point.
     */
    const std::shared_ptr<EvalPoint> _frameCenter;

    /**
     The Mads mesh can be available if this is called during a Search method. If not, it is set to \c nullptr. When available, trials points can be projected on it.
     */
    const std::shared_ptr<MeshBase> _madsMesh;

    std::shared_ptr<SGTELIB::TrainingSet>   _trainingSet; ///<
    std::shared_ptr<SGTELIB::Surrogate>     _model;

public:
    /// Constructor
    /**
     \param parentStep       The parent of this step -- \b IN.
     \param frameCenter    The frame center -- \b IN.
     \param k                              The iteration number -- \b IN.
     \param madsMesh            Mads Mesh for trial point projection (can be null) -- \b IN.
     */
    explicit QuadModelIteration(const Step *parentStep,
                                const std::shared_ptr<EvalPoint> &frameCenter,
                                const size_t k,
                                std::shared_ptr<MeshBase> madsMesh)
      : Iteration(parentStep, k) ,
        _frameCenter(frameCenter),
        _madsMesh(madsMesh)
    {
        init();
    }


    /// \brief Destructor
    /// When iteration is done, Flush prints output queue.
    virtual ~QuadModelIteration()
    {
        reset();
    }

    /// Reset the model and the training set.
    void reset();

    /// Access to the quadratic model
    const std::shared_ptr<SGTELIB::Surrogate> getModel() const { return _model;}

    /// Access to the training set
    const std::shared_ptr<SGTELIB::TrainingSet> getTrainingSet() const { return _trainingSet; }

    /// Reimplement to have access to the frame center (can be undefined)
    const std::shared_ptr<EvalPoint> getFrameCenter() const override { return _frameCenter ; }

    /// Reimplement to have access to the mesh (can be null)
    const std::shared_ptr<MeshBase> getMesh() const override { return _madsMesh; }


protected:

    /// Manage  quad model update
    /**
     Use the cache to determine a quad model.
     */
    virtual void startImp() override;

    virtual bool runImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_ITERATION__

#ifndef __NOMAD400_QUAD_MODEL_UPDATE__
#define __NOMAD400_QUAD_MODEL_UPDATE__

#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

class QuadModelUpdate : public Step
{
private:
    OutputLevel _displayLevel;

    const Point * _frameCenter;
    ArrayOfDouble _radiuses;

public:
    explicit QuadModelUpdate(const Step* parentStep)
      : Step(parentStep),
        _displayLevel(OutputLevel::LEVEL_INFO)
    {
        init();
    }

    virtual ~QuadModelUpdate();

private:
    void init();

    /**
     No start task is required
     */
    virtual void startImp() override {}

    /// Implementation of the run task.
    /**
     Update the SGTELIB::TrainingSet and SGTELIB::Surrogate contained in the QuadModelIteration ancestor:
     - Get relevant points in cache around current frame center.
     - Add points to training set and build new model.
     - Assess if model is ready. Update its bounds.
     \return \c true if model is ready \c false otherwise.
     */
    virtual bool runImp() override;

    /**
     No end task is required
     */
    virtual void    endImp() override {}

    bool isValidForUpdate(const EvalPoint& evalPoint) const; ///< Helper function for cache find.

    bool isValidForIncludeInModel(const EvalPoint& evalPoint) const; ///< Helper function for cache find.
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_UPDATE__

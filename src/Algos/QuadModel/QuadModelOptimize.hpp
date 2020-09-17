#ifndef __NOMAD400_QUAD_MODEL_OPTIMIZE__
#define __NOMAD400_QUAD_MODEL_OPTIMIZE__

#include "../../Algos/Step.hpp"
#include "../../Algos/QuadModel/QuadModelIterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to create trial points by performing quadratic model optimization using Mads
/**
 - Start, run and end tasks are performed.
 - Start: the quadratic model optimization problem is setup and solved by calling startImp. Call ::generateTrialPoints.
 - Run: trial (oracle) points are evaluated with EvalType::BB. Set the stop reason.
 - End: Remove from cache EvalType::SGTE only cache points.
 */
class QuadModelOptimize : public Step, public QuadModelIterationUtils
{
private:

    OutputLevel                         _displayLevel;

    ArrayOfDouble _modelLowerBound; ///> Lower bound: min of trainingSet points
    ArrayOfDouble _modelUpperBound; ///> Upper bound: max of trainingSet points
    Point         _modelFixedVar;   ///> Fixed variables: fixed variables detected from trainingSet


    const std::shared_ptr<PbParameters> _refPbParams; ///< Reference to the original problem parameters.

    std::shared_ptr<RunParameters>      _optRunParams; ///< run parameters for model optimization
    std::shared_ptr<PbParameters>       _optPbParams; ///< pb parameters for model optimization


public:
    /// Constructor
    /* Parent must explicitely be a (pointer to a) QuadModelAlgo.
     * Run parameters will be recomputed for model optimization.
     */
    explicit QuadModelOptimize(const Step* parentStep,
                               const std::shared_ptr<PbParameters>               refPbParams )
      : Step(parentStep),
      QuadModelIterationUtils (parentStep),
        _displayLevel(OutputLevel::LEVEL_INFO),
        _modelLowerBound(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelUpperBound(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelFixedVar(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _refPbParams(refPbParams),
        _optRunParams(nullptr),
        _optPbParams(nullptr)
    {
        init();
    }

    /// Generate new points to evaluate
    /**
     - Setup the evaluator control parameters.
     - Manage display of sub-optimization.
     - Setup evaluator (EvalType::SGTE) and success type identification function.
     - Setup the bounds and fixed variables from the trainingSet of the quadratic model.
     - Setup run and pb parameters for Mads
     - Perform start, run and end tasks on Mads.
     - best feasible and best infeasible (if available) are inserted as trial points.
     */
    void generateTrialPoints() override;


private:
    void init();

    virtual void startImp() override; ///< The quadratic model optimization problem is setup and solved by calling startImp. Calls ::generateTrialPoints.
    virtual bool runImp() override; ///< Trial (oracle) points are evaluated with EvalType::BB. Set the stop reason.
    virtual void endImp() override; ///< Remove from cache EvalType::SGTE only cache points.

    // Helpers
    void setupRunParameters();
    void setupPbParameters();
    void setModelBoundsAndFixedVar();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_OPTIMIZE__

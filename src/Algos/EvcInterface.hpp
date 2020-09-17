
#ifndef __NOMAD400_EVCINTERFACE__
#define __NOMAD400_EVCINTERFACE__

#include "../Algos/Step.hpp"
#include "../Eval/EvaluatorControl.hpp"

#include "../nomad_nsbegin.hpp"


/// Class interface with EvaluatorControl, used by an Algorithm step through IterationUtils
/**
 \todo Complete documentation
 */
class EvcInterface
{
private:
    const Step* _step;      ///< Step that uses the EvaluatorControl
    Point _fixedVariable;   ///< Full dimension point including fixed variables

    static std::shared_ptr<EvaluatorControl> _evaluatorControl; ///< Static EvaluatorControl

public:
    /// Constructor
    /**
     \param step            The step using this EvcInterface
     */
    explicit EvcInterface(Step* step )
      : _step(step )
    {
        init();
    }

    /*---------*/
    /* Get/Set */
    /*---------*/

    static const std::shared_ptr<EvaluatorControl> getEvaluatorControl()
    {
        return _evaluatorControl;
    }

    /**
     If the evaluatorControl attribute is NULL, throws an exception.
     */
    static void setEvaluatorControl(const std::shared_ptr<EvaluatorControl>& evaluatorControl);

    /// Interface for EvaluatorControl::setBarrier.
    /**
     Transform from subBarrier to fullBarrier, that is transform points in sub-space to full-space.
     Set the barrier to a full space barrier.
     */
    void setBarrier(const std::shared_ptr<Barrier>& subBarrier);

    /// Get all evaluated points
    /**
     \note              _evaluatedPoints is cleared
     */
    std::vector<EvalPoint> retrieveAllEvaluatedPoints();


    /*---------------*/
    /* Other Methods */
    /*-------------- */

    // This method may be used by MegaIteration, or by a SearchMethod or by Poll
    /**
     *  For each point, look if it is in the cache.
     *  If it is, count a cache hit.
     *  If not, convert it to an EvalQueuePoint and add it to EvaluatorControl's Queue.

     \param trialPoints The trial points -- \b IN/OUT.
     \param useMesh     Flag to use mesh or not -- \b IN.
     */
    void keepPointsThatNeedEval(const EvalPointSet &trialPoints, bool useMesh = true);


    /**
     When points are generated and added to queue, we can start evaluation.
     */
    SuccessType startEvaluation();

    /// Evaluate a single point.
    /**
     Useful for X0. \n
     This method will convert a point from subspace to full space
     before calling EvaluatorControl's method of the same name.

     \param evalPoint   The poin to evaluate -- \b IN/OUT.
     \param hMax        The max infeasibility for keeping points in barrier -- \b IN.
     \return            \c true if evaluation worked (evalOk), \c false otherwise.
    */
    bool evalSinglePoint(EvalPoint &evalPoint, const Double &hMax = INF);

private:
    /// Helper for constructor
    void init();

    /// Helper for init
    /**
     Utility that throws an exception when not verified.
     */
    void verifyStepNotNull();

    /// Helper for init and setEvaluatorControl
    /**
     Utility that throws an exception when not verified.
     */
    static void verifyEvaluatorControlNotNull();

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVCINTERFACE__

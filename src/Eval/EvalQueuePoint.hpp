/**
 \file   EvalQueuePoint.hpp
 \brief  Point specificly designed for EvalQueue.
 \author Viviane Rochon Montplaisir
 \date   October 2018
 \see    EvalQueuePoint.cpp
 */

#ifndef __NOMAD400_EVALQUEUEPOINT__
#define __NOMAD400_EVALQUEUEPOINT__

#include "../Eval/EvalPoint.hpp"

#include <functional>   // For std::function

#include "../nomad_nsbegin.hpp"


/**
 *  Elements of the evaluation queue:
 * - evalPoint: Point to eval and its Eval. It is what goes in the cache.
 * - success: result of the comparison of evalPoint's eval with EvaluatorControl's barrier
 */
class EvalQueuePoint : public EvalPoint
{
private:
    const EvalType  _evalType;          ///< Must this point be evaluated by BB or by SGTE evaluator
    SuccessType     _success;           ///< Result of the comparison of evalPoint's eval with barrier
    bool            _relativeSuccess;   ///< Did better than the previous evaluation

    std::string     _comment;           ///< Algo comment to be printed out in DISPLAY_STATS
    std::string     _genStep;           ///< Generating step, also for stats

    ArrayOfDouble _meshSize; ///< Remenbers size of mesh that created point.
    ArrayOfDouble _frameSize; ///< Remenbers size of frame that created point.

    // TODO Support Mesh Index
    //ArrayOfDouble _meshIndex();

    size_t _k; ///< The number of the iteration that generated this point. For sorting purposes.
    int _threadAlgoNum; ///< Remembers the number of the main thread that created point.

public:

    /// Constructor
    /**
     \param evalPoint       The point to eval and its evaluation. It is what goes in the cache.-- \b IN.
     */
    explicit EvalQueuePoint(const EvalPoint& evalPoint, const EvalType& evalType)
      : EvalPoint(evalPoint),
        _evalType(evalType),
        _success(SuccessType::NOT_EVALUATED),
        _relativeSuccess(false),
        _comment(""),
        _genStep(""),
        _meshSize(),
        _frameSize(),
        _k(0),
        _threadAlgoNum(0)
    {}

    const EvalType& getEvalType() const { return _evalType; }

    void setSuccess(const SuccessType success) { _success = success; }
    const SuccessType& getSuccess() const { return _success; }

    void setRelativeSuccess(const bool relativeSuccess) { _relativeSuccess = relativeSuccess; }
    bool getRelativeSuccess() const { return _relativeSuccess; }

    void setComment(const std::string& comment) { _comment = comment; }
    const std::string& getComment() const { return _comment; }

    void setGenStep(const std::string& genStep) { _genStep = genStep; }
    const std::string& getGenStep() const { return _genStep; }

    void setMeshSize(const ArrayOfDouble& meshSize) { _meshSize = meshSize; };
    const ArrayOfDouble& getMeshSize() const { return _meshSize; }
    void setFrameSize(const ArrayOfDouble& frameSize) { _frameSize = frameSize; };
    const ArrayOfDouble& getFrameSize() const { return _frameSize; }

    void setK(const size_t k) { _k = k; };
    size_t getK() const { return _k; }

    void setThreadAlgo(const int threadAlgoNum) { _threadAlgoNum = threadAlgoNum; };
    int getThreadAlgo() const { return _threadAlgoNum; }
};

/// Smart pointer to EvalQueuePoint
typedef std::shared_ptr<EvalQueuePoint> EvalQueuePointPtr;

/// Block of EvalQueuePointPtrs
typedef std::vector<EvalQueuePointPtr> BlockForEval;


// Class for comparison using a direction.
class OrderByDirection
{
private:
    static Direction _lastSuccessfulDir;

public:
    // Get/Set
    static void setLastSuccessfulDir(const Direction& dir) { _lastSuccessfulDir = dir; }
    static Direction& getLastSuccessfulDir() { return _lastSuccessfulDir; }

    static bool comp(EvalQueuePointPtr& point1, EvalQueuePointPtr& point2);
};


/// Class to compare priority of two EvalQueuePoint
class ComparePriority
{
private:
    static std::function<bool(EvalQueuePointPtr &p1, EvalQueuePointPtr &p2)> _comp; ///< Comparison function between two eval queue points

public:
    /// Constructor
    ComparePriority() {}

    ///  Function call operator
    /**
     \param p1  First eval queue point -- \b IN.
     \param p2  Second eval queue point -- \b IN.
     \return    \c true if p1 has a lower priority than p2. \c false otherwise.
     */
    bool operator()(EvalQueuePointPtr& p1, EvalQueuePointPtr& p2);

    static void setComp(std::function<bool(EvalQueuePointPtr &p1, EvalQueuePointPtr &p2)> compFunction) { _comp = compFunction; }

    /** Currently only compares iteration number (k).
     \todo Return a basic default like tag comparison or lexicographical order.
     \param p1  First eval queue point -- \b IN.
     \param p2  Second eval queue point -- \b IN.
     \return    \c true if p1 has a lower priority than p2. \c false otherwise.
     */
    static bool basicDefaultComp(EvalQueuePointPtr& p1, EvalQueuePointPtr& p2);

    /** Compare angle with last successful direction.
     \todo figure out where the "last successful direction" comes from: is it on the points themselves, or from another structure? See unit test for EvaluatorControl. Detail that info here.
     \param p1  First eval queue point -- \b IN.
     \param p2  Second eval queue point -- \b IN.
     \return    \c true if p1 is closer in angle to last successful direction than p2. \c false otherwise.
     */
    static bool lastSuccessfulDirComp(EvalQueuePointPtr& p1, EvalQueuePointPtr& p2);
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVALQUEUEPOINT__



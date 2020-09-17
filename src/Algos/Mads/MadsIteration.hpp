#ifndef __NOMAD400_MADSITERATION__
#define __NOMAD400_MADSITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Algos/MeshBase.hpp"
#include "../../Eval/EvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for MADS iteration
/**
 A MADS iteration consists of a Search step followed by a Poll step depending on the stop reasons and successes.
 */
class MadsIteration: public Iteration
{
private:
    const std::shared_ptr<EvalPoint> _frameCenter; ///< Center around which the points are generated
    const std::shared_ptr<MeshBase>  _mesh;        ///< Mesh on which the points are
    SuccessType                      _success;     ///< Success type of this iteration

#ifdef TIME_STATS
    /// Time counters
    static double       _iterTime;          ///< Total time spent running this class
    static double       _searchTime;        ///< Total time spent running searches
    static double       _searchEvalTime;    ///< Total time spent evaluating search points
    static double       _pollTime;          ///< Total time spent running polls
    static double       _pollEvalTime;      ///< Total time spent evaluating poll points
    double              _iterStartTime;     ///< Time at which the start method was called
#endif // TIME_STATS

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param frameCenter        Frame center of this iteration -- \b IN.
     \param k                  The iteration number -- \b IN.
     \param mesh               The mesh of the iteration -- \b IN.
     */
    explicit MadsIteration(const Step *parentStep,
                           const std::shared_ptr<EvalPoint>& frameCenter,
                           const size_t k,
                           const std::shared_ptr<MeshBase> mesh)
      : Iteration(parentStep, k),
        _frameCenter(frameCenter),
        _mesh(mesh),
        _success(SuccessType::NOT_EVALUATED)
#ifdef TIME_STATS
        ,_iterStartTime(0.0)
#endif // TIME_STATS
    {
        init();
    }


    // Gets/Sets

    /**
     The Mads algorithm iteration possesses a frame center, unlike the base iteration that has none.
     \remark Used by Step::getIterationFrameCenter() to pass the frame center whenever needed
     */
    const std::shared_ptr<EvalPoint> getFrameCenter() const override { return _frameCenter; }

    /**
     The Mads algorithm iteration possesses a mesh, unlike the base iteration that has none.
     \remark Used by Step::getIterationMesh() to pass the mesh whenever needed
     */
    const std::shared_ptr<MeshBase> getMesh() const override { return _mesh; }

    /// Return current SuccessType
    const SuccessType& getSuccessType() const { return _success; }

    /// Set SuccessType member
    void setSuccessType(const SuccessType& success) { _success = success; }

#ifdef TIME_STATS
    /// Time stats
    static double getIterTime()         { return _iterTime; }
    static double getSearchTime()       { return _searchTime; }
    static double getSearchEvalTime()   { return _searchEvalTime; }
    static double getPollTime()         { return _pollTime; }
    static double getPollEvalTime()     { return _pollEvalTime; }
#endif // TIME_STATS

    /*---------------------*/
    /* Other class methods */
    /*---------------------*/

    /// Is this the main iteration of the current MegaIteration?
    bool isMainIteration() const override;


private:
    /// Helper for constructor
    void init();

    virtual void startImp() override;

    /// Implementation of the run tasks of MADS algorithm.
    /**
     Run a MADS iteration: a Search step followed by a Poll step depending on the stop reasons and successes.
     */
    virtual bool runImp() override;

#ifdef TIME_STATS
    virtual void endImp() override;
#endif // TIME_STATS
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_MADSITERATION__

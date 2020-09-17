/**
 \file   Evaluator.hpp
 \brief  Evaluation of blackbox functions.
 \author Viviane Rochon Montplaisir
 \date   September 2017
 \see    Evaluator.cpp
 */

#ifndef __NOMAD400_EVALUATOR__
#define __NOMAD400_EVALUATOR__

#include "../Eval/BBOutput.hpp"
#include "../Eval/EvalPoint.hpp"
#include "../Param/EvalParameters.hpp"
#include "../Type/EvalType.hpp"

#include "../nomad_nsbegin.hpp"

/// Enum for the type of Evaluator.
enum class EvalXDefined
{
    EVAL_BLOCK_DEFINED_BY_USER, ///< User redefined eval_block() in library mode; Default value
    EVAL_X_DEFINED_BY_USER,     ///< User redefined eval_x() in library mode
    USE_BB_EVAL                 ///< Neither eval_x() nor eval_block() were redefined by library mode. An external executable is provided.
};


class Barrier;

/// Class for the evaluator
/**
 * Evaluation of a point can be done by calling an external executable
 * (provided in BB_EXE parameter) or by redefining the evaluation function
 * Evaluator::eval_x(). /n
 * To evaluate a block of points, the user must redefine Evaluator::eval_block() or make sure the external executable can evaluate all the provided points. /n
 *
 */
class Evaluator
{
protected:
    std::shared_ptr<EvalParameters> _evalParams; ///< The parameters controlling the behavior of the evaluator

private:
    std::vector<std::string> _tmpFiles; ///< One file per thread.

    /// Did the user redefine eval_x() for single point, or should we use BB_EXE ?
    mutable EvalXDefined _evalXDefined;

    /** If we are using SGTE, it means EvalPoint's surrogate evaluation needs to be updated.
     *  If we are using BB, the blackbox evaluation is updated.
     */
    const EvalType _evalType;

public:

    /// Constructor
    /**
     \param evalParams      The parameters to control the behavior of the evaluator
     \param evalType        Which type of Eval will be updated by this Evaluator:                        blackbox (BB) or surrogate (SGTE)
     \param nbThreads       The number of threads the evaluators will run on. Used to initialize temporary files.
     \param evalXDefined    Flag.
     */
    explicit Evaluator(const std::shared_ptr<EvalParameters> &evalParams,
                       const EvalType evalType = EvalType::BB,
                       int nbThreads = 0,
                       const EvalXDefined evalXDefined = EvalXDefined::EVAL_BLOCK_DEFINED_BY_USER);

    /// Destructor.
    /**
     Needed to remove temporary files.
     */
    virtual ~Evaluator();

    /*---------*/
    /* Get/Set */
    /*---------*/
    std::shared_ptr<EvalParameters> getEvalParams() const
    {
        return _evalParams;
    }

    const EvalType& getEvalType() const { return _evalType; }

    /*---------------*/
    /* Other methods */
    /*---------------*/

    /// Evaluate the blackbox functions at a given trial point.
    /**
     * - May be user-defined.
     * - Default implementation is to use executable defined by parameter BB_EXE.

     \param x           The point to evaluate -- \b IN/OUT.
     \param countEval   Indicates if the evaluation has to be counted or not -- \b OUT.
     \param hMax        Maximum h acceptable for constraint violation -- b IN.
     \return            \c true if the evaluation succeded, \c false otherwise.
     */
    virtual bool eval_x(EvalPoint &x,
                        const Double& hMax,
                        bool &countEval) const;

    /// Evaluate the blackbox functions for a block of trial points.
    /**
     * - May be user-defined.
     * - Default implementation is to use executable defined by parameter BB_EXE.

     \param block       The block of points to evaluate -- \b IN/OUT.
     \param hMax        Maximum h acceptable for constraint violation -- \b IN.
     \param countEval   Indicates if the evaluation has to be counted or not -- \b OUT.
     \return            \c true if the evaluation succeded, \c false otherwise.
     */
    virtual std::vector<bool> eval_block(Block &block,
                                         const Double &hMax,
                                         std::vector<bool> &countEval) const;

private:
    /// Helper for eval_block()
    virtual std::vector<bool> evalXBBExe(Block &block,
                                         const Double &hMax,
                                         std::vector<bool> &countEval) const;
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVALUATOR__

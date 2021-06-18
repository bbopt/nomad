/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/
/**
 \file   EvalQueuePoint.hpp
 \brief  Point specificly designed for EvalQueue.
 \author Viviane Rochon Montplaisir
 \date   October 2018
 \see    EvalQueuePoint.cpp
 */

#ifndef __NOMAD_4_0_EVALQUEUEPOINT__
#define __NOMAD_4_0_EVALQUEUEPOINT__

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
    const EvalType  _evalType;          ///< EvalType of the evaluator that must evaluate this point (BB, MODEL, SURROGATE)
    SuccessType     _success;           ///< Result of the comparison of evalPoint's eval with barrier
    bool            _relativeSuccess;   ///< Did better than the previous evaluation

    ArrayOfDouble   _meshSize; ///< Remenbers size of mesh that created point.
    ArrayOfDouble   _frameSize; ///< Remenbers size of frame that created point.

    // Support Mesh Index (issue (feature) #381)

    size_t          _k; ///< The number of the iteration that generated this point. For sorting purposes.

public:

    /// Constructor
    /**
     \param evalPoint       The point to eval and its evaluation. It is what goes in the cache.-- \b IN.
     \param evalType        The type of evaluation (BB, MODEL,...).-- \b IN.
     */
    explicit EvalQueuePoint(const EvalPoint& evalPoint, const EvalType& evalType)
      : EvalPoint(evalPoint),
        _evalType(evalType),
        _success(SuccessType::NOT_EVALUATED),
        _relativeSuccess(false),
        _meshSize(),
        _frameSize(),
        _k(0)
    {}

    const EvalType& getEvalType() const { return _evalType; }

    void setSuccess(const SuccessType success) { _success = success; }
    const SuccessType& getSuccess() const { return _success; }

    void setRelativeSuccess(const bool relativeSuccess) { _relativeSuccess = relativeSuccess; }
    bool getRelativeSuccess() const { return _relativeSuccess; }

    void setMeshSize(const ArrayOfDouble& meshSize) { _meshSize = meshSize; };
    const ArrayOfDouble& getMeshSize() const { return _meshSize; }
    void setFrameSize(const ArrayOfDouble& frameSize) { _frameSize = frameSize; };
    const ArrayOfDouble& getFrameSize() const { return _frameSize; }

    void setK(const size_t k) { _k = k; };
    size_t getK() const { return _k; }
};

/// Smart pointer to EvalQueuePoint
typedef std::shared_ptr<EvalQueuePoint> EvalQueuePointPtr;

/// Block of EvalQueuePointPtrs
typedef std::vector<EvalQueuePointPtr> BlockForEval;


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_EVALQUEUEPOINT__



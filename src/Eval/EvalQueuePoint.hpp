/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#ifndef __NOMAD400_EVALQUEUEPOINT__
#define __NOMAD400_EVALQUEUEPOINT__

#include "../Eval/EvalPoint.hpp"

#include <functional>   // For std::function

#include "../nomad_nsbegin.hpp"


/**
 *  Elements of the evalutaion queue:
 * - evalPoint: Point to eval and its Eval. It is what goes in the cache.
 * - success: result of the comparison of evalPoint's eval with EvaluatorControl's barrier
 */
class EvalQueuePoint : public EvalPoint
{
private:

    SuccessType     _success;           ///< Result of the comparison of evalPoint's eval with barrier
    bool            _relativeSuccess;   ///< Did better than the previous evaluation

    std::string     _comment;           ///< Algo comment to be printed out in DISPLAY_STATS
    std::string     _genStep;           ///< Generating step, also for stats

    ArrayOfDouble _meshSize; ///< Remenbers size of mesh that created point.
    ArrayOfDouble _frameSize; ///< Remenbers size of frame that created point.
    
    // TODO Support Mesh Index
    //ArrayOfDouble _meshIndex();

    size_t _k; ///< The number of the iteration that generated this point. For sorting purposes.

public:
    
    /// Constructor
    /**
     \param evalPoint       The point to eval and its evaluation. It is what goes in the cache.-- \b IN.
     */
    explicit EvalQueuePoint(const EvalPoint &evalPoint)
      : EvalPoint(evalPoint),
        _success(SuccessType::NOT_EVALUATED),
        _relativeSuccess(false),
        _comment(""),
        _genStep(""),
        _meshSize(),
        _frameSize(),
        _k(0)
    {}

    
    
    void setSuccess(const SuccessType success) { _success = success; }
    SuccessType getSuccess() const { return _success; }
    
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
};

/// Smart pointer to EvalQueuePoint
typedef std::shared_ptr<EvalQueuePoint> EvalQueuePointPtr;

/// Block of EvalQueuePointPtrs
typedef std::vector<EvalQueuePointPtr> BlockForEval;

/// Class to compare priority of two EvalQueuePoint
class ComparePriority
{
private:
    std::function<bool(EvalQueuePointPtr &p1, EvalQueuePointPtr &p2)> _comp; ///< Comparison function between two eval queue points

public:
    /// Constructor
    /**
     \param comp Comparison function -- \b IN.
     */
    ComparePriority(const std::function<bool(EvalQueuePointPtr &p1, EvalQueuePointPtr &p2)> comp =  std::function<bool(EvalQueuePointPtr &p1, EvalQueuePointPtr &p2)>(defaultComp))
      : _comp(comp)
    {}

    ///  Function call operator
    /**
     \param p1  First eval queue point -- \b IN.
     \param p2  Second eval queue point -- \b IN.
     \return    \c true if p1 has a lower priority than p2. \c false otherwise.
     */
    bool operator()(EvalQueuePointPtr& p1, EvalQueuePointPtr& p2);


    /// Currently always returns false.
    /**
     \param p1  First eval queue point -- \b IN.
     \param p2  Second eval queue point -- \b IN.
     \return    \c true if p1 has a lower priority than p2. \c false otherwise.
     */
    static bool defaultComp(EvalQueuePointPtr& p1, EvalQueuePointPtr& p2);
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVALQUEUEPOINT__



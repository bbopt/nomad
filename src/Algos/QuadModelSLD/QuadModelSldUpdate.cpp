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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldAlgo.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldIteration.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldUpdate.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"

NOMAD::QuadModelSldUpdate::~QuadModelSldUpdate()
{
}

void NOMAD::QuadModelSldUpdate::init()
{
    setStepType(NOMAD::StepType::UPDATE);
    verifyParentNotNull();
    
    // Clear old model info from cache.
    // Very important because we need to update model info with new bb evaluations info.
    NOMAD::CacheBase::getInstance()->clearModelEval(NOMAD::getThreadNum());
    
    _flagUseTrialPointsToDefineBox = ( _trialPoints.size() > 0);
    
}


std::string NOMAD::QuadModelSldUpdate::getName() const
{
    if (_flagUseTrialPointsToDefineBox)
    {
        return NOMAD::stepTypeToString(_stepType);
    }
    else
    {
        return getAlgoName() + NOMAD::stepTypeToString(_stepType);
    }
}



// Update the QuadModelSld contained in the
// ancestor.
//
// 1- Get relevant points in cache, around current center (can be the frame center or not).
// 2- Add points to build new model.
// 3- Assess if model is ready. Update its bounds.
//
// Note: Update uses blackbox values
bool NOMAD::QuadModelSldUpdate::runImp()
{
    std::string s;  // Used for output
    bool updateSuccess = false;

    const NOMAD::QuadModelSldIteration * iter = getParentOfType<NOMAD::QuadModelSldIteration*>();

    if (nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Update must have a Iteration among its ancestors.");
    }
    
    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    auto model=iter->getModel();

    // Go through cache points
    OUTPUT_INFO_START
    AddOutputInfo("Review of the cache");
    OUTPUT_INFO_END
    NOMAD::Double v;

    //
    // 1- Get relevant points in cache, around current frame centers.
    //
    std::vector<NOMAD::EvalPoint> evalPointList;
    NOMAD::Double boxFactor;
    if (NOMAD::EvcInterface::getEvaluatorControl()->getUseCache())
    {
        
        // Select valid points that are close enough to model center.
        // Compute distances that must not be violated for each variable (box size).
        // For Quad Model Search: Get a Delta (frame size) if available and multiply by factor
        if (_flagUseTrialPointsToDefineBox)
        {
        
            throw NOMAD::Exception(__FILE__,__LINE__,"Quad Model SLD Update not intended for sort.");
        }
        else if (nullptr != iter->getMesh())
        {
            if (nullptr == iter->getModelCenter() )
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"Update must have a model center.");
            }
            
            // Model center is the frame center from the QuadModelIteration
            // Box size is the current frame size
            _modelCenter = *iter->getModelCenter()->getX();
            _boxSize = iter->getMesh()->getDeltaFrameSize();

            // Multiply by box factor parameter
            boxFactor = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_SEARCH_BOX_FACTOR");
            _boxSize *= boxFactor;
            
            OUTPUT_INFO_START
            s = "Mesh size: " + iter->getMesh()->getdeltaMeshSize().display();
            AddOutputInfo(s);
            s = "Frame size: " + iter->getMesh()->getDeltaFrameSize().display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
            
        }
        else
        {
           if (nullptr == iter->getModelCenter() )
           {
               throw NOMAD::Exception(__FILE__,__LINE__,"Update must have a model center.");
           }
            _boxSize = NOMAD::ArrayOfDouble(n,NOMAD::INF);
            _modelCenter = *iter->getModelCenter()->getX();
            OUTPUT_INFO_START
            AddOutputInfo("Box size set to infinity. All evaluated points in cache will be used for update.");
            OUTPUT_INFO_END
        }
        
        if ( ! _modelCenter.isComplete() || ! _boxSize.isComplete() )
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Update must have complete frame center and box size.");
        }

        OUTPUT_INFO_START
        s = "Box size: " + _boxSize.display();
        AddOutputInfo(s);
        s = "Model center: " + _modelCenter.display();
        AddOutputInfo(s);
        OUTPUT_INFO_END
        
        // Get valid points: notably, they have a BB evaluation.
        // Use CacheInterface to ensure the points are converted to subspace
        NOMAD::CacheInterface cacheInterface(this);
        auto crit = [&](const NOMAD::EvalPoint& evalPoint){return this->isValidForIncludeInModel(evalPoint);};
        cacheInterface.find(crit, evalPointList, true /*find in subspace*/);
        
        // If the number of points is less than n, enlarge the box size and repeat
        if (evalPointList.size()<n)
        {
            evalPointList.clear();
            _boxSize *= 2;
            cacheInterface.find(crit, evalPointList, true /*find in subspace*/);
        }
    }

    // Minimum and maximum number of valid points to build a model
    const size_t minNbPoints = _runParams->getAttributeValue<size_t>("SGTELIB_MIN_POINTS_FOR_MODEL");
    const size_t maxNbPoints = _runParams->getAttributeValue<size_t>("SGTELIB_MAX_POINTS_FOR_MODEL");

    size_t nbValidPoints = evalPointList.size();

    // Keep the maxNbPoints points closest to the frame center
    if (nbValidPoints > maxNbPoints)
    {
        // Sort the eval points list based on distance to model center
        std::sort(evalPointList.begin(), evalPointList.end(),
                  [&](const NOMAD::EvalPoint& x, const NOMAD::EvalPoint &y)
                        {
                            return NOMAD::Point::dist(*x.getX(),_modelCenter) < NOMAD::Point::dist(*y.getX(),_modelCenter) ;
                        }
                 );
        // Keep only the first maxNbPoints elements of list
        evalPointList.resize(maxNbPoints);
                 
        OUTPUT_INFO_START
        s = "QuadModel found " + std::to_string(nbValidPoints);
        s += " valid points to build model. Keep only ";
        s += std::to_string(maxNbPoints);
        AddOutputInfo(s);
        OUTPUT_INFO_END
        nbValidPoints = evalPointList.size();
    }

    if (nbValidPoints < 2 || nbValidPoints < minNbPoints)
    {
        
        // If no points available, it is impossible to build a model, a stop reason is set.
        if (!_flagUseTrialPointsToDefineBox)
        {
            auto stopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
            stopReason->set(NOMAD::ModelStopType::NOT_ENOUGH_POINTS );
        }

        OUTPUT_INFO_START
        AddOutputInfo("QuadModel has not enough points to build model");
        OUTPUT_INFO_END
        return false;
    }
    OUTPUT_INFO_START
    s = "Found " + std::to_string(nbValidPoints);
    s += " point";
    if (nbValidPoints > 1)
    {
        s += "s";
    }
    s += " in cache of size " + std::to_string(NOMAD::CacheBase::getInstance()->size());
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_INFO);
    OUTPUT_INFO_END

    // Set the training set
    // --------------------
    model->setY(evalPointList);
    
    
    // define scaling and identify fixed variables:
    // ---------------
    // The min box around the interpolation set Y
    //   is scaled to [-r;r] with r=MODEL_RADIUS_FACTOR.
    model->define_scaling ( boxFactor );
    
    
    if ( nbValidPoints <= static_cast<size_t>(model->get_nfree()) )
    {
        // If not enough points are available (wrt nfree), it is impossible to build a model, a stop reason is set.
        if (!_flagUseTrialPointsToDefineBox)
        {
            auto stopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
            stopReason->set(NOMAD::ModelStopType::NOT_ENOUGH_POINTS );
        }

        OUTPUT_INFO_START
        AddOutputInfo("QuadModel has not enough points to build model");
        OUTPUT_INFO_END
        return false;
    }
    
#ifdef DEBUG
    model->display_Y ( "scaled interpolation set Ys" );
#endif
    
    // construct model:
    // ----------------
    const bool quad_use_WP        = false;
    const double SVD_EPS      = 1e-13;      ///< Epsilon for SVD
    const int    SVD_MAX_MPN  = 1500;       ///< Matrix maximal size (\c m+n )
    model->construct ( quad_use_WP , SVD_EPS , SVD_MAX_MPN , static_cast<int>(maxNbPoints) );


    //
    // 3- Assess if model is ready. Update its bounds.
    //
    // Check if the model is ready
    if ( model->check() )
    {
        updateSuccess = true;
    }
    else
    {
        updateSuccess = false;
    }

    OUTPUT_INFO_START
    s = "New nb of points in model: " + std::to_string(model->get_nY());
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_INFO);
    s = "Ready: " + NOMAD::boolToString(model->check());
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_INFO);
    OUTPUT_INFO_END

    return updateSuccess;
}


bool NOMAD::QuadModelSldUpdate::isValidForIncludeInModel(const NOMAD::EvalPoint& evalPoint) const
{

    if(! isValidForUpdate(evalPoint))
        return false;
    
    NOMAD::ArrayOfDouble diff = (*evalPoint.getX() - _modelCenter);

    diff *= 2.0; // Comparison with half of the box size. But instead we multiply the diff by two.
    
    return diff.abs() <= _boxSize ;

}


bool NOMAD::QuadModelSldUpdate::isValidForUpdate(const NOMAD::EvalPoint& evalPoint) const
{
    // Verify that the point is valid
    // - Not a NaN
    // - Not a fail
    // - All outputs defined
    // - Blackbox OBJ available (Not MODEL)
    bool validPoint = true;
    NOMAD::ArrayOfDouble bbo;

    auto evalType = NOMAD::EvcInterface::getEvaluatorControl()->getEvalType();
    auto eval = evalPoint.getEval(evalType);
    if (NOMAD::EvalType::BB != evalType)
    {
        validPoint = false;
    }
    else if (nullptr == eval)
    {
        // Eval must be available
        validPoint = false;
    }
    else
    {
        // Note: it could be discussed if points that have h > hMax should still be used
        // to build the model (Nomad 3). We validate them to comply with Nomad 3.
        // If f or h greater than MODEL_MAX_OUTPUT (default=1E10) the point is not valid (same as Nomad 3)
        if (   ! eval->isBBOutputComplete()
            || !(NOMAD::EvalStatusType::EVAL_OK == eval->getEvalStatus())
            || !eval->getF(NOMAD::ComputeType::STANDARD).isDefined()
            || !eval->getH(NOMAD::ComputeType::STANDARD).isDefined()
            || eval->getF(NOMAD::ComputeType::STANDARD) > NOMAD::MODEL_MAX_OUTPUT
            || eval->getH(NOMAD::ComputeType::STANDARD) > NOMAD::MODEL_MAX_OUTPUT)
        {
            validPoint = false;
        }
    }

    return validPoint;
}


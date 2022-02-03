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
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelUpdate.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"

NOMAD::QuadModelUpdate::~QuadModelUpdate()
{
}

void NOMAD::QuadModelUpdate::init()
{
    setStepType(NOMAD::StepType::UPDATE);
    verifyParentNotNull();
    
    // Clear old model info from cache.
    // Very important because we need to update model info with new bb evaluations info.
    NOMAD::CacheBase::getInstance()->clearModelEval(NOMAD::getThreadNum());
    
    _flagUseTrialPointsToDefineBox = ( _trialPoints.size() > 0);
    
    _flagUseScaledModel = (_scalingDirections.size() > 0);
    
}


std::string NOMAD::QuadModelUpdate::getName() const
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


// Update the SGTELIB::TrainingSet and SGTELIB::Surrogate contained in the
// ancestor QuadModel (modelAlgo).
//
// 1- Get relevant points in cache, around current center (can be the frame center or not).
// 2- Add points to training set, and build new model.
// 3- Assess if model is ready. Update its bounds.
//
// Note: Update uses blackbox values
bool NOMAD::QuadModelUpdate::runImp()
{
    std::string s;  // Used for output
    bool updateSuccess = false;

    const NOMAD::QuadModelIteration * iter = getParentOfType<NOMAD::QuadModelIteration*>();

    if (nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Update must have a Iteration among its ancestors.");
    }
    
    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    const auto bbot = NOMAD::QuadModelAlgo::getBBOutputType();
    size_t nbConstraints = NOMAD::getNbConstraints(bbot);
    size_t nbModels = nbConstraints+1; /// Constraint plus a single objective

    // row_X and row_Z are one-liner matrices to stock current point.
    SGTELIB::Matrix row_X("row_X", 1, static_cast<int>(n));
    SGTELIB::Matrix row_Z("row_Z", 1, static_cast<int>(nbModels));

    // Matrices to stock all points to send to the model
    auto add_X = std::make_shared<SGTELIB::Matrix>("add_X", 0, static_cast<int>(n));
    auto add_Z = std::make_shared<SGTELIB::Matrix>("add_Z", 0, static_cast<int>(nbModels));

    auto model=iter->getModel();
    auto trainingSet=iter->getTrainingSet();


    // Go through cache points
    OUTPUT_INFO_START
    AddOutputInfo("Review of the cache");
    OUTPUT_INFO_END
    NOMAD::Double v;

    //
    // 1- Get relevant points in cache, around current frame centers.
    //
    std::vector<NOMAD::EvalPoint> evalPointList;
    if (NOMAD::EvcInterface::getEvaluatorControl()->getUseCache())
    {
        
        // Select valid points that are close enough to model center.
        // Compute distances that must not be violated for each variable (box size).
        // For Quad Model Search: Get a Delta (frame size) if available and multiply by factor
        // For Quad Model Sort: use min and max of the points to sort
        if (_flagUseTrialPointsToDefineBox)
        {
            // Case of Quad Model used for SORT
            
            NOMAD::ArrayOfDouble minVal(n,INF), maxVal(n,M_INF);
            for (const auto & pt: _trialPoints)
            {
                for (size_t i=0; i < n ; i++)
                {
                    minVal[i] = min(minVal[i], (*pt.getX())[i]);
                    maxVal[i] = max(maxVal[i], (*pt.getX())[i]);
                }
            }
            _modelCenter.reset(n);
            _boxSize.reset(n);
            
            NOMAD::Double boxFactor = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_SORT_BOX_FACTOR");
            
            // Model center is obtained by averaging min and max of the trial points.
            // Multiply box size by box factor parameter
            for (size_t i=0; i < n ; i++)
            {
                _modelCenter[i] = ( minVal[i] + maxVal[i] ) / 2.0;
                _boxSize[i] = ( maxVal[i] - minVal[i] ) * boxFactor;
            }
        }
        else if (nullptr != iter->getMesh())
        {
            
            // Case of Quad Model used for SEARCH
            
            if (nullptr == iter->getModelCenter() )
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"Update must have a model center.");
            }
            
            // Model center is the frame center from the QuadModelIteration
            // Box size is the current frame size
            _modelCenter = *iter->getModelCenter()->getX();
            _boxSize = iter->getMesh()->getDeltaFrameSize();

            NOMAD::Double boxFactor = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_SEARCH_BOX_FACTOR");
            
            // Multiply by box factor parameter
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
            OUTPUT_INFO_START
            s = "Enlarge box size to get more points: " + _boxSize.display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
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

    if (nbValidPoints < 2 || nbValidPoints < minNbPoints || nbValidPoints < n)
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

    OUTPUT_DEBUG_START
    for (const auto & dir : _scalingDirections)
    {
        s = "Scaling direction: ( " + dir.display() + ")" ;
        AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_DEBUG);
    }
    OUTPUT_DEBUG_END
    
    
    for (auto evalPoint : evalPointList)
    {
        NOMAD::Point x = *evalPoint.getX();
                
        bool success = true;
        if (_flagUseScaledModel)
        {
            success = scalingByDirections(x);
        }
        if (! success)
        {
            s = "xNew = " + evalPoint.display() + " cannot be scaled";
            AddOutputInfo(s, OutputLevel::LEVEL_DEBUG);
            return false;
        }
        else
        {
            s = "xNew = " + x.display();
            AddOutputInfo(s, OutputLevel::LEVEL_DEBUG);
        }

        for (size_t j = 0; j < n; j++)
        {
            // X
            row_X.set(0, static_cast<int>(j), x[j].todouble());
        }
        add_X->add_rows(row_X);

        // Objective
        // Update uses blackbox values
        if (!evalPoint.getF(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD).isDefined())
        {
            s = "Error: In QuadModelUpdate, this point should have a BB Eval: ";
            s += evalPoint.displayAll();
            throw NOMAD::Exception(__FILE__,__LINE__,s);
        }
        row_Z.set(0, 0, evalPoint.getF(NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD).todouble()); // 1st column: constraint model

        NOMAD::ArrayOfDouble bbo = evalPoint.getEval(NOMAD::EvalType::BB)->getBBOutput().getBBOAsArrayOfDouble();
        // Constraints
        int k = 1;
        for (size_t j = 0; j < bbo.size(); j++)
        {
            if (NOMAD::isConstraint(bbot[j]))
            {
                row_Z.set(0, k, bbo[j].todouble());
                k++;
            }
        }
        add_Z->add_rows(row_Z);

    }   // Done going through valid points

    //
    // 2- Add points to training set, and build new model.
    //
    if (add_X->get_nb_rows() > 0)
    {
        // Build the model
        OUTPUT_INFO_START
        AddOutputInfo("Add points to training set...", _displayLevel);
        OUTPUT_INFO_END

        trainingSet->partial_reset_and_add_points(*add_X, *add_Z);

        OUTPUT_INFO_START
        AddOutputInfo("OK.", _displayLevel);
        AddOutputInfo("Build model from training set...", _displayLevel);
        OUTPUT_INFO_END

        if (model->build())
        {
            OUTPUT_INFO_START
            AddOutputInfo("OK.", _displayLevel);
            OUTPUT_INFO_END
        }
        else
        {
            OUTPUT_INFO_START
            AddOutputInfo("Cannot build model.", _displayLevel);
            OUTPUT_INFO_END
        }
    }

    //
    // 3- Assess if model is ready. Update its bounds.
    //
    // Check if the model is ready
    if ( model->is_ready() )
    {
        updateSuccess = true;
    }
    else
    {
        updateSuccess = false;
    }

    OUTPUT_INFO_START
    s = "New nb of points: " + std::to_string(trainingSet->get_nb_points());
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_INFO);
    s = "Ready: " + NOMAD::boolToString(model->is_ready());
    AddOutputInfo(s, NOMAD::OutputLevel::LEVEL_INFO);
    OUTPUT_INFO_END

    return updateSuccess;
}


bool NOMAD::QuadModelUpdate::isValidForIncludeInModel(const NOMAD::EvalPoint& evalPoint) const
{

    if(! isValidForUpdate(evalPoint))
        return false;
    
    NOMAD::ArrayOfDouble diff = (*evalPoint.getX() - _modelCenter);
    
    diff *= 2.0; // Comparison with half of the box size. But instead we multiply the diff by two.
    
    return diff.abs() <= _boxSize ;

}


bool NOMAD::QuadModelUpdate::isValidForUpdate(const NOMAD::EvalPoint& evalPoint) const
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


bool NOMAD::QuadModelUpdate::scalingByDirections ( NOMAD::Point & x )
{
    
    if (_scalingDirections.size() <= 0)
    {
        // For now, we expect to have scaling directions to do scaling. Maybe we could just return.
        throw NOMAD::Exception(__FILE__,__LINE__,"Scaling directions not provided");
        
    }
    if (! _modelCenter.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Defining the scaling requires a model center");
    }
    
    size_t n = _modelCenter.size();
    if (_scalingDirections.size() != n)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Number of scaling directions must be n");
    }
    

    
    // Scale with rotation based on direction and center (see paper Reducing the number of function evaluations in Mesh Adaptive Direct Search algorithms, Audet, Ianni, LeDigabel, Tribes, 2014
    // T(y)=(D')^-1*(center-x)/delta_m/(1-epsilon) - epsilon*1/(1-epsilon)
    // (D')^-1=(D')^T/normCol^2
    NOMAD::Point temp(n,0.0);
    std::vector<NOMAD::Direction>::const_iterator itDir=_scalingDirections.begin();
    for ( size_t i = 0 ; i < n ; ++i )
    {
        temp[i]=(_modelCenter[i].todouble()-x[i].todouble())/(1.0-epsilon);
        
        x[i]=0.0;
    }
    size_t j=0;
    for (itDir=_scalingDirections.begin(); itDir != _scalingDirections.end(); itDir++,j++)
    {
        double normCol2= (*itDir).squaredL2Norm().todouble();
        if (normCol2 == 0)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Norm of a scaling direction cannot be null");
        }
        for ( size_t i = 0 ; i < n ; ++i )
        {
            x[j]+=temp[i].todouble()*(*itDir)[i].todouble()/normCol2;
        }
        x[j]-=epsilonDelta;
    }

    
    return true;
}

void NOMAD::QuadModelUpdate::unscalingByDirections( NOMAD::EvalPoint & x)
{
    if (_scalingDirections.size() <= 0)
    {
        // For now, we expect to have scaling directions to do unscaling. Maybe we could just return.
        throw NOMAD::Exception(__FILE__,__LINE__,"Scaling directions not provided");
        
    }
    if (! _modelCenter.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Defining the scaling requires a model center");
    }
    
    size_t n = _modelCenter.size();
    
    
    // UnScale with rotation see paper Reducing the number of function evaluations in Mesh Adaptive Direct Search algorithms, Audet, Ianni, LeDigabel, Tribes, 2014
    //T^−1(x) = center + Dp ((ε−1)x−ε1)
    NOMAD::Point temp(n,0.0);
    for ( size_t i = 0 ; i < n ; ++i )
    {
        temp[i]=(x[i]*(epsilon-1.0)-epsilon);
        x[i]=0.0;
    }
    std::vector<NOMAD::Direction>::const_iterator itDir;
    int j=0;
    for (itDir=_scalingDirections.begin(); itDir != _scalingDirections.end(); itDir++,j++)
    {
        for (size_t i=0 ; i< n ; i++)
        {
            x[i]+=temp[j]*(*itDir)[i];
        }
    }
    for ( size_t i = 0 ; i < n ; ++i )
    {
        x[i]+=_modelCenter[i];
    }

}

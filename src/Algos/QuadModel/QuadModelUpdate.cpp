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

}


std::string NOMAD::QuadModelUpdate::getName() const
{
    return getAlgoName() + NOMAD::stepTypeToString(_stepType);
}


// Update the SGTELIB::TrainingSet and SGTELIB::Surrogate contained in the
// ancestor QuadModel (modelAlgo).
//
// 1- Get relevant points in cache, around current frame center.
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
    AddOutputInfo("Review of the cache", _displayLevel);
    NOMAD::Double v;

    _radiuses = NOMAD::ArrayOfDouble(n,INF);
    //
    // 1- Get relevant points in cache, around current frame centers.
    //
    std::vector<NOMAD::EvalPoint> evalPointList;
    if (NOMAD::EvcInterface::getEvaluatorControl()->getUseCache())
    {
        _frameCenter = iter->getFrameCenter()->getX();
        NOMAD::EvalPoint frameCenterEvalPoint;
        NOMAD::CacheInterface cacheInterface(this);
        if (0 == cacheInterface.find(*_frameCenter, frameCenterEvalPoint))
        {
            // In the context of using the cache, ensure the frame center is in it.
            cacheInterface.smartInsert(*iter->getFrameCenter());
        }

        // Select valid points that are close enough to frame centers
        // Compute distances that must not be violated for each variable.
        // Get a Delta (frame size) if available.

        if (nullptr != iter->getMesh())
        {
            _radiuses = iter->getMesh()->getDeltaFrameSize();

            // Multiply by radius parameter
            auto radiusFactor = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_RADIUS_FACTOR");
            _radiuses *= radiusFactor;
        }
        OUTPUT_DEBUG_START
        if (nullptr != iter->getMesh())
        {
            s = "Mesh size: " + iter->getMesh()->getdeltaMeshSize().display();
            AddOutputInfo(s);
            s = "Frame size: " + iter->getMesh()->getDeltaFrameSize().display();
            AddOutputInfo(s);
        }
        s = "Radiuses: " + _radiuses.display();
        AddOutputInfo(s);
        OUTPUT_DEBUG_END

        // Get valid points: notably, they have a BB evaluation.
        // Use CacheInterface to ensure the points are converted to subspace
        auto crit = [&](const NOMAD::EvalPoint& evalPoint){return this->isValidForIncludeInModel(evalPoint);};
        cacheInterface.find(crit, evalPointList, true /*find in subspace*/);
    }

    // Minimum and maximum number of valid points to build a model
    const size_t minNbPoints = _runParams->getAttributeValue<size_t>("SGTELIB_MIN_POINTS_FOR_MODEL");
    // not using SGTELIB_MAX_POINTS_FOR_MODEL, see justification below.
    //const size_t maxNbPoints = _runParams->getAttributeValue<size_t>("SGTELIB_MAX_POINTS_FOR_MODEL");


    size_t nbValidPoints = evalPointList.size();

    /*
    // This is not a safe way to select a subset of points.
    // First, we should ensure that the frame center is part of that subset.
    // Second, how to select the other points? The points closest to
    // the frame center seem like an interesting option.
    if (nbValidPoints > maxNbPoints)
    {

        OUTPUT_INFO_START
        s = "QuadModel found " + std::to_string(nbValidPoints);
        s += " valid points to build model. Keep only ";
        s += std::to_string(maxNbPoints);
        AddOutputInfo(s);
        OUTPUT_INFO_END

        // First, sort evalPointList by dominance
        std::sort(evalPointList.begin(), evalPointList.end());
        // Next, get first maxNbPoints points
        std::vector<NOMAD::EvalPoint> evalPointListShort;
        for (size_t i = 0; i < maxNbPoints; i++)
        {
            evalPointListShort.push_back(evalPointList[i]);
        }
        evalPointList = evalPointListShort;
        nbValidPoints = evalPointList.size();
    }
    */

    if (nbValidPoints < minNbPoints)
    {
        // If no points available, it is impossible to build a model.
        auto stopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        stopReason->set(NOMAD::ModelStopType::NOT_ENOUGH_POINTS );

        OUTPUT_INFO_START
        AddOutputInfo("QuadModel has not enough points to build model");
        OUTPUT_INFO_END
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

    for (auto evalPoint : evalPointList)
    {
        NOMAD::Point x = *evalPoint.getX();
        s = "xNew = " + evalPoint.display();
        AddOutputInfo(s, _displayLevel);

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
        AddOutputInfo("Add points to training set...", _displayLevel);

        trainingSet->partial_reset_and_add_points(*add_X, *add_Z);

        AddOutputInfo("OK.", _displayLevel);
        AddOutputInfo("Build model from training set...", _displayLevel);

        model->build();
        AddOutputInfo("OK.", _displayLevel);
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
    if ( ! _frameCenter->NOMAD::ArrayOfDouble::isComplete() || ! _radiuses.isComplete() )
    {
        std::cerr << "QuadModelUpdate : isValidForIncludeInModel : frameCenter or radiuses not defined  " << std::endl;
    }

    return (isValidForUpdate(evalPoint) && NOMAD::ArrayOfDouble(*evalPoint.getX() - *_frameCenter).abs() <= _radiuses ) ;

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


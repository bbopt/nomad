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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelEvaluator.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelMegaIteration.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelUpdate.hpp"
#include "../../Type/SgtelibModelFeasibilityType.hpp"
#include "../../Type/SgtelibModelFormulationType.hpp"
#include "../../../ext/sgtelib/src/Surrogate.hpp"

NOMAD::SgtelibModelUpdate::~SgtelibModelUpdate()
{
}


void NOMAD::SgtelibModelUpdate::init()
{
    _name = getAlgoName() + "Update";
    verifyParentNotNull();
}


void NOMAD::SgtelibModelUpdate::startImp()
{
    auto modelDisplay = _runParams->getAttributeValue<std::string>("SGTELIB_MODEL_DISPLAY");
    _displayLevel = (std::string::npos != modelDisplay.find("U"))
                        ? NOMAD::OutputLevel::LEVEL_INFO
                        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

}


// Update the SGTELIB::TrainingSet and SGTELIB::Surrogate contained in the
// ancestor SgtelibModel (modelAlgo).
//
// 1- Get relevant points in cache, around current frame centers.
// 2- Add points to training set, and build new model.
// 3- Assess if model is ready. Update its bounds.
//
// Note: Update uses blackbox values
bool NOMAD::SgtelibModelUpdate::runImp()
{
    bool updateDone = true;
    std::string s;  // Used for output

    auto modelFeasibility = _runParams->getAttributeValue<NOMAD::SgtelibModelFeasibilityType>("SGTELIB_MODEL_FEASIBILITY");
    auto modelFormulation = _runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION");

    const NOMAD::SgtelibModel* modelAlgoConst = dynamic_cast<const NOMAD::SgtelibModel*>(getParentOfType<NOMAD::SgtelibModel*>());
    auto modelAlgo = const_cast<NOMAD::SgtelibModel*>(modelAlgoConst);
    if (nullptr == modelAlgo)
    {
        s = "Error: In SgtelibModelUpdate, need a SgtelibModel parent.";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    const auto bbot = NOMAD::SgtelibModel::getBBOutputType();
    size_t nbConstraints = NOMAD::getNbConstraints(bbot);
    size_t nbModels = NOMAD::SgtelibModel::getNbModels(modelFeasibility, nbConstraints);

    if (NOMAD::SgtelibModelFormulationType::EXTERN == modelFormulation)
    {
        // Extern SGTE. Early out
        AddOutputInfo("FORMULATION: EXTERN.", _displayLevel);
        return updateDone;
    }

    // row_X and row_Z are one-liner matrices to stock current point.
    SGTELIB::Matrix row_X("row_X", 1, static_cast<int>(n));
    SGTELIB::Matrix row_Z("row_Z", 1, static_cast<int>(nbModels));

    // Matrices to stock all points to send to the model
    auto add_X = std::make_shared<SGTELIB::Matrix>("add_X", 0, static_cast<int>(n));
    auto add_Z = std::make_shared<SGTELIB::Matrix>("add_Z", 0, static_cast<int>(nbModels));

    // Go through cache points
    AddOutputInfo("Review of the cache", _displayLevel);
    int k = 0;
    NOMAD::Double v;

    //
    // 1- Get relevant points in cache, around current frame centers.
    //
    std::vector<NOMAD::EvalPoint> evalPointList;
    // Get valid points: notably, they have a BB evaluation.
    // Use CacheInterface to ensure the points are converted to subspace
    NOMAD::CacheInterface cacheInterface(this);
    cacheInterface.find(validForUpdate, evalPointList);

    // Minimum and maximum number of valid points to build a model
    const size_t minNbPoints = _runParams->getAttributeValue<size_t>("SGTELIB_MIN_POINTS_FOR_MODEL");
    const size_t maxNbPoints = _runParams->getAttributeValue<size_t>("SGTELIB_MAX_POINTS_FOR_MODEL");

    // Select valid points that are close enough to frame centers
    // Compute distances that must not be violated for each variable.
    // Get a Delta (frame size) if available.
    NOMAD::ArrayOfDouble radius(n, 1.0);
    if (nullptr != modelAlgo->getMesh())
    {
        radius = modelAlgo->getMesh()->getDeltaFrameSize();
    }
    // Multiply by radius parameter
    auto radiusFactor = _runParams->getAttributeValue<NOMAD::Double>("MODEL_RADIUS_FACTOR");
    radius *= radiusFactor;

    // Get all frame centers
    auto megaIter = dynamic_cast<const NOMAD::SgtelibModelMegaIteration*>(getParentOfType<NOMAD::SgtelibModelMegaIteration*>());
    auto allCenters = megaIter->getBarrier()->getAllPoints();
    size_t nbCenters = allCenters.size();
    std::vector<NOMAD::EvalPoint> evalPointListWithinRadius;
    for (auto evalPoint : evalPointList)
    {
        bool pointAdded = false;
        for (size_t centerIndex = 0;
             centerIndex < nbCenters && !pointAdded;
             centerIndex++)
        {
            auto center = allCenters[centerIndex];
            auto distances = NOMAD::Point::vectorize(*center->getX(), evalPoint);
            distances = distances.abs();
            if (distances <= radius)
            {
                // This evalPoint is within radius. Keep it.
                evalPointListWithinRadius.push_back(evalPoint);
                pointAdded = true;
            }
        }
    }
    evalPointList = evalPointListWithinRadius;
    size_t nbValidPoints = evalPointList.size();

    if (nbValidPoints > maxNbPoints)
    {
        s = "SgtelibModel found " + std::to_string(nbValidPoints);
        s += " valid points to build model. Keep only ";
        s += std::to_string(maxNbPoints);
        AddOutputInfo(s);
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

    if (nbValidPoints < minNbPoints)
    {
        // If no points available, it is impossible to build a model.
        auto sgteStopReason = NOMAD::AlgoStopReasons<NOMAD::SgtelibModelStopType>::get(_stopReasons);
        sgteStopReason->set(NOMAD::SgtelibModelStopType::NO_POINTS);

        AddOutputError("Sgtelib Model has not enough points to build model");
    }
    s = "Found " + std::to_string(nbValidPoints);
    s += " point";
    if (nbValidPoints > 1)
    {
        s += "s";
    }
    s += " in cache of size " + std::to_string(NOMAD::CacheBase::getInstance()->size());
    AddOutputInfo(s, _displayLevel);

    for (auto evalPoint : evalPointList)
    {
        NOMAD::Point x = *evalPoint.getX();
        s = "xNew = " + x.display();
        AddOutputInfo(s, _displayLevel);

        for (size_t j = 0; j < n; j++)
        {
            // X
            row_X.set(0, static_cast<int>(j), x[j].todouble());
        }
        add_X->add_rows(row_X);

        // Objective
        // Update uses blackbox values
        row_Z.set(0, 0, evalPoint.getF(NOMAD::EvalType::BB).todouble()); // 1st column: constraint model

        NOMAD::ArrayOfDouble bbo = evalPoint.getEval(NOMAD::EvalType::BB)->getBBOutput().getBBOAsArrayOfDouble();
        // Constraints
        switch (modelFeasibility)
        {
            case NOMAD::SgtelibModelFeasibilityType::C:
                k = 1;
                for (size_t j = 0; j < bbo.size(); j++)
                {
                    if (NOMAD::isConstraint(bbot[j]))
                    {
                        row_Z.set(0, k, bbo[j].todouble());
                        k++;
                    }
                }
                break;

            case NOMAD::SgtelibModelFeasibilityType::H:
                NOMAD::SgtelibModelEvaluator::evalH(bbo, bbot, v);
                row_Z.set(0, 1, v.todouble()); // 2nd column: constraint model
                break;

            case NOMAD::SgtelibModelFeasibilityType::B:
                row_Z.set(0, 1, evalPoint.isFeasible(NOMAD::EvalType::BB)); // 2nd column: constraint model
                break;

            case NOMAD::SgtelibModelFeasibilityType::M:
                v = -NOMAD::INF;
                for (size_t j = 0; j < bbot.size(); j++)
                {
                    if (NOMAD::isConstraint(bbot[j]))
                    {
                        v = max(v, bbo[j]);
                    }
                }
                row_Z.set(0, 1, v.todouble()); // 2nd column: constraint model
                break;

            case NOMAD::SgtelibModelFeasibilityType::UNDEFINED:
            default:
                AddOutputInfo("UNDEFINED", _displayLevel);
                break;
        } // end switch

        add_Z->add_rows(row_Z);

        if (evalPoint.isFeasible(NOMAD::EvalType::BB))
        {
            if (!modelAlgo->getFoundFeasible())
            {
                AddOutputInfo(" (feasible)", _displayLevel);
                modelAlgo->setFoundFeasible(true);
            }
        }
    }   // Done going through valid points

    //
    // 2- Add points to training set, and build new model.
    //
    std::shared_ptr<SGTELIB::TrainingSet> trainingSet = modelAlgo->getTrainingSet();
    std::shared_ptr<SGTELIB::Surrogate> model = modelAlgo->getModel();

    s = "Current nb of points: " + std::to_string(trainingSet->get_nb_points());
    AddOutputInfo(s, _displayLevel);

    if (add_X->get_nb_rows() > 0)
    {
        // Build the model
        AddOutputInfo("Add points...", _displayLevel);

        trainingSet->add_points(*add_X, *add_Z);

        AddOutputInfo("OK", _displayLevel);
        AddOutputInfo("Build model...", _displayLevel);

        model->build();
        AddOutputInfo("OK.", _displayLevel);
    }

    //
    // 3- Assess if model is ready. Update its bounds.
    //
    // Check if the model is ready
    bool ready = model->is_ready();
    modelAlgo->setReady(ready);

    s = "New nb of points: " + std::to_string(trainingSet->get_nb_points());
    AddOutputInfo(s, _displayLevel);
    s = "Ready: " + NOMAD::boolToString(ready);
    AddOutputInfo(s, _displayLevel);

    // Update the bounds of the model
    modelAlgo->setModelBounds(add_X);

    return updateDone;
}


void NOMAD::SgtelibModelUpdate::endImp()
{
}


bool NOMAD::SgtelibModelUpdate::validForUpdate(const NOMAD::EvalPoint& evalPoint)
{
    // Verify that the point is valid
    // - Not a NaN
    // - Not a fail
    // - All outputs defined
    // - Blackbox OBJ available (Not sgte)
    bool validPoint = true;
    NOMAD::ArrayOfDouble bbo;

    auto eval = evalPoint.getEval(NOMAD::EvalType::BB);
    if (nullptr == eval)
    {
        // Blackbox Eval must be available
        validPoint = false;
    }
    else
    {
        bbo = eval->getBBOutput().getBBOAsArrayOfDouble();

        if (   ! bbo.isComplete()
            || ! (NOMAD::EvalStatusType::EVAL_OK == eval->getEvalStatus())
            || ! eval->getF().isDefined() )
        {
            validPoint = false;
        }
    }

    return validPoint;
}

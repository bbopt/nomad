/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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

#include "../../Algos/Mads/MadsMegaIteration.hpp"

#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelEvaluator.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelInitialization.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelMegaIteration.hpp"

#include "../../../ext/sgtelib/src/Surrogate_Factory.hpp"
//
// Reference: File Sgtelib_Model_Manager.cpp in NOMAD 3.9.1
// Author: Bastien Talgorn

void NOMAD::SgtelibModel::init()
{
    setName("SgtelibModel");
    verifyParentNotNull();

    auto modelFormulation = _runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION");
    auto modelFeasibility = _runParams->getAttributeValue<NOMAD::SgtelibModelFeasibilityType>("SGTELIB_MODEL_FEASIBILITY");
    auto modelDefinition = _runParams->getAttributeValue<NOMAD::ArrayOfString>("SGTELIB_MODEL_DEFINITION");

    if (NOMAD::SgtelibModelFormulationType::EXTERN == modelFormulation)
    {
        // Extern SGTE. Early out.
        return;
    }

    // Check
    if (   (NOMAD::SgtelibModelFormulationType::FS == modelFormulation) 
        || (NOMAD::SgtelibModelFormulationType::EIS == modelFormulation) )
    {
        if (NOMAD::SgtelibModelFeasibilityType::C != modelFeasibility)
        {
            std::cerr << "ERROR : Formulations FS and EIS can only be used with FeasibilityMethod C" << std::endl;
            throw SGTELIB::Exception(__FILE__, __LINE__, "SgtelibModel: SGTELIB_MODEL_FEASIBILITY not valid");
        }
    }

    // Count the number of constraints
    const auto bbot = NOMAD::SgtelibModel::getBBOutputType();
    size_t nbConstraints = NOMAD::getNbConstraints(bbot);
    _nbModels = getNbModels(modelFeasibility, nbConstraints);


    // Init the TrainingSet
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    SGTELIB::Matrix empty_X("empty_X", 0, static_cast<int>(n));
    SGTELIB::Matrix empty_Z("empty_Z", 0, static_cast<int>(_nbModels));
    _trainingSet = std::make_shared<SGTELIB::TrainingSet>(empty_X, empty_Z);

    // Build the Sgtelib Model
    _model = std::shared_ptr<SGTELIB::Surrogate>(SGTELIB::Surrogate_Factory(*_trainingSet, modelDefinition.display()));

    // Instanciate Sgte initialization class
    _initialization = std::make_unique<NOMAD::SgtelibModelInitialization>(this);
}


/*-------------------------*/
/*       Destructor        */
/*-------------------------*/
NOMAD::SgtelibModel::~SgtelibModel()
{
    reset();
}


/*-------------------------*/
/*           reset         */
/*-------------------------*/
void NOMAD::SgtelibModel::reset()
{
    if (nullptr != _model)
    {
        _model.reset();
    }

    if (nullptr != _trainingSet)
    {
        _trainingSet.reset();
    }

    _ready = false;
}


/*-------------------------------------------------------------*/
/*                      isReady                                */
/* Return a boolean to know if the SgtelibModel model has been */
/* build and can be called for predictions.                    */
/*-------------------------------------------------------------*/
bool NOMAD::SgtelibModel::isReady() const
{
    bool retReady = _ready;

    if (!retReady)
    {
        auto modelFormulation = _runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION");
        if (NOMAD::SgtelibModelFormulationType::EXTERN == modelFormulation) 
        {
            // Extern SGTE.
            _ready = true;
            retReady = true;
        }
    }

    if (!retReady)
    {
        if (!_trainingSet)
        {
            throw NOMAD::Exception ( __FILE__, __LINE__ ,
                                    "SgtelibModel::isReady : no training set");
        }

        if ( _trainingSet->is_ready())
        {
            const int pvar = _trainingSet->get_pvar();
            _ready = ( _model->is_ready() ) && ( pvar > 10 );
            retReady = _ready;
        }
    }

    return retReady;
}


/*-------------------------*/
/*           info          */
/*-------------------------*/
void NOMAD::SgtelibModel::info()
{
    std::cout << "  #===================================================== #" << std::endl;
    std::cout << "SgtelibModel::info" << std::endl ;
    std::cout << "SgtelibModel : " << this << std::endl;
    std::cout << "Model : " << _model << std::endl;

    std::cout << "Cache size : " << NOMAD::CacheBase::getInstance()->size() << std::endl;
    std::cout << "Found feasible : " << _foundFeasible << std::endl;

    // Display of the model's bounds.
    std::cout << "Model Bounds, lower bounds : ( " << _modelLowerBound.display() << " ";

    std::cout << ") , upper bounds : ( " << _modelUpperBound.display() << " )" << std::endl;

    std::cout << "Model Extended Bounds, lower bounds : ( " << getExtendedLowerBound().display() << " ";

    std::cout << ") , upper bounds : ( " << getExtendedUpperBound() << " )" << std::endl;

    if ( _ready )
    {
        std::cout << "sgtelibModel model is ready" << std::endl;
    }
    else
    {
        std::cout << "sgtelibModel model is NOT ready" << std::endl;
    }


    std::cout << "  #===================================================== #" << std::endl;
}


// Used during SgtelibModelUpdate step.
// Update the bounds of the model.
void NOMAD::SgtelibModel::setModelBounds(std::shared_ptr<SGTELIB::Matrix> X)
{
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    if (n != (size_t)X->get_nb_cols())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "SgtelibModel::setModelBounds() dimensions do not match");
    }

    int nbDim = X->get_nb_cols();
    int nbPoints = X->get_nb_rows();

    // Build model bounds
    NOMAD::Double lb;
    NOMAD::Double ub;

    for (int j = 0; j < nbDim; j++)
    {
        lb = _modelLowerBound[j];
        ub = _modelUpperBound[j];
        for (int p = 0; p < nbPoints; p++)
        {
            lb = NOMAD::min(lb, NOMAD::Double(X->get(p,j)));
            ub = NOMAD::max(ub, NOMAD::Double(X->get(p,j)));
        }
        _modelLowerBound[j] = lb;
        _modelUpperBound[j] = ub;
    }


} // end setModelBounds


/*------------------------------------------------------------------------*/
/*                          Extended Bounds                               */
/*------------------------------------------------------------------------*/
NOMAD::ArrayOfDouble NOMAD::SgtelibModel::getExtendedLowerBound() const
{
    auto extLowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");

    for (size_t i = 0; i < extLowerBound.size(); i++)
    {
        if (!extLowerBound[i].isDefined())
        {
            extLowerBound[i] = _modelLowerBound[i]
                               - max(Double(10.0), _modelUpperBound[i] - _modelLowerBound[i]);
        }
    }

    return extLowerBound;
} // end getExtendedLowerBound


NOMAD::ArrayOfDouble NOMAD::SgtelibModel::getExtendedUpperBound() const
{
    auto extUpperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");

    for (size_t i = 0; i < extUpperBound.size(); i++)
    {
        if (!extUpperBound[i].isDefined())
        {
            extUpperBound[i] = _modelUpperBound[i]
                               + max(Double(10.0), _modelUpperBound[i] - _modelUpperBound[i]);
        }
    }

    return extUpperBound;
} // end getExtendedUpperBound


// Start is executed when SgtelibModel is used as an algorithm on its own.
void NOMAD::SgtelibModel::startImp()
{

    // Manages initialization among other things.
    NOMAD::Algorithm::startImp();
    
    
    // Comment to appear at the end of stats lines
    NOMAD::MainStep::setAlgoComment("(SgtelibModel)");


    // Setup EvalPoint success computation to be based on sgte rather than bb.
    NOMAD::ComputeSuccessType::setComputeSuccessTypeFunction(
                                NOMAD::ComputeSuccessType::computeSuccessTypeSgte);

    // There is no upper step, so barrier is not inherited from an Algorithm Ancestor.
    // Barrier was computed in the Initialization step.
    // This barrier is in subspace.
    // X0s are found relative to BB, not SGTE
    _barrierForX0s = _initialization->getBarrier();

}


bool NOMAD::SgtelibModel::runImp()
{
    size_t k = 0;   // Iteration number

    if (!_termination->terminate(k))
    {
        // Barrier constructor automatically finds the best points in the cache.
        // This barrier is not the same as the _barrierForX0s member, which
        // is used for model optimization.
        // This barrier is used for MegaIteration management.
        std::vector<NOMAD::EvalPoint> evalPointList;
        auto barrier = _initialization->getBarrier();
        if (nullptr == barrier)
        {
            auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
            barrier = std::make_shared<NOMAD::Barrier>(hMax, getSubFixedVariable(),
                                                       NOMAD::EvalType::BB);
        }
        NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::NOT_EVALUATED;

        // TODO fix this
        /*
        if (nullptr != _megaIteration)
        {
            // Case hot restart
            k       = _megaIteration->getK();
            barrier = _megaIteration->getBarrier();
            megaIterSuccess = _megaIteration->getSuccessType();
        }
        */

        while (!_termination->terminate(k))
        {
            // Create an MegaIteration: manage multiple iterations around
            // different frame centers at the same time.
            NOMAD::SgtelibModelMegaIteration megaIteration(this, k, barrier, megaIterSuccess);
            megaIteration.start();
            megaIteration.run();
            megaIteration.end();

            // Remember these values to construct the next MegaIteration.
            k       = megaIteration.getK();
            barrier = megaIteration.getBarrier();
            megaIterSuccess = megaIteration.NOMAD::MegaIteration::getSuccessType();

            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }
        }


        // member _megaIteration is used for hot restart (read and write)
        // Update it here.
        _megaIteration = std::make_shared<NOMAD::SgtelibModelMegaIteration>(this, k, barrier, megaIterSuccess);
    }


    _termination->start();
    _termination->run();
    _termination->end();

    NOMAD::OutputQueue::Flush();

    return true;
}


void NOMAD::SgtelibModel::endImp()
{
    // Remove any remaining points from eval queue.
    EvcInterface::getEvaluatorControl()->clearQueue();

    // Reset success computation function
    NOMAD::ComputeSuccessType::setComputeSuccessTypeFunction(
                                            NOMAD::ComputeSuccessType::defaultComputeSuccessType);

    NOMAD::MainStep::resetPreviousAlgoComment();
    NOMAD::Algorithm::endImp();
}


/*------------------------------------------------------------------------*/
/*          Compute which formulation must be used in the evalX           */
/*------------------------------------------------------------------------*/
const NOMAD::SgtelibModelFormulationType NOMAD::SgtelibModel::getFormulation() const
{
    auto formulation = _runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION");
    if ( (formulation != NOMAD::SgtelibModelFormulationType::EXTERN) && ( ! _ready) )
    {
        formulation = NOMAD::SgtelibModelFormulationType::D;
    }

    return formulation;
}

// Compute "norm" of deltaMesh
NOMAD::Double NOMAD::SgtelibModel::getDeltaMNorm() const
{
    NOMAD::Double deltaMNorm;

    if (nullptr != _mesh)
    {
        auto deltaMesh = _mesh->getdeltaMeshSize();
        NOMAD::Double squaredNorm = 0;
        for (size_t i = 0; i < deltaMesh.size(); i++)
        {
            squaredNorm += deltaMesh[i] * deltaMesh[i];
        }
        deltaMNorm = squaredNorm.sqrt();
    }

    return deltaMNorm;
}


/*------------------------------------------------------------------------*/
/*          Check that h & f are defined, and if not, correct it          */
/*------------------------------------------------------------------------*/
void NOMAD::SgtelibModel::checkHF(NOMAD::EvalPoint& x) const
{
    NOMAD::Double f = x.getF(NOMAD::EvalType::SGTE);
    NOMAD::Double h = x.getH(NOMAD::EvalType::SGTE);

    if (!f.isDefined())
    {
        // Should not happen
        // was:
        //f = x.get_bb_outputs().get_coord(_p.get_index_obj().front());
        AddOutputWarning("Warning: SgtelibModel::checkHF(): f is not defined and needs to be recomputed.");
    }

    if (!h.isDefined())
    {
        // Also should not happen
        const auto bbo  = x.getEval(NOMAD::EvalType::SGTE)->getBBOutput().getBBOAsArrayOfDouble();
        const auto bbot = NOMAD::SgtelibModel::getBBOutputType();
        NOMAD::SgtelibModelEvaluator::evalH(bbo, bbot, h);
        AddOutputWarning("Warning: SgtelibModel::checkHF(): h is not defined and needs to be recomputed.");
    }

    if ( !f.isDefined() || !h.isDefined() )
    {
        f = NOMAD::INF;
        h = NOMAD::INF;
    }

    x.setF(f, NOMAD::EvalType::SGTE);
    x.setH(h, NOMAD::EvalType::SGTE);
}


/*------------------------------------------------------------------------*/
/*         get fmin from the training set                                 */
/*------------------------------------------------------------------------*/
NOMAD::Double NOMAD::SgtelibModel::getFMin() const
{
    NOMAD::Double fMin;

    if (_trainingSet->is_ready())
    {
        std::cout << "(getFMin : training set is ready:) " << _trainingSet->get_nb_points() << ")" << std::endl;
        fMin = _trainingSet->get_f_min();
    }
    else
    {
        std::cout << "(getFMin : training set is not ready) " << std::endl;
        // Do nothing: fMin remains undefined.
    }

    return fMin;
}//


size_t NOMAD::SgtelibModel::getNbModels(const NOMAD::SgtelibModelFeasibilityType modelFeasibility,
                                        const size_t nbConstraints)
{
    size_t nbModels = -1;

    switch(modelFeasibility)
    {
        case NOMAD::SgtelibModelFeasibilityType::C:
            nbModels = 1 + nbConstraints;
            break;
        case NOMAD::SgtelibModelFeasibilityType::H:
        case NOMAD::SgtelibModelFeasibilityType::B:
        case NOMAD::SgtelibModelFeasibilityType::M:
            nbModels = 2;
            break;
        case NOMAD::SgtelibModelFeasibilityType::UNDEFINED:
            std::cerr << "UNDEFINED SGTELIB_MODEL_FEASIBILITY" << std::endl;
            break;
    }

    return nbModels;
}


std::vector<NOMAD::EvalPoint> NOMAD::SgtelibModel::getX0s() const
{
    std::vector<NOMAD::EvalPoint> x0s;

    if (nullptr != _barrierForX0s)
    {
        x0s = _barrierForX0s->getAllPoints();
    }

    return x0s;
}


// To be used outside of SgtelibModel, e.g., in SgtelibSearchMethod.
NOMAD::EvalPointSet NOMAD::SgtelibModel::createOraclePoints()
{
    // As long as we are managing points using their SGTE evaluation, 
    // setup EvalPoint success computation to be based on SGTE rather than BB.
    // Setting the ComputeSuccessType function ensures that at all steps,
    // we compare oranges with oranges (i.e., SGTE with SGTE).
    NOMAD::ComputeSuccessType::setComputeSuccessTypeFunction(
                                NOMAD::ComputeSuccessType::computeSuccessTypeSgte);

    // Create one MegaIteration. It will not be run. It is used to
    // generate oracle points.
    // For this reason, k and success are irrelevant.
    // Barrier points are used to create Iterations. A sub Mads will be run
    // on every Iteration, using sgte evaluation.
    NOMAD::SgtelibModelMegaIteration megaIteration(this, 0, _barrierForX0s, NOMAD::SuccessType::NOT_EVALUATED);
    megaIteration.generateTrialPoints();

    NOMAD::OutputQueue::Flush();

    // Reset success computation function
    NOMAD::ComputeSuccessType::setComputeSuccessTypeFunction(
                                NOMAD::ComputeSuccessType::defaultComputeSuccessType);

    // The returned EvalPoints are not evaluated by the blackbox, but they will be soon.
    return megaIteration.getTrialPoints();
}

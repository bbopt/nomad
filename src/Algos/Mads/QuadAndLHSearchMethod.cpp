/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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


#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/QuadAndLHSearchMethod.hpp"
#include "../../Algos/Mads/QuadSearchMethod.hpp"
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Math/LHS.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::QuadAndLHSearchMethod::init()
{
    
    setStepType(NOMAD::StepType::SEARCH_METHOD_QUAD_MODEL);
    verifyParentNotNull();

    const auto parentSearch = getParentStep()->getParentOfType<NOMAD::QuadAndLHSearchMethod*>(false);

    
    
    // For some testing, it is possible that _runParams is null or evaluator control is null
    setEnabled((nullptr == parentSearch)
               && (nullptr !=_runParams)
               && _runParams->getAttributeValue<bool>("QUAD_MODEL_AND_LH_SEARCH")
               &&  (nullptr != EvcInterface::getEvaluatorControl()));
#ifndef USE_SGTELIB
    if (isEnabled())
    {
        OUTPUT_INFO_START
        AddOutputInfo(getName() + " cannot be performed because NOMAD is compiled without sgtelib library");
        OUTPUT_INFO_END
        setEnabled(false);
    }
#endif

#ifdef USE_SGTELIB
    // Check that there is exactly one objective
    if (isEnabled())
    {
        auto nbObj = NOMAD::Algorithm::getNbObj();
        if (0 == nbObj)
        {
            OUTPUT_INFO_START
            AddOutputInfo(getName() + " not performed when there is no objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }
        else if (nbObj > 1)
        {
            OUTPUT_INFO_START
            AddOutputInfo(getName() + " not performed on multi-objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }

        auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");
        _displayLevel = modelDisplay.empty()
                            ? NOMAD::OutputLevel::LEVEL_DEBUGDEBUG
                            : NOMAD::OutputLevel::LEVEL_INFO;
    }
#endif
}


bool NOMAD::QuadAndLHSearchMethod::runImp()
{
    if ( isEnabled() )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"QuadAndLHSearchMethod: not fully implemented.");
    }
    
    bool foundBetter = false;
#ifdef USE_SGTELIB
    // The trial points are generated for a feasible frame center and an infeasible one.
    if ( ! _stopReasons->checkTerminate() )
    {
        
        auto nbConsecutiveFail = getSuccessStats().getStatsNbConsecutiveFail();
        
        NOMAD::EvalPointSet trialPts;
        
        if(nbConsecutiveFail > 3) 
        {
    
            auto madsIteration = getParentOfType<MadsIteration*>();
            
            // MegaIteration's barrier member is already in sub dimension.
            auto bestX = madsIteration->getMegaIterationBarrier()->getCurrentIncumbentFeas();
            
            if (nullptr == bestX)
            {
                bestX = madsIteration->getMegaIterationBarrier()->getCurrentIncumbentInf();
            }
            if (nullptr == bestX)
            {
                return false;
            }
            
            auto n = bestX->size();
            
            auto frameSize = madsIteration->getMesh()->getDeltaFrameSize();
            frameSize *= 100;
            
            // LH sampling around best feas and best infeasible.
            NOMAD::ArrayOfDouble lowerBound = *(bestX->getX()) - frameSize ;
            auto upperBound = *(bestX->getX()) + frameSize ;

            NOMAD::Double scaleFactor = sqrt(-log(NOMAD::DEFAULT_EPSILON));
            // Apply Latin Hypercube algorithm (provide frameCenter, deltaFrameSize, and scaleFactor for updating bounds)
            NOMAD::LHS lhs(n, 3*n, lowerBound, upperBound, *bestX, frameSize, scaleFactor);
            auto pointVector = lhs.Sample();
            
            for (const auto& point: pointVector )
            {
                NOMAD::EvalPoint evalPoint(point);
                evalPoint.setPointFrom(std::make_shared<NOMAD::EvalPoint>(*bestX), NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this)); // !!! Point from is a copy of frame center
                evalPoint.addGenStep(getStepType());
                insertTrialPoint(evalPoint);
            }
            if ( ! _stopReasons->checkTerminate() )
            {
                foundBetter = evalTrialPoints(this);
            }
            
        }
        
        // QuadSearchMethod
        auto quadSearch             = std::make_shared<NOMAD::QuadSearchMethod>(this);
        quadSearch->generateTrialPoints();
        
        trialPts = quadSearch->getTrialPoints();

        for (auto evalPoint: trialPts )
        {
            evalPoint.addGenStep(getStepType());
            insertTrialPoint(evalPoint);
        }

        if ( ! _stopReasons->checkTerminate() )
        {
            foundBetter = evalTrialPoints(this);
        }

        // From IterationUtils. Update megaIteration barrier.
        postProcessing();
    }
#endif
    return foundBetter;
    
}
void NOMAD::QuadAndLHSearchMethod::generateTrialPointsFinal()
{

    throw NOMAD::Exception(__FILE__,__LINE__,"QuadAndLHSearchMethod: not fully implemented.");
}


 // end generateTrialPoints

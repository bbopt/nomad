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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoIteration.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoRandom.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Math/RNG.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::TemplateAlgoRandom::init()
{
    setStepType(NOMAD::StepType::ALGORITHM_RANDOM);
    
    
    // Let's get the frame size. The frame size may not be available when no mesh is available. The frame center (best incumbent) may not be available.
    const NOMAD::TemplateAlgoIteration * iter = getParentOfType<NOMAD::TemplateAlgoIteration*>();

    if(nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"An iteration is required.");
    }
    
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    _boxSize = NOMAD::ArrayOfDouble(n,1);
    if (nullptr != iter->getMesh())
    {
        _boxSize = iter->getMesh()->getDeltaFrameSize();
    }

    verifyParentNotNull();
}



void NOMAD::TemplateAlgoRandom::startImp()
{
    
    // Create EvalPoints
    generateTrialPoints();

    if (_iterAncestor->getMesh())
    {
        if (_projectOnMesh && !verifyPointsAreOnMesh(getName()))
        {
            OUTPUT_INFO_START
            AddOutputInfo("At least one trial point is not on mesh. May need investigation if this happens too often.");
            OUTPUT_INFO_END
        }
    }
}


bool NOMAD::TemplateAlgoRandom::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }

    // From IterationUtils. Update megaIteration barrier.
    postProcessing();

    return foundBetter;
}


void NOMAD::TemplateAlgoRandom::generateTrialPointsImp()
{
    
    OUTPUT_INFO_START
    AddOutputInfo("Generate point for " + getName() );
    OUTPUT_INFO_END

    // Clear the trial points.
    clearTrialPoints();
    
    const NOMAD::TemplateAlgoIteration * iter = getParentOfType<NOMAD::TemplateAlgoIteration*>();
    if(nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"An iteration is required.");
    }
    _center = iter->getFrameCenter();
    
    if(nullptr == _center)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"A frame center is required.");
    }
    
    // The pb params handle only variables (fixed variables are not considered)
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    
    auto lb = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto ub = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    auto k = _runParams->getAttributeValue<size_t>("RANDOM_ALGO_DUMMY_FACTOR");
    

    // Creation of points
    // Sample randomly between -(j+1) and j+1
    for (size_t j = 0; j < n*k; j++)
    {
        NOMAD::EvalPoint xt = *_center;
        for (size_t i = 0; i < n; i++)
        {
            xt[i] += RNG::rand(-(j+1.0)*_boxSize[i].todouble()/2.0,(j+1.0)*_boxSize[i].todouble()/2.0);
        }
        
        NOMAD::EvalPointPtr pointFrom = nullptr;
        auto barrier = getMegaIterationBarrier();
        if (nullptr != barrier)
        {
            pointFrom = std::make_shared<NOMAD::EvalPoint>(*(barrier->getFirstPoint())); // Make a copy of Eval Point from the barrier
            xt.setPointFrom(pointFrom, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        }
        
        if (snapPointToBoundsAndProjectOnMesh(xt, lb, ub))
        {
            xt.addGenStep(getStepType());
            bool inserted = insertTrialPoint(xt);
            
            OUTPUT_INFO_START
            std::string s = "xt:";
            s += (inserted) ? " " : " not inserted: ";
            s += xt.display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
    }
    
    if (_iterAncestor->getMesh())
    {
        if (!verifyPointsAreOnMesh(getName()))
        {
            OUTPUT_INFO_START
            AddOutputInfo("xt is not on mesh. May need investigation if this happens too often.");
            OUTPUT_INFO_END
        }
    }
}


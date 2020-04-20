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

#include "../../Algos/NelderMead/NMShrink.hpp"
#include "../../Algos/NelderMead/NMIteration.hpp"
#include "../../Algos/NelderMead/NMUpdate.hpp"


void NOMAD::NMShrink::init()
{
    _name = getAlgoName() + "Shrink";
    
    _currentStepType = NOMAD::NMStepType::SHRINK;
    
    _gamma = _runParams->getAttributeValue<NOMAD::Double>("NM_GAMMA");

    if ( _gamma <= 0.0 || _gamma > 1 )
        throw NOMAD::Exception(__FILE__,__LINE__,"Gamma value not compatible with shrink");
    
    verifyParentNotNull();

}



void NOMAD::NMShrink::startImp()
{
    
    // Update main barrier.
    NOMAD::NMUpdate update( this );
    update.start();
    update.run();
    update.end();

    // Create EvalPoints
    generateTrialPoints();
    
}


bool NOMAD::NMShrink::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }
    if ( getNbEvalPointsThatNeededEval() == 0 )
        setStopReason( );
    
    return foundBetter;
}


void NOMAD::NMShrink::endImp()
{    
    // From IterationUtils
    postProcessing(getEvalType());
    
}


void NOMAD::NMShrink::generateTrialPoints ()
{

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    size_t YSize = _nmY->size();
    
    OUTPUT_INFO_START
    AddOutputInfo("Shrink simplex with " + _name +" (gamma=" + _gamma.tostring() +") with " + std::to_string(YSize) + " points.");
    OUTPUT_INFO_END

    
    // Shrink simplex Y
    std::set<NOMAD::EvalPoint>::const_iterator it = _nmY->begin() ;
    const NOMAD::EvalPoint & y0 = (*it);
    int i=0;
    for ( ; it !=_nmY->end(); ++it, ++i )
    {
        OUTPUT_INFO_START
        AddOutputInfo("y" + std::to_string(i) + ": " + (*it).display() );
        OUTPUT_INFO_END
        
        NOMAD::Point yi(n,0);
        for (size_t k = 0 ; k < n ; ++k )
        {
            yi[k] = (*it)[k];
        }
        NOMAD::Point y(n,0);
        for (size_t k = 0 ; k < n ; ++k )
        {
            y[k] = y0[k] + _gamma*(yi[k]-y0[k]);
        }
        
        
        // Shrink should not generate points outside the bounds
        // Shrink is not used with mesh
        // ----> No need to use snapPointToBoundsAndProjectOnMesh
        
        // New EvalPoint to be evaluated.
        // Add it to the list.
        bool inserted = insertTrialPoint(NOMAD::EvalPoint(y));
        
        OUTPUT_INFO_START
        std::string s = "xr:";
        s += (inserted) ? " " : " not inserted: ";
        s += y.display();
        AddOutputInfo(s);
        OUTPUT_INFO_END
        
        // Test for too small shrink (not the first point)
        if (i > 0 && y == yi )
        {
            OUTPUT_INFO_START
            AddOutputInfo("Shrink point to close to simplex point.");
            OUTPUT_INFO_END
            setStopReason( );
            clearTrialPoints();
            return;
        }
    }
}


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

#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsMegaIteration.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::DMultiMadsMegaIteration::init()
{
    setStepType(NOMAD::StepType::MEGA_ITERATION);

    
    // Create a DMultiMads iteration.
    // Set current best point during the Update step: use xFeas OR xInf if XFeas is not available.
    // During MegaIteration Run step we use a single iteration object with several start, run, end for the various iterations of the algorithm.
    
    _dMultiMadsIteration = std::make_shared<NOMAD::DMultiMadsIteration>(this,
                                                                        nullptr, /* the best point will be updated before start */
                                                                        0, /* the counter will be updated at start */
                                                                        _initialMesh);
    
}


void NOMAD::DMultiMadsMegaIteration::startImp()
{
    
    if ( ! _stopReasons->checkTerminate() )
    {
    
        OUTPUT_DEBUG_START
        auto frameCenter = _dMultiMadsIteration->getFrameCenter();
        AddOutputDebug("Previous frame center: " + (frameCenter ? frameCenter->display() : "NULL"));
        OUTPUT_DEBUG_END
    }
}


bool NOMAD::DMultiMadsMegaIteration::runImp()
{
    bool successful = false;
    std::string s;

    if ( _stopReasons->checkTerminate() )
    {
        OUTPUT_DEBUG_START
        s = getName() + ": stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
        return false;
    }

    if ( _dMultiMadsIteration == nullptr )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "No iteration to run");
    }

    if (! _stopReasons->checkTerminate())
    {
        _dMultiMadsIteration->start();

        successful = _dMultiMadsIteration->run(); // Is this iteration successful

        _dMultiMadsIteration->end();

        if (successful)
        {
            OUTPUT_DEBUG_START
            s = getName() + ": new success " + NOMAD::enumStr(getSuccessType());
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }
            
        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }
    }
    OUTPUT_DEBUG_START
    // Display MegaIteration's stop reason
    AddOutputDebug(getName() + " stop reason set to: " + _stopReasons->getStopReasonAsString());
    OUTPUT_DEBUG_END

    // MegaIteration is a success if either a better xFeas or
    // a dominating or partial success for xInf was found.
    // See Algorithm 12.2 from DFBO.

    // return true if we have a partial or full success.
    return successful;
}


void NOMAD::DMultiMadsMegaIteration::display( std::ostream& os ) const
{
    NOMAD::MegaIteration::display(os);
}


void NOMAD::DMultiMadsMegaIteration::read(  std::istream& is )
{
    
    throw NOMAD::Exception(__FILE__,__LINE__,"DMultiMadsMegationIteration is not yet available.");
//    // Set up structures to gather member info
//    size_t k=0;
//    // Read line by line
//    std::string name;
//    while (is >> name && is.good() && !is.eof())
//    {
//        if ("ITERATION_COUNT" == name)
//        {
//            is >> k;
//        }
//        else if ("BARRIER" == name)
//        {
//            if (nullptr != _barrier)
//            {
//                is >> *_barrier;
//            }
//            else
//            {
//                std::string err = "Error: Reading a Barrier onto a NULL pointer";
//                std::cerr << err;
//            }
//        }
//        else
//        {
//            for (size_t i = 0; i < name.size(); i++)
//            {
//                is.unget();
//            }
//            break;
//        }
//    }
//
//    setK(k);
}




std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::DMultiMadsMegaIteration& megaIteration)
{
    megaIteration.display ( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::DMultiMadsMegaIteration& megaIteration)
{

    megaIteration.read( is );
    return is;

}

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

#include <sstream>

#include "../../Algos/Mads/MadsMegaIteration.hpp"

#include "../../Algos/NelderMead/NMMegaIteration.hpp"
#include "../../Algos/NelderMead/NMInitializeSimplex.hpp"

#include "../../Algos/EvcInterface.hpp"

void NOMAD::NMMegaIteration::init()
{
    _name = getAlgoName() + NOMAD::MegaIteration::getName();
}

void NOMAD::NMMegaIteration::startImp()
{
    // Create a Nelder Mead iteration for a frame center.
    // Use xFeas or xInf if XFeas is not available.
    // During NM we use a single iteration object with several start, run, end for the various iterations of the algorithm.


    if ( ! _stopReasons->checkTerminate() )
    {
        // MegaIteration's barrier member is already in sub dimension.
        auto bestXFeas = _barrier->getFirstXFeas();
        auto bestXInf  = _barrier->getFirstXInf();

        // Note: getParentOfType with argument "false" gets over the "Algorithm" parents.
        // Here, we are looking for a MadsMegaIteration which would be ancestor of
        // the NM (Algorithm) parent.
        auto madsMegaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>(getParentOfType<NOMAD::MadsMegaIteration*>(false));
        std::shared_ptr<NOMAD::MeshBase> mesh = nullptr;

        if ( madsMegaIter != nullptr )
        {
            mesh = madsMegaIter->getMesh();
        }

        if (nullptr != bestXFeas)
        {
            _nmIteration = std::make_shared<NOMAD::NMIteration>(this,
                                    std::make_shared<NOMAD::EvalPoint>(*bestXFeas),
                                    0, /*counter at 0 for start */
                                    mesh);
        }
        else if (nullptr != bestXInf)
        {
            _nmIteration = std::make_shared<NOMAD::NMIteration>(this,
                                    std::make_shared<NOMAD::EvalPoint>(*bestXInf),
                                    0, /*counter at 0 for start */
                                    mesh);
        }

        auto frameCenter = _nmIteration->getFrameCenter();
        AddOutputDebug("Frame center: " + frameCenter->display());
        auto previousFrameCenter = frameCenter->getPointFrom();
        AddOutputDebug("Previous frame center: " + (previousFrameCenter ? previousFrameCenter->display() : "NULL"));
    }
}


bool NOMAD::NMMegaIteration::runImp()
{
    bool successful = false;
    std::string s;
    
    if ( _stopReasons->checkTerminate() )
    {
        s = _name + ": stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
        return false;
    }
    
    if ( _nmIteration == nullptr )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "No iteration to run");
    }
    
    const size_t maxIter = _runParams->getAttributeValue<size_t>("MAX_ITERATION_PER_MEGAITERATION");
    size_t nbMegaIter = 0;
    while ( ! _stopReasons->checkTerminate() && nbMegaIter < maxIter )
    {
        _nmIteration->start();
        
        bool iterSuccessful = _nmIteration->run();          // Is this iteration successful
        successful = iterSuccessful || successful;  // Is the whole MegaIteration successful
        
        _nmIteration->end();
        
        if (iterSuccessful)
        {
            s = _name + ": new success " + NOMAD::enumStr(getSuccessType());
            AddOutputDebug(s);
        }
        
        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }
        
        nbMegaIter++;
    }
    // Display MegaIteration's stop reason
    AddOutputDebug(_name + " stop reason set to: " + _stopReasons->getStopReasonAsString());
    
    // MegaIteration is a success if either a better xFeas or
    // a dominating or partial success for xInf was found.
    // See Algorithm 12.2 from DFBO.
    
    // return true if we have a partial or full success.
    return successful;
}


void NOMAD::NMMegaIteration::display( std::ostream& os ) const
{
// TODO display simplex
//    os << "MAIN_MESH " << std::endl;
//    os << *_mainMesh ;
    NOMAD::MegaIteration::display(os);
}


void NOMAD::NMMegaIteration::read(  std::istream& is )
{
    // TODO read simplex
    // Set up structures to gather member info
    size_t k=0;
    // Read line by line
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
//        if ("MAIN_MESH" == name)
//        {
//            if (nullptr != _mainMesh)
//            {
//                is >> *_mainMesh;
//            }
//            else
//            {
//                std::string err = "Error: Reading a mesh onto a NULL pointer";
//                std::cerr << err;
//            }
//        }
//        else
        if ("ITERATION_COUNT" == name)
        {
            is >> k;
        }
        else if ("BARRIER" == name)
        {
            if (nullptr != _barrier)
            {
                is >> *_barrier;
            }
            else
            {
                std::string err = "Error: Reading a Barrier onto a NULL pointer";
                std::cerr << err;
            }
        }
        else
        {
            for (size_t i = 0; i < name.size(); i++)
            {
                is.unget();
            }
            break;
        }
    }

    setK(k);
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::NMMegaIteration& megaIteration)
{
    megaIteration.display ( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::NMMegaIteration& megaIteration)
{

    megaIteration.read( is );
    return is;

}

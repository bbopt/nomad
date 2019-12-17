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

#include "../Algos/MegaIteration.hpp"
#include "../Algos/EvcInterface.hpp"


void NOMAD::MegaIteration::init()
{
    _name = "MegaIteration " + std::to_string(_k);
    verifyParentNotNull();
}


size_t NOMAD::MegaIteration::getNextK() const
{
    size_t nextK;
    if (_nbIterRun > 0)
    {
        // Next MegaIteration number will start where the Iterations ended.
        nextK = _k + _nbIterRun;
    }
    else
    {
        // Must increment MegaIteration number, cannot have
        // twice the same number.
        nextK = _k + 1;
    }

    return nextK;
}


void NOMAD::MegaIteration::endImp()
{
    if (_runParams->getAttributeValue<bool>("USER_CALLS_ENABLED"))
    {
        bool stop = false;
        runCallback(NOMAD::CallbackType::MEGA_ITERATION_END, *this, stop);
        if (!_stopReasons->checkTerminate() && stop)
        {
            _stopReasons->set(NOMAD::BaseStopType::USER_STOPPED);
        }
    }

    _iterList.clear();
}


void NOMAD::MegaIteration::computeMaxXFeasXInf(size_t &maxXFeas, size_t &maxXInf)
{
    const size_t maxIter = _runParams->getAttributeValue<size_t>("MAX_ITERATION_PER_MEGAITERATION");
    const size_t maxXFeas0 = maxXFeas;
    const size_t maxXInf0 = maxXInf;

    // If maxXFeas + maxXInf does not exceed maxIter, do nothing.
    if (maxXFeas + maxXInf > maxIter)
    {
        if (maxXFeas <= maxIter / 2)
        {
            // Use all xFeas, and the remaining in xInf.
            maxXInf = maxIter - maxXFeas;
        }
        else if (maxXInf < maxIter / 2)
        {
            // Use all xInf, and the remaining in xFeas.
            maxXFeas = maxIter - maxXInf;
        }
        else
        {
            // Both the number of xFeas and xInf is over.
            // Use half of maxIter for xFeas and half for xInf.
            maxXInf = maxIter / 2;
            maxXFeas = maxIter - maxXInf;
        }
        if (maxXFeas + maxXInf > maxIter)
        {
            // This case should not happen and should be debugged.
            std::cerr << "Warning: Bad computation in computeMaxXFeasXInf. maxIter = " << maxIter << " maxXFeas = " << maxXFeas << " (was " << maxXFeas0 << ") maxXInf = " << maxXInf << " (was " << maxXInf0 << ")" << std::endl;
        }
    }
}


void NOMAD::MegaIteration::read(std::istream& is)
{
    // Set up structures to gather member info
    size_t k = 0;

    // Read line by line
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
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
                std::cerr << err << std::endl;
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


void NOMAD::MegaIteration::display(  std::ostream& os ) const
{
    os << "ITERATION_COUNT " << _k << std::endl;
    os << "BARRIER " << std::endl;
    os << *_barrier;
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::MegaIteration& megaIteration)
{

    megaIteration.display( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::MegaIteration& megaIteration)
{
    megaIteration.read(is);
    return is;
}

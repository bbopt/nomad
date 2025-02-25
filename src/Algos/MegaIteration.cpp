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

#include "../Algos/MegaIteration.hpp"
#include "../Algos/EvcInterface.hpp"

// Constructor
NOMAD::MegaIteration::MegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<BarrierBase> barrier,
                              SuccessType success)
  : Step(parentStep),
    _barrier(barrier),
    _k(k),
    _megaIterationSuccess(success)
{
    if (nullptr == _barrier)
    {
        throw NOMAD::StepException(__FILE__, __LINE__, "MegaIteration constructor: barrier must not be NULL.", this);
    }

    init();
}

void NOMAD::MegaIteration::startImp()
{
    if (_runParams->getAttributeValue<bool>("USER_CALLS_ENABLED"))
    {
        bool stop = false;
        runCallback(NOMAD::CallbackType::MEGA_ITERATION_START, *this, stop);
        if (!_stopReasons->checkTerminate() && stop)
        {
            _stopReasons->set(NOMAD::BaseStopType::USER_GLOBAL_STOP);
        }
    }
}


void NOMAD::MegaIteration::init()
{
    setStepType(NOMAD::StepType::MEGA_ITERATION);
    verifyParentNotNull();
}


std::string NOMAD::MegaIteration::getName() const
{
    return getAlgoName() + NOMAD::stepTypeToString(_stepType) + " " + std::to_string(_k);
}


void NOMAD::MegaIteration::endImp()
{
    if (_runParams->getAttributeValue<bool>("USER_CALLS_ENABLED"))
    {
        // Run callback and set stop reason if overall stop is requested
        bool stop = false;
        runCallback(NOMAD::CallbackType::MEGA_ITERATION_END, *this, stop);
        if (!_stopReasons->checkTerminate() && stop)
        {
            _stopReasons->set(NOMAD::BaseStopType::USER_GLOBAL_STOP);
        }
        
        
        // Reset user iteration stop reason
        if (_stopReasons->testIf(NOMAD::IterStopType::USER_ITER_STOP))
        {
            _stopReasons->set(NOMAD::IterStopType::STARTED);
        }
    }

    // If last megaIteration, update hmax and barrier incumbents to ensure consistent display
    // (if stopping criteria reached during search for ex, hmax and incumbents have not been updated if it was not a full success)
    if(_stopReasons->checkTerminate())
    {
            // Update of hmax and feasible/infeasible incumbents
            bool barrierModified = false;
            std::vector<NOMAD::EvalPoint> evalPointList;  // eval point list empty just to call updateWithPoints
            if(_barrier!= nullptr)
            {
                barrierModified = _barrier->updateWithPoints(
                                    evalPointList,
                                    false /* not used by progressive barrier */,
                                    true /* update incumbents and hmax*/ );
            }
    }

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
            std::cout << "Warning: Bad computation in computeMaxXFeasXInf. maxIter = " << maxIter << " maxXFeas = " << maxXFeas << " (was " << maxXFeas0 << ") maxXInf = " << maxXInf << " (was " << maxXInf0 << ")" << std::endl;
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
                throw NOMAD::Exception(__FILE__,__LINE__, err);
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

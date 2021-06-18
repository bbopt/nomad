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

#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelMegaIteration.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::QuadModelMegaIteration::init()
{
    setStepType(NOMAD::StepType::MEGA_ITERATION);
}


NOMAD::QuadModelMegaIteration::~QuadModelMegaIteration()
{
    // Clear model info from cache.
    // Very important so we don't have false info in a later MegaIteration.
    NOMAD::CacheBase::getInstance()->clearModelEval(NOMAD::getThreadNum());
}


void NOMAD::QuadModelMegaIteration::startImp()
{
    // Create an iteration for a frame center.
    // Use xFeas or xInf if XFeas is not available.
    // Use a single iteration object with several start, run, end for the various iterations of the algorithm.

    // See issue (feature) #384

    if ( ! _stopReasons->checkTerminate() )
    {
        // MegaIteration's barrier member is already in sub dimension.
        auto bestXFeas = _barrier->getFirstXFeas();
        auto bestXInf  = _barrier->getFirstXInf();

        if (nullptr != bestXFeas)
        {
            auto sqmIteration = std::make_shared<NOMAD::QuadModelIteration>(
                                            this,
                                            bestXFeas,
                                            0,/*counter at 0 for start */
                                            nullptr);
            _iterList.push_back(sqmIteration);

        }
        else if (nullptr != bestXInf)
        {
            auto sqmIteration = std::make_shared<NOMAD::QuadModelIteration>(
                                            this,
                                            bestXInf,
                                            0,  /*counter at 0 for start */
                                            nullptr);
            _iterList.push_back(sqmIteration);
        }

        size_t nbIter = _iterList.size();

        AddOutputInfo(getName() + " has " + NOMAD::itos(nbIter) + " iteration" + ((nbIter > 1)? "s" : "") + ".");

        AddOutputDebug("Iterations generated:");
        for (size_t i = 0; i < nbIter; i++)
        {
            auto sqmIteration = _iterList[i];
            if ( sqmIteration == nullptr )
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "Invalid shared pointer");
            }

            AddOutputDebug( _iterList[i]->getName());
            // Ensure we get frame center from a QuadModelIteration.
            auto frameCenter = sqmIteration->getFrameCenter();
            AddOutputDebug("Frame center: " + frameCenter->display());
            auto previousFrameCenter = frameCenter->getPointFrom();
            AddOutputDebug("Previous frame center: " + (previousFrameCenter ? previousFrameCenter->display() : "NULL"));

            if (nullptr != sqmIteration->getMesh())
            {
                NOMAD::ArrayOfDouble meshSize  = sqmIteration->getMesh()->getdeltaMeshSize();
                NOMAD::ArrayOfDouble frameSize = sqmIteration->getMesh()->getDeltaFrameSize();

                AddOutputDebug("Mesh size:  " + meshSize.display());
                AddOutputDebug("Frame size: " + frameSize.display());
            }

            NOMAD::OutputQueue::Flush();
        }
    }
}


bool NOMAD::QuadModelMegaIteration::runImp()
{
    bool successful = false;
    std::string s;

    if (_iterList.empty())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "No iterations to run");
    }


    for (size_t i = 0; i < _iterList.size(); i++)
    {

        auto sqmIteration = _iterList[i];
        if ( sqmIteration == nullptr )
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "No iteration to run");
        }

        if (!_stopReasons->checkTerminate())
        {
            sqmIteration->start();

            bool iterSuccessful = sqmIteration->run();          // Is this iteration successful
            successful = iterSuccessful || successful;  // Is the whole MegaIteration successful

            sqmIteration->end();

            if (iterSuccessful)
            {
                s = getName() + ": new success " + NOMAD::enumStr(getSuccessType());
                AddOutputDebug(s);
            }

            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }
        }
    }
    // Display MegaIteration's stop reason
    AddOutputDebug(getName() + " stop reason set to: " + _stopReasons->getStopReasonAsString());


    // MegaIteration is a success if either a better xFeas or
    // a dominating or partial success for xInf was found.
    // See Algorithm 12.2 from DFBO.
    // return true if we have a partial or full success.
    return successful;
}

void NOMAD::QuadModelMegaIteration::endImp()
{
    // Clear model info from cache.
    // Very important so we don't have false info in a later MegaIteration.
    NOMAD::CacheBase::getInstance()->clearModelEval(NOMAD::getThreadNum());
    NOMAD::MegaIteration::endImp();
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::QuadModelMegaIteration& megaIteration)
{
    megaIteration.display ( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::QuadModelMegaIteration& megaIteration)
{

    megaIteration.read( is );
    return is;

}

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
#ifndef __NOMAD400_POLL__
#define __NOMAD400_POLL__

#include <set>

#include "../../Algos/Mads/MadsIterationUtils.hpp"

#include "../../Eval/EvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the poll (step 3) of MADS algorithm.
/**
 Generate the trial points (Poll::startImp), launch evaluation (Poll::runImp) and postprocecssing (Poll::endImp).
 */
class Poll: public Step , public MadsIterationUtils
{


public:
    /// Constructor
    /**
     \param parentStep The parent of this poll step
     */
    explicit Poll(const Step* parentStep)
      : Step(parentStep),
        MadsIterationUtils(parentStep)
    {
        init();
    }
    virtual ~Poll() {}

    /// Generate new points to evaluate
    /**
     The trial points are obtained by:
        - adding poll directions (Poll::setPollDirections) to the poll center (frame center).
        - snaping points (and directions) to bounds.
        - projecting points on mesh.
     */
    void generateTrialPoints() override ;

    
private:
    /// Helper for constructor
    void init();

    /// Implementation for start tasks for MADS poll.
    /**
     Call to generate trial points and test for mesh precision
     */
    virtual void    startImp() override ;
    
    /// Implementation for run tasks for MADS poll.
    /**
     Start trial points evaluation.
     \return Flag \c true if found better solution \c false otherwise.
     */
    virtual bool    runImp() override;
    
    /// Implementation for end tasks for MADS poll.
    /**
     Call the IterationUtils::postProcessing of the points.
     */
    virtual void    endImp() override ;
    
    /*------------------------------*/
    /* Private methods used by poll */
    /*------------------------------*/
    
    /// Helper for Poll::generateTrialPoints
    /**
     - A single direction on unit n-sphere is computed (Poll::computeDirOnUnitSphere).
     - This direction is transformed into 2n directions on a unit n-sphere using the householder transformation.
     - The 2n directions are scaled and projected on the mesh.
     
     \param directions  The directions obtained for this poll -- \b OUT.
     */
    void setPollDirections(std::list<Direction> &directions) const;

    /// Helper for Poll::generateTrialPoints
    void householder(const Direction &dir, bool completeTo2n, Direction ** H) const;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_POLL__

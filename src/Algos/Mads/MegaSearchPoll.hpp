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
#ifndef __NOMAD_4_0_MEGASEARCHPOLL__
#define __NOMAD_4_0_MEGASEARCHPOLL__

#include "../../Algos/IterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the mega search and poll of MADS
/**
 Calling the start function generates search and poll trial points at the same time before starting evaluation.
 Calling the run function starts the evaluaions.
 The postprocessing is performed when calling the end funcion.
 */
class MegaSearchPoll: public Step, public IterationUtils
{
public:
    /// Constructor
    /**
     \param parentStep The parent of this step
     */
    explicit MegaSearchPoll(const Step* parentStep )
      : Step( parentStep ),
        IterationUtils( parentStep )
    {
        init();
    }

    // Destructor
    virtual ~MegaSearchPoll()
    {
    }

private:

    /// Generate the trial poins for the search and poll steps.
    /**
     Call MegaSearchPoll::generateTrialPoints.
     */
    virtual void    startImp() override;

    ///Start evaluations
    virtual bool    runImp() override;

    /**
     Call for postprocessing: computation of a new hMax and update of the barrier.
     */
    virtual void    endImp() override ;

    void init();

    /// Generate new points to evaluate
    /**
     The trial points are produced using poll and search. The duplicates are removed and they are merged all together.
     */
    void generateTrialPoints() override ;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_MEGASEARCHPOLL__

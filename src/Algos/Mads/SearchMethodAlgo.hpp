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
#ifndef __NOMAD_4_0_SEARCHMETHODALGO__
#define __NOMAD_4_0_SEARCHMETHODALGO__

#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for generic search method of MADS. Run by Search.
class SearchMethodAlgo: public SearchMethodBase
{


public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SearchMethodAlgo( const Step* parentStep )
    : SearchMethodBase( parentStep ) {}

    /**
         An empty (disabled) startImp is required for a search method that launches an iterative algorithm during the run.

     */
    void startImp() override {}

    /**
     - Pure virtual function.
     - This function must be implemented for algo based search methods that can perform several iterations.
     - The derived runImp implementation of this function launches the sequence of start, run and end on an algorithm.
     - When running the algorithm, evaluations must be performed.
     - Calling this function implies that ::generateTrialPointsImp is not called.
     - This function is used only when the option to generate all points before evaluation is disabled, that is the ::generateTrialPointsImp is not called.
     */
    virtual bool runImp() override = 0 ;

    /**
     - Pure virtual function.
     - This function must be implemented for algo based search methods that can perform a single iteration for generating points.
     - This function is used only when the option to generate all points before evaluation is enabled, that is the runImp is not called.
     - Evaluations are automatically performed when running SearchMethodSimple::runImp.
     */
    virtual void generateTrialPointsImp() override = 0 ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_SEARCHMETHODALGO__


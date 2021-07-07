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
#ifndef __NOMAD_4_0_VNSSEARCHMETHOD__
#define __NOMAD_4_0_VNSSEARCHMETHOD__

#include "../../Algos/Mads/SearchMethodAlgo.hpp"
#include "../../Algos/VNSMads/VNS.hpp"

#include "../../nomad_nsbegin.hpp"

/// Implementation of VNS Mads search method
class VNSSearchMethod final: public SearchMethodAlgo
{
private:
    OutputLevel _displayLevel;

    DLL_ALGO_API static Point   _refFrameCenter;    ///< The reference frame center. If frame center same as reference, do not perform search
    DLL_ALGO_API static size_t  _bbEvalByVNS;       ///< Counter of VNS evals, used to trigger VNS.
    
    DLL_ALGO_API static size_t  _nbVNSSearchRuns;   ///< Counter of VNS runs (for display).
    
/*----------------------------------------------------------------------*/


public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit VNSSearchMethod(const Step* parentStep)
      : SearchMethodAlgo(parentStep),
        _displayLevel(OutputLevel::LEVEL_NORMAL)
    {
        init();
    }

    /// Reset the reference frame center and the number of evals done by VNS so far.
    static void reset();
    
private:
    void init();

    bool runImp() override;

    ///Generate new points (no evaluation)
    /**
     \copydoc SearchMethod::generateTrialPointsImp \n
     This function is used only when a VNS MADS search with
     the option to generate all points before evaluation. It performs a single
     Mads iteration (search and poll) around the best incumbent points in the Barrier.
     */
    void generateTrialPointsImp() override;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_VNSSEARCHMETHOD__


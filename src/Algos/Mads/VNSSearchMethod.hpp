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
#ifndef __NOMAD_4_4_VNSSEARCHMETHOD__
#define __NOMAD_4_4_VNSSEARCHMETHOD__

#include <string>

#include "../../Algos/Mads/SearchMethodAlgo.hpp"
#include "../../Algos/VNSMads/VNS.hpp"

#include "../../nomad_nsbegin.hpp"

/// Implementation of VNS Mads search method
class VNSSearchMethod final: public SearchMethodAlgo
{
private:
    OutputLevel _displayLevel;

    Point   _refFrameCenter;    ///< The reference frame center for the last call. If frame center same as reference, do not perform search.
        
    double _trigger; ///< Evaluation ratio (vns evals vs all evals) to trigger vns search

    bool _useSurrogate; ///< Flag for surrogate or bb use.
    
    /**
        The algorithm used by the search method.
     */
    std::unique_ptr<VNS> _vnsAlgo;
    
    /**
        VNS has its own stop reasons
     */
    std::shared_ptr<NOMAD::AlgoStopReasons<NOMAD::VNSStopType>>        _vnsStopReasons;

    
/*----------------------------------------------------------------------*/


public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit VNSSearchMethod(const Step* parentStep)
      : SearchMethodAlgo(parentStep),
        _displayLevel(OutputLevel::LEVEL_NORMAL),
        _useSurrogate(false),
        _vnsAlgo(nullptr)
    {
        init();
    }

private:
    void init();
    
    bool runImp() override;

    ///Generate new points (no evaluation)
    /*! \copydoc SearchMethodBase::generateTrialPointsFinal() /n
     * This function is used only when a VNS MADS search with
     * the option to generate all points before evaluation. It performs a single
     * Mads iteration (search and poll) around the best incumbent points in the Barrier.
     */
    void generateTrialPointsFinal() override;


    
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_VNSSEARCHMETHOD__


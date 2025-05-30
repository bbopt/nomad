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
#ifndef __NOMAD_4_5_DMULTIMADSQUADMODSEARCHMETHOD__
#define __NOMAD_4_5_DMULTIMADSQUADMODSEARCHMETHOD__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/DMultiMads/DMultiMadsBarrier.hpp"
#include "../../Algos/Mads/SearchMethodSimple.hpp"



#include "../../nomad_nsbegin.hpp"

/// Class to perform a Search method for DMultiMads using Quad model optimization algorithm.
/**
 Quadratic model optimization works only for single objective. We need to modify the multi objective problem
 into a single objective problem before launching quadratic model optimization.
 Once quad model optim is done the DMultiMads barrier
 must be updated with the new points.
 The regular QuadModelSearchMethod is disabled when running DMultiMads.
 */
class DMultiMadsQuadModSearchMethod final : public SearchMethodSimple
{
private:
    
    ComputeType _ref_compute_type;
    
    bool _flagPriorCombineObjsForModel;

    bool _use_dom_strategy;
    
public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit DMultiMadsQuadModSearchMethod(const Step* parentStep)
      : SearchMethodSimple(parentStep)
    {
        init();
    }

private:

    /// Helper for constructor.
    /**
     Test if the quad model search is enabled or not. 
     */
    void init();

    /// Generate new points (no evaluation)
    /**
     \copydoc SearchMethodAlgo::generateTrialPointsFinal 
     
     Perform one quad model optimization to produce trial points.
     */
     virtual void generateTrialPointsFinal() override;
    
    // Helpers for generateTrialPoints

    // MultiMads strategy
    void prepareSingleObjectiveRun(const NOMAD::ArrayOfDouble& ref);
    void prepareMultiMadsRun(const NOMAD::ArrayOfDouble& ref);
    void runMultiMadsStrategy();

    // DoM strategy
    void runDoMStrategy();

    NOMAD::ArrayOfDouble computeReferencePoint(const NOMAD::DMultiMadsBarrier& barrier) const;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_DMULTIMADSQUADMODSEARCHMETHOD__


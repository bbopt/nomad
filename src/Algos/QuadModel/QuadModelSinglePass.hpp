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
#ifndef __NOMAD_4_0_QUAD_MODEL_SINGLE_PASS__
#define __NOMAD_4_0_QUAD_MODEL_SINGLE_PASS__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelIterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/**
 Final class to generate points for single pass on a given frame center. Used as Search method.
 The start, run and end tasks are empty. No evaluations are performed.
 The QuadModelSinglePass::generateTrialPoints function manages the creation process. The sample set to build the quad model is created by calling QuadModelIteration::startImp(). The points are not projected on mesh (done in SearchMethodBase).
 */
class QuadModelSinglePass final: public QuadModelIteration, public QuadModelIterationUtils
{

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param frameCenter    The "center" around which to construct the training set  -- \b IN.
     \param madsMesh           The Mads mesh for constructing the training set (no snap on mesh is performed) -- \b IN.
     */
    explicit QuadModelSinglePass(const Step* parentStep,
                                 const std::shared_ptr<EvalPoint>& frameCenter,
                                 const std::shared_ptr<MeshBase>& madsMesh)
      : QuadModelIteration(parentStep, frameCenter, 0, madsMesh),
        QuadModelIterationUtils(parentStep)
    {
        _stopReasons = std::make_shared<AlgoStopReasons<ModelStopType>>();
    }
    // No Destructor needed - keep defaults.

    /// Implementation of start task. Nothing to do.
    void startImp() override {}

    /// Implementation of run task. Nothing to do.
    bool runImp() override { return  false;}

    /// Implementation of run task. Nothing to do.
    void endImp() override {}

    /// Generate trial points
    /**
     - Update the quadratic model.
     - Optimize the quadratic model problem.
     - Insert the best feasible and best infeasible (if available) as trial points.
     */
    void generateTrialPoints() override;

};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_NMALLREFLECTIVE__

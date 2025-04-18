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
#ifndef __NOMAD_4_5_QPSOLVERALGOSINGLEPASS__
#define __NOMAD_4_5_QPSOLVERALGOSINGLEPASS__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverAlgoIteration.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelIterationUtils.hpp"
#include "../../Algos/IterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/**
 Class to generate points with a single pass of QP solver generation (no BB evaluation is performed).
 The points are projected on mesh in the search method.
 */
class QPSolverAlgoSinglePass: public QuadModelIteration, public QuadModelIterationUtils
{
private:
    
    bool _flagUseScaledModel;   ///< The model can be scaled between [0,1]^n (rotation facilitates setting the bounds for optimization)
    
    const std::vector<Direction> & _scalingDirections;
public:
public:
    /// Constructor
    /**
     \param parentStep                  The parent step of this step -- \b IN.
     \param frameCenter                The "center" around which to construct the training set  -- \b IN.
     \param madsMesh                       The Mads mesh for constructing the training set (no snap on mesh is performed) -- \b IN.
     \param scalingDirections   The directions used for scaling and bounding the model (optional) -- \b IN.
     */
    explicit QPSolverAlgoSinglePass(const Step* parentStep,
                                 const EvalPointPtr frameCenter,
                                 const MeshBasePtr madsMesh,
                                 const std::vector<Direction> & scalingDirections )
      : QuadModelIteration(parentStep, frameCenter, 0, madsMesh, {} /* no trial points */),
        QuadModelIterationUtils(parentStep),
        _scalingDirections(scalingDirections)
    {
        _stopReasons = std::make_shared<AlgoStopReasons<ModelStopType>>();
        _flagUseScaledModel = (_scalingDirections.size() > 0);
        
    }
    
    // No special destructor needed - keep defaults.


private:

    /// Implementation of start tasks.
    /**
     - call the default Iteration::startImp
     - create the model.
     - generate trial points using the model
     - verify that trial points are on mesh.
     */
    void startImp() override {}

    /// Implementation of run task. Nothing to do.
    bool runImp() override { return  false;}

    /// Implementation of run task. Nothing to do.
    void endImp() override {}

    /// Generation of trial points uses the random generator
    void generateTrialPointsImp() override;
    

    
};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_QPSOLVERALGOSINGLEPASS__

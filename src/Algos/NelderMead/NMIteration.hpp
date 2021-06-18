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
#ifndef __NOMAD_4_0_NMITERATION__
#define __NOMAD_4_0_NMITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Algos/MeshBase.hpp"
#include "../../Algos/NelderMead/NMSimplexEvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for Nelder Mead (NM) iterations
/**
 The start function for a Nelder Mead iteration creates the initial simplex using a frame ("simplex") center and the points in cache. \n
 The run function of this class iterates between the different reflective step
 (REFLECT, EXPAND, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION) and the SHRINK step
 if it is required. The function also updates the type of success to pass to the NMMegaIteration (if it exists) and manages the stop reason.
 */
class NMIteration: public Iteration
{
private:
    /// Helper for constructor
    void init ();

    /**
     The simplex is shared among Nelder Mead components. The simplex is initially built by NMInitializeSimplex::run
     */
    std::shared_ptr<NMSimplexEvalPointSet> _nmY ;

    /**
     The simplex "center" at the creation of this iteration.
     The initial simplex is built around this point.
     The frame center of MADS is used when Nelder Mead is used for a search method of MADS.
     */
    const std::shared_ptr<EvalPoint> _simplexCenter;

    /**
     The Mads mesh can be available if Nelder Mead is used as a Search method. If not, it is set to \c nullptr. When available, trials can be projected on it.
     */
    const std::shared_ptr<MeshBase> _madsMesh;

    SuccessType _bestSuccess; ///< The best success obtained during the Nelder Mead iterations.

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param frameCenter        The frame center -- \b IN.
     \param k                  The iteration number -- \b IN.
     \param madsMesh           Mads Mesh for trial point projection (can be null) -- \b IN.
     */
    explicit NMIteration(const Step *parentStep,
                         const std::shared_ptr<EvalPoint> &frameCenter,
                         const size_t k,
                         std::shared_ptr<MeshBase> madsMesh)
      : Iteration(parentStep, k),
        _simplexCenter(frameCenter),
        _madsMesh(madsMesh)
    {
        init();

        // Create an empty simplex to be shared among Nelder Mead components
        _nmY = std::make_shared<NMSimplexEvalPointSet>();
    }


    // Get/Set

    const std::shared_ptr<EvalPoint> getSimplexCenter() const { return _simplexCenter ; }

    const std::shared_ptr<MeshBase> getMesh() const override { return _madsMesh; }

    const std::shared_ptr<NMSimplexEvalPointSet> getY( void ) const { return _nmY; }

protected:

    /// Implementation of run task.
    /**
     Sequential run of Nelder Mead steps among INITIAL, ( REFLECT, EXPANSION, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION ), SHRINK.
     Set the stop reason and also updates the Nelder Mead mega iteration success.
     */
    virtual bool runImp() override ;

    /// Implementation of start task.
    /**
     Update the barrier and create the initial simplex if it is empty.
     */
    virtual void startImp() override;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_NMITERATION__

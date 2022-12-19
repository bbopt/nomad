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
#ifndef __NOMAD_4_3_CSSITERATION__
#define __NOMAD_4_3_CSSITERATION__

#include "../../Algos/CoordinateSearch/CSPoll.hpp"
#include "../../Algos/Iteration.hpp"
#include "../../Algos/MeshBase.hpp"
#include "../../Eval/EvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for CS iteration
/**
 Unlike a Mads iteration that possesses a Search and Poll step, a CS iteration consists of a Poll step.  A stop reasons and a success type are identified.
  Note: it would be possible to add the same search step as in Mads.
 */
class CSIteration: public Iteration
{
private:
    const MeshBasePtr  _mesh;        ///< Mesh on which the points are

    std::unique_ptr<CSPoll> _csPoll;

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param k                  The iteration number -- \b IN.
     \param mesh               The mesh of the iteration -- \b IN.
     */
    explicit CSIteration(const Step *parentStep,
                           const size_t k,
                           const MeshBasePtr mesh)
      : Iteration(parentStep, k),
        _mesh(mesh)
    {
        init();
    }



    // Gets/Sets

    /**
     The CS algorithm iteration possesses a mesh, unlike the base iteration that has none.
     \remark Used by Step::getIterationMesh() to pass the mesh whenever needed
     */
    const MeshBasePtr getMesh() const override { return _mesh; }

    /*---------------------*/
    /* Other class methods */
    /*---------------------*/


private:
    /// Helper for constructor
    void init();

    virtual void startImp() override;

    /// Implementation of the run tasks of CS algorithm.
    /**
     Run a CS  iteration: a sinle CS Poll step is performed.
     */
    virtual bool runImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_3_CSITERATION__

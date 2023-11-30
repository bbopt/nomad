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

#ifndef __NOMAD_4_4_CSMEGAITERATION__
#define __NOMAD_4_4_CSMEGAITERATION__

#include "../../Algos/CoordinateSearch/CSIteration.hpp"
#include "../../Algos/MegaIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the mega iterations of CS.
/**
Manager for CS iterations.
 Steps:
 - Generate a lot of points over multiple meshes, using a single Poll strategy.
 - Evaluate points
 - Post-processing
*/
class CSMegaIteration: public MegaIteration
{
protected:

    /**
     Main mesh that holds the mesh size and frame size  that we would use in the standard CS algorithm
     */
    MeshBasePtr _mainMesh;
    
    std::unique_ptr<CSIteration> _csIteration;

    void init();

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling -- \b IN.
     \param mesh            Mesh on which other Iteration meshes are based -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit CSMegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<BarrierBase> barrier,
                              MeshBasePtr mesh,
                              SuccessType success)
      : MegaIteration(parentStep, k, barrier, success),
        _mainMesh(mesh),
        _csIteration(nullptr)
    {
        init();
    }
    virtual ~CSMegaIteration() {}

    /// Implementation of the start tasks for CS mega iteration.
    /**
     Creates a CSIteration for each frame center and each desired mesh size.
     Use all xFeas and xInf available.
     For now, not using other frame centers.
     */
    virtual void startImp() override ;

    /// Implementation of the run tasks for CS mega iteration.
    /**
     Manages the generation of points: either all poll and search points are generated all together before starting evaluation using the MegaSearchPoll or they are generated using a CSIteration with search and poll separately. A run parameter controls the behavior.
     */
    virtual bool runImp() override;


    const MeshBasePtr getMesh() const override         { return _mainMesh; }
    void setMesh(const MeshBasePtr &mesh)      { _mainMesh = mesh; }

    void read(  std::istream& is ) override;
    void display(  std::ostream& os ) const override ;

};

/**
 Display useful values so that a new MegaIteration could be constructed using these values.
 */
std::ostream& operator<<(std::ostream& os, const CSMegaIteration& megaIteration);

/// Get an MegaIteration values from a stream
std::istream& operator>>(std::istream& is, CSMegaIteration& megaIteration);

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_CSMEGAITERATION__


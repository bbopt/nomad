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
#ifndef __NOMAD_4_0_MADSMEGAITERATION__
#define __NOMAD_4_0_MADSMEGAITERATION__


#include "../../Algos/MegaIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the mega iterations of MADS.
/**
Manager for Mads iterations.
 Steps:
 - Generate a lot of points over multiple meshes, using different Search and Poll strategies.
 - Evaluate points
 - Post-processing

 \note As an hypothesis, the time load is taken by the evaluation,
  which is parallelized over all evaluations simultaneously.
  The iteration generation, including trial points generation,
  has little time load, so they do not need to be parallelized.
  It is also preferable to keep parallelization to the only place where
  it matters the most to avoid errors.
  There is no parallelization at the algorithmic level.
  Algorithms are run in main thread(s) only; Secundary threads are available for evaluations.
*/
class MadsMegaIteration: public MegaIteration
{
protected:

    /**
     Main mesh that holds the mesh size and frame size that we would use in the standard MADS algorithm or other Mesh-based algorithm.
     */
    std::shared_ptr<MeshBase> _mainMesh;

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
    explicit MadsMegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<Barrier> barrier,
                              std::shared_ptr<MeshBase> mesh,
                              SuccessType success)
      : MegaIteration(parentStep, k, barrier, success),
        _mainMesh(mesh)
    {
        init();
    }
    virtual ~MadsMegaIteration() {}

    NOMAD::ArrayOfPoint suggest() override;

    void observe(const std::vector<NOMAD::EvalPoint>& evalPointList) override;

    /// Implementation of the start tasks for MADS mega iteration.
    /**
     Creates a MadsIteration for each frame center and each desired mesh size.
     Use all xFeas and xInf available.
     For now, not using other frame centers.
     */
    virtual void startImp() override ;

    /// Implementation of the run tasks for MADS mega iteration.
    /**
     Manages the generation of points: either all poll and search points are generated all together before starting evaluation using the MegaSearchPoll or they are generated using a MadsIteration with search and poll separately. A run parameter controls the behavior.
     */
    virtual bool runImp() override;


    const std::shared_ptr<MeshBase> getMesh() const          { return _mainMesh; }
    void setMesh(const std::shared_ptr<MeshBase> &mesh)      { _mainMesh = mesh; }

    void read(  std::istream& is ) override;
    void display(  std::ostream& os ) const override ;

};

/**
 Display useful values so that a new MegaIteration could be constructed using these values.
 */
std::ostream& operator<<(std::ostream& os, const MadsMegaIteration& megaIteration);

/// Get an MegaIteration values from a stream
std::istream& operator>>(std::istream& is, MadsMegaIteration& megaIteration);

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_MADSMEGAITERATION__

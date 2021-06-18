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
#ifndef __NOMAD_4_0_PSDMADSMEGAITERATION__
#define __NOMAD_4_0_PSDMADSMEGAITERATION__

#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the iterations of PSD MADS

class PSDMadsMegaIteration: public MadsMegaIteration
{
private:
    std::shared_ptr<Mads>   _madsOnSubPb; ///< Mads on a subproblem or pollster.
    const Point             _x0; ///< Best point found yet by all subproblems, used as x0
    const Point             _fixedVariable; ///< Fixed variable defining the subproblem

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling -- \b IN.
     \param mesh            Mesh on which other Iteration meshes are based -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit PSDMadsMegaIteration(const Step* parentStep,
                                  const size_t k,
                                  const std::shared_ptr<Barrier>& barrier,
                                  const std::shared_ptr<MeshBase>& mesh,
                                  const SuccessType& success,
                                  const Point& x0,
                                  const Point& fixedVariable)
      : MadsMegaIteration(parentStep, k, barrier, mesh, success),
        _madsOnSubPb(nullptr),
        _x0(x0),
        _fixedVariable(fixedVariable)
    {
    }

    virtual ~PSDMadsMegaIteration()
    {
        destroy();
    }

    ///Get / Set
    const std::shared_ptr<Mads>& getMads() const { return _madsOnSubPb; }
    const Point& getFixedVariable() const { return _fixedVariable; }

    /// Implementation of the start tasks for PSD-MADS mega iteration.
    /**
     Creates a MadsIteration for each frame center and each desired mesh size.
     Use all xFeas and xInf available.
     For now, not using other frame centers.
     */
    virtual void startImp() override ;

    /// Implementation of the run tasks for PSD-MADS mega iteration.
    /**
     Manages the generation of points: either all poll and search points are generated all together before starting evaluation using the MegaSearchPoll or they are generated using a MadsIteration with search and poll separately. A run parameter controls the behavior.
     */
    virtual bool runImp() override;

private:
    void setupSubproblemParams(std::shared_ptr<PbParameters> &subProblemPbParams,
                               std::shared_ptr<RunParameters> &subProblemRunParams,
                               const bool isPollster);

    void destroy();
};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_PSDMADSMEGAITERATION__

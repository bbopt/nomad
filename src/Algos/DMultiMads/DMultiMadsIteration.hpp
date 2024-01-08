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
#ifndef __NOMAD_4_4_DMULTIMADSITERATION__
#define __NOMAD_4_4_DMULTIMADSITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Algos/DMultiMads/DMultiMadsUpdate.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Algos/Mads/Search.hpp"
#include "../../Eval/MeshBase.hpp"


#include "../../nomad_nsbegin.hpp"

/// Class for template algo iterations
/**
Iteration manages the update of the center point (best feasible or best infeasible) around which the trial points are generated.
Trial points generation is performed by random sampling. Evaluation is done. Iteration success is passed to mega iteration.
 */
class DMultiMadsIteration: public Iteration
{
private:
    /// Helper for constructor
    void init ();

    /**
     The current frame center of this iteration. This frame center is either the primary or secondary frame center used by poll and search methods.
     The frame center stores a mesh used by poll and search methods.
     Important note: primary and secondary frame centers will use the single mesh associated to the frame center.
     */
    EvalPointPtr _frameCenter;
    
    std::unique_ptr<NOMAD::DMultiMadsUpdate> _DMultiMadsAlgoUpdate;
    
    std::unique_ptr<Poll> _poll;
    std::unique_ptr<Search> _search;
    
    MeshBasePtr  _mesh;        ///< Mesh of the frame center on which the points are generated
    
    SuccessType _previousSuccess ; ///<  Step success is initialized before calling update. This stores the success of previous iteration.
    
public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param frameCenter        The frame center -- \b IN.
     \param k                  The iteration number -- \b IN.
     */
    explicit DMultiMadsIteration(const Step *parentStep,
                                 const EvalPointPtr frameCenter,
                                 const size_t k,
                                 const MeshBasePtr initialMesh)
      : Iteration(parentStep, k),
        _frameCenter(frameCenter),
        _DMultiMadsAlgoUpdate(nullptr),
        _poll(nullptr),
        _search(nullptr),
        _mesh(initialMesh),
        _previousSuccess(SuccessType::UNDEFINED)
    {
        init();
    }


    // Get/Set
    const EvalPointPtr getFrameCenter() const { return _frameCenter ; }
    void setFrameCenter(EvalPointPtr frameCenter) { _frameCenter = frameCenter ;}
    
    /**
     The DMultiMads algorithm iteration possesses a mesh. This mesh changes during iterations. For the first iteration, a mesh obtained during algo initialization is used. After that, it is the frame center mesh. It is updated during start.
     \remark Used by Step::getIterationMesh() to pass the mesh whenever needed
     */
    const MeshBasePtr getMesh() const override { return _mesh; }
    
    SuccessType getPreviousSuccessType() const { return _previousSuccess ; }
    
protected:

    /// Implementation of run task.
    /**
     Evaluate trial point(s).
     Set the stop reason and also updates the Nelder Mead mega iteration success.
     */
    virtual bool runImp() override ;

    /// Implementation of start task.
    /**
     Update the barrier and create the trial point randomly from the best incumbent.
     */
    virtual void startImp() override;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_DMULTIMADSITERATION__

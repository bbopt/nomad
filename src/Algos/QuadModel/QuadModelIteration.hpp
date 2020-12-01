/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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
#ifndef __NOMAD400_QUAD_MODEL_ITERATION__
#define __NOMAD400_QUAD_MODEL_ITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Algos/MeshBase.hpp"
#include "../../Eval/EvalPoint.hpp"
#include "../../../ext/sgtelib/src/Surrogate.hpp"
#include "../../../ext/sgtelib/src/TrainingSet.hpp"

#include "../../nomad_nsbegin.hpp"

///
class QuadModelIteration: public Iteration
{
private:

    void init();

    /**
     - The center point of the model.
     - Cache points used to build the model are taken around this point.
     */
    const std::shared_ptr<EvalPoint> _frameCenter;

    /**
     The Mads mesh can be available if this is called during a Search method. If not, it is set to \c nullptr. When available, trials points can be projected on it.
     */
    const std::shared_ptr<MeshBase> _madsMesh;

    std::shared_ptr<SGTELIB::TrainingSet>   _trainingSet; ///<
    std::shared_ptr<SGTELIB::Surrogate>     _model;

public:
    /// Constructor
    /**
     \param parentStep       The parent of this step -- \b IN.
     \param frameCenter    The frame center -- \b IN.
     \param k                              The iteration number -- \b IN.
     \param madsMesh            Mads Mesh for trial point projection (can be null) -- \b IN.
     */
    explicit QuadModelIteration(const Step *parentStep,
                                const std::shared_ptr<EvalPoint> &frameCenter,
                                const size_t k,
                                std::shared_ptr<MeshBase> madsMesh)
      : Iteration(parentStep, k) ,
        _frameCenter(frameCenter),
        _madsMesh(madsMesh)
    {
        init();
    }


    /// \brief Destructor
    /// When iteration is done, Flush prints output queue.
    virtual ~QuadModelIteration()
    {
        reset();
    }

    /// Reset the model and the training set.
    void reset();

    /// Access to the quadratic model
    const std::shared_ptr<SGTELIB::Surrogate> getModel() const { return _model;}

    /// Access to the training set
    const std::shared_ptr<SGTELIB::TrainingSet> getTrainingSet() const { return _trainingSet; }

    /// Reimplement to have access to the frame center (can be undefined)
    const std::shared_ptr<EvalPoint> getFrameCenter() const override { return _frameCenter ; }

    /// Reimplement to have access to the mesh (can be null)
    const std::shared_ptr<MeshBase> getMesh() const override { return _madsMesh; }


protected:

    /// Manage  quad model update
    /**
     Use the cache to determine a quad model.
     */
    virtual void startImp() override;

    virtual bool runImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_ITERATION__

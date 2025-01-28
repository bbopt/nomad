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
#ifndef __NOMAD_4_5_QUAD_MODEL_ITERATION__
#define __NOMAD_4_5_QUAD_MODEL_ITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Eval/EvalPoint.hpp"
#include "../../Eval/MeshBase.hpp"
#include "../../../ext/sgtelib/src/Surrogate.hpp"
#include "../../../ext/sgtelib/src/TrainingSet.hpp"

#include "../../nomad_nsbegin.hpp"

///
class QuadModelIteration: public Iteration
{
private:

    void init();


    /**
     - The reference center point
     - Can be null, so trial points are used to define model center.
    */
    const EvalPointPtr _refCenter;


    /**
     The trial points use to create the radius's to select the training set when building the model
     */
    const EvalPointSet & _trialPoints;
    
    /**
     For multi-objective we can combine all objectives into a single objective for the model
     */
    bool _flagPriorCombineObjsForModel;
    
    /**
     The Mads mesh can be available if this is called during a Search method. If not, it is set to \c nullptr. When available, trials points can be projected on it.
     */
    const MeshBasePtr _madsMesh;

    bool _useForSortingTrialPoints;
    
protected:
    std::shared_ptr<SGTELIB::TrainingSet>   _trainingSet; ///<
    std::shared_ptr<SGTELIB::Surrogate>     _model;

    
public:
    /// Constructor
    /**
     \param parentStep      The parent of this step -- \b IN.
     \param center     The frame center (cen be null, so model center is defined with trial points) -- \b IN.
     \param k               The iteration number -- \b IN.
     \param madsMesh        Mads Mesh for trial point projection (can be null) -- \b IN.
     \param trialPoints   Trial points used to define the selection box  (can be empty, so box is defined with mesh)  -- \b IN.
     */
    explicit QuadModelIteration(const Step *parentStep,
                                const EvalPointPtr center,
                                const size_t k = 0,
                                const MeshBasePtr madsMesh = nullptr,
                                const EvalPointSet & trialPoints = emptyEvalPointSet,
                                bool flagPriorCombineObjsForModel = false)
      : Iteration(parentStep, k) ,
        _refCenter(center),
        _madsMesh(madsMesh),
        _useForSortingTrialPoints(false),
        _trialPoints(trialPoints),
        _flagPriorCombineObjsForModel(flagPriorCombineObjsForModel)
    {
        init();
    }


    /// \brief Destructor
    /// When iteration is done, Flush prints output queue.
    virtual ~QuadModelIteration()
    {
        /// Reset the model and the training set.
        if (nullptr != _model)
        {
            _model.reset();
        }
        
        if (nullptr != _trainingSet)
        {
            _trainingSet.reset();
        }
    }


    /// Access to the quadratic model
    const std::shared_ptr<SGTELIB::Surrogate> getModel() const { return _model;}
    
    /// Access to the training set
    const std::shared_ptr<SGTELIB::TrainingSet> getTrainingSet() const { return _trainingSet; }

    /// Access to the frame center (can be undefined)
    const EvalPointPtr getRefCenter() const { return _refCenter ; }

    /// Reimplement to have access to the mesh (can be null)
    const MeshBasePtr getMesh() const override { return _madsMesh; }
    
    /// Reimplement the access to the name. If the class is used for sorting trial points we get name without algo name.
    std::string getName() const override;

protected:

    /// Manage  quad model update
    /**
     Use the cache to determine a quad model.
     */
    virtual void startImp() override;

    virtual bool runImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_QUAD_MODEL_ITERATION__

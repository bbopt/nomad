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
#ifndef __NOMAD_4_0_QUAD_MODEL_UPDATE__
#define __NOMAD_4_0_QUAD_MODEL_UPDATE__

#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

class QuadModelUpdate : public Step
{
private:
    OutputLevel _displayLevel;

    const Point * _frameCenter;
    ArrayOfDouble _radiuses;

public:
    explicit QuadModelUpdate(const Step* parentStep)
      : Step(parentStep),
        _displayLevel(OutputLevel::LEVEL_INFO)
    {
        init();
    }

    virtual ~QuadModelUpdate();

    std::string getName() const override;

private:
    void init();

    /**
     No start task is required
     */
    virtual void startImp() override {}

    /// Implementation of the run task.
    /**
     Update the SGTELIB::TrainingSet and SGTELIB::Surrogate contained in the QuadModelIteration ancestor:
     - Get relevant points in cache around current frame center.
     - Add points to training set and build new model.
     - Assess if model is ready. Update its bounds.
     \return \c true if model is ready \c false otherwise.
     */
    virtual bool runImp() override;

    /**
     No end task is required
     */
    virtual void    endImp() override {}

    bool isValidForUpdate(const EvalPoint& evalPoint) const; ///< Helper function for cache find.

    bool isValidForIncludeInModel(const EvalPoint& evalPoint) const; ///< Helper function for cache find.

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_QUAD_MODEL_UPDATE__

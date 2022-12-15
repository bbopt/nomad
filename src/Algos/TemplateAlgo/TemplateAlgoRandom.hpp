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
#ifndef __NOMAD_4_3_TEMPLATEALGORANDOM__
#define __NOMAD_4_3_TEMPLATEALGORANDOM__

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Step.hpp"
#include "../../Eval/EvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to perform random generation of trial points.
/**
  
 */
class TemplateAlgoRandom: public Step, public IterationUtils
{
private:
    
    NOMAD::EvalPointPtr _center; // The point around which we sample randomly
    
    NOMAD::ArrayOfDouble _boxSize; // The sample is made within boxes of size {1..k}*Delta
    
    
public:
    /// Constructor
    /**
     \param parentStep The parent of this step
     */
    explicit TemplateAlgoRandom(const Step* parentStep)
      : Step(parentStep),
        IterationUtils(parentStep),
        _center(nullptr)
    {
        init();
    }
    virtual ~TemplateAlgoRandom() {}

    /// Generate new points to evaluate
    /**
     A new point is obtained using the simplex. xt = yc + delta*d. Delta is the multiplicative factor that can increase or decrease with success. \n

     The point is snapped to bounds and projected on the mesh.
     */
    void generateTrialPointsImp() override;


private:

    /**
     The delta parameter used to create the trial point is different. The possible delta parameters are obtained from _runParams. The validity of the parameters are checked. \n
     The flag to perform a standalone template random algo optimization is also set.
     */
    void init();

    /// Implementation of the start task.
    /**
     Call TemplateAlgoRandom::generateTrialPoints and update the trial points.
     */
    virtual void    startImp() override ;

    /// Implementation of the run task.
    /**
     Evaluate the trial point and store it locally. Call IterationUtils::postProcessing.

     \return \c true if a better point is obtained \c false otherwise.
     */
    virtual bool    runImp() override ;

    /// No end task is required
    virtual void    endImp() override {}

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_3_TEMPLATEALGORANDOM__

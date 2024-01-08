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
#ifndef __NOMAD_4_4_TEMPLATEALGOSINGLEPASS__
#define __NOMAD_4_4_TEMPLATEALGOSINGLEPASS__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoIteration.hpp"
#include "../../Algos/IterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/**
 Class to generate points for single pass of random generation (no iteration).
 The points are projected on mesh.
 */
class TemplateAlgoSinglePass: public TemplateAlgoIteration, public IterationUtils
{
public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param center             The frame center around which is done the random sampling  -- \b IN.
     */
    explicit TemplateAlgoSinglePass(const Step* parentStep,
                                    const EvalPointPtr center)
      : TemplateAlgoIteration(parentStep, center, 0 ),
        IterationUtils(parentStep)
    {
        init();
    }
    
    // No special destructor needed - keep defaults.


private:

    /// Implementation of start tasks.
    /**
     - call the default Iteration::startImp
     - create the initial simplex if it is empty.
     - generate trial points using
     - verify that trial points are on mesh.
     */
    void startImp() override ;

    /// Implementation of run task. Nothing to do.
    bool runImp() override { return  false;}

    /// Implementation of run task. Nothing to do.
    void endImp() override {}

    /// Generation of trial points uses the random generator
    void generateTrialPointsImp() override;
    
    void init();

protected:
    // Override default Step::getName that seeks the algo name. Here we have a single Iteration and the class does not derive from an Algo.
    std::string getName() const override;
    
};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_TEMPLATEALGOSINGLEPASS__

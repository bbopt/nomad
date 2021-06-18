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
#ifndef __NOMAD_4_0_NMSHRINK__
#define __NOMAD_4_0_NMSHRINK__

#include "../../Algos/NelderMead/NMIterationUtils.hpp"
#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for shrink step for NM algorithm.
/**
 The SHRINK step is executed when NM is called as a standalone algorithm and when all reflective steps (REFLECT, EXPAND, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION) fail to obtain improvement. \n
 The class manages the creation and evaluation of the shrunk simplex (start, run, end).
 In the start function, before shrinking the simplex, an update of the main barrier must be performed by calling NMUpdate start, run and end. The shrunk simplex is obtained in start task. \n
 Evaluation of the new simplex is performed in run.
 */
class NMShrink: public Step , public NMIterationUtils
{
private:
    Double _gamma;
public:
    /// Constructor
    /**
     \param parentStep The parent of this NM step
     */
    explicit NMShrink(const Step* parentStep )
      : Step( parentStep ) ,
        NMIterationUtils ( parentStep )
    {
        init();
    }
    virtual ~NMShrink() {}

    /// Implementation of the start tasks for simplex shrink.
    /**
     - update the barrier
     - call NMShrink::generateTrialPoints
     */
    virtual void    startImp() override ;

    /// Implementation of the run task for simplex shrink.
    /**
     Evaluate the trial points.
     */
    virtual bool    runImp() override ;

    /// Implementation of the end task for simplex shrink.
     /**
      Call default IterationUtils::postProcessing.
      */
    virtual void    endImp() override ;

    /// Generate new points to evaluate
    /**
     The new shrunk simplex is obtained with the formula y[k] = y0[k] + _gamma*(yi[k]-y0[k]). Where y0 is the frame center and yi are elements of the previous simplex. Gamma must be a parameter in ]0;1]
     */
    void generateTrialPoints() override;


private:

    /// Helper for constructor
    void init();

    /*---------------------------------*/
    /* Private methods used by Shrink */
    /*---------------------------------*/


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_NMSHRINK__

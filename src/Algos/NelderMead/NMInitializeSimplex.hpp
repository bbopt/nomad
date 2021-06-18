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
#ifndef __NOMAD_4_0_NMINITIALIZESIMPLEX__
#define __NOMAD_4_0_NMINITIALIZESIMPLEX__

#include "../../Algos/NelderMead/NMIterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for simplex initialization of NM algorithm.
/**
 The simplex is obtained by calling run (done in NMIteration).
 */
class NMInitializeSimplex: public Step, public NMIterationUtils
{
public:
    /// Constructor
    /**
     \param parentStep The parent of this NM step
     */
    explicit NMInitializeSimplex(const Step* parentStep)
      : Step(parentStep) ,
        NMIterationUtils(parentStep)
    {
        init();
    }
    virtual ~NMInitializeSimplex() {}

    /// No new points are generated
    void generateTrialPoints() override {}

private:
    /// Helper for constructor
    void init();

    /// No start task is required.
    virtual void    startImp() override {}

    /// Implementation of the run task for simplex initialization.
    /**
     Calls createSimplex if required.
     */
    virtual bool    runImp() override ;

    /// No start task is required
    virtual void    endImp() override {}

    /// Helper for run
    /**
     From the Cache or from the Barrier, a set of points within a radius of the current frame center is considered before creating the initial simplex. The radius depends on a two given parameters and the frame size.
     The initial simplex is obtained by adding points to obtain dim=n+1 and simplex being affinely independant. We use the rank of DZ=[(y1-y0 (y2-y0) ...(yk-y0)] to decide if a point yk increases the rank or not.
     The characteristics of the initial simplex (volumes and diameter) are updated.
     */
    bool createSimplex();


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_NMINITIALIZESIMPLEX__

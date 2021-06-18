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
#ifndef __NOMAD_4_0_TERMINATION__
#define __NOMAD_4_0_TERMINATION__

#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

///  Class for termination of an algorithm.
/**
 The terminate function checks for termination criterions such as MAX_ITERATIONS, MAX_TIME, STOP_IF_FEASIBLE and set the stop reason.
 */
class Termination: public Step
{
public:
    /// Constructor
    explicit Termination(const Step* parentStep,
                         const std::shared_ptr<RunParameters>& runParams = nullptr,
                         const std::shared_ptr<PbParameters>& pbParams = nullptr)
      : Step(parentStep, runParams, pbParams)
    {
        init();
    }

    /// Destructor
    virtual ~Termination() {}

    /**
     The terminate function is called when algorithm are performing iterations during a run. At each iteration, we test if a stop criterion is reached.
     */
    virtual bool terminate(size_t iteration);

    virtual void    startImp() override; ///< Will update the step name

    /// Implementation for run task of algorithm Termination.
    /**
     \return \c true is a stop reason requires termination of an algorithm, \c false otherwise.
     */
    virtual bool    runImp()   override;

    /// Implementation for end tasks of algorithm Termination.
    /**
     Upon completing an algorithm run, this end function is called to display termination info.
     */
    virtual void    endImp()   override;

private:

    /// Helper for constructor
    void init();
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_TERMINATION__

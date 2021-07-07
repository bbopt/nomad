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
#ifndef __NOMAD_4_0_NM__
#define __NOMAD_4_0_NM__


#include "../../Algos/Algorithm.hpp"
#include "../../Algos/AlgoStopReasons.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class implementing Nelder Mead non-linear simplex algorithm for constrained problems.
/**
 See the NM-Mads paper: https://link.springer.com/article/10.1007/s10589-018-0016-0 for details.
 */
class NM: public Algorithm
{
public:
    /// Constructor
    /**
     \param parentStep          The parent of this Step -- \b IN.
     \param stopReasons         The stop reasons for NM -- \b IN.
     \param runParams           The run parameters that control NM -- \b IN.
     \param pbParams            The problem parameters that control NM -- \b IN.
     */
    explicit NM(const Step* parentStep,
                std::shared_ptr<AlgoStopReasons<NMStopType>> stopReasons,
                const std::shared_ptr<RunParameters>& runParams,
                const std::shared_ptr<PbParameters>& pbParams)
      : Algorithm(parentStep, stopReasons, runParams, pbParams)
    {
        init();
    }

    /// Destructor
    virtual ~NM() {}

    virtual void readInformationForHotRestart() override ;

private:
    /// Helper for constructor
    void init();

    /// Implementation for run tasks.
    /**
     - Algorithm execution for single-objective.
     - Loop on NMMegaIteration (start, run, end) until a stop reason to terminate is obtained.
     - Update the succes type
     - Perform Termination tasks (start, run, end)
     - Update the SearchMethod success type with best success found.
     \return \c true
     */
    virtual bool runImp() override;

    /// Implementation for start tasks.
    /**
     - Set the stop reason to STARTED
     - Reset sub-algorithm counter
     - Perform Initialization tasks (start, run, end)
     */
    virtual void startImp() override;

    virtual void endImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_NM__

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
#ifndef __NOMAD_4_0_SEARCH__
#define __NOMAD_4_0_SEARCH__

#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to manage the SearchMethods used by MADS algorithm during its search step.
class Search final : public Step , public IterationUtils
{
private:
    std::vector<std::shared_ptr<SearchMethodBase>> _searchMethods;
#ifdef TIME_STATS
    DLL_ALGO_API static std::vector<double> _searchTime;        ///< Total time spent running each search
    DLL_ALGO_API static std::vector<double> _searchEvalTime;    ///< Total time spent evaluating search points
#endif // TIME_STATS

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit Search(const Step* parentStep )
      : Step( parentStep ),
        IterationUtils( parentStep ),
        _searchMethods()
    {
        init();
    }

    virtual ~Search() {}

    /**
     - Generate new points to evaluate. Use all enabled search methods.
     - To be used only when parameter MEGA_SEARCH_POLL is true.
     */
    void generateTrialPoints() override;

#ifdef TIME_STATS
    /// Time stats
    static std::vector<double> getSearchTime()       { return _searchTime; }
    static std::vector<double> getSearchEvalTime()   { return _searchEvalTime; }
#endif // TIME_STATS

private:

    void init();

    /// Implementation of the start task.
    /**
     Just perform a sanity check on MEGA_SEARCH_POLL that must be false.
     */
    virtual void startImp() override;

    /// The implementation of run tasks.
    /**
      Perform start+run+end for all search methods in the vector _searchMethods.
     */
    virtual bool runImp() override;

    /// Implementation of the end tasks
    /**
      If a sub optimization is used during search we probably set a stop reason to terminate. The parent optimization must go on. The stop reason is set to started if sub optimization reached its evaluation budget.
     */
    virtual void endImp() override;

    /**
     Identify if there is at least one search enabled. If there are none, do not print Search step at all.
     */
    bool isEnabled() const;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_SEARCH__


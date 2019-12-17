/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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
#ifndef __NOMAD400_SEARCH__
#define __NOMAD400_SEARCH__

#include "../../Algos/Mads/SearchMethod.hpp"

#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to manage the SearchMethods used by MADS algorithm during its search step.
class Search final : public Step , public MadsIterationUtils
{
private:
    std::vector<std::shared_ptr<SearchMethod>> _searchMethods;

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit Search(const Step* parentStep )
      : Step( parentStep ),
        MadsIterationUtils( parentStep ),
        _searchMethods()
    {
        init();
    }

    virtual ~Search() {}

    /**
     Generate new points to evaluate. Use all enabled search methods.
     */
    void generateTrialPoints() override ;
    
private:

    void init();

    /// Implementation of the start task.
    /**
     Sanity check on GENERATE_ALL_POINTS_BEFORE_EVAL that must be false.
     */
    virtual void startImp() override;
    
    /// The implementation of run tasks.
    /**
      Perform start+run+end for all search methods in the vector _searchMethods.
     */
    virtual bool runImp() override ;
    
    /// Implementation of the end tasks
    /**
      If a sub optimization is used during search we probably set a stop reason to terminate. The parent optimization must go on. The stop reason is set to started if sub optimization reached its evaluation budget.
     */
    virtual void endImp() override ;
    
    /**
     Identify if there is at least one search enabled. If there are none, do not print Search step at all.
     */
    bool isEnabled() const;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SEARCH__


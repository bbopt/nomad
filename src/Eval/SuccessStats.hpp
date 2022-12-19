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
#ifndef __NOMAD_4_3_STEPSUCCESSSTATS__
#define __NOMAD_4_3_STEPSUCCESSSTATS__

#include <map>

#include "../Type/StepType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"

#include "../nomad_platform.hpp"
#include "../nomad_nsbegin.hpp"

/// Class to manage success type stats of step
/**
  Stats for success type of steps.
  Some stats are propagated from step to parent step (nbSuccessAndFail) but the number of consecutive successes and fails of a step is not propagated.
 */
class DLL_EVAL_API SuccessStats
{
protected:
    
    // This success type is propagated from step to parent step.
    std::map<std::pair<StepType,SuccessType>, size_t> _nbSuccessAndFail;
    
    // Success: PARTIAL_SUCCESS and FULL_SUCCESS
    // Fail: NO_TRIALS, UNSUCCESSFUL
    // UNDEFINED success type is not recorded
    // This is not propagated from step to parent step. No need to have a map with StepType
    size_t _nbConsecutiveSuccess;
    size_t _nbConsecutiveFail;
    
public:
    /// Constructor
    explicit SuccessStats():
    _nbConsecutiveSuccess(0),
    _nbConsecutiveFail(0)
    {}
       
    /// Update the stats for a given success type and step.
    void updateStats(SuccessType successType, StepType stepType , size_t val=1);
    
    /// Update the stats with a given success stats
    void updateStats(const SuccessStats & evalStats);
    
    /// Access to consecutive successes and fails
    size_t getStatsNbConsecutiveSuccess() const { return _nbConsecutiveSuccess ; }
    size_t getStatsNbConsecutiveFail() const { return _nbConsecutiveFail ; }
    
    
    // size_t getNbSuccessAndFail(SuccessType successType, StepType stepType) const;
     
    // Reset stats that are passed to parents
    void resetCurrentStats() { _nbSuccessAndFail.clear(); } ///< Reset the current stats for all registered step types. Resets are performed at start of step.
    
    bool hasStatsForPropagation() const { return !_nbSuccessAndFail.empty();} 
    
    std::string display() const;
    
private:
    
    const std::map<std::pair<StepType,SuccessType>, size_t> & getStatsMapSuccessAndFail() const { return _nbSuccessAndFail ;}
    
    void setNbConsecutiveSuccessAndFail(SuccessType successType, size_t val);
    
    /// Update the stats for a given success type and step.
    void updateSuccessAndFailStats(SuccessType successType, StepType stepType , size_t val);
    
    
};
///   Display method stats.
std::ostream& operator<<(std::ostream& os, const SuccessStats& stats);


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_3_BASEEVALSTATS__

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
/**
\file   SuccessStats.cpp
\brief  Manage step success stats.
\author Christophe Tribes
\date   June 2022
*/

#include "../Eval/SuccessStats.hpp"

// This function is used for incrementing the SuccessStats of a step
void NOMAD::SuccessStats::updateStats(SuccessType successType, StepType stepType ,size_t val )
{
    
    // UNDEFINED is not used for stats
    if ( NOMAD::SuccessType::UNDEFINED == successType )
    {
        return;
    }
 
    setNbConsecutiveSuccessAndFail(successType,val);
    
    updateSuccessAndFailStats(successType, stepType, val);
}

void NOMAD::SuccessStats::updateSuccessAndFailStats(SuccessType successType, StepType stepType ,size_t val )
{
    
    auto p = std::make_pair(stepType,successType);
    auto iter = _nbSuccessAndFail.find(p);
    
    // First time insert
    if (iter == _nbSuccessAndFail.end())
    {
        _nbSuccessAndFail.insert(std::make_pair(p, val));
    }
    else
    {
        _nbSuccessAndFail.at(p) += val;
    }
}


// This function is used for propagation.
void NOMAD::SuccessStats::updateStats(const SuccessStats & evalStats )
{
    // We may have more stats
    auto statsMap = evalStats.getStatsMapSuccessAndFail();
    
    for (auto it=statsMap.begin(); it!=statsMap.end(); ++it )
    {
        auto p = it->first;
        auto val = it->second;
        updateSuccessAndFailStats(p.second, p.first,val);
    }
    
}

void NOMAD::SuccessStats::setNbConsecutiveSuccessAndFail(SuccessType successType, size_t val)
{
  
    
    if (successType >= NOMAD::SuccessType::PARTIAL_SUCCESS)
    {
        // Increment success count for this step and reset fail count
        _nbConsecutiveSuccess +=val;
        _nbConsecutiveFail = 0;
    }
    else
    {
        // Increment fail count for this step and reset success count
        _nbConsecutiveFail += val;
        _nbConsecutiveSuccess = 0;
    }
}

//size_t NOMAD::SuccessStats::getNbSuccessAndFail(SuccessType successType, StepType stepType) const
//{
//    auto iter = _nbSuccessAndFail.find(std::make_pair(stepType,successType));
//    
//    if (iter == _nbSuccessAndFail.end())
//    {
//        return 0;
//    }
//    else
//    {
//        return iter->second;
//    }
//}


/*---------------*/
/*    display    */
/*---------------*/
std::string NOMAD::SuccessStats::display() const
{
    std::string s;
    throw NOMAD::Exception(__FILE__, __LINE__,"Not yet implemented ");
    return s;
}

std::ostream& operator<<(std::ostream& os, NOMAD::SuccessStats& stats)
{
    
    std::ostringstream oss;
    oss << stats.display();
    return os;
}

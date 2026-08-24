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

#include "../../Algos/CatMads/CatCustomOrder.hpp"
#include "../../Output/OutputQueue.hpp"

NOMAD::CatCustomOrder::CatCustomOrder(CustomOrderCompWrapper compWrapper, NOMAD::BBOutputTypeList bbot)
    : NOMAD::OrderByEval({NOMAD::EvalType::CAT_MODEL, NOMAD::defaultFHComputeTypeS}),
      _compWrapper(compWrapper), _bbot(bbot)
{
    setName("CatCustomOrder");
}


// Comparison function
bool NOMAD::CatCustomOrder::comp(NOMAD::EvalQueuePointPtr& p1, NOMAD::EvalQueuePointPtr& p2) const {

    // From function comp() in ComparePriority.cpp:
    // sorting from less interesting to most interesting point
    // return true if point1 is less interesting than point2

    // For now, use mixed model for ordering the points in the union of the Cat Poll (categorical var) and regular ortho Poll (quantitative var)

    auto eval1 = p1->getEval(NOMAD::EvalType::CAT_MODEL);
    auto eval2 = p2->getEval(NOMAD::EvalType::CAT_MODEL);
    if (nullptr == eval2)
    {
        return false;
    }
    else if (nullptr == eval1)
    {
        return true;
    }

    return NOMAD::OrderByEval::comp(p1, p2); 
    
}



// Complete trial points information
void NOMAD::CatCustomOrder::completeTrialPointsInformation(const NOMAD::Step *step, NOMAD::EvalPointSet & trialPoints)
{

    if (trialPoints.empty())
    {
        OUTPUT_INFO_START
        std::string s = "Empty trial points no need to complete information";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
        OUTPUT_INFO_END
        return;
    }
        
        
    // Convert trial points to send to CallbackCatMadsComp function
    std::vector<std::vector<double>> convTrialPoints;
    for (const auto &trialPoint: trialPoints)
    {
        std::vector<double> tp;
        for (size_t i=0; i < trialPoint.size() ; i++)
        {
            tp.push_back(trialPoint[i].todouble());
        }
        convTrialPoints.push_back(tp);
    }

    
    std::vector<std::vector<double>> trialPointsModelEval = _compWrapper(convTrialPoints);
    
    if (trialPointsModelEval.empty())
    {
        OUTPUT_INFO_START
        std::string s = "No model available to complete information on trial points.";
        NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
        OUTPUT_INFO_END
        return;
    }

    // Lambda function to convert a vector<double>  to a string
    auto vectorToString = [](const std::vector<double>& vec) -> std::string
    {
        std::ostringstream oss; // Create a string stream
        for (size_t i = 0; i < vec.size(); ++i) {
            oss << vec[i]; // Add the value
            if (i != vec.size() - 1) oss << " "; // Add space between values
        }
        return oss.str(); // Return the string
    };
    
    // Insert CatMads model values in the copies of eval points
    NOMAD::EvalPointSet copyTrialPoints;
    int i = 0;
    for (const auto & evalPoint : trialPoints)
    {
        NOMAD::EvalPoint copyEvalPoint = evalPoint;
        
        copyEvalPoint.setBBO(vectorToString(trialPointsModelEval[i]), _bbot, NOMAD::EvalType::CAT_MODEL, true);
        copyEvalPoint.setEvalStatus(NOMAD::EvalStatusType::EVAL_OK, NOMAD::EvalType::CAT_MODEL);
        copyTrialPoints.insert(copyEvalPoint);
        i++;
    }
    
    // Replace original with the copy
    trialPoints.clear();
    trialPoints = copyTrialPoints;
    
    return;
    
}


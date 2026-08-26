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


#ifndef __NOMAD_4_6_CATALGOUTILS__
#define __NOMAD_4_6_CATALGOUTILS__

#include "../../Algos/Mads/Mads.hpp"

#include "../../nomad_nsbegin.hpp"

/// The Cat-Algo utils class to manage callbacks.
class DLL_ALGO_API CatAlgoUtils
{
    
protected:
    
    ListOfVariableGroup _catVG;
    
    // Callback functions
    std::function<bool()> _buildCatModel = [](){return false;};
    std::function<std::vector<Direction>(EvalPointPtr)> _catPollFree= [](EvalPointPtr ep){return std::vector<Direction>();};
    std::function<std::vector<EvalPoint>(EvalPointPtr, EvalPointPtr)> _catCategoricalSearch= [](EvalPointPtr feasEp, EvalPointPtr infeasEp){return std::vector<EvalPoint>();};
    
public:
    /// Constructor

    explicit CatAlgoUtils()
    {
    }
    
    
    bool setBuildModelCallback(std::function<bool()> buildCatModel);
    bool buildCatModel() const ;

    // For categorical free poll
    bool setPollFreeCallback(std::function<std::vector<Direction>(EvalPointPtr)> catPollFree);

    std::vector<Direction> runCatPollFree(EvalPointPtr frameCenter) const;

    // For categorical search
    bool setCategoricalSearchCallback(std::function<std::vector<EvalPoint>(EvalPointPtr, EvalPointPtr)> catCategoricalSearch);
    
    std::vector<EvalPoint> runCatCategoricalSearch(EvalPointPtr feasFrameCenter, EvalPointPtr infeasFrameCenter) const ;
    
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_6_CATALGOUTILS__

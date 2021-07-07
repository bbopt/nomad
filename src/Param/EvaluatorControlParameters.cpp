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

#include "../Param/EvaluatorControlParameters.hpp"
#include "../Type/EvalSortType.hpp"


/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::EvaluatorControlParameters::init()
{
    _typeName = "EvaluatorControl";

    try
    {
        #include "../Attribute/evaluatorControlAttributesDefinition.hpp"
        registerAttributes(_definition);

        // Note: we cannot call checkAndComply() here, the default values
        // are not valid, for instance DIMENSION, X0, etc.

    }
    catch (NOMAD::Exception& e)
    {
        std::string errorMsg = "Attribute registration failed: ";
        errorMsg += e.what();
        throw NOMAD::Exception(__FILE__,__LINE__, errorMsg);
    }

}

/*----------------------------------------*/
/*            check the parameters        */
/*----------------------------------------*/
void NOMAD::EvaluatorControlParameters::checkAndComply(
                        const std::shared_ptr<NOMAD::EvaluatorControlGlobalParameters>& evaluatorControlGlobalParams,
                        const std::shared_ptr<NOMAD::RunParameters>& runParams)
{
    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    // When runParameters are provided, update internal parameter SUBPROBLEM_MAX_BB_EVAL.
    if (nullptr != runParams)
    {
        auto psdMadsOpt = runParams->getAttributeValue<bool>("PSD_MADS_OPTIMIZATION");
        auto ssdMadsOpt = runParams->getAttributeValue<bool>("SSD_MADS_OPTIMIZATION");
        if (psdMadsOpt)
        {
            setAttributeValue("SUBPROBLEM_MAX_BB_EVAL", getAttributeValueProtected<size_t>("PSD_MADS_SUBPROBLEM_MAX_BB_EVAL", false));
        }
        else if (ssdMadsOpt)
        {
            setAttributeValue("SUBPROBLEM_MAX_BB_EVAL", getAttributeValueProtected<size_t>("SSD_MADS_SUBPROBLEM_MAX_BB_EVAL", false));
        }
        else
        {
            setAttributeValue("SUBPROBLEM_MAX_BB_EVAL", NOMAD::INF_SIZE_T);
        }
    }

    if (nullptr != evaluatorControlGlobalParams)
    {
        if (evaluatorControlGlobalParams->toBeChecked())
        {
            evaluatorControlGlobalParams->checkAndComply();
        }
        auto maxSurrogateEval = evaluatorControlGlobalParams->getAttributeValue<size_t>("MAX_SURROGATE_EVAL_OPTIMIZATION");
        bool isSurrogateOptimization = getAttributeValueProtected<bool>("EVAL_SURROGATE_OPTIMIZATION", false);
        if (isSurrogateOptimization)
        {
            // If this is a surrogate optimization, either it has to
            // have a cost relative to bb evaluations (for MAX_EVAL to have effect - it could be set if all variables are granular),
            // or it has to have a maximum number of surrogate evaluations.
            auto surrogateCost = evaluatorControlGlobalParams->getAttributeValue<size_t>("EVAL_SURROGATE_COST");
            if (NOMAD::INF_SIZE_T == surrogateCost && NOMAD::INF_SIZE_T == maxSurrogateEval)
            {
                throw NOMAD::Exception(__FILE__, __LINE__,
                    "Parameter MAX_SURROGATE_EVAL_OPTIMIZATION or EVAL_SURROGATE_COST must be non-infinite when EVAL_SURROGATE_OPTIMIZATION is used.");
            }
            if (evaluatorControlGlobalParams->getAttributeValue<size_t>("MAX_BB_EVAL") < NOMAD::INF_SIZE_T)
            {
                throw NOMAD::Exception(__FILE__, __LINE__,
                    "Parameter MAX_BB_EVAL should not be set when EVAL_SURROGATE_OPTIMIZATION is used. Use MAX_SURROGATE_EVAL_OPTIMIZATION instead.");
            }
            if (NOMAD::EvalSortType::SURROGATE == evaluatorControlGlobalParams->getAttributeValue<NOMAD::EvalSortType>("EVAL_QUEUE_SORT"))
            {
                throw NOMAD::InvalidParameter(__FILE__, __LINE__, "Parameter EVAL_QUEUE_SORT cannot be SURROGATE when EVAL_SURROGATE_OPTIMIZATION is set");
            }
        }
        else
        {
            if (maxSurrogateEval < NOMAD::INF_SIZE_T)
            {
                throw NOMAD::InvalidParameter(__FILE__,__LINE__, "Parameter MAX_SURROGATE_EVAL_OPTIMIZATION should be set only when EVAL_SURROGATE_OPTIMIZATION is used.");
            }
        }

    }

    _toBeChecked = false;

}
// End checkAndComply()





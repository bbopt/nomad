/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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

#include "../Param/EvaluatorControlGlobalParameters.hpp"


/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::EvaluatorControlGlobalParameters::init()
{
    _typeName = "EvaluatorControl";

    try
    {
        #include "../Attribute/evaluatorControlGlobalAttributesDefinition.hpp"
        registerAttributes(_definition);

        // Note: we cannot call checkAndComply() here, the default values
        // are not valid, for instance DIMENSION, X0, etc.

    }
    catch (NOMAD::Exception & e)
    {
        std::string errorMsg = "Attribute registration failed: ";
        errorMsg += e.what();
        throw NOMAD::Exception(__FILE__,__LINE__, errorMsg);
    }

}

/*----------------------------------------*/
/*            check the parameters        */
/*----------------------------------------*/
void NOMAD::EvaluatorControlGlobalParameters::checkAndComply(const std::shared_ptr<PbParameters> pbParams )
{
    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    if (getAttributeValueProtected<size_t>("MAX_BB_EVAL",false) <= 0)
    {
        setAttributeValue("MAX_BB_EVAL", NOMAD::INF_SIZE_T);
    }

    if (getAttributeValueProtected<size_t>("MAX_EVAL",false) <= 0)
    {
        setAttributeValue("MAX_EVAL", NOMAD::INF_SIZE_T);
    }

    // Safety net using MAX_EVAL when variables are all granular -> setting a value for MAX_EVAL prevent Nomad from circling around the solution forever when all neighbors have been visited.
    if (pbParams != nullptr && getAttributeValueProtected<size_t>("MAX_EVAL",false) == NOMAD::INF_SIZE_T )
    {
        auto granularity = pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("GRANULARITY");
        bool allGranular = true;

        size_t dimension = granularity.size();
        for (size_t i = 0; i < dimension ; i++)
        {
            if (0.0 == granularity[i])
            {
                allGranular = false;
                break;
            }
        }

        if (allGranular)
        {
            // First guess, set the limit to 100*3^n evaluations, with a maximum of 1000000.
            size_t new_max_eval = 1000000; // default value.
            if (dimension < 9)
            {
                new_max_eval = 100 * (size_t)pow(3, dimension);
            }

            // MAX_EVAL must be >> MAX_BB_EVAL.
            // If MAX_BB_EVAL is set, set MAX_EVAL to at least 10 times MAX_BB_EVAL. Don't exceed INF_SIZE_T.

            auto max_bb_eval = getAttributeValueProtected<size_t>("MAX_BB_EVAL",false);

            if (max_bb_eval < NOMAD::INF_SIZE_T)
            {
                if (max_bb_eval < NOMAD::INF_SIZE_T / 10)
                {
                    if (new_max_eval < 10 * max_bb_eval)
                    {
                        new_max_eval = 10 * max_bb_eval;
                    }
                }
                else
                {
                    new_max_eval = NOMAD::INF_SIZE_T;
                }
            }

            std::cerr << "All variables are granular. MAX_EVAL is set to " << new_max_eval << " to prevent algorithm from circling around best solution indefinetely" << std::endl;

            setAttributeValue("MAX_EVAL",new_max_eval);
        }
    }

    auto sgteBlockSize = getAttributeValueProtected<size_t>("SGTE_MAX_BLOCK_SIZE", false);
    if (0 == sgteBlockSize)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter SGTE_MAX_BLOCK_SIZE must be positive");
    }
    else if (sgteBlockSize > getAttributeValueProtected<size_t>("MAX_SGTE_EVAL", false))
    {
        setAttributeValue("SGTE_MAX_BLOCK_SIZE", getAttributeValueProtected<size_t>("MAX_SGTE_EVAL", false));
    }

    _toBeChecked = false;

}
// End checkAndComply()





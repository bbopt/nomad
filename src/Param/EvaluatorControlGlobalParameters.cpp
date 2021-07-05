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

#include "../Param/EvaluatorControlGlobalParameters.hpp"
#ifdef WINDOWS
#include <tchar.h>
#include <windows.h>
#endif


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
void NOMAD::EvaluatorControlGlobalParameters::checkAndComply(const std::shared_ptr<PbParameters> pbParams)
{
    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    if (isAttributeDefaultValue<std::string>("TMP_DIR"))
    {
#ifdef WINDOWS
	TCHAR tempPathBuffer[MAX_PATH];
	GetTempPath(MAX_PATH, tempPathBuffer);
        setAttributeValue("TMP_DIR", std::string(tempPathBuffer));
#else
        setAttributeValue("TMP_DIR", std::string("/tmp/"));
#endif
    }

    if (0 == getAttributeValueProtected<size_t>("MAX_BB_EVAL",false))
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter MAX_BB_EVAL must be positive");
    }

    if (0 == getAttributeValueProtected<size_t>("MAX_EVAL",false))
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter MAX_EVAL must be positive");
    }

    if (0 == getAttributeValueProtected<size_t>("MAX_SURROGATE_EVAL_OPTIMIZATION",false))
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter MAX_SURROGATE_EVAL_OPTIMIZATION must be positive");
    }

    // Safety net using MAX_EVAL when variables are all granular -> setting a value for MAX_EVAL prevent Nomad from circling around the solution forever when all neighbors have been visited.
    if (pbParams != nullptr && getAttributeValueProtected<size_t>("MAX_EVAL",false) == NOMAD::INF_SIZE_T)
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

            std::cerr << "All variables are granular. MAX_EVAL is set to " << new_max_eval << " to prevent algorithm from circling around best solution indefinitely" << std::endl;

            setAttributeValue("MAX_EVAL",new_max_eval);
        }
    }

    // Set value of MODEL_MAX_EVAL (internal) to QUAD_MODEL_MAX_EVAL or SGTELIB_MODEL_MAX_EVAL.
    auto quadModelMaxEval = getAttributeValueProtected<size_t>("QUAD_MODEL_MAX_EVAL", false);
    auto sgtelibModelMaxEval = getAttributeValueProtected<size_t>("SGTELIB_MODEL_MAX_EVAL", false);
    bool quadModelChanged = isSetByUser("QUAD_MODEL_MAX_EVAL") || !isAttributeDefaultValue<size_t>("QUAD_MODEL_MAX_EVAL");
    bool sgtelibModelChanged = isSetByUser("SGTELIB_MODEL_MAX_EVAL") || !isAttributeDefaultValue<size_t>("SGTELIB_MODEL_MAX_EVAL");
    auto modelMaxEval = quadModelMaxEval;
    if (quadModelChanged && sgtelibModelChanged && quadModelMaxEval != sgtelibModelMaxEval)
    {
        // See issue #505
        std::cerr << "Warning: Currently not supported: QUAD_MODEL_MAX_EVAL (";
        std::cerr << quadModelMaxEval << ") different than SGTELIB_MODEL_MAX_EVAL (";
        std::cerr << sgtelibModelMaxEval << "). Using only the value of QUAD_MODEL_MAX_EVAL." << std::endl;
        setAttributeValue("MODEL_MAX_EVAL", quadModelMaxEval);
    }
    else if (quadModelChanged)
    {
        setAttributeValue("MODEL_MAX_EVAL", quadModelMaxEval);
    }
    else
    {
        setAttributeValue("MODEL_MAX_EVAL", sgtelibModelMaxEval);
        modelMaxEval = sgtelibModelMaxEval;
    }

    // Set value of MODEL_MAX_BLOCK_SIZE (internal) to QUAD_MODEL_MAX_BLOCK_SIZE or SGTELIB_MODEL_MAX_BLOCK_SIZE.
    auto quadModelBlockSize = getAttributeValueProtected<size_t>("QUAD_MODEL_MAX_BLOCK_SIZE", false);
    auto sgtelibModelBlockSize = getAttributeValueProtected<size_t>("SGTELIB_MODEL_MAX_BLOCK_SIZE", false);
    quadModelChanged = isSetByUser("QUAD_MODEL_MAX_BLOCK_SIZE") || !isAttributeDefaultValue<size_t>("QUAD_MODEL_MAX_BLOCK_SIZE");
    sgtelibModelChanged = isSetByUser("SGTELIB_MODEL_MAX_BLOCK_SIZE") || !isAttributeDefaultValue<size_t>("SGTELIB_MODEL_MAX_BLOCK_SIZE");

    if (0 == quadModelBlockSize)
    {
        throw NOMAD::InvalidParameter(__FILE__, __LINE__, "Parameter QUAD_MODEL_MAX_BLOCK_SIZE must be positive");
    }
    else if (quadModelBlockSize > modelMaxEval)
    {
        setAttributeValue("QUAD_MODEL_MAX_BLOCK_SIZE", modelMaxEval);
        quadModelBlockSize = modelMaxEval;
        quadModelChanged = !isAttributeDefaultValue<size_t>("QUAD_MODEL_MAX_BLOCK_SIZE");
    }

    if (0 == sgtelibModelBlockSize)
    {
        throw NOMAD::InvalidParameter(__FILE__, __LINE__, "Parameter SGTELIB_MODEL_MAX_BLOCK_SIZE must be positive");
    }
    else if (sgtelibModelBlockSize > modelMaxEval)
    {
        setAttributeValue("SGTELIB_MODEL_MAX_BLOCK_SIZE", modelMaxEval);
        sgtelibModelBlockSize = modelMaxEval;
        sgtelibModelChanged = !isAttributeDefaultValue<size_t>("SGTELIB_MODEL_MAX_BLOCK_SIZE");
    }

    if (quadModelChanged && sgtelibModelChanged && quadModelBlockSize != sgtelibModelBlockSize)
    {
        // See issue #505
        std::cerr << "Warning: Currently not supported: QUAD_MODEL_MAX_BLOCK_SIZE (";
        std::cerr << quadModelBlockSize << ") different than SGTELIB_MODEL_MAX_BLOCK_SIZE (";
        std::cerr << sgtelibModelBlockSize << "). Using only the value of QUAD_MODEL_MAX_BLOCK_SIZE." << std::endl;
        setAttributeValue("MODEL_MAX_BLOCK_SIZE", quadModelBlockSize);
    }
    else if (quadModelChanged)
    {
        setAttributeValue("MODEL_MAX_BLOCK_SIZE", quadModelBlockSize);
    }
    else
    {
        setAttributeValue("MODEL_MAX_BLOCK_SIZE", sgtelibModelBlockSize);
    }

    if (0 == getAttributeValueProtected<size_t>("EVAL_SURROGATE_COST", false))
    {
        throw NOMAD::InvalidParameter(__FILE__, __LINE__, "Parameter EVAL_SURROGATE_COST must be positive");
    }

    _toBeChecked = false;

}
// End checkAndComply()





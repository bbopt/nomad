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

#include "../Param/EvalParameters.hpp"
#include "../Type/BBInputType.hpp"
#include "../Type/BBOutputType.hpp"
#include "../Util/fileutils.hpp"

/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::EvalParameters::init()
{
    _typeName = "Eval";

    try
    {
        #include "../Attribute/evalAttributesDefinition.hpp"
        registerAttributes( _definition );

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
void NOMAD::EvalParameters::checkAndComply(const std::shared_ptr<NOMAD::RunParameters>& runParams,
                                           const std::shared_ptr<NOMAD::PbParameters>& pbParams)
{
    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    /*--------------------------*/
    /* BB_EXE and SURROGATE_EXE */
    /*--------------------------*/
    updateExeParam(runParams, "BB_EXE");
    updateExeParam(runParams, "SURROGATE_EXE");

    /*----------------*/
    /* BB_OUTPUT_TYPE */
    /*----------------*/
    // The default value is empty: set a single OBJ
    auto bbOType = getAttributeValueProtected<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE", false);
    if ( bbOType.size() == 0 )
    {
        bbOType.push_back(NOMAD::BBOutputType::OBJ);
        setAttributeValue("BB_OUTPUT_TYPE", bbOType );
    }

    /*---------------------------*/
    /* BB_EVAL_FORMAT (internal) */
    /*---------------------------*/
    
    // Copy of POINT_FORMAT
    auto pointFormat = pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("POINT_FORMAT");
    setAttributeValue("BB_EVAL_FORMAT", pointFormat);

    _toBeChecked = false;

}
// End checkAndComply()


void NOMAD::EvalParameters::updateExeParam(const std::shared_ptr<NOMAD::RunParameters>& runParams, const std::string& paramName)
{
    // Update BB_EXE / SURROGATE_EXE if it was set by user:
    // - Set full path
    // - Remove '$' indicating a global call (e.g. python, perl)
    // - Verify file is executable
    if (isSetByUser(paramName))
    {
        auto exe = getAttributeValueProtected<std::string>(paramName, false);

        bool localExe = true;
        auto problemDir = runParams->getAttributeValue<std::string>("PROBLEM_DIR");

        if ('$' == exe[0])
        {
            // When the '$' character is put in first
            // position of a string, it is considered
            // as global and no path will be added.
            localExe = false;
        }

        // Convert arguments; add path as needed.
        auto exeAsArray = NOMAD::ArrayOfString(exe);
        exe.clear();
        for (size_t i = 0; i < exeAsArray.size(); i++)
        {
            std::string word = exeAsArray[i];
            if (i > 0)
            {
                exe += " ";
            }

            if ('$' == word[0])
            {
                exe += word.substr(1, word.size());
            }
            else
            {
                // word is relative to problem directory.
                NOMAD::completeFileName(word, problemDir);
                exe += word;
            }
        }

        setAttributeValue(paramName, exe);
        exeAsArray = NOMAD::ArrayOfString(exe);

        if (localExe && !exe.empty() && !checkExeFile(exeAsArray[0]))
        {
            throw NOMAD::Exception(__FILE__, __LINE__, paramName + " needs to be an executable file: " + exeAsArray[0]);
        }
    }
}




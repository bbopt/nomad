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

#include "../Param/EvalParameters.hpp"
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
    catch ( NOMAD::Exception & e)
    {
        std::string errorMsg = "Attribute registration failed: ";
        errorMsg += e.what();
        throw NOMAD::Exception(__FILE__,__LINE__, errorMsg);
    }
}

/*----------------------------------------*/
/*            check the parameters        */
/*----------------------------------------*/
void NOMAD::EvalParameters::checkAndComply( const std::shared_ptr<NOMAD::RunParameters> & runParams )
{
    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    // Update BB_EXE if it was set by user:
    // - Set full path
    // - Remove '$' indicating a global call (e.g. python, perl)
    // - Verify file is executable
    if (isSetByUser("BB_EXE"))
    {
        auto bbExe = getAttributeValueProtected<std::string>("BB_EXE", false);
        if (bbExe.empty())
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "BB_EXE is not defined");
        }

        bool localExe = true;
        if ('$' == bbExe[0])
        {
            // When the '$' character is put in first
            // position of a string, it is considered
            // as global and no path will be added.
            bbExe = bbExe.substr(1, bbExe.size());
            setAttributeValue("BB_EXE", bbExe);
            localExe = false;
        }
        else
        {
            auto problemDir = runParams->getAttributeValue<std::string>("PROBLEM_DIR");
            // bbExe is relative to problem directory.
            completeFileName(bbExe, problemDir);
            setAttributeValue("BB_EXE", bbExe);
        }

        if (localExe && !checkExeFile(bbExe))
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "BB_EXE needs to be an executable file: " + bbExe);
        }
    }

    // The default value is empty: set a single OBJ
    auto bbOType = getAttributeValueProtected<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE", false);
    if ( bbOType.size() == 0 )
    {
        bbOType.push_back(NOMAD::BBOutputType::OBJ);
        setAttributeValue("BB_OUTPUT_TYPE", bbOType );
    }


    _toBeChecked = false;

}
// End checkAndComply()





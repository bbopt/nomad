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

#include "../Param/AllParameters.hpp"
#include "../Util/fileutils.hpp"

// Do we need to call checkAndComply() ?
bool NOMAD::AllParameters::toBeChecked() const
{
    bool check =    (!_pbParams          || _pbParams->toBeChecked())
                 || (!_evalParams        || _evalParams->toBeChecked())
                 || (!_evaluatorControlGlobalParams || _evaluatorControlGlobalParams->toBeChecked())
                 || (!_evaluatorControlParams || _evaluatorControlParams->toBeChecked())
                 || (!_runParams         || _runParams->toBeChecked())
                 || (!_cacheParams       || _cacheParams->toBeChecked())
                 || (!_dispParams        || _dispParams->toBeChecked());
    return check;
}



void NOMAD::AllParameters::resetToDefaultValues() noexcept
{
    _runParams->resetToDefaultValues();
    _pbParams->resetToDefaultValues();
    _cacheParams->resetToDefaultValues();
    _dispParams->resetToDefaultValues();
    _evalParams->resetToDefaultValues();
    _evaluatorControlGlobalParams->resetToDefaultValues();
    _evaluatorControlParams->resetToDefaultValues();
}

/*----------------------------------------------------------------*/
/*          read a parameters file and interpret attributes       */
/*----------------------------------------------------------------*/
void NOMAD::AllParameters::read(const std::string &paramFile, bool overwrite , bool resetAllEntries )
{

    // Read all entries
    NOMAD::Parameters::readParamFileAndSetEntries(paramFile, overwrite ,resetAllEntries );

    // Read entries to detect deprecated entries explicitly set
    _deprecatedParams->readAndDetectExplicitlySet();

    // Read entries and set attribute values for each type of parameters
    _runParams->readEntries();
    _pbParams->readEntries(false, NOMAD::dirname(paramFile));
    _evalParams->readEntries();
    _evaluatorControlGlobalParams->readEntries();
    _evaluatorControlParams->readEntries();
    _cacheParams->readEntries();
    _dispParams->readEntries();


}


/*----------------------------------------------------------------*/
/*   read a parameter line, interpret it, and update the links    */
/*----------------------------------------------------------------*/
void NOMAD::AllParameters::readParamLine(const std::string &line)
{
    // Note: Here we have the line, but we do not know which parameter type
    // it applies to. For now, the workaround is to create a temporary
    // ParameterEntry and to extract the parameter name from it.
    // This could be cleaner, since we will create another ParameterEntry
    // in readParamLine().
    auto pe = std::make_unique<NOMAD::ParameterEntry>(line);
    std::string name = pe->getName();

    bool overwrite = true;
    if ( _cacheParams->isRegisteredAttribute(name) )
    {
        _cacheParams->readParamLine(line, overwrite);
    }
    else if ( _dispParams->isRegisteredAttribute(name) )
    {
        _dispParams->readParamLine(line, overwrite);
    }
    else if (_evalParams->isRegisteredAttribute(name))
    {
        _evalParams->readParamLine(line, overwrite);
    }
    else if ( _evaluatorControlParams->isRegisteredAttribute(name) )
    {
        _evaluatorControlParams->readParamLine(line, overwrite);
    }
    else if ( _evaluatorControlGlobalParams->isRegisteredAttribute(name) )
    {
        _evaluatorControlGlobalParams->readParamLine(line, overwrite);
    }
    else if ( _pbParams->isRegisteredAttribute(name) )
    {
        _pbParams->readParamLine(line, overwrite);
    }
    else if ( _runParams->isRegisteredAttribute(name) )
    {
        _runParams->readParamLine(line, overwrite);
    }
    else
    {
        std::string err = "Unknown parameter: " + name;
        std::cerr << err << std::endl;
    }

}

void NOMAD::AllParameters::eraseAllEntries()
{
    NOMAD::Parameters::eraseAllEntries();
}


bool NOMAD::AllParameters::isAlgoCompatible(const NOMAD::AllParameters& allP_tmp) const
{
    return _pbParams->isAlgoCompatible(allP_tmp.getPbParams().get()) &&
           _dispParams->isAlgoCompatible(allP_tmp.getDispParams().get()) &&
           _runParams->isAlgoCompatible(allP_tmp.getRunParams().get()) &&
           _evaluatorControlGlobalParams->isAlgoCompatible(allP_tmp.getEvaluatorControlGlobalParams().get()) &&
           _evaluatorControlParams->isAlgoCompatible(allP_tmp.getEvaluatorControlParams().get()) &&
           _evalParams->isAlgoCompatible(allP_tmp.getEvalParams().get()) ;
}


/*----------------------------------------------------------------*/
/*          read a parameters file and interpret attributes       */
/*----------------------------------------------------------------*/
void NOMAD::AllParameters::displayHelp(const std::string &helpSubject , bool devHelp , std::ostream & os )
{

    std::ostringstream ossBasic,ossAdvanced;

    _evalParams->displayHelp(helpSubject, devHelp, ossBasic , ossAdvanced);
    _evaluatorControlGlobalParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _evaluatorControlParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _runParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _pbParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _cacheParams->displayHelp(helpSubject,devHelp,ossBasic , ossAdvanced);
    _dispParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);

    if ( !devHelp )
    {
        if (ossBasic.str().empty() && ossAdvanced.str().empty())
        {
            os << "No help found for " << helpSubject << std::endl << std::endl;
        }

        if (!ossBasic.str().empty())
        {
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << "-------------------------------- BASIC PARAMETERS -----------------------------" << std::endl;
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << std::endl << ossBasic.str() << std::endl << std::endl;
        }

        if (!ossAdvanced.str().empty())
        {
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << "------------------------------ ADVANCED PARAMETERS ----------------------------" << std::endl;
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << std::endl << ossAdvanced.str() << std::endl << std::endl;
        }
    }
    else
    {
        if ( ossBasic.str().empty() )
        {
            os << "No help found for " << helpSubject << std::endl << std::endl;
        }
        else
        {
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << "----------------------------- DEVELOPER PARAMETERS ----------------------------" << std::endl;
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << ossBasic.str() << std::endl << std::endl;
        }
    }
}

/*----------------------------------------*/
/*            check the parameters        */
/*----------------------------------------*/
void NOMAD::AllParameters::checkAndComply()
{
    if (!toBeChecked())
    {
        // Early out
        return;
    }

    _pbParams->checkAndComply();
    _evaluatorControlGlobalParams->checkAndComply(_pbParams);
    _runParams->checkAndComply(_evaluatorControlGlobalParams, _pbParams);
    _evaluatorControlParams->checkAndComply(_evaluatorControlGlobalParams, _runParams);
    _evalParams->checkAndComply(_runParams, _pbParams);
    _cacheParams->checkAndComply(_runParams);
    _dispParams->checkAndComply(_runParams, _pbParams);

}
// End checkAndComply()


std::string NOMAD::AllParameters::getSetAttributeAsString() const
{
    std::string str = _runParams->getSetAttributeAsString()
           + _pbParams->getSetAttributeAsString()
           + _evalParams->getSetAttributeAsString()
           + _evaluatorControlGlobalParams->getSetAttributeAsString()
           + _evaluatorControlParams->getSetAttributeAsString()
           + _cacheParams->getSetAttributeAsString()
           + _dispParams->getSetAttributeAsString() ;
    if ( str.empty() )
        str = "All attributes have default value";
    return str;
}


void NOMAD::AllParameters::display(std::ostream & os, bool flagHelp )
{
    if (toBeChecked())
    {
        std::cerr << "Warning: AllParameters::display(): Parameters are not checked." << std::endl;
    }

    os << "----- RUN PARAMETERS -----" << std::endl;
    _runParams->display(os,flagHelp);

    os << "----- PROBLEM PARAMETERS -----" << std::endl;
    _pbParams->display(os,flagHelp);

    os << "----- EVAL PARAMETERS -----" << std::endl;
    _evalParams->display(os,flagHelp);

    os << "----- EVALUATOR CONTROL PARAMETERS (GLOBAL) -----" << std::endl;
    _evaluatorControlGlobalParams->display(os,flagHelp);

    os << "----- EVALUATOR CONTROL PARAMETERS (BY MAIN THREAD)-----" << std::endl;
    _evaluatorControlParams->display(os,flagHelp);

    os << "----- CACHE PARAMETERS -----" << std::endl;
    _cacheParams->display(os,flagHelp);

    os << "----- DISPLAY PARAMETERS -----" << std::endl;
    _dispParams->display(os,flagHelp);
}



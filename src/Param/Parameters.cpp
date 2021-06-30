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

#include "../Param/Parameters.hpp"
#include "../Type/BBInputType.hpp"
#include "../Type/BBOutputType.hpp"
#include "../Type/DirectionType.hpp"
#include "../Type/EvalSortType.hpp"
#include "../Type/EvalType.hpp"
#include "../Type/LHSearchType.hpp"
#include "../Type/SgtelibModelFeasibilityType.hpp"
#include "../Type/SgtelibModelFormulationType.hpp"
#include "../Util/fileutils.hpp"
#ifdef _WIN32
#include <io.h>    // For _access
#define access _access
#define R_OK 04
#endif

// The deep copy of parameters. Used only by derived object that implemented the copy constructor and copy assignment
void NOMAD::Parameters::copyParameters(const Parameters& params)
{

    for(const auto &att : params.getAttributes() )
    {
        auto paramName = att->getName();
        auto paramType = getAttributeType(paramName);


        // BOOL
        if (paramType == typeid(bool).name())
        {
            auto value = params.getAttributeValue<bool>(paramName);
            setAttributeValue(paramName, value);
        }
        // SIZE_T
        else if (paramType == typeid(size_t).name())
        {
            auto value = params.getAttributeValue<size_t>(paramName);
            setAttributeValue(paramName, value);
        }
        // INTEGER
        else if (paramType == typeid(int).name())
        {
            auto value = params.getAttributeValue<int>(paramName);
            setAttributeValue(paramName, value);
        }
        // STRING
        else if ( paramType == typeid(std::string).name() )
        {
            auto value = params.getAttributeValue<std::string>(paramName);
            setAttributeValue(paramName , value );
        }
        // ARRAYOFSTRING
        else if ( paramType == typeid(NOMAD::ArrayOfString).name() )
        {
            auto value = params.getAttributeValue<NOMAD::ArrayOfString>(paramName);
            setAttributeValue(paramName, value);
        }
        // BBInputTypeList
        else if (paramType == typeid(NOMAD::BBInputTypeList).name())
        {
            auto value = params.getAttributeValue<NOMAD::BBInputTypeList>(paramName);
            setAttributeValue(paramName, value );
        }
        // BBOutputTypeList
        else if (paramType == typeid(NOMAD::BBOutputTypeList).name())
        {
            auto value = params.getAttributeValue<NOMAD::BBOutputTypeList>(paramName);
            setAttributeValue(paramName, value );
        }
        // LHSearchType
        else if (paramType == typeid(NOMAD::LHSearchType).name())
        {
            auto value = params.getAttributeValue<NOMAD::LHSearchType>(paramName);
            setAttributeValue(paramName, value );
        }
        // SgtelibModelFeasibilityType
        else if (paramType == typeid(NOMAD::SgtelibModelFeasibilityType).name())
        {
            auto value = params.getAttributeValue<NOMAD::SgtelibModelFeasibilityType>(paramName);
            setAttributeValue(paramName, value );
        }
        // SgtelibModelFormulationType
        else if (paramType == typeid(NOMAD::SgtelibModelFormulationType).name())
        {
            auto value = params.getAttributeValue<NOMAD::SgtelibModelFormulationType>(paramName);
            setAttributeValue(paramName, value );
        }
        // EvalType
        else if (paramType == typeid(NOMAD::EvalType).name())
        {
            auto value = params.getAttributeValue<NOMAD::EvalType>(paramName);
            setAttributeValue(paramName, value );
        }
        // EvalSortType
        else if (paramType == typeid(NOMAD::EvalSortType).name())
        {
            auto value = params.getAttributeValue<NOMAD::EvalSortType>(paramName);
            setAttributeValue(paramName, value );
        }
        // DirectionTypeList
        else if (paramType == typeid(NOMAD::DirectionTypeList).name())
        {
            auto value = params.getAttributeValue<NOMAD::DirectionTypeList>(paramName);
            setAttributeValue(paramName, value );
        }
        // DOUBLE
        else if (paramType == typeid(NOMAD::Double).name())
        {
            auto value = params.getAttributeValue<NOMAD::Double>(paramName);
            setAttributeValue(paramName, value );
        }
        // ARRAYOFDOUBLE
        else if (paramType == typeid(NOMAD::ArrayOfDouble).name() )
        {
            auto value = params.getAttributeValue<NOMAD::ArrayOfDouble>(paramName);
            setAttributeValue(paramName, value);
        }
        // POINT
        else if (paramType == typeid(NOMAD::Point).name())
        {
            auto value = params.getAttributeValue<NOMAD::Point>(paramName);
            setAttributeValue(paramName, value);
        }
        // ARRAY OF POINTS
        else if (paramType == typeid(NOMAD::ArrayOfPoint).name())
        {
            auto value = params.getAttributeValue<NOMAD::ArrayOfPoint>(paramName);
            setAttributeValue(paramName, value);
        }
        // LIST OF VARIABLE GROUP
        else if (paramType == typeid(NOMAD::ListOfVariableGroup).name())
        {
            auto value = params.getAttributeValue<NOMAD::ListOfVariableGroup>(paramName);
            setAttributeValue(paramName, value);
        }
        else
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Cannot copy a Parameter for a type not defined");
        }
    }
}


// Initialize static members
std::map<std::string,std::string> NOMAD::Parameters::_typeOfAttributes;
NOMAD::ParameterEntries NOMAD::Parameters::_paramEntries;

const NOMAD::AttributeSet NOMAD::Parameters::getAttributes() const
{
    return _attributes;
}


// Do we need to call checkAndComply() ?
bool NOMAD::Parameters::toBeChecked() const
{
    return _toBeChecked;
}


std::vector<std::string> NOMAD::Parameters::getAttributeNames() const
{
    std::vector<std::string> names;
    for (auto att : _attributes)
    {
        names.push_back(att->getName());
    }

    return names;
}


// All registered attributes are reset to their default value
void NOMAD::Parameters::resetToDefaultValues() noexcept
{
    for_each(_attributes.begin(), _attributes.end(), [](SPtrAtt att)
    {
        att->resetToDefaultValue();
    });

    _toBeChecked = true;
}


// Reset this particular attribute to its default value
void NOMAD::Parameters::resetToDefaultValue(const std::string& paramName)
{
    // Get attribute from which to get value
    SPtrAtt att = getAttribute(paramName);

    if (nullptr == att)
    {
        // At this point, Verify att is non-null.
        std::string err = "resetToDefaultValue: attribute " + paramName + " does not exist";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    att->resetToDefaultValue();

    _toBeChecked = true;
}


bool NOMAD::Parameters::isAlgoCompatible(const NOMAD::Parameters *p)
{
    bool isCompatible = true;
    const bool debugAlgoCompatibility = false; // Set to true for debugging purposes
    std::string sdebug;

    // Loop on all registered attributes.
    // Check if values are the same for attributes that are marked
    for(const auto &att : getAttributes())
    {
        if ( att->isForAlgoCompatibilityCheck() )
        {
            auto paramName = att->getName();
            auto paramType = getAttributeType(paramName);

            // The Protected version of getAttribute value is used to get the value even if the checkAndComply test has not been performed.

            // BOOL
            if ( paramType == typeid(bool).name() )
            {
                if ( getAttributeValueProtected<bool>(paramName,false) != p->getAttributeValueProtected<bool>(paramName,false) )
                {
                    sdebug += std::to_string(getAttributeValueProtected<bool>(paramName,false)) + "\n";
                    sdebug += std::to_string(p->getAttributeValueProtected<bool>(paramName,false));
                    isCompatible = false;
                }
            }
            // SIZE_T
            else if ( paramType == typeid(size_t).name() )
            {
                if ( getAttributeValueProtected<size_t>(paramName,false) != p->getAttributeValueProtected<size_t>(paramName,false) )
                {
                    sdebug += std::to_string(getAttributeValueProtected<size_t>(paramName,false)) + "\n";
                    sdebug += std::to_string(p->getAttributeValueProtected<size_t>(paramName,false));
                    isCompatible = false;
                }

            }
            // INTEGER
            else if ( paramType == typeid(int).name() )
            {
                if ( getAttributeValueProtected<int>(paramName,false) != p->getAttributeValueProtected<int>(paramName,false) )
                {
                    sdebug += std::to_string(getAttributeValueProtected<int>(paramName,false)) + "\n";
                    sdebug += std::to_string(p->getAttributeValueProtected<int>(paramName,false));
                    isCompatible = false;
                }
            }
            // NOMAD::DOUBLE
            else if ( paramType == typeid(NOMAD::Double).name() )
            {
                if ( getAttributeValueProtected<NOMAD::Double>(paramName,false) != p->getAttributeValueProtected<NOMAD::Double>(paramName,false) )
                {
                    sdebug += getAttributeValueProtected<NOMAD::Double>(paramName,false).tostring() + "\n";
                    sdebug += p->getAttributeValueProtected<NOMAD::Double>(paramName,false).tostring();
                    isCompatible = false;
                }
            }
            // double
            else if ( paramType == typeid(double).name() )
            {
                if ( getAttributeValueProtected<double>(paramName,false) != p->getAttributeValueProtected<double>(paramName,false) )
                {
                    sdebug += std::to_string(getAttributeValueProtected<double>(paramName,false)) + "\n";
                    sdebug += std::to_string(p->getAttributeValueProtected<double>(paramName,false));
                    isCompatible = false;
                }
            }
            // STRING
            else if ( paramType == typeid(std::string).name() )
            {
                if ( getAttributeValueProtected<std::string>(paramName,false).compare( p->getAttributeValueProtected<std::string>(paramName,false)) != 0 )
                {
                    sdebug += getAttributeValueProtected<std::string>(paramName,false) + "\n";
                    sdebug += p->getAttributeValueProtected<std::string>(paramName,false);
                    isCompatible = false;
                }
            }
            // ARRAY OF STRING
            else if ( paramType == typeid(NOMAD::ArrayOfString).name() )
            {
                if ( getAttributeValueProtected<NOMAD::ArrayOfString>(paramName,false) != p->getAttributeValueProtected<NOMAD::ArrayOfString>(paramName,false) )
                {
                    sdebug += getAttributeValueProtected<NOMAD::ArrayOfString>(paramName,false).display() + "\n";
                    sdebug += p->getAttributeValueProtected<NOMAD::ArrayOfString>(paramName,false).display();
                    isCompatible = false;
                }
            }
            // ARRAY OF DOUBLE
            else if ( paramType == typeid(NOMAD::ArrayOfDouble).name() )
            {
                if ( getAttributeValueProtected<NOMAD::ArrayOfDouble>(paramName,false) != p->getAttributeValueProtected<NOMAD::ArrayOfDouble>(paramName,false) )
                {
                    sdebug += getAttributeValueProtected<NOMAD::ArrayOfDouble>(paramName,false).display() + "\n";
                    sdebug += p->getAttributeValueProtected<NOMAD::ArrayOfDouble>(paramName,false).display();
                    isCompatible = false;
                }
            }
            // LIST OF VARIABLE GROUP
            else if ( paramType == typeid(NOMAD::ListOfVariableGroup).name() )
            {
                // For the comparison, if the groups of variables are not in the same order
                // the two instances are not algo compatible (the runs will be different!)
                auto lvg = getAttributeValueProtected<NOMAD::ListOfVariableGroup>(paramName,false);
                auto plvg = p->getAttributeValueProtected<NOMAD::ListOfVariableGroup>(paramName,false);

                if (lvg.size() != plvg.size() || lvg != plvg )
                {
                    std::ostringstream sds,sds2;
                    sds << lvg;
                    sds2 << plvg;
                    sdebug += sds.str() + "\n";
                    sdebug += sds2.str() ;
                    isCompatible = false;
                }
            }
            // Typesdebug defined in the Type directory
            // BBInputType
            else if ( paramType == typeid(NOMAD::BBInputType).name() )
            {
                if ( getAttributeValueProtected<NOMAD::BBInputType>(paramName,false) != p->getAttributeValueProtected<NOMAD::BBInputType>(paramName,false) )
                {
                    std::ostringstream oss;
                    oss << getAttributeValueProtected<NOMAD::BBInputType>(paramName,false) << std::endl;
                    oss << p->getAttributeValueProtected<NOMAD::BBInputType>(paramName,false);
                    sdebug = oss.str();
                    isCompatible = false;
                }
            }
            // BBOutputType
            else if ( paramType == typeid(NOMAD::BBOutputType).name() )
            {
                if ( getAttributeValueProtected<NOMAD::BBOutputType>(paramName,false) != p->getAttributeValueProtected<NOMAD::BBOutputType>(paramName,false) )
                {
                    std::ostringstream oss;
                    oss << getAttributeValueProtected<NOMAD::BBOutputType>(paramName,false) << std::endl;
                    oss << p->getAttributeValueProtected<NOMAD::BBOutputType>(paramName,false);
                    sdebug = oss.str();
                    isCompatible = false;
                }
            }
            // EvalType
            else if ( paramType == typeid(NOMAD::EvalType).name() )
            {
                if ( getAttributeValueProtected<NOMAD::EvalType>(paramName,false) != p->getAttributeValueProtected<NOMAD::EvalType>(paramName,false) )
                {
                    sdebug += NOMAD::evalTypeToString(getAttributeValueProtected<NOMAD::EvalType>(paramName,false)) + "\n";
                    sdebug += NOMAD::evalTypeToString(p->getAttributeValueProtected<NOMAD::EvalType>(paramName,false));
                    isCompatible = false;
                }
            }
            // EvalSortType
            else if ( paramType == typeid(NOMAD::EvalSortType).name() )
            {
                if ( getAttributeValueProtected<NOMAD::EvalSortType>(paramName,false) != p->getAttributeValueProtected<NOMAD::EvalSortType>(paramName,false) )
                {
                    sdebug += NOMAD::evalSortTypeToString(getAttributeValueProtected<NOMAD::EvalSortType>(paramName,false)) + "\n";
                    sdebug += NOMAD::evalSortTypeToString(p->getAttributeValueProtected<NOMAD::EvalSortType>(paramName,false));
                    isCompatible = false;
                }
            }
            // DirectionTypeList
            else if ( paramType == typeid(NOMAD::DirectionTypeList).name() )
            {
                if ( getAttributeValueProtected<NOMAD::DirectionTypeList>(paramName,false) != p->getAttributeValueProtected<NOMAD::DirectionTypeList>(paramName,false) )
                {
                    sdebug += NOMAD::directionTypeListToString(getAttributeValueProtected<NOMAD::DirectionTypeList>(paramName,false)) + "\n";
                    sdebug += NOMAD::directionTypeListToString(p->getAttributeValueProtected<NOMAD::DirectionTypeList>(paramName,false));
                    isCompatible = false;
                }
            }
            // LHSearchType
            else if ( paramType == typeid(NOMAD::LHSearchType).name() )
            {
                if ( getAttributeValueProtected<NOMAD::LHSearchType>(paramName,false) != p->getAttributeValueProtected<NOMAD::LHSearchType>(paramName,false) )
                {
                    std::ostringstream oss;
                    oss << getAttributeValueProtected<NOMAD::LHSearchType>(paramName,false) << std::endl;
                    oss << p->getAttributeValueProtected<NOMAD::LHSearchType>(paramName,false);
                    sdebug = oss.str();
                    isCompatible = false;
                }
            }
            // SgtelibModelFeasibilityType
            else if ( paramType == typeid(NOMAD::SgtelibModelFeasibilityType).name() )
            {
                if ( getAttributeValueProtected<NOMAD::SgtelibModelFeasibilityType>(paramName,false) != p->getAttributeValueProtected<NOMAD::SgtelibModelFeasibilityType>(paramName,false) )
                {
                    sdebug += NOMAD::SgtelibModelFeasibilityTypeToString(getAttributeValueProtected<NOMAD::SgtelibModelFeasibilityType>(paramName,false)) + "\n";
                    sdebug += NOMAD::SgtelibModelFeasibilityTypeToString(p->getAttributeValueProtected<NOMAD::SgtelibModelFeasibilityType>(paramName,false));
                    isCompatible = false;
                }
            }
            // SgtelibModelFormulationType
            else if ( paramType == typeid(NOMAD::SgtelibModelFormulationType).name() )
            {
                if ( getAttributeValueProtected<NOMAD::SgtelibModelFormulationType>(paramName,false) != p->getAttributeValueProtected<NOMAD::SgtelibModelFormulationType>(paramName,false) )
                {
                    sdebug += NOMAD::SgtelibModelFormulationTypeToString(getAttributeValueProtected<NOMAD::SgtelibModelFormulationType>(paramName,false)) + "\n";
                    sdebug += NOMAD::SgtelibModelFormulationTypeToString(p->getAttributeValueProtected<NOMAD::SgtelibModelFormulationType>(paramName,false));
                    isCompatible = false;
                }
            }
            else
            {
                std::string err = "Error: Cannot test the type " + paramType + " for compatibility for parameter " + paramName;
                throw NOMAD::Exception(__FILE__, __LINE__, err);
            }

            if (!isCompatible)
            {
                // Debug info about which parameter names and values
                // were not compatible.
                if (debugAlgoCompatibility)
                {
                    // Prepend info to string sdebug
                    sdebug = "Parameter values are not compatible: Parameter name: " + paramName + "; parameter type: " + paramType + "; parameter values:\n" + sdebug;
                    std::cerr << sdebug << std::endl;
                }
                break;
            }
        }
    }

    return isCompatible;
}


/*----------------------------------------------------------------*/
/*          read a parameters file and interpret attributes       */
/*----------------------------------------------------------------*/
void NOMAD::Parameters::readParamFileAndSetEntries(const std::string &paramFile,
                                                   bool overwrite,
                                                   bool resetAllEntries          )
{
    // Warning:
    // This method is static, so it cannot set _toBeChecked on the
    // Parameters object. It will be set by readEntries().

    // Open the parameters file:
    std::string err = "Could not open parameters file \'" + paramFile + "\'";
    std::ifstream fin;
    if ( access ( paramFile.c_str() , R_OK ) == 0 )
    {
        fin.open ( paramFile.c_str() );
        if ( !fin.fail() )
        {
            err.clear();
        }
    }
    if ( !err.empty() )
    {
        fin.close();
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    // the file is read: fill the set '_paramEntries' of ParameterEntry:
    std::string line;
    int line_number = 0;

    if ( resetAllEntries )
        eraseAllEntries();

    while (fin.good() && !fin.eof())
    {
        line.clear();

        getline(fin, line);
        line_number++;

        if (!fin.fail() && !line.empty())
        {
            readParamLine(line, paramFile, line_number, overwrite);
        }
    }

    // done with the file
    fin.close();
}


void NOMAD::Parameters::readParamLine(const std::string &line, bool overwrite)
{
    // parameters will have to be checked:
    _toBeChecked = true;

    // Read a single parameter from standard input.
    std::string paramFile = "Standard Input";
    int line_number = 0;

    // entries will be a set of 1 entry. Later the method read(ParameterEntries) is re-used.
    // NOMAD::ParameterEntries entries;

    readParamLine(line, paramFile, line_number, overwrite);

    try
    {
        readEntries();
    }
    catch (NOMAD::Exception &e)
    {
        // Show exceptions thrown when reading entries from standard input,
        // but do not re-throw them.
        std::cerr << "Warning: " << e.what() << std::endl;
    }
}


void NOMAD::Parameters::readParamLine(const std::string &line,
                                      const std::string &paramFile,
                                      const int line_number,
                                      bool overwrite)
{
    auto pe = std::make_shared<NOMAD::ParameterEntry>(line);
    if (nullptr == pe)
    {
        std::string err = "readParamLine: Error: Could not create parameter entry for parameter " + pe->getName();
        throw NOMAD::Exception(paramFile, line_number, err);
    }
    pe->setParamFile(paramFile);
    pe->setLine(line_number);

    if (pe->isOk())
    {
        if (overwrite)
        {
            std::shared_ptr<NOMAD::ParameterEntry> refPe = _paramEntries.find(pe->getName());
            // Erase all parameters with this name.
            if (nullptr != refPe)
            {
                _paramEntries.erase(refPe);
            }
        }
        _paramEntries.insert(pe);
    }
    else
    {
        if (pe->getName() != "" && pe->getNbValues() == 0)
        {
            std::string err = "Invalid parameter: " + pe->getName();
            // If reading a file (positive line number), throw an exception.
            // Else, only show a warning.
            if (line_number > 0)
            {
                throw NOMAD::Exception(paramFile, line_number, err);
            }
            else
            {
                std::cerr << "Warning: " << err << std::endl;
            }
        }
    }

}


/*----------------------------------------*/
/*          read parameter entries        */
/*----------------------------------------*/
void NOMAD::Parameters::readEntries(const bool overwrite, std::string problemDir)
{
    if (problemDir.empty())
    {
        problemDir = std::string(".") + NOMAD::DIR_SEP;
    }

    // parameters will have to be checked:
    _toBeChecked = true;

    // interpret and set the entries using SET methods:
    size_t sz;
    int i;
    NOMAD::Double d;
    bool flag;
    std::shared_ptr<NOMAD::ParameterEntry> pe;
    std::string err;

    // First set DIMENSION from entries
    std::string paramName = "DIMENSION";
    auto attDim = getAttribute(paramName);

    if ( nullptr != attDim )
    {
        pe = _paramEntries.find(paramName);
        if ( pe )
        {
            if ( !pe->isUnique() )
            {
                err = "Multiple entries of parameter: " + paramName + " at line #" + std::to_string(pe->getLine()) + "! ";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }

            checkFormatSizeT(pe, sz);
            setAttributeValue("DIMENSION", sz);
            pe->setHasBeenInterpreted();
        }
    }
    /*
    else
    {
        err = "Missing mandatory DIMENSION parameter! ";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
    */

    // Loop on all registered attributes.
    // Set the attribute value from entries
    // according to its type
    for(const auto &att : getAttributes())
    {
        paramName = att->getName();
        auto paramType = getAttributeType(paramName);
        bool paramFromUniqueEntry = att->getParamFromUniqueEntry();
        bool internal = att->isInternal();

        // Dimension already set
        if (paramName == "DIMENSION")
        {
            continue;
        }

        pe = _paramEntries.find(paramName);
        while (pe)
        {
            // Test if internal
            if (internal)
            {
                err = "Internal parameter: " + paramName + " at line #" + std::to_string(pe->getLine()) + " cannot be set in a parameter file! ";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }

            // Test if multiple entries
            if (paramFromUniqueEntry && !pe->isUnique() && !overwrite)
            {
                err = "Multiple entries of parameter: " + paramName + " at line #" + std::to_string(pe->getLine()) + "! ";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }

            // BOOL
            if (paramType == typeid(bool).name())
            {
                checkFormatBool(pe);
                flag = NOMAD::stringToBool( *(pe->getValues().begin()) );
                setAttributeValue(paramName, flag);
            }
            // SIZE_T
            else if (paramType == typeid(size_t).name())
            {
                checkFormatSizeT(pe, sz);
                setAttributeValue(paramName, sz);
            }
            // INTEGER
            else if (paramType == typeid(int).name())
            {
                checkFormatInt(pe, i);
                setAttributeValue(paramName, i);
            }
            // STRING
            else if ( paramType == typeid(std::string).name() )
            {
                checkFormatString(pe);
                setAttributeValue(paramName , *(pe->getValues().begin()) );
            }
            // ARRAY OF STRING
            else if ( paramType == typeid(NOMAD::ArrayOfString).name() )
            {
                checkFormatArrayOfString(pe);
                NOMAD::ArrayOfString aos(pe->getAllValues());
                setAttributeValue(paramName, aos);
            }
            // BBInputTypeList
            else if (paramType == typeid(NOMAD::BBInputTypeList).name())
            {
                checkFormat1(pe);
                setAttributeValue(paramName, NOMAD::stringToBBInputTypeList(pe->getAllValues()));
            }
            // BBOutputTypeList
            else if (paramType == typeid(NOMAD::BBOutputTypeList).name())
            {
                checkFormat1(pe);
                setAttributeValue(paramName, NOMAD::stringToBBOutputTypeList(pe->getAllValues()));
            }
            // LHSearchType
            else if (paramType == typeid(NOMAD::LHSearchType).name())
            {
                checkFormatNbEntries(pe, 2);
                checkFormatAllSizeT(pe);
                std::string allValues = pe->getAllValues();
                auto lhSearch = NOMAD::LHSearchType(allValues);
                setAttributeValue(paramName, lhSearch);
            }
            // SgtelibModelFeasibilityType
            else if (paramType == typeid(NOMAD::SgtelibModelFeasibilityType).name())
            {
                checkFormat1(pe);
                setAttributeValue(paramName, NOMAD::stringToSgtelibModelFeasibilityType(pe->getAllValues()));
            }
            // SgtelibModelFormulationType
            else if (paramType == typeid(NOMAD::SgtelibModelFormulationType).name())
            {
                checkFormat1(pe);
                setAttributeValue(paramName, NOMAD::stringToSgtelibModelFormulationType(pe->getAllValues()));
            }
            // EvalType
            else if (paramType == typeid(NOMAD::EvalType).name())
            {
                checkFormat1(pe);
                setAttributeValue(paramName, NOMAD::stringToEvalType(pe->getAllValues()));
            }
            // EvalSortType
            else if (paramType == typeid(NOMAD::EvalSortType).name())
            {
                checkFormat1(pe);
                setAttributeValue(paramName, NOMAD::stringToEvalSortType(pe->getAllValues()));
            }
            // DirectionTypeList
            else if (paramType == typeid(NOMAD::DirectionTypeList).name())
            {
                checkFormat1(pe);

                auto dirTypeList = getAttributeValueProtected<NOMAD::DirectionTypeList>(paramName, false);
                if (isAttributeDefaultValue<NOMAD::DirectionTypeList>(paramName))
                {
                    dirTypeList.clear();
                }
                auto newDirType = NOMAD::stringToDirectionType(pe->getAllValues());
                if (dirTypeList.end() == std::find(dirTypeList.begin(), dirTypeList.end(), newDirType))
                {
                    dirTypeList.push_back(newDirType);
                }

                setAttributeValue(paramName, dirTypeList);
            }
            // DOUBLE
            else if (paramType == typeid(NOMAD::Double).name())
            {
                checkFormatDouble(pe, d);
                setAttributeValue(paramName , d);
            }
            // ARRAYOFDOUBLE; POINT
            else if (paramType == typeid(NOMAD::ArrayOfDouble).name()
                     || paramType == typeid(NOMAD::Point).name())
            {
                if ( isAttributeDefaultValue<size_t>("DIMENSION") )
                {
                    throw NOMAD::Exception(__FILE__,__LINE__,"Dimension must be set!");
                }

                // Dimension is needed for ArrayOfDouble because
                // several formats of the provided values are allowed
                const size_t n = getAttributeValueProtected<size_t>("DIMENSION", false);

                NOMAD::ArrayOfDouble aod(n);
                readValuesAsArray(*pe, aod);

                if ( paramType == typeid(NOMAD::Point).name() )
                {
                    //Convert aod to Point before setting parameter
                    NOMAD::Point point(n);
                    for (size_t index = 0; index < n; index++)
                    {
                        point[index] = aod[index];
                    }
                    setAttributeValue(paramName, point);
                }
                else
                {
                    setAttributeValue(paramName, aod);
                }
            }
            // Vector of POINTS
            else if (paramType == typeid(NOMAD::ArrayOfPoint).name())
            {
                if ( isAttributeDefaultValue<size_t>("DIMENSION") )
                {
                    throw NOMAD::Exception(__FILE__,__LINE__,"Dimension must be set!");
                }
                const size_t n = getAttributeValueProtected<size_t>("DIMENSION", false);

                auto aopRef = getAttributeValueProtected<NOMAD::ArrayOfPoint>(paramName, false);
                NOMAD::ArrayOfPoint aop = aopRef;

                if (1 == pe->getValues().size())
                {
                    // Consider we have a file and read points in this file.
                    std::string pointFile = *pe->getValues().begin();
                    NOMAD::completeFileName(pointFile, problemDir);
                    auto aopNew = readPointValuesFromFile(pointFile);
                    for (auto newPoint: aopNew)
                    {
                        aop.push_back(newPoint);
                    }
                }
                else
                {
                    NOMAD::Point updatePoint(n);
                    size_t pointIndex = readValuesForArrayOfPoint(*pe, updatePoint);
                    // If the aop at this index is already set, update it.
                    // Else, create it.
                    if (pointIndex < aop.size())
                    {
                        auto pointToUpdate = aop[pointIndex];
                        for (size_t index = 0; index < n; index++)
                        {
                            if (updatePoint[index].isDefined())
                            {
                                pointToUpdate[index] = updatePoint[index];
                            }
                        }
                        aop[pointIndex] = pointToUpdate;
                    }
                    else
                    {
                        aop.resize(pointIndex+1);
                        aop[pointIndex] = updatePoint;
                    }
                }

                setAttributeValue(paramName, aop);
            }
            // List of VARIABLE_GROUP
            else if (paramType == typeid(NOMAD::ListOfVariableGroup).name())
            {
                if ( isAttributeDefaultValue<size_t>("DIMENSION") )
                {
                    throw NOMAD::Exception(__FILE__,__LINE__,"Dimension must be set!");
                }
                const size_t n = getAttributeValueProtected<size_t>("DIMENSION", false);

                auto lvg = getAttributeValueProtected<NOMAD::ListOfVariableGroup>(paramName, false);

                NOMAD::VariableGroup aVariableGroup;
                size_t nbIndex = readValuesForVariableGroup(*pe, aVariableGroup);
                // If the aop at this index is already set, update it.
                // Else, create it.
                if (nbIndex < n)
                {
                    lvg.push_back(aVariableGroup);
                }
                else
                {
                    err = "Number of indices for VARIABLE_GROUP must be smaller than DIMENSION (" + std::to_string(n) + ").";
                    throw NOMAD::Exception(__FILE__,__LINE__, err);
                }

                setAttributeValue(paramName, lvg);
            }
            else
            {
                err = "Parameter " + paramName + " has been registered but its type cannot be read";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            pe->setHasBeenInterpreted();
            // Get next parameter entry with this name, if multiple entries are permitted
            pe = pe->getNext();
        }
    }
}


void NOMAD::Parameters::checkFormat1(const std::shared_ptr<NOMAD::ParameterEntry> pe) const
{
    if (pe->getNbValues() < 1)
    {
        std::string err = "Invalid format for parameter: ";
        err += pe->getName() + " at line " + std::to_string(pe->getLine());
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
}


void NOMAD::Parameters::checkFormatNbEntries(const std::shared_ptr<NOMAD::ParameterEntry> pe, const size_t nbEntries) const
{
    if (pe->getNbValues() != nbEntries)
    {
        std::string err = "Parameter ";
        err += pe->getName();
        err += " expects exactly " + NOMAD::itos(nbEntries);
        err += " values, at line " + std::to_string(pe->getLine());
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
}


void NOMAD::Parameters::checkFormatBool(const std::shared_ptr<NOMAD::ParameterEntry> pe) const
{
    if (pe->getNbValues() != 1 )
    {
        std::string err = "Invalid format for bool parameter: ";
        err += pe->getName() + " at line " + std::to_string(pe->getLine());
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


void NOMAD::Parameters::checkFormatSizeT(const std::shared_ptr<NOMAD::ParameterEntry> pe, size_t &sz) const
{
    int i = -1;
    if (pe->getNbValues() != 1
        || !NOMAD::atoi ( *(pe->getValues().begin()), i)
        || i < 0)
    {
        std::string err = "Invalid format for size_t parameter: ";
        err += pe->getName() + " at line " + std::to_string(pe->getLine());
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
    sz = i;
}


void NOMAD::Parameters::checkFormatAllSizeT(const std::shared_ptr<NOMAD::ParameterEntry> pe) const
{
    int i;

    for (auto value : pe->getValues())
    {
        if (!NOMAD::atoi(value, i) || i < 0)
        {
            std::string err = "Invalid format for size_t parameter: ";
            err += pe->getName() + " at line " + std::to_string(pe->getLine());
            throw NOMAD::Exception(__FILE__,__LINE__, err);
        }
    }
}


void NOMAD::Parameters::checkFormatInt(const std::shared_ptr<NOMAD::ParameterEntry> pe, int &i) const
{
    if (pe->getNbValues() != 1 || ! NOMAD::atoi ( *(pe->getValues().begin()), i ))
    {
        std::string err = "Invalid format for integer parameter: ";
        err += pe->getName() + " at line " + std::to_string(pe->getLine());
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
}


void NOMAD::Parameters::checkFormatString(const std::shared_ptr<NOMAD::ParameterEntry> pe) const
{
    if ( pe->getNbValues() != 1 )
    {
        std::string err = "Invalid format for string parameter: ";
        err += pe->getName() + " at line " + std::to_string(pe->getLine());
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
}


void NOMAD::Parameters::checkFormatArrayOfString(const std::shared_ptr<NOMAD::ParameterEntry> pe) const
{
    // Do nothing
}


void NOMAD::Parameters::checkFormatDouble(const std::shared_ptr<NOMAD::ParameterEntry> pe, NOMAD::Double &d) const
{
    if (pe->getNbValues() != 1 || !d.atof ( *(pe->getValues().begin())))
    {
        std::string err = "Invalid format for double parameter: ";
        err += pe->getName() + " at line " + std::to_string(pe->getLine());
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
}

void NOMAD::Parameters::checkInfo( void ) const
{
    // Test if Info has been provided
    for_each(_attributes.begin(),_attributes.end(), [](SPtrAtt att){
        if ( att->hasEmptyInfo() )
        {
            std::string err = "Check: empty info (Short info and/or Help info) for attribute " + att->getName() + "!";
            throw NOMAD::Exception(__FILE__,__LINE__, err );
        }
    });
}


void NOMAD::Parameters::readValuesAsArray(const NOMAD::ParameterEntry &pe,
                                          NOMAD::ArrayOfDouble &array)
{
    // Convert list of strings to ArrayOfStrings...
    const std::list<std::string> peValues = pe.getValues();
    NOMAD::ArrayOfString arrayOfStrings;
    std::list<std::string>::const_iterator it;
    for (it = peValues.begin(); it != peValues.end(); ++it)
    {
        arrayOfStrings.add((*it));
    }

    array.readValuesAsArray(arrayOfStrings);
}


size_t NOMAD::Parameters::readValuesForArrayOfPoint(const NOMAD::ParameterEntry &pe,
                                                    NOMAD::Point &point)
{
    size_t index = 0;

    // Convert list of strings to ArrayOfStrings...
    const std::list<std::string> peValues = pe.getValues();
    NOMAD::ArrayOfString arrayOfStrings;
    std::list<std::string>::const_iterator it;
    for (it = peValues.begin(); it != peValues.end(); ++it)
    {
        arrayOfStrings.add((*it));
    }

    // Verify if there is an index. If so, save it, and remove it from the
    // arrayOfStrings.
    std::string firstElem = arrayOfStrings[0];
    NOMAD::Double d;
    d.atof(firstElem);
    if (d.isInteger())
    {
        index = size_t(d.todouble());
        arrayOfStrings.erase(0);
    }
    //
    // Convert arrayOfStrings to point.
    point.readValuesAsArray(arrayOfStrings);

    return index;
}


NOMAD::ArrayOfPoint NOMAD::Parameters::readPointValuesFromFile(const std::string& pointFile)
{
    if (!NOMAD::checkReadFile(pointFile))
    {
        std::string err = "File does not exist or is not readable: " + pointFile;
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    const size_t n = getAttributeValueProtected<size_t>("DIMENSION", false);
    NOMAD::ArrayOfPoint aop;
    NOMAD::Point point(n);  // Empty point used to pass the dimension
    aop.push_back(point);
    NOMAD::read<NOMAD::ArrayOfPoint>(aop, pointFile);  // Calls ArrayOfPoint::operator>>

    return aop;
}


size_t NOMAD::Parameters::readValuesForVariableGroup(const NOMAD::ParameterEntry &pe,
                                                     NOMAD::VariableGroup &vg )
{
    int i;

    std::list<std::string>::const_iterator it , end;
    std::pair<NOMAD::VariableGroup::iterator,bool> ret;
    // just one variable index (can be '*' or a range of indices 'i-j'):
    if ( pe.getNbValues() == 1 )
    {
        int j, k;
        it = pe.getValues().begin();
        if ( !NOMAD::stringToIndexRange ( *it , i , j ) )
        {
            std::string err = "Invalid format for index range: ";
            err += pe.getName() + " at line " + std::to_string(pe.getLine());
            throw NOMAD::Exception(__FILE__,__LINE__, err);
        }

        for ( k = i ; k <= j ; k++ )
        {
            ret = vg.insert(k);
            if (!ret.second)
            {
                std::string err = "Invalid index. Duplicate index not allowed: ";
                err += pe.getName() + " at line " + std::to_string(pe.getLine());
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
        }
    }

    // list of variable indexes:
    else
    {
        end = pe.getValues().end();
        for ( it = pe.getValues().begin() ; it != end ; ++it )
        {
            size_t ist = (size_t)i;
            if ( !NOMAD::atost ( *it , ist ) )
            {
                    std::string err = "Invalid format for index list: ";
                    err += pe.getName() + " at line " + std::to_string(pe.getLine());
                    throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            i = (int)ist;
            ret = vg.insert(i);
            if (!ret.second)
            {
                std::string err = "Invalid index. Duplicate index not allowed: ";
                err += pe.getName() + " at line " + std::to_string(pe.getLine());
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
        }
    }

    return vg.size();
}

// Get an attribute by its name (works lower/upper case)
NOMAD::SPtrAtt NOMAD::Parameters::getAttribute(std::string name) const
{
    NOMAD::toupper(name);

    auto it =   find_if(_attributes.begin(), _attributes.end(),
                        [name](NOMAD::SPtrAtt const& Att)
                        {
                            return name.compare(Att->getName()) == 0;
                        });

    if (it != _attributes.end())
    {
        return *it;
    }
    else
    {
        return nullptr;
    }
}


bool NOMAD::Parameters::isSetByUser(const std::string& paramName) const
{
    return (nullptr != _paramEntries.find(paramName));
}


void NOMAD::Parameters::displayHelp(const std::string & helpSubject , bool devHelp, std::ostringstream & ossBasic, std::ostringstream & ossAdvanced)
{
    // Search is performed on touppered strings

    // Display help as Basic, Advanced or Developer
    // Separate Basic and Advanced into sections
    std::ostringstream oss;
    for(const auto &att: _attributes)
    {
        oss.str("");
        oss.clear();
        std::string paramName(att->getName());
        std::string helpInfo(att->getHelpInfo());
        std::string keywords(att->getKeywords());
        NOMAD::toupper(paramName);
        NOMAD::toupper(helpInfo);
        NOMAD::toupper(keywords);

        if (   helpSubject.empty()
            || paramName.find(helpSubject) != std::string::npos
            || keywords.find(helpSubject)  != std::string::npos
            || helpInfo.find(helpSubject)  != std::string::npos )
        {
            if ( !devHelp || (keywords.find("DEVELOPER") != std::string::npos) )
            {
                std::string typeOfHelp = (devHelp) ? "(Developer)":"(Basic)";
                typeOfHelp = ( keywords.find("ADVANCED")!= std::string::npos )? "(Advanced)":typeOfHelp;

                oss << att->getName() << " {" ;
                oss << att->getHelpInfo() << std::endl;
                oss << "}" << std::endl;

                if ( typeOfHelp == "(Basic)" || typeOfHelp == "(Developer)" )
                {
                    ossBasic << oss.str() << std::endl;
                }
                else
                {
                    ossAdvanced << oss.str() << std::endl;
                }
            }
        }
    }
}


void NOMAD::Parameters::display(std::ostream & os, bool helpInfo)
{
    if ( !helpInfo && toBeChecked())
    {
        std::cerr << "Warning: Parameters::display(): Parameters are not checked." << std::endl;
    }

    for(const auto &att: _attributes)
    {
        if ( helpInfo )
            os << att->getHelpInfo() << std::endl;
        else
            os << *att << std::endl;
    }

}


void NOMAD::Parameters::registerAttributes( const std::vector<NOMAD::AttributeDefinition> & attributeDef )
{


    // Developers can provide new attributes type -> add "else if" below.
    // Developers MUST ALSO ADD the corresponding type in isAlgoCompatible

    for ( auto &att : attributeDef )
    {
        bool algoCompatibilityCheck = NOMAD::stringToBool(att._algoCompatibilityCheck);
        bool restartAttribute = NOMAD::stringToBool(att._restartAttribute);
        bool uniqueEntry = NOMAD::stringToBool(att._uniqueEntry);
        if ( att._type == "int" )
        {
            int defVal;
            if ( ! NOMAD::atoi( att._defaultValue.c_str() ,defVal ) )
            {
                std::string err = "Invalid int attribute definition: ";
                err +=  att._name + ". Infinity is set with INF or +INF or -INF.";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            registerAttribute<int>(att._name, defVal,
                            algoCompatibilityCheck, restartAttribute, uniqueEntry,
                            att._shortInfo, att._helpInfo, att._keywords);
        }
        else if ( att._type == "size_t" )
        {
            size_t defVal;
            if ( ! NOMAD::atost( att._defaultValue.c_str() ,defVal ) )
            {
                std::string err = "Invalid size_t attribute definition: ";
                err +=  att._name + ". Infinity is set with INF or +INF.";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            registerAttribute<size_t>(  att._name, defVal, algoCompatibilityCheck , restartAttribute, uniqueEntry,
                                        att._shortInfo, att._helpInfo, att._keywords);
        }
        else if ( att._type== "Double" || att._type== "NOMAD::Double" )
        {
            NOMAD::Double defVal;
            if ( defVal.atof( att._defaultValue) )
            {
                registerAttribute<NOMAD::Double>(att._name, defVal,
                            algoCompatibilityCheck, restartAttribute, uniqueEntry,
                            att._shortInfo, att._helpInfo, att._keywords);
            }
            else
            {
                std::string err = "Invalid format for NOMAD::Double attribute definition: ";
                err +=  att._name  + ". Infinity is set with NOMAD::INF or INF.";;
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
        }
        else if (att._type== "bool" )
        {
            // Triggers an exception if not a proper bool
            bool flag = NOMAD::stringToBool(  att._defaultValue );
            registerAttribute<bool>(att._name,
                                    flag,
                                    algoCompatibilityCheck,
                                    restartAttribute,
                                    uniqueEntry,
                                    att._shortInfo, att._helpInfo, att._keywords);
        }
        else if ( att._type== "string" || att._type== "std::string" )
        {
            registerAttribute<std::string>(att._name,
                                           att._defaultValue,
                                           algoCompatibilityCheck,
                                           restartAttribute,
                                           uniqueEntry,
                                           att._shortInfo, att._helpInfo, att._keywords);
        }
        else if ( att._type== "NOMAD::ArrayOfDouble" || att._type== "ArrayOfDouble"   )
        {
            if (  att._defaultValue != "-" &&  att._defaultValue != "N/A" )
            {
                std::string err = "Invalid attribute definition: ";
                err +=  att._name + " (" + att._defaultValue + "). An ArrayOfDouble must have an undefined default value \"-\" or \"N/A\" ";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            registerAttribute<NOMAD::ArrayOfDouble>(att._name,
                                                    NOMAD::ArrayOfDouble(),
                                                    algoCompatibilityCheck,
                                                    restartAttribute, uniqueEntry,
                                                    att._shortInfo, att._helpInfo, att._keywords);
        }
        else if ( att._type == "NOMAD::LHSearchType" || att._type == "LHSearchType" )
        {
            if (  att._defaultValue != "-" &&  att._defaultValue != "N/A" )
            {
                registerAttribute<NOMAD::LHSearchType>(att._name,
                                                       NOMAD::LHSearchType(att._defaultValue),
                                                       algoCompatibilityCheck,
                                                       restartAttribute, uniqueEntry,
                                                       att._shortInfo, att._helpInfo,
                                                       att._keywords);
            }
            else
            {
                registerAttribute<NOMAD::LHSearchType>(att._name,                                                                 NOMAD::LHSearchType(),
                                                       algoCompatibilityCheck,
                                                       restartAttribute, uniqueEntry,
                                                       att._shortInfo, att._helpInfo,
                                                       att._keywords);
            }

        }
        else if ( att._type== "NOMAD::Point" || att._type== "Point" )
        {
            // The default for this type of attribute must be undefined
            if (  att._defaultValue != "-" &&  att._defaultValue != "N/A" )
            {
                std::string err = "Invalid attribute definition: ";
                err +=  att._name + ". A Point must have an undefined default value \"-\" or \"N/A\" ";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            registerAttribute<NOMAD::Point>(att._name,
                                            NOMAD::Point(),
                                            algoCompatibilityCheck,
                                            restartAttribute,
                                            uniqueEntry,
                                            att._shortInfo, att._helpInfo, att._keywords);
        }
        // ARRAYOFPOINTS
        else if (att._type== "NOMAD::ArrayOfPoint" || att._type == "ArrayOfPoint")
        {
            if (  att._defaultValue != "-" &&  att._defaultValue != "N/A" )
            {
                std::string err = "Invalid attribute definition: ";
                err +=  att._name + " (" + att._defaultValue + "). A ArrayOfPoint must have an undefined default value \"-\" or \"N/A\" ";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            registerAttribute<NOMAD::ArrayOfPoint>(att._name,
                                                   NOMAD::ArrayOfPoint(),
                                                   algoCompatibilityCheck,
                                                   restartAttribute,
                                                   uniqueEntry,
                                                   att._shortInfo, att._helpInfo, att._keywords);
        }
        // ARRAYOFSTRING
        else if ( att._type== "NOMAD::ArrayOfString" || att._type=="ArrayOfString" )
        {
            registerAttribute( att._name,
                              NOMAD::ArrayOfString(att._defaultValue),
                              algoCompatibilityCheck,
                              restartAttribute,
                              uniqueEntry,
                              att._shortInfo, att._helpInfo, att._keywords);
        }
        // BBInputTypeList
        else if ( att._type== "NOMAD::BBInputTypeList" || att._type== "BBInputTypeList" )
        {
            registerAttribute( att._name,
                              NOMAD::stringToBBInputTypeList(att._defaultValue),
                              algoCompatibilityCheck, restartAttribute, uniqueEntry,
                              att._shortInfo, att._helpInfo, att._keywords);
        }
        // BBOutputTypeList
        else if ( att._type== "NOMAD::BBOutputTypeList" || att._type== "BBOutputTypeList" )
        {
            registerAttribute( att._name,
                              NOMAD::stringToBBOutputTypeList(att._defaultValue),
                              algoCompatibilityCheck, restartAttribute, uniqueEntry,
                              att._shortInfo, att._helpInfo, att._keywords);
        }
        // SgtelibModelFeasibilityType
        else if (   att._type== "NOMAD::SgtelibModelFeasibilityType"
                 || att._type== "SgtelibModelFeasibilityType")
        {
            registerAttribute( att._name,
                              NOMAD::stringToSgtelibModelFeasibilityType(att._defaultValue),
                              algoCompatibilityCheck, restartAttribute,
                              uniqueEntry, att._shortInfo , att._helpInfo, att._keywords );
        }
        // SgtelibModelFormulationType
        else if (   att._type== "NOMAD::SgtelibModelFormulationType"
                 || att._type== "SgtelibModelFormulationType")
        {
            registerAttribute( att._name,
                              NOMAD::stringToSgtelibModelFormulationType(att._defaultValue),
                              algoCompatibilityCheck , restartAttribute, uniqueEntry,
                              att._shortInfo , att._helpInfo, att._keywords );
        }
        // EvalType
        else if (   att._type== "NOMAD::EvalType"
                 || att._type== "EvalType")
        {
            registerAttribute( att._name,
                              NOMAD::stringToEvalType(att._defaultValue),
                              algoCompatibilityCheck, restartAttribute, uniqueEntry,
                              att._shortInfo , att._helpInfo, att._keywords );
        }
        // EvalSortType
        else if (   att._type== "NOMAD::EvalSortType"
                 || att._type== "EvalSortType")
        {
            registerAttribute( att._name,
                              NOMAD::stringToEvalSortType(att._defaultValue),
                              algoCompatibilityCheck, restartAttribute, uniqueEntry,
                              att._shortInfo , att._helpInfo, att._keywords );
        }
        // DirectionTypeList
        else if (   att._type== "NOMAD::DirectionTypeList"
                 || att._type== "DirectionTypeList")
        {
            NOMAD::DirectionTypeList dirTypeList { NOMAD::stringToDirectionType(att._defaultValue) };
            registerAttribute( att._name,
                              dirTypeList,
                              algoCompatibilityCheck, restartAttribute, uniqueEntry,
                              att._shortInfo , att._helpInfo, att._keywords );
        }
        // ListOfVariableGroup
        else if (att._type== "NOMAD::ListOfVariableGroup" || att._type == "ListOfVariableGroup")
        {
            if (  att._defaultValue != "-" &&  att._defaultValue != "N/A" )
            {
                std::string err = "Invalid attribute definition: ";
                err +=  att._name + " (" + att._defaultValue + "). A ListOfVariableGroup must have an undefined default value \"-\" or \"N/A\" ";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
            registerAttribute<NOMAD::ListOfVariableGroup>(att._name,
                                                   NOMAD::ListOfVariableGroup(),
                                                   algoCompatibilityCheck,
                                                   restartAttribute,
                                                   uniqueEntry,
                                                   att._shortInfo, att._helpInfo, att._keywords);
        }
        // Unrecognized type
        else
        {
            std::string err = "Unknown attribute type definition: ";
            err +=  att._name + " type: " + att._type ;
            throw NOMAD::Exception(__FILE__,__LINE__, err);
        }

    }
}

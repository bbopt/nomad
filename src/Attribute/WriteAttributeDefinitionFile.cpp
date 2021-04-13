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
// This executable converts the user-friendly txt files that describe the
// parameters, to computer-readable header files.
// When a parameter needs to be added or edited, only the txt files need
// to be modified.

#include <algorithm>    // for for_each
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

// Registered attribute definition names
const std::string attributeDefinitionNames[13] = { "deprecatedAttributesDefinition",
    "displayAttributesDefinition",
    "evalAttributesDefinition",
    "cacheAttributesDefinition",
    "evaluatorControlAttributesDefinition",
    "evaluatorControlGlobalAttributesDefinition",
    "pbAttributesDefinition",
    "runAttributesDefinition",
    "runAttributesDefinitionLH",
    "runAttributesDefinitionNM",
    "runAttributesDefinitionPSDSSD",
    "runAttributesDefinitionQuadModel",
    "runAttributesDefinitionSgtelibModel"
};

/// \brief Registered attribute flags and default value
/**
The order of the registered attribute flags is important for attribute definition in the parameter classes.
Default value for the registered flag is given here.
 */
const std::map<std::string, bool> attributeFlagsRegistered {
            {"RESTART_ATTRIBUTE", false},
            {"ALGO_COMPATIBILITY_CHECK", false},
            {"UNIQUE_ENTRY", true} };

/// \brief Toupper utility
void  toUpperCase(std::string& str)
{
    for_each(str.begin(), str.end(), [](char& in){ in = std::toupper(in); });
}

/// \brief Exception utility
class Exception : public std::exception
{

private:

    mutable std::string _what;  ///< Error message.
    size_t                 _line;  ///< Line number at which the exception is thrown.

public:

    Exception ( size_t line , const std::string & msg )
    : _what ( msg  ) ,
    _line ( line )  {}

    /// Destructor.
    virtual ~Exception ( void ) {}

    /// Access to the error message.
    /**
     \return A string with the error message.
     */
    const char * what ( void ) const noexcept { return _what.c_str(); }
    size_t getLineNumber ( void ) const noexcept { return _line; }
};


/// \brief Utility to interpret plural words
void duplicateParPlurals(std::string & blockOfStrings )
{
    // Pad with blanks
    blockOfStrings.insert(0," ");
    blockOfStrings.append(" ");

    size_t posBegin = 0;

    // Loop on ()
    while ( true )
    {
        size_t posPar1 = blockOfStrings.find( "(" , posBegin );
        if ( posPar1 != std::string::npos )
        {
            size_t posBlank = blockOfStrings.rfind( " " , posPar1 );
            size_t posPar2 = blockOfStrings.find( ")" ,posBegin );
            if ( posPar2 != std::string::npos && posBlank < posPar2 )
            {
                // Word to duplicate
                std::string word = blockOfStrings.substr(posBlank+1,posPar1-posBlank-1)+" ";

                // Remove parenthesis
                blockOfStrings.erase(posPar2,1);
                blockOfStrings.erase(posPar1,1);

                // Insert word without parenthesis
                blockOfStrings.insert(posBlank+1,word);

                posBegin = posBlank;
            }
            else
                break;
        }
        else
            break;

    }
}


/// \brief Utility to read a block of lines
std::string readBlockOfLines(std::ifstream & fin,
                             size_t & lineNumber,
                             const std::string & blockDelimiterOpen,
                             const std::string & blockDelimiterClose)

{
    std::string line;

    std::getline(fin, line);
    lineNumber++;


    if ( fin.fail() )
        throw ( Exception(lineNumber,"Cannot read file") );

    // The block must be provided between delimiters
    if (line.find( blockDelimiterOpen ) != 0)
    {
        std::string errMsg = " Cannot find opening delimiter ";
        errMsg = errMsg.append(blockDelimiterOpen) ;
        throw ( Exception(lineNumber,errMsg) );
    }
    std::string blockOfLines = line + " \\n ";

    // Take the content of line and put it in string
    // with end of line
    while ( line.empty() || line.find( blockDelimiterClose ) == std::string::npos )
    {
        std::getline(fin, line);
        lineNumber++;

        if ( fin.fail() )
        {
            std::string errMsg = " Reached end of file without finding end delimiter ";
            errMsg = errMsg.append(blockDelimiterClose) ;
            throw ( Exception(lineNumber,errMsg) );
        }

        blockOfLines += line + " \\n ";

    }

    // Remove block delimiters
    blockOfLines.erase(0,blockDelimiterOpen.size() );
    size_t posEnd = blockOfLines.find( blockDelimiterClose );
    blockOfLines.erase(posEnd);

    return blockOfLines;
}


int main(int argc, char *argv[])
{
    // Get directory where the input .txt files are.
    std::string inputDir, outputDir;
    if (argc > 1)
    {
        std::string s = argv[1];
        // Ensure directory ends with '/'.
        if (s[s.size()-1] != '/')
        {
            s += '/';
        }
        inputDir    = s;
        outputDir   = s;
    }

    std::ifstream fin;
    std::ofstream fout;
    std::string errMsg ;

    std::istringstream iss;
    std::ostringstream oss;

    std::map<std::string,bool> attributeFlags;

    for ( auto attDefName : attributeDefinitionNames )
    {
        std::string attDefFile      = inputDir + attDefName + ".txt";
        std::string attDefHeader    = outputDir + attDefName + ".hpp";

        // Clear input stream
        iss.str("");
        iss.clear();

        // Clear output stream
        oss.str("");
        oss.clear();
        std::string UpperAttDefName = attDefName ;
        for (auto & c: UpperAttDefName) c = toupper(c);

        oss << "//////////// THIS FILE MUST BE CREATED BY EXECUTING WriteAttributeDefinitionFile ////////////" << std::endl;
        oss << "//////////// DO NOT MODIFY THIS FILE MANUALLY ///////////////////////////////////////////////" << std::endl << std::endl;

        oss << "#ifndef __NOMAD400_"<< UpperAttDefName << "__" << std::endl;
        oss << "#define __NOMAD400_"<< UpperAttDefName << "__" << std::endl << std::endl;
        oss << "_definition = {" ;

        bool flagInFile;
        try
        {
            flagInFile=true;

            // Default error message for problem opening file
            errMsg = " Failed to open file " + attDefFile + ".";

            size_t lineNumber=0;
            fin.open ( attDefFile );

            if ( !fin.fail() )
            {
                errMsg.clear();
            }
            if ( !errMsg.empty() )
            {
                throw ( Exception( 0,errMsg ) );
            }

            std::string line;
            bool firstAttribute = true;

            while ( fin.good() && !fin.eof() )
            {

                // Read attribute name
                if ( firstAttribute )
                {
                    getline(fin, line);
                    lineNumber++;
                }

                // the previous getline may have reached end of file
                if ( fin.eof() )
                {
                    continue;
                }

                if ( fin.fail() )
                {
                    errMsg = " Failed to read file " + attDefFile + ".";
                    throw ( Exception(lineNumber,errMsg) );
                }

                if ( line.size() > 0 && line.find_first_of("#")== 0 )
                {
                    if (!firstAttribute)
                    {
                        getline(fin, line);
                        lineNumber++;
                    }
                    continue;
                }

                // After the first attribute we need comma separation between attributes
                if ( ! firstAttribute )
                    oss << "," << std::endl;
                else
                    oss << std::endl;

                std::string attributeName;
                iss.clear();
                iss.str(line);
                iss >> attributeName;

                if ( attributeName.size()==0)
                {
                    errMsg = " The first word after a comment must be an attribute name. \n An empty line is not valid.";
                    throw ( Exception(lineNumber,errMsg) );
                }
                oss << "{ \"" << attributeName << "\", ";

                // Read attribute type
                getline(fin, line);
                lineNumber++;

                if ( fin.fail() || line.empty() )
                {
                    errMsg = "Failed to read file " + attDefFile + ".";
                    throw ( Exception(lineNumber,errMsg) );
                }

                std::string attributeType;
                iss.clear();
                iss.str(line);

                iss >> attributeType ;

                if ( attributeType.size()==0 )
                {
                    errMsg = " The next word after an attribute name must be the attribute type. \n An empty line is not valid.";
                    throw ( Exception(lineNumber,errMsg) );
                }
                oss << " \"" << attributeType << "\", " ;

                // Read attribute default value (can be undefined)
                getline(fin, line);
                lineNumber++;

                if ( fin.fail() || line.empty() )
                {
                    errMsg = " Failed to read file " + attDefFile + ".";
                    throw ( Exception(lineNumber,errMsg) );
                }

                std::string attributeDefault= line;


                // Special case: attributeType is std::string or NOMAD::ArrayOfString, and attributeDefault
                // is set to a string meaning no default: set attributeDefault to empty string.
                if ("std::string" == attributeType || "NOMAD::ArrayOfString" == attributeType)
                {
                    if ("\"\"" == attributeDefault || "-" == attributeDefault || "N/A" == attributeDefault)
                    {
                        attributeDefault.clear();
                    }
                }

                oss << " \"" << attributeDefault << "\", ";

                // Read attribute short info
                // The block delimiter string must be "\(" and "\)" in the text
                // file but the string passed to the function must be preceeded
                // by the escape symbol "\"
                // We use raw string literals (R) to avoid escaping of any character. Anything between the delimiters becomes part of the string.
                std::string attributeShortInfo = readBlockOfLines( fin, lineNumber ,R"f(\()f" , R"f(\))f" );
                oss << " \"" << attributeShortInfo << "\", " ;


                // Read attribute help info
                // The block delimiter string must be "\(" and "\)" in the text
                // file but the string passed to the function must be preceeded
                // by the escape symbol "\"
                // We use raw string literals (R) to avoid escaping of any character. Anything between the delimiters becomes part of the string.
                std::string attributeHelpInfo = readBlockOfLines( fin , lineNumber , R"f(\()f", R"f(\))f");


                //Put attributeDefault value as a special line of HelpInfo
                bool flagEmpty = attributeDefault.empty() || attributeDefault=="\"\"";
                bool flagNoDefault = attributeDefault=="-";
                attributeHelpInfo += ". ";
                if (flagNoDefault)
                {
                    attributeHelpInfo += "No default value.";
                }
                else if (flagEmpty)
                {
                    attributeHelpInfo += "Default: Empty string.";
                }
                else
                {
                    attributeHelpInfo += "Default: " + attributeDefault;
                }
                attributeHelpInfo += "\\n\\n";

                oss << " \"" << attributeHelpInfo << "\", ";

                // Read keywords
                // The block delimiter string must be "\(" and "\)" in the text
                // file but the string passed to the function must be preceeded
                // by the escape symbol "\"
                // We use raw string literals (R) to avoid escaping of any character. Anything between the delimiters becomes part of the string.
                std::string attributeKeywords = readBlockOfLines( fin , lineNumber , R"f(\()f", R"f(\))f");
                duplicateParPlurals(attributeKeywords);
                oss << " \"" << attributeKeywords << "\" ";


                // Read optional flags for attribute
                getline(fin, line);
                lineNumber++;
                std::string attributeFlagName,attributeFlagValue;

                // Set the flag value to default
                attributeFlags = attributeFlagsRegistered ;

                while ( line.size() > 0 && line.find_first_of("#") != 0 )
                {
                    iss.clear();
                    iss.str(line);
                    iss >> attributeFlagName >> attributeFlagValue;

                    if ( attributeFlagName.size()==0 || attributeFlagValue.size()==0 )
                    {
                        errMsg = " An attribute flag name and value (bool) is expected (ex. RESTART_ATTRIBUTE yes). \n An empty line or a missing flag value (bool) is not valid.";
                        throw ( Exception(lineNumber,errMsg) );
                    }

                    // Toupper the attribute flag name
                    toUpperCase( attributeFlagName );

                    if ( attributeFlags.find(attributeFlagName) != attributeFlags.end() )
                    {
                        // The flag name is registered --> update the flag value

                        toUpperCase( attributeFlagValue );
                        if ( attributeFlagValue.compare("YES")==0 || attributeFlagValue.compare("TRUE")==0 )
                            attributeFlags[attributeFlagName] = true;
                        else if ( attributeFlagValue.compare("NO")==0 || attributeFlagValue.compare("FALSE")==0 )
                            attributeFlags[attributeFlagName] = false;
                        else
                        {
                            errMsg = " An attribute flag name and value (bool) is expected (ex. RESTART_ATTRIBUTE yes). \n The value must be YES/TRUE or NO/FALSE.";
                            throw ( Exception(lineNumber,errMsg) );
                        }
                    }
                    else
                    {
                        errMsg = " An attribute flag name and value (bool) is expected (ex. RESTART_ATTRIBUTE yes).\n  This flag name: " + attributeFlagName + " is not registered.";
                        throw ( Exception(lineNumber,errMsg) );
                    }

                    getline(fin, line);
                    lineNumber++;
                }

                // Write the flag in the registered order
                for ( auto attFlag : attributeFlags )
                {
                    oss << " , \"" << ((attFlag.second) ? "true":"false") << "\"";
                }
                oss << " }";

                getline(fin, line);
                lineNumber++;

                firstAttribute = false;

            }
            oss << " };";
            oss << std::endl << std::endl << "#endif" << std::endl;


            // input file is closed:
            fin.close();

            // write header
            // Default error message for problem opening file
            fout.open( attDefHeader );

            if ( fout.fail() )
            {
                errMsg = " Failed to open the header file for writing attribute definition.";
                throw ( Exception(0,errMsg) );
            }
            fout << oss.str() ;
            fout.close();

            flagInFile = false;

        }
        catch ( Exception & e)
        {
            fin.close();
            std::cerr << "ERROR: Problem handling file for attribute definition: " << ((flagInFile) ? attDefFile:attDefHeader) << " at line " <<  e.getLineNumber() << ". \n" << e.what() << std::endl << std::endl;
            return -1;
        }
    }
    return 0;
}


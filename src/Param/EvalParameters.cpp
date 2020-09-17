
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





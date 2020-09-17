
#include "../Param/CacheParameters.hpp"
#include "../Util/fileutils.hpp"

/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::CacheParameters::init()
{
    _typeName = "Cache";

    try
    {
        #include "../Attribute/cacheAttributesDefinition.hpp"
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
void NOMAD::CacheParameters::checkAndComply( std::shared_ptr<NOMAD::RunParameters> runParams )
{
    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    /*----------------------------*/
    /* Update cache file names    */
    /*----------------------------*/
    auto problemDir = runParams->getAttributeValue<std::string>("PROBLEM_DIR",false);
    std::string cacheFileName = getAttributeValueProtected<std::string>("CACHE_FILE",false);
    if (!cacheFileName.empty())
    {
        NOMAD::completeFileName(cacheFileName, problemDir);
        setAttributeValue("CACHE_FILE", cacheFileName);
    }

    // Hot restart needs a cache file
    bool hotRestartRead = runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES", false);
    bool hotRestartWrite = runParams->getAttributeValue<bool>("HOT_RESTART_WRITE_FILES", false);
    if (hotRestartRead || hotRestartWrite)
    {
        if (cacheFileName.empty())
        {
            cacheFileName = "cache.txt";
            std::cerr << "Warning: " << ((hotRestartWrite) ? "HOT_RESTART_WRITE_FILES" : "HOT_RESTART_READ_FILES") << " is set. CACHE_FILE set to \"" << cacheFileName << "\"" << std::endl;

            NOMAD::completeFileName(cacheFileName, problemDir);
            setAttributeValue("CACHE_FILE", cacheFileName);
        }
    }

    _toBeChecked = false;

}
// End checkAndComply()



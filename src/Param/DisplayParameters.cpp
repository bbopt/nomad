
#include "../Param/DisplayParameters.hpp"
#include "../Util/fileutils.hpp"

/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::DisplayParameters::init()
{
    _typeName = "Display";

    try
    {
        #include "../Attribute/displayAttributesDefinition.hpp"
        registerAttributes( _definition );

        // Note: we cannot call checkAndComply() here, the default values
        // are not valid, for instance DIMENSION, X0, etc.

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
void NOMAD::DisplayParameters::checkAndComply(
                    const std::shared_ptr<NOMAD::RunParameters> &runParams,
                    const std::shared_ptr<NOMAD::PbParameters> &pbParams)
{

    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    // Pb params must be checked before accessing its value
    size_t n = pbParams->getAttributeValue<size_t>("DIMENSION");
    if (n == 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DIMENSION must be positive" );
    }


    // SOL_FORMAT (internal)
    auto solFormat = getAttributeValueProtected<NOMAD::ArrayOfDouble>("SOL_FORMAT",false);
    if ( !solFormat.isDefined() )
    {
        solFormat.reset(n, NOMAD::DISPLAY_PRECISION_STD);
        setAttributeValue("SOL_FORMAT", solFormat);
    }

    if ( ! pbParams->isAttributeDefaultValue<NOMAD::ArrayOfDouble>("GRANULARITY") )
    {

        // Update SOL_FORMAT.
        // We want to remember the number of decimals in
        // GRANULARITY arguments, to be used later for formatting.
        //
        // Default is DISPLAY_PRECISION_STD.

        auto newSolFormat = setFormatFromGranularity( pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("GRANULARITY") );
        setAttributeValue("SOL_FORMAT", newSolFormat);
    }

    // The default value is empty: set a basic display stats: BBE OBJ
    auto displayStats = getAttributeValueProtected<NOMAD::ArrayOfString>("DISPLAY_STATS",false);
    if ( displayStats.size() == 0 )
    {
        NOMAD::ArrayOfString aos("BBE OBJ");
        setAttributeValue("DISPLAY_STATS", aos );
    }


    /*------------------------------------------------------*/
    /* Stats file                                           */
    /*------------------------------------------------------*/

    auto statsFileParam = getAttributeValueProtected<NOMAD::ArrayOfString>("STATS_FILE",false) ;
    std::string statsFileName;
    if (statsFileParam.size() > 0)
    {
        statsFileName = statsFileParam[0];
        if (statsFileParam.size() == 1)
        {
            // Default stats: BBE OBJ.
            statsFileParam.add("BBE");
            statsFileParam.add("OBJ");
        }
    }



    // Update stats file name
    auto addSeedToFileNames = runParams->getAttributeValue<bool>("ADD_SEED_TO_FILE_NAMES");
    auto problemDir = runParams->getAttributeValue<std::string>("PROBLEM_DIR");
    if (!statsFileName.empty())
    {
        auto seed = runParams->getAttributeValue<int>("SEED");
        NOMAD::completeFileName(statsFileName, problemDir, addSeedToFileNames, seed);
        statsFileParam.replace(0, statsFileName);
        setAttributeValue("STATS_FILE", statsFileParam);
    }

    _toBeChecked = false;

}
// End checkAndComply()




NOMAD::ArrayOfDouble NOMAD::DisplayParameters::setFormatFromGranularity( const NOMAD::ArrayOfDouble & aod )
{
    size_t n = aod.size();
    NOMAD::ArrayOfDouble solFormat(n, NOMAD::DISPLAY_PRECISION_STD);

    // Use GRANULARITY as an ArrayOfDouble.
    size_t nbDecimals;
    for ( size_t i=0 ; i < n ; i++ )
    {
        if ( aod[i] > 0 )
        {
            nbDecimals = aod[i].nbDecimals( );
            solFormat.set(i, nbDecimals);
        }
    }
    return solFormat;

}





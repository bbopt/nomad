/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#include "../Eval/Eval.hpp"
#include "../Math/RNG.hpp"
#include "../Param/RunParameters.hpp"

#include "../nomad_version.hpp"


// Static members initialization
bool NOMAD::RunParameters::_warningUnknownParamShown = false;

/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::RunParameters::init()
{
    _typeName = "Run";

    try
    {
        #include "../Attribute/runAttributesDefinition.hpp"
        registerAttributes( _definition );
    
        // Registered attributes using defined keywords (not in preprocessed special header file)
        registerAttribute<NOMAD::Double>("EPSILON", NOMAD::DEFAULT_EPSILON, false,
            " NOMAD precision for comparison of values ",
            " \n \n . NOMAD precision for comparison of values \n "," advanced precision double(s) value(s) number(s) ");
        registerAttribute<std::string>("UNDEF_STR", NOMAD::DEFAULT_UNDEF_STR, false,
            "String for undefined values ",
            " \n \n String for undefined values \n "," advanced undef(ined) value(s) string(s) ");
        registerAttribute<std::string>("NOMAD_VERSION", NOMAD_VERSION_NUMBER, true,
            " NOMAD version number (for runner) ",
            " \n \n . NOMAD version number (optional) \n  . If not compatible with current version will trigger exception \n " ," advanced nomad version(s) release(s) revision(s) " );
        registerAttribute<std::string>("INF_STR", NOMAD::DEFAULT_INF_STR, false,
            "String for infinite values", 
            " \n \n . String for infinite values \n "," advanced string(s) inf(inite) value(s) ");
        registerAttribute<std::string>("PROBLEM_DIR", std::string(".") + NOMAD::DIR_SEP , false,
            "Problem directory " , 
            "\n \n . Problem directory \n . To complete \n "," problem dir(ectory) folder(s) ");
        // Note: we cannot call checkAndComply() here, the default values
        // are not valid.
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
void NOMAD::RunParameters::checkAndComply(
        const std::shared_ptr<NOMAD::EvaluatorControlParameters>& evaluatorControlParams,
        const std::shared_ptr<NOMAD::PbParameters>& pbParams)
{
    std::string err;
    checkInfo();
    
    if (!toBeChecked())
    {
        // Early out
        return;
    }
    
    // check the non-interpreted parameters:
    //const std::shared_ptr<NOMAD::ParameterEntry> pe = getNonInterpretedParamEntry();
    std::vector<std::shared_ptr<NOMAD::ParameterEntry>> allNonInterp = getAllNonInterpretedParamEntries();
    if (allNonInterp.size() > 0)
    {
        err = "Unrecognized parameters in file " + allNonInterp[0]->getParamFile() + ":\n";
        for (auto pe : allNonInterp)
        {
            err += "line " + std::to_string(pe->getLine());
            err += ": Unrecognized parameter: " + pe->getName() + "\n";
        }
        if (getAttributeValueProtected<bool>("REJECT_UNKNOWN_PARAMETERS", false))
        {
            throw NOMAD::Exception(__FILE__,__LINE__, err);
        }
        else
        {
            if (!_warningUnknownParamShown)
            {
                std::cerr << "Warning: " << err << "Ignoring unknown parameters." << std::endl;
            }
        }
    }
    
    auto problemDir = getAttributeValueProtected<std::string>("PROBLEM_DIR", false);
    
    auto seed = getAttributeValueProtected<int>("SEED" ,false);
    if ( seed < 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: SEED must be non-negative" );
    }
    
    // CT todo modify to manage retro-compatiblity
    auto version_number = getAttributeValueProtected<std::string>("NOMAD_VERSION", false);
    if ( version_number != NOMAD_VERSION_NUMBER )
    {
       throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: VERSION_NUMBER is not compatible" );
    }
    
    auto anisotropyFactor = getAttributeValueProtected<NOMAD::Double>("ANISOTROPY_FACTOR", false);
    if ( anisotropyFactor <= 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: ANISOTROPY_FACTOR must be positive" );
    }
    
    setStaticParameters();
    
    /*---------------------------*/
    /* Sgtelib Search parameters */
    /*---------------------------*/
    const size_t bigDim = 50;
    size_t n = pbParams->getAttributeValue<size_t>("DIMENSION");
#ifdef USE_SGTELIB
    bool showDisableWarn = true;
#endif
    // If dimension is too large, disable models.
    if (n >= bigDim)
    {
        setAttributeValue("DISABLE", std::string("MODELS"));
        std::cerr << "Warning: Dimension is higher than " << bigDim << ". Models are disabled." << std::endl;
#ifdef USE_SGTELIB
        showDisableWarn = false;
#endif
    }
#ifdef USE_SGTELIB
    // If models are disabled, set SGTELIB_SEARCH to false.
    auto disableModels = getAttributeValueProtected<std::string>("DISABLE", false);
    if ("MODELS" == disableModels)
    {
        if (getAttributeValueProtected<bool>("SGTELIB_SEARCH", false))
        {
            if (showDisableWarn)
            {
                std::cerr << "Warning: Models are disabled. SGTELIB_SEARCH set to false." << std::endl;
            }
            setAttributeValue("SGTELIB_SEARCH", false);
        }
    }
#endif

    // Set default value, if the parameter is not set.
    // Default value: TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM KERNEL_COEF OPTIM RIDGE 0 METRIC AOECV
    auto sgtelibModelDefinition = getAttributeValueProtected<NOMAD::ArrayOfString>("SGTELIB_MODEL_DEFINITION",false);
    if (sgtelibModelDefinition.size() == 0)
    {
        //NOMAD::ArrayOfString aos("TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM KERNEL_COEF OPTIM RIDGE 0 METRIC AOECV");
        NOMAD::ArrayOfString aos("TYPE PRS DEGREE 2");
        setAttributeValue("SGTELIB_MODEL_DEFINITION", aos );
    }

    auto sgtelibModelTrials = getAttributeValueProtected<size_t>("SGTELIB_MODEL_TRIALS", false);
    if (0 == sgtelibModelTrials)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter SGTELIB_MODEL_TRIALS must be positive");
    }

    /*--------------------------------*/
    /* Parallelism related parameters */
    /*--------------------------------*/
    // Ensure we can get value for BB_MAX_BLOCK_SIZE without throwing an exception.
    if (evaluatorControlParams->toBeChecked())
    {
        evaluatorControlParams->checkAndComply();
    }
    auto blockSize = evaluatorControlParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
    if (0 == blockSize)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter BB_MAX_BLOCK_SIZE must be positive");
    }
#ifdef _OPENMP
    else if (blockSize > 1)
    {
        if (getAttributeValueProtected<int>("NB_THREADS_OPENMP", false) != 1)
        {
            std::cerr << "Warning: Parallelism management: BB_MAX_BLOCK_SIZE is ";
            std::cerr << "larger than 1 (value is " << blockSize << "). ";
            std::cerr << "Setting parameter NB_THREADS_OPENMP to 1." << std::endl;
        }
        setAttributeValue("NB_THREADS_OPENMP", 1);
    }
#endif

#ifndef USE_SGTELIB
    // Look for SgtelibModel parameters that are set but cannot be used
    if (getAttributeValueProtected<bool>("SGTELIB_MODEL_EVAL", false))
    {
        err = "Sgtelib Model sampling cannot be used. Either set parameter ";
        err += "SGTELIB_MODEL_EVAL to false, or recompile NOMAD using option USE_SGTELIB=1.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
    if (getAttributeValueProtected<bool>("SGTELIB_SEARCH", false))
    {
        err = "Warning: Parameter SGTELIB_SEARCH is set to true, but ";
        err += "Sgtelib Model sampling cannot be used. To be able to use Sgtelib Model ";
        err += "search method, NOMAD must be recompiled using option USE_SGTELIB=1.";
        std::cerr << err << std::endl;
    }
#endif

    // Algorithm parameters: use an algorithm other than MADS.
    // They are mutually-exclusive.
    bool useAlgoLH = (getAttributeValueProtected<size_t>("LH_EVAL", false) > 0);
    bool useAlgoNM = getAttributeValueProtected<bool>("NM_OPTIMIZATION", false);
    bool useAlgoSgtelibModel = getAttributeValueProtected<bool>("SGTELIB_MODEL_EVAL", false);
    int totalAlgoSet = (int)useAlgoLH + (int)useAlgoNM + (int)useAlgoSgtelibModel;
    if (totalAlgoSet >= 2)
    {
        err = "Multiple parameters for algorithms are set. ";
        err += "These parameters are mutually exclusive:";
        if (useAlgoLH)
        {
            err += " LH_EVAL";
        }
        if (useAlgoNM)
        {
            err += " NM_OPTIMIZATION";
        }
        if (useAlgoSgtelibModel)
        {
            err += " SGTELIB_MODEL_EVAL";
        }
        err += ". Please review parameters settings and choose only one algorithm.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
    
    /*------------------------*/
    /* Hot restart parameters */
    /* Update file names      */
    /*------------------------*/
    std::string hotRestartFileName  = getAttributeValueProtected<std::string>("HOT_RESTART_FILE",false);
    if (!hotRestartFileName.empty())
    {
        NOMAD::completeFileName(hotRestartFileName, problemDir);
        setAttributeValue("HOT_RESTART_FILE", hotRestartFileName);
    }
    

    auto hMax = getAttributeValueProtected<NOMAD::Double>("H_MAX_0", false);
    if (hMax <= 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: H_MAX_0 must be positive");
    }

    _warningUnknownParamShown = true;

    _toBeChecked = false;
}
// End checkAndComply()


void NOMAD::RunParameters::setStaticParameters()
{
    // Sub-method of checkAndComply() to set static variables of some classes.
    NOMAD::RNG::setSeed ( getAttributeValueProtected<int>("SEED",false) );
    NOMAD::Double::setEpsilon ( getAttributeValueProtected<NOMAD::Double>("EPSILON",false).todouble() );
    NOMAD::Double::setUndefStr ( getAttributeValueProtected<std::string>("UNDEF_STR",false) );
    NOMAD::Double::setInfStr ( getAttributeValueProtected<std::string>("INF_STR",false) );
    
    // Reset parameter values from these static values, to ensure coherence.
    // This is bad because we have twice the same value for some parameters.
    setAttributeValue ( "SEED", NOMAD::RNG::getSeed() );
    setAttributeValue ( "EPSILON", NOMAD::Double(NOMAD::Double::getEpsilon()) );
    setAttributeValue ( "UNDEF_STR", NOMAD::Double::getUndefStr() );
    setAttributeValue ( "INF_STR", NOMAD::Double::getInfStr() );
}


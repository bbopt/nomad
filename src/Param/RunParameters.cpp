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

#include "../Math/RNG.hpp"  // for setSeed()
#include "../Param/RunParameters.hpp"
#include "../Type/DirectionType.hpp"
#include "../Util/fileutils.hpp"

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

        #include "../Attribute/runAttributesDefinitionLH.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionNM.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionPSDSSD.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionQuadModel.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionSgtelibModel.hpp"
        registerAttributes( _definition );
        
        #include "../Attribute/runAttributesDefinitionVNS.hpp"
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
void NOMAD::RunParameters::checkAndComply(
        const std::shared_ptr<NOMAD::EvaluatorControlGlobalParameters>& evaluatorControlGlobalParams,
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

    auto version_number = getAttributeValueProtected<std::string>("NOMAD_VERSION", false);
    if ( version_number != NOMAD_VERSION_NUMBER )
    {
        err = "Parameters check: NOMAD_VERSION is not compatible with registered version " + std::string(NOMAD_VERSION_NUMBER) + "\n";
       throw NOMAD::Exception(__FILE__,__LINE__, err );
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
    // If dimension is too large, disable models.
    if (n >= bigDim)
    {
        if (   getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH", false)
            || getAttributeValueProtected<bool>("SGTELIB_MODEL_SEARCH", false))
        {
            setAttributeValue("QUAD_MODEL_SEARCH", false);
            setAttributeValue("SGTELIB_MODEL_SEARCH", false);
            std::cerr << "Warning: Dimension " << n << " is greater than (or equal to) " << bigDim << ". Models are disabled." << std::endl;
        }
    }

    // Set default value, if the parameter is not set.
    // Default value: TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM KERNEL_COEF OPTIM RIDGE 0 METRIC AOECV
    auto sgtelibModelDefinition = getAttributeValueProtected<NOMAD::ArrayOfString>("SGTELIB_MODEL_DEFINITION",false);
    if (sgtelibModelDefinition.size() == 0)
    {
        //NOMAD::ArrayOfString aos("TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM KERNEL_COEF OPTIM RIDGE 0 METRIC AOECV");
        NOMAD::ArrayOfString aos("TYPE PRS DEGREE 2");
        setAttributeValue("SGTELIB_MODEL_DEFINITION", aos );
    }

    auto sgtelibModelSearchTrials = getAttributeValueProtected<size_t>("SGTELIB_MODEL_SEARCH_TRIALS", false);
    if (0 == sgtelibModelSearchTrials)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter SGTELIB_MODEL_SEARCH_TRIALS must be positive");
    }

    // Misc.
    auto frameCenterUseCache = getAttributeValueProtected<bool>("FRAME_CENTER_USE_CACHE", false);
    if (!frameCenterUseCache)
    {
        // Void parameter MAX_ITERATION_PER_MEGAITERATION
        setAttributeValue("MAX_ITERATION_PER_MEGAITERATION", INF_SIZE_T);
    }


    /*--------------------------------*/
    /* Parallelism related parameters */
    /*--------------------------------*/
    // Ensure we can get value for BB_MAX_BLOCK_SIZE without throwing an exception.
    if (nullptr != evaluatorControlGlobalParams)
    {
        if (evaluatorControlGlobalParams->toBeChecked())
        {
            evaluatorControlGlobalParams->checkAndComply();
        }
        auto bbBlockSize = evaluatorControlGlobalParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
        if (0 == bbBlockSize)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Parameter BB_MAX_BLOCK_SIZE must be positive");
        }
    }

#ifndef USE_SGTELIB
    // Look for SgtelibModel parameters that are set but cannot be used
    if (getAttributeValueProtected<bool>("SGTELIB_MODEL_EVAL", false))
    {
        err = "Sgtelib Model sampling cannot be used. Either set parameter ";
        err += "SGTELIB_MODEL_EVAL to false, or recompile NOMAD using option USE_SGTELIB=1.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
    if (getAttributeValueProtected<bool>("SGTELIB_MODEL_SEARCH", false))
    {
        err = "Warning: Parameter SGTELIB_MODEL_SEARCH is set to true, but ";
        err += "Sgtelib Model sampling cannot be used. To be able to use Sgtelib Model ";
        err += "search method, NOMAD must be recompiled using option USE_SGTELIB=1.";
        std::cerr << err << std::endl;
    }
    if (getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH", false))
    {
        err = "Warning: Parameter QUAD_MODEL_SEARCH is set to true, but ";
        err += "Quad Model sampling cannot be used. To be able to use Quad Model Search ";
        err += "search method, NOMAD must be recompiled using option USE_SGTELIB=1.";
        std::cerr << err << std::endl;
    }
#endif

    // Update secondary poll direction based on primary poll direction
    // If DIRECTION_TYPE contains ORTHO, do nothing.
    // Else, set DIRECTION_TYPE_SECONDARY_POLL to SINGLE.
    // This is not exactly the behavior of NOMAD 3, but it is close enough.
    bool orthoInDirTypes = false;
    auto primaryDirTypes = getAttributeValueProtected<NOMAD::DirectionTypeList>("DIRECTION_TYPE", false);
    for (auto primaryDirType : primaryDirTypes)
    {
        if (   NOMAD::DirectionType::ORTHO_2N == primaryDirType
            || NOMAD::DirectionType::ORTHO_NP1_NEG == primaryDirType
            || NOMAD::DirectionType::ORTHO_NP1_QUAD == primaryDirType)
        {
            orthoInDirTypes = true;
            break;
        }
    }
    if (!orthoInDirTypes)
    {
        std::vector<NOMAD::DirectionType> dirTypes;
        dirTypes.push_back(NOMAD::DirectionType::SINGLE);
        setAttributeValue("DIRECTION_TYPE_SECONDARY_POLL", dirTypes);
    }

    // Precisions on MEGA_SEARCH_POLL
    if (getAttributeValueProtected<bool>("MEGA_SEARCH_POLL", false))
    {
        // MEGA_SEARCH_POLL does not support ORTHO_NP1_NEG and ORTHO_NP1_QUAD.
        for (auto dirType : primaryDirTypes)
        {
            if (   NOMAD::DirectionType::ORTHO_NP1_NEG == dirType
                || NOMAD::DirectionType::ORTHO_NP1_QUAD == dirType)
            {
                err = "Parameters check: Direction type " + NOMAD::directionTypeToString(dirType) + " is not supported with MEGA_SEARCH_POLL";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
        }
        for (auto dirType : getAttributeValueProtected<NOMAD::DirectionTypeList>("DIRECTION_TYPE_SECONDARY_POLL", false))
        {
            if (   NOMAD::DirectionType::ORTHO_NP1_NEG == dirType
                || NOMAD::DirectionType::ORTHO_NP1_QUAD == dirType)
            {
                err = "Parameters check: Direction type " + NOMAD::directionTypeToString(dirType) + " is not supported with MEGA_SEARCH_POLL";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
        }
    }

    // PSD-Mads and SSD-Mads parameters
    bool useAlgoPSDMads = getAttributeValueProtected<bool>("PSD_MADS_OPTIMIZATION", false);

#ifndef _OPENMP
    if (useAlgoPSDMads)
    {
        err = "Error: PSD_MADS_OPTIMIZATION can only be used when OpenMP is available. If that is not the case, use SSD_MADS_OPTIMIZATION.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
#endif // _OPENMP

    bool useAlgoSSDMads = getAttributeValueProtected<bool>("SSD_MADS_OPTIMIZATION", false);
    if (useAlgoPSDMads || useAlgoSSDMads)
    {
        std::string nbVarParamName = (useAlgoPSDMads ? "PSD_MADS_NB_VAR_IN_SUBPROBLEM" : "SSD_MADS_NB_VAR_IN_SUBPROBLEM");
        const size_t nbVariablesInSubproblem = getAttributeValueProtected<size_t>(nbVarParamName, false);
        if (0 == nbVariablesInSubproblem || nbVariablesInSubproblem > n)
        {
            err = "Parameter " + nbVarParamName + " must be between 1 and " + NOMAD::itos(n);
            err += ". Value provided: " + NOMAD::itos(nbVariablesInSubproblem);
            throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
        }

        size_t nbMadsSubproblem = getAttributeValueProtected<size_t>((useAlgoPSDMads ? "PSD_MADS_NB_SUBPROBLEM" : "SSD_MADS_NB_SUBPROBLEM"), false);
        // Distribute all the variables between subproblems of reduced dimension.
        bool nbMadsSubproblemSetByUser = true;
        if (nbMadsSubproblem == INF_SIZE_T)
        {
            nbMadsSubproblem = (size_t)std::round(n/nbVariablesInSubproblem)+1; // Add an additional mads for the pollster
            nbMadsSubproblemSetByUser = false;
        }
        if (useAlgoPSDMads)
        {
            // Cannot have more subproblems than the number of threads
            size_t nbThreads = (size_t)getAttributeValueProtected<int>("NB_THREADS_OPENMP", false);
            if (nbMadsSubproblem > nbThreads)
            {
                nbMadsSubproblem = nbThreads;
                if (nbMadsSubproblemSetByUser)
                {
                    // Warn the user
                    std::cerr << "Warning: parameter PSD_MADS_NB_SUBPROBLEM reset to number of available threads (" << nbThreads << ")" <<  std::endl;
                }
            }
        }
        setAttributeValue((useAlgoPSDMads ? "PSD_MADS_NB_SUBPROBLEM" : "SSD_MADS_NB_SUBPROBLEM"), nbMadsSubproblem);

        // Check parameter for coverage
        if (useAlgoPSDMads)
        {
            std::string covParamName = "PSD_MADS_SUBPROBLEM_PERCENT_COVERAGE";
            auto coverage = getAttributeValueProtected<NOMAD::Double>(covParamName, false);
            if (coverage < 0.0 || coverage > 100.0)
            {
                err = "Parameter " + covParamName + " must be between 0.0 and 100.0";
                throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
            }
        }
    }

    // Algorithm parameters: use an algorithm other than MADS.
    // They are mutually-exclusive.
    bool useAlgoLH = (getAttributeValueProtected<size_t>("LH_EVAL", false) > 0);
    bool useAlgoNM = getAttributeValueProtected<bool>("NM_OPTIMIZATION", false);
    bool useAlgoQuadOpt = getAttributeValueProtected<bool>("QUAD_MODEL_OPTIMIZATION", false);
    bool useAlgoSgtelibModel = getAttributeValueProtected<bool>("SGTELIB_MODEL_EVAL", false);
    int totalAlgoSet = (int)useAlgoLH + (int)useAlgoNM + (int)useAlgoQuadOpt
                       + (int)useAlgoPSDMads + (int)useAlgoSgtelibModel + (int)useAlgoSSDMads;

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
        if (useAlgoPSDMads)
        {
            err += " PSD_MADS_OPTIMIZATION";
        }
        if (useAlgoQuadOpt)
        {
            err += " QUAD_MODEL_OPTIMIZATION";
        }
        if (useAlgoSgtelibModel)
        {
            err += " SGTELIB_MODEL_EVAL";
        }
        if (useAlgoSSDMads)
        {
            err += " SSD_MADS_OPTIMIZATION";
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

    /*------------------------*/
    /* Speculative search     */
    /*  parameters            */
    /*------------------------*/
    auto baseFactor = getAttributeValueProtected<NOMAD::Double>("SPECULATIVE_SEARCH_BASE_FACTOR", false);
    if (baseFactor <= 1)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: SPECULATIVE_SEARCH_BASE_FACTOR must be strictly greater than 1");
    }

    /*------------------------*/
    /* Quad model search      */
    /*  parameters            */
    /*------------------------*/
    auto reductionFactor = getAttributeValueProtected<NOMAD::Double>("QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR", false);
    if (reductionFactor <= 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR must be strictly greater than 0");
    }

    _warningUnknownParamShown = true;

    _toBeChecked = false;
}
// End checkAndComply()


void NOMAD::RunParameters::setStaticParameters()
{
    // Sub-method of checkAndComply() to set static variables of some classes.
    int currentRNGSeed = NOMAD::RNG::getSeed();
    int seedToSet = getAttributeValueProtected<int>("SEED",false);
    // If the seed has changed we call setSeed which reset the private RNG seed to its default and runs the rand() function to set the seed (seed and private seed are different).
    // If the seed is the same as before we do nothing.
    if (currentRNGSeed != seedToSet)
    {
        NOMAD::RNG::setSeed ( seedToSet );
    }
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

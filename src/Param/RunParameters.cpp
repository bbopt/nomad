/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
#include "../Type/BBOutputType.hpp"
#include "../Util/fileutils.hpp"

#include <thread>

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

        #include "../Attribute/runAttributesDefinitionDMulti.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionIBEX.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionLH.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionCS.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionNM.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionPSDSSD.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionCOOP.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionQPSolver.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionQuadModel.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionSgtelibModel.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionVNS.hpp"
        registerAttributes( _definition );

        #include "../Attribute/runAttributesDefinitionDisco.hpp"
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
    std::vector<std::shared_ptr<NOMAD::ParameterEntry>> allNonInterp = getAllNonInterpretedParamEntries();
    if (!allNonInterp.empty())
    {
        err = "Unrecognized parameters in file " + allNonInterp[0]->getParamFile() + ":\n";
        for (const auto& pe : allNonInterp)
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
                std::cout << "Warning: " << err << "Ignoring unknown parameters." << std::endl;
            }
        }
    }

    auto problemDir = getAttributeValueProtected<std::string>("PROBLEM_DIR", false);

    auto seed = getAttributeValueProtected<int>("SEED" ,false);
    if ( seed < 0 && seed != -1)
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: SEED must be non-negative (-1 is ok)" );
    }

    auto version_number = getAttributeValueProtected<std::string>("NOMAD_VERSION", false);
    if ( version_number != NOMAD_VERSION_NUMBER )
    {
        err = "Parameters check: NOMAD_VERSION is not compatible with registered version " + std::string(NOMAD_VERSION_NUMBER) + "\n";
       throw NOMAD::Exception(__FILE__,__LINE__, err );
    }

    auto anisotropyFactor = getAttributeValueProtected<NOMAD::Double>("ANISOTROPY_FACTOR", false);
    if ( anisotropyFactor < 0)
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
            || getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH_SIMPLE_MADS", false)
            || getAttributeValueProtected<bool>("SGTELIB_MODEL_SEARCH", false))
        {
            setAttributeValue("QUAD_MODEL_SEARCH", false);
            setAttributeValue("QUAD_MODEL_SEARCH_SIMPLE_MADS", false);
            setAttributeValue("SGTELIB_MODEL_SEARCH", false);
            std::cout << "Warning: Dimension " << n << " is greater than (or equal to) " << bigDim << ". Models are disabled." << std::endl;
        }
    }

    // Set default value, if the parameter is not set.
    // Default value: TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM KERNEL_COEF OPTIM RIDGE 0 METRIC AOECV
    auto sgtelibModelDefinition = getAttributeValueProtected<NOMAD::ArrayOfString>("SGTELIB_MODEL_DEFINITION",false);
    if (sgtelibModelDefinition.empty())
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
        std::cout << err << std::endl;
    }
    if (getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH", false) || getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH_SIMPLE_MADS", false))
    {
        err = "Warning: Parameter QUAD_MODEL_SEARCH or QUAD_MODEL_SEARCH_SIMPLE_MADS is set to true, but ";
        err += "Quad Model sampling cannot be used. To be able to use Quad Model Search ";
        err += "search method, NOMAD must be recompiled using option USE_SGTELIB=1.";
        std::cout << err << std::endl;
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

    // Test for extra trial points. Only valid for a unique Ortho 2n direction type.
    auto extraTrialPointsAddUp = getAttributeValueProtected<size_t>("TRIAL_POINT_MAX_ADD_UP",false);
    if ( extraTrialPointsAddUp > 0 && ( primaryDirTypes.size() > 1 || primaryDirTypes[0] != NOMAD::DirectionType::ORTHO_2N ))
    {
        err = "TRIAL_POINT_MAX_ADD_UP can only be used with a single ORTHO 2N direction type";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }

    // Test for CS
    bool useAlgoCS = getAttributeValueProtected<bool>("CS_OPTIMIZATION", false) ;
    if (primaryDirTypes.size() > 1 || !useAlgoCS)
    {
        if (std::count(primaryDirTypes.begin(), primaryDirTypes.end(),NOMAD::DirectionType::CS) >= 1)
        {
            err = "CS Direction_Type can only be used by CS Optimization. It cannot be combined with any other direction type";
            throw NOMAD::Exception(__FILE__,__LINE__, err);
        }
    }

    // Test for DMultiMads
    bool useAlgoDMultiMads = getAttributeValueProtected<bool>("DMULTIMADS_OPTIMIZATION", false) ;
    if (useAlgoDMultiMads)
    {
        for (auto primaryDirType : primaryDirTypes)
        {
            if (  NOMAD::DirectionType::ORTHO_NP1_QUAD == primaryDirType)
            {
                err = "Direction_Type for DMultiMads Optimization do not support ORTHO_NP1_QUAD. To deactivate, set DIRECTION_TYPE ORTHO 2N or DIRECTION_TYPE ORTHO N+1 NEG.";
                throw NOMAD::Exception(__FILE__,__LINE__, err);
            }
        }
        if (getAttributeValueProtected<bool>("MEGA_SEARCH_POLL", false))
        {
            err = "DMultiMads does not support the MEGA_SEARCH_POLL option. To deactivate, set MEGA_SEARCH_POLL to FALSE";
            throw NOMAD::Exception(__FILE__,__LINE__, err);
        }
    }

    // Precisions on MEGA_SEARCH_POLL
    if (getAttributeValueProtected<bool>("MEGA_SEARCH_POLL", false))
    {
        for (auto dirType : getAttributeValueProtected<NOMAD::DirectionTypeList>("DIRECTION_TYPE", false))
        {
            if (NOMAD::DirectionType::USER_FREE_POLL == dirType)
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

    // Coordinate Search algorithm
    if (useAlgoCS)
    {
        if (primaryDirTypes.size() > 1)
        {
            err = "CS Optimization can only use one direction type CS. It cannot be combined with any other direction type";
            throw NOMAD::Exception(__FILE__,__LINE__, err);
        }

        // Get the default direction type
        auto defaultPrimaryDirTypes = getAttributeValueProtected<NOMAD::DirectionTypeList>("DIRECTION_TYPE", false, true);
        if (primaryDirTypes[0] != NOMAD::DirectionType::CS)
        {
            if ( primaryDirTypes[0] != defaultPrimaryDirTypes[0] )
            {
                // Warn the user if he has changed the direction type. Otherwise, change silently to CS.
                std::cout << "Warning: parameter DIRECTION_TYPE reset to CS because CS_OPTIMIZATION is enabled" <<  std::endl;
            }
            setAttributeValue("DIRECTION_TYPE", std::vector<NOMAD::DirectionType> {NOMAD::DirectionType::CS});
        }
    }

    // DiscoMads algorithm

    bool useAlgoDiscoMads = getAttributeValueProtected<bool>("DISCO_MADS_OPTIMIZATION", false);
    if (useAlgoDiscoMads)
    {
        // If DiscoMads used to reveal hidden constraints...
        bool hiddenConstraints = getAttributeValueProtected<bool>("DISCO_MADS_HID_CONST", false);
        if (hiddenConstraints)
        {
            // In this case, detection radius and limit rate are useless, detection radius should be put to 0
            std::cout << "Warning: DiscoMads is used to reveal hidden constraints, so the detection radius and limit rate will not be considered. Detection radius is forced to 0." <<  std::endl;
            setAttributeValue<NOMAD::Double>("DISCO_MADS_DETECTION_RADIUS",0.0);

            // Check high value chosen for return of OBJ and PB constraints of failed evaluations
            const NOMAD::Double hiddConstOutputValue = getAttributeValueProtected<NOMAD::Double>("DISCO_MADS_HID_CONST_OUTPUT_VALUE", false);

                 // Ensure that points leading to failed eval are not used to build models
            if(hiddConstOutputValue<MODEL_MAX_OUTPUT){
                throw NOMAD::Exception(__FILE__,__LINE__,"The high value ("+hiddConstOutputValue.tostring()+") attributed to objective function and PB constraints for failed evaluations should be more than MODEL_MAX_OUTPUT ("+std::to_string(MODEL_MAX_OUTPUT)+") to ensure these points are not used to construct models.");
            }

            if(hiddConstOutputValue>=NOMAD::INF){
                throw NOMAD::Exception(__FILE__,__LINE__,"The high value ("+hiddConstOutputValue.tostring()+") attributed to objective function and PB constraints for failed evaluations should be less than NOMAD::INF. ");
            }
        }

        // If DiscoMads used to reveal discontinuities...
        auto detectionRadius = getAttributeValueProtected<NOMAD::Double>("DISCO_MADS_DETECTION_RADIUS", false);
            // Detection Radius
        if (detectionRadius < 0)
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DISCO_MADS_DETECTION_RADIUS must be positive" );
        }

        // Limit rate for detection
        auto limitRate = getAttributeValueProtected<NOMAD::Double>("DISCO_MADS_LIMIT_RATE", false);
        if (limitRate <= 0)
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DISCO_MADS_LIMIT_RATE must be strictly positive" );
        }

        // Exclusion Radius
        auto exclusionRadius = getAttributeValueProtected<NOMAD::Double>("DISCO_MADS_EXCLUSION_RADIUS", false);
        if (exclusionRadius <= 0)
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DISCO_MADS_EXCLUSION_RADIUS must be strictly positive" );
        }

        // Revealing poll radius
        auto revealingRadius = getAttributeValueProtected<NOMAD::Double>("DISCO_MADS_REVEALING_POLL_RADIUS", false);
        if (revealingRadius <= exclusionRadius+detectionRadius)
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DISCO_MADS_REVEALING_POLL_RADIUS must be strictly greater than DISCO_MADS_DETECTION_RADIUS + DISCO_MADS_EXCLUSION_RADIUS (use for instance 1.01*(DISCO_MADS_DETECTION_RADIUS + DISCO_MADS_EXCLUSION_RADIUS))" );
        }

        // Number of points for revealing poll
        const size_t revealingPointsNb = getAttributeValueProtected<size_t>("DISCO_MADS_REVEALING_POLL_NB_POINTS", false);
        if (revealingPointsNb==0)
        {
            // Warn the user
            std::cout << "Warning: the revealing poll is disabled as DISCO_MADS_REVEALING_POLL_NB_POINTS is null. This should only be used for testing as a strictly positive value is required for the convergence analysis." <<  std::endl;
        }

        if (revealingPointsNb<0) // probably useless as size_t negative values are already checked in Parameters::checkFormatSizeT
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DISCO_MADS_REVEALING_POLL_NB_POINTS must be strictly positive" );
        }

        //--- Compatibility with other options

        // Check compatibility with block of evaluations
        if (nullptr != evaluatorControlGlobalParams)
        {
            if (evaluatorControlGlobalParams->toBeChecked())
            {
                evaluatorControlGlobalParams->checkAndComply();
            }
            auto bbBlockSize = evaluatorControlGlobalParams->getAttributeValue<size_t>("BB_MAX_BLOCK_SIZE");
            if(bbBlockSize>1)
            {
                throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DiscoMads cannot be used with BB_MAX_BLOCK_SIZE > 1." );
            }
        }

        // Check here the initial frame size
        auto initialFrameSize = pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE");
        bool warning = false;
        for(size_t i=0;i<initialFrameSize.size();i++)
        {
            if(initialFrameSize[i]>1.1*revealingRadius)
            {
                warning=true;
                break;
            }
        }
        if(warning)
        {   // Warn the user
            std::cout << "Warning: INITIAL_FRAME_SIZE > 1.1*REVEALING_POLL_RADIUS,  this may lead to poor convergence. Choose an INITIAL_FRAME_SIZE similar to REVEALING_POLL_RADIUS may help." <<  std::endl;
        }

        // Use with openMP
        int nb_threads = evaluatorControlGlobalParams->getAttributeValue<int>("NB_THREADS_PARALLEL_EVAL");
        if(nb_threads>1)
        {
            std::cerr << "Warning: NB_THREADS_PARALLEL_EVAL>1. DiscoMads should not return any errors but it was not extensively validated with OpenMP. Prefer run on one thread if you want to stick to the theory." <<  std::endl;
        }

        // Use with quad models search
        bool quadModelSearch = getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH",false) || getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH_SIMPLE_MADS",false);
        if(quadModelSearch)
        {
            std::cerr << "Warning: it is currently not recommended to activate QUAD_MODEL_SEARCH or QUAD_MODEL_SEARCH_SIMPLE_MADS with DiscoMads as it may be much slower." <<  std::endl;
        }

        // DiscoMads is currently not compatible with megaSearchPoll (because revealingPoll is not seen by the megaSearchPoll)
        bool megaSearchPoll = getAttributeValueProtected<bool>("MEGA_SEARCH_POLL",false);
        if(megaSearchPoll)
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DiscoMads is not compatible with MEGA_SEARCH_POLL." );
        }

        // Variable group during revealing poll
        auto varGroups = pbParams->getAttributeValue<NOMAD::ListOfVariableGroup>("VARIABLE_GROUP",false);
        if (!varGroups.empty())
        {
            throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: DiscoMads is not compatible with VARIABLE_GROUP. This breaks the density properties of the revealing poll." );
        }
    }


    // COOP-Mads
    bool useAlgoCoopMads = getAttributeValueProtected<bool>("COOP_MADS_OPTIMIZATION", false);
#ifndef _OPENMP
    if (useAlgoCoopMads)
    {
        err = "Error: COOP_MADS_OPTIMIZATION can only be used when OpenMP is available. Please rebuild Nomad with OpenMP enabled.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
#else
    size_t nbCoopMads = getAttributeValueProtected<size_t>("COOP_MADS_NB_PROBLEM", false);
    // Test COOP-Mads nb problem
    if (nbCoopMads <= 1)
    {
        err = "Error: COOP-Mads requires to have more than one problem to solve in parallel. COOP_MADS_NB_PROBLEM must be greater than 1.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
#endif // _OPENMP


    // PSD-Mads
    bool useAlgoPSDMads = getAttributeValueProtected<bool>("PSD_MADS_OPTIMIZATION", false);

#ifndef _OPENMP
    if (useAlgoPSDMads)
    {
        err = "Error: PSD_MADS_OPTIMIZATION can only be used when OpenMP is available.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
#else
    size_t nbPSDMads = getAttributeValueProtected<size_t>("PSD_MADS_NB_SUBPROBLEM", false);
    // Test PSD-Mads nb problem
    if (nbPSDMads <= 1)
    {
        err = "Error: PSD-Mads requires to have more than one sub-problem to solve in parallel. PSD_MADS_NB_SUBPROBLEM must be greater than 1.";
        throw NOMAD::Exception(__FILE__,__LINE__, err);
    }
#endif // _OPENMP

#ifdef _OPENMP
    if (useAlgoPSDMads)
    {
        std::string nbVarParamName = "PSD_MADS_NB_VAR_IN_SUBPROBLEM";
        const size_t nbVariablesInSubproblem = getAttributeValueProtected<size_t>(nbVarParamName, false);
        if (0 == nbVariablesInSubproblem || nbVariablesInSubproblem > n)
        {
            err = "Parameter " + nbVarParamName + " must be between 1 and " + NOMAD::itos(n);
            err += ". Value provided: " + NOMAD::itos(nbVariablesInSubproblem);
            throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
        }

        size_t nbMadsSubproblem = getAttributeValueProtected<size_t>("PSD_MADS_NB_SUBPROBLEM", false);
        // Distribute all the variables between subproblems of reduced dimension.
        if (nbMadsSubproblem == INF_SIZE_T)
        {
            nbMadsSubproblem = (size_t)std::round(n/nbVariablesInSubproblem)+1; // Add a mads for the pollster
            setAttributeValue("PSD_MADS_NB_SUBPROBLEM", nbMadsSubproblem);
        }

        if (useAlgoPSDMads)
        {
            int nbThreadsHard = static_cast<int>(std::thread::hardware_concurrency());
            if (nbMadsSubproblem > nbThreadsHard)
            {
                std::string s = "Warning: PSD_MADS_NB_SUBPROBLEM exceeds the number of threads registered for this hardware: ";
                s += NOMAD::itos(nbThreadsHard);
                s += ". If this is true, it is not efficient. Let's continue anyway.";
                std::cout << s << std::endl;
            }
        }

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
#endif

    // Test quad model search regular or simple mads. Cannot be both
    if ( getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH",false) && getAttributeValueProtected<bool>("QUAD_MODEL_SEARCH_SIMPLE_MADS",false))
    {
        throw NOMAD::InvalidParameter(__FILE__,__LINE__,"Quad model search using simple mads is incompatible with the regular quad model search. Please deactivate: QUAD_MODEL_SEARCH no.");
    }

    // Algorithm parameters: use an algorithm other than MADS.
    // They are mutually-exclusive.
    bool useAlgoLH = (getAttributeValueProtected<size_t>("LH_EVAL", false) > 0);
    bool useAlgoNM = getAttributeValueProtected<bool>("NM_OPTIMIZATION", false);
    bool useAlgoQuadOpt = getAttributeValueProtected<bool>("QUAD_MODEL_OPTIMIZATION", false);
    bool useAlgoSgtelibModel = getAttributeValueProtected<bool>("SGTELIB_MODEL_EVAL", false);

    int totalAlgoSet = (int)useAlgoLH + (int)useAlgoCS +(int)useAlgoNM + (int)useAlgoQuadOpt
                       + (int)useAlgoPSDMads + (int)useAlgoSgtelibModel + (int)useAlgoDMultiMads + (int)useAlgoDiscoMads;

    if (totalAlgoSet >= 2)
    {
        err = "Multiple parameters for algorithms are set. ";
        err += "These parameters are mutually exclusive:";
        if (useAlgoCS)
        {
            err += " CS_OPTIMIZATION";
        }
        if (useAlgoDMultiMads)
        {
            err += " DMULTIMADS_OPTIMIZATION";
        }
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
        if (useAlgoDiscoMads)
        {
            err += " DISCO_MADS_OPTIMIZATION";
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

    auto projectOnMesh = getAttributeValueProtected<bool>("SEARCH_METHOD_MESH_PROJECTION", false);
    auto granularity = pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("GRANULARITY");
    if (!projectOnMesh && granularity != NOMAD::ArrayOfDouble(n,0.0))
    {
        throw NOMAD::Exception(__FILE__,__LINE__, "Parameters check: GRANULARITY is defined and search method mesh projection is disabled. Mesh projection is required to maintain granularity of variable.");
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
        // With the alternative way of seeding the RNG we set xdef to s.
        // This cannot be done with s=0. Exception is triggered.
        bool rngAltSeeding = getAttributeValueProtected<bool>("RNG_ALT_SEEDING",false);
        if (rngAltSeeding)
        {
            NOMAD::RNG::setSeedForXDef(seedToSet);
        }
        else
        {
            NOMAD::RNG::setSeed (seedToSet);
        }
    }
    NOMAD::Double::setEpsilon ( getAttributeValueProtected<NOMAD::Double>("EPSILON",false).todouble() );
    NOMAD::Double::setHMin ( getAttributeValueProtected<NOMAD::Double>("H_MIN",false).todouble() );
    NOMAD::Double::setUndefStr ( getAttributeValueProtected<std::string>("UNDEF_STR",false) );
    NOMAD::Double::setInfStr ( getAttributeValueProtected<std::string>("INF_STR",false) );

    // Reset parameter values from these static values, to ensure coherence.
    // This is bad because we have twice the same value for some parameters.
    setAttributeValue ( "SEED", NOMAD::RNG::getSeed() );
    setAttributeValue ( "EPSILON", NOMAD::Double(NOMAD::Double::getEpsilon()) );
    setAttributeValue ( "H_MIN", NOMAD::Double(NOMAD::Double::getHMin()) );
    setAttributeValue ( "UNDEF_STR", NOMAD::Double::getUndefStr() );
    setAttributeValue ( "INF_STR", NOMAD::Double::getInfStr() );
}

// These set methods are for CatMads (not yet available) and are not used in the current version of NOMAD.
// Probably not needed once a CatMads algorithm is available.
bool NOMAD::RunParameters::setMapDirTypeToVG(const std::shared_ptr<NOMAD::PbParameters>& pbParams, std::map<NOMAD::DirectionType,NOMAD::ListOfVariableGroup> & mapDirTypeToVG)
{
    if (_toBeChecked)
    {
        std::string errorMsg = "Cannot set map between direction type and variable group before checkAndComply is done";
        throw NOMAD::Exception(__FILE__,__LINE__, errorMsg);
    }

    auto listVG = pbParams->getAttributeValue<ListOfVariableGroup>("VARIABLE_GROUP");

    // Check that provided map is consistent with variables groups defined in pb
    for (const auto &vgBase: listVG)
    {
        bool vgBaseFound = false;
        for(const auto & elMap: mapDirTypeToVG )
        {
            for(const auto & vgMap: elMap.second)
            {
                if ( vgMap == vgBase )
                {
                    vgBaseFound = true;
                    break;
                }
            }
            if (vgBaseFound)
            {
                break;
            }
        }
        if (!vgBaseFound)
        {
            return false;
        }
    }
    _mapDirTypeToVG = mapDirTypeToVG;
    return true;
}

bool NOMAD::RunParameters::setListFixVGForQuadModelSearch(const std::shared_ptr<NOMAD::PbParameters>& pbParams, const NOMAD::ListOfVariableGroup & listFixVG)
{
    if (_toBeChecked)
    {
        std::string errorMsg = "Cannot set fixed variable group for QMS before checkAndComply is done";
        throw NOMAD::Exception(__FILE__,__LINE__, errorMsg);
    }

    auto listVG = pbParams->getAttributeValue<ListOfVariableGroup>("VARIABLE_GROUP");

    if (listVG.empty() && !listFixVG.empty())
    {
        std::string errorMsg = "Cannot set fixed variable group for QMS if no variable group is defined";
        throw NOMAD::Exception(__FILE__,__LINE__, errorMsg);
    }

    // Check that list of VG is consistent with variables groups defined in pb
    for(const auto & vg: listFixVG )
    {
        bool vgFound = false;
        for (const auto &vgBase: listVG)
        {
            if ( vgBase == vg )
            {
                vgFound = true;
                break;
            }
        }
        if (!vgFound)
        {
            return false;
        }
    }


    _fixVGForQMS = listFixVG;
    return true;
}

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

#include <iomanip>  // For std::setprecision
#include "../Math/RNG.hpp"
#include "../Type/BBInputType.hpp"
#include "../Param/PbParameters.hpp"

/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::PbParameters::init()
{
    
    _typeName = "Problem";
    
    try {
        #include "../Attribute/pbAttributesDefinition.hpp"
        registerAttributes( _definition ); // Registering attributes must be done for each instance of PbParameters
        
        // Note: we cannot call checkAndComply() here, the default values
        // are not valid, for instance DIMENSION, X0, etc.
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
void NOMAD::PbParameters::checkAndComply( )
{
    std::string err;

    checkInfo();
    
    if (!toBeChecked())
    {
        // Early out
        return;
    }
    
    // DIMENSION is the most important parameter
    size_t n = getAttributeValueProtected<size_t>("DIMENSION",false);
    if (n == 0)
    {
        throw NOMAD::InvalidParameter(__FILE__,__LINE__, "Parameters check: DIMENSION must be positive" );
    }
    
    auto lb = getAttributeValueProtected<NOMAD::ArrayOfDouble>("LOWER_BOUND",false);
    if (lb.size() != n)
    {
        if (lb.size() > 0)
        {
            err = "Warning: Parameter LOWER_BOUND resized from ";
            err += std::to_string(lb.size()) + " to " + std::to_string(n);
            err += ". Values may be lost.";
            std::cerr << err << std::endl;
        }

        lb.resize(n);
        setAttributeValue("LOWER_BOUND", lb);
    }
    
    auto ub = getAttributeValueProtected<NOMAD::ArrayOfDouble>("UPPER_BOUND",false);
    if (ub.size() != n)
    {
        if (ub.size() > 0)
        {
            err = "Warning: Parameter UPPER_BOUND resized from ";
            err += std::to_string(ub.size()) + " to " + std::to_string(n);
            err += ". Values may be lost.";
            std::cerr << err << std::endl;
        }

        ub.resize(n);
        setAttributeValue("UPPER_BOUND", ub);
    }
    
    auto x0s = getAttributeValueProtected<NOMAD::ArrayOfPoint>("X0",false);
    if (x0s.empty())
    {
        // When Initialization::eval_x0s() is called, the best
        // points in the cache will be used.
        // For now, set a x0 point "to be defined".
        NOMAD::Point x0(n);
        for (size_t i = 0 ; i < n ; i++)
        {
            x0[i].setToBeDefined();
        }
        x0s.push_back(x0);
        setAttributeValue("X0", x0s);
    }
    
    NOMAD::ArrayOfPoint x0SResized;
    for (size_t x0index = 0; x0index < x0s.size(); x0index++)
    {
        auto x0 = x0s[x0index];
        if ( x0.size() != n )
        {
            if (x0.size() > n)
            {
                err = "Warning: Parameter X0 resized from ";
                err += std::to_string(x0.size()) + " to " + std::to_string(n);
                err += ". Values may be lost or reset.";
                std::cerr << err << std::endl;
                x0.resize(n, x0[0]);
            }
            else
            {
                // X0 empty
                err = "Error: X0 " + x0.display() + " has dimension ";
                err += std::to_string(x0.size()) + " which is different from ";
                err += "problem dimension " + std::to_string(n);
                throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
            }
        }
        x0SResized.push_back(x0);
    }
    setAttributeValue("X0", x0SResized);
    
    for ( size_t i = 0 ; i < n ; ++i )
    {
        if ( lb[i].isDefined() && ub[i].isDefined() )
        {
            if ( lb[i] > ub[i] )
            {
                err = "Check: LOWER_BOUND is greater than UPPER_BOUND at index ";
                err += std::to_string(i);
                err += ". Lower bound = " + lb[i].tostring();
                err += "; Upper bound = " + ub[i].tostring();
                throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
            }
            
            if ( lb[i] == ub[i] )
            {
                err = "Check: LOWER_BOUND is equal to UPPER_BOUND at index ";
                err += std::to_string(i);
                err += ". Value = " + lb[i].tostring();
                throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
            }
        }
    }
    
    // Check for deprecated ***_POLL_SIZE parameters --> convert into ***_FRAME_SIZE
    auto iPS = getAttributeValueProtected<NOMAD::ArrayOfDouble>("INITIAL_POLL_SIZE",false);
    auto iFS = getAttributeValueProtected<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE",false);
    if ( iFS.isDefined() && iPS.isDefined() )
    {
        std::string s = "Check: INITIAL_POLL_SIZE (deprecated) and INITIAL_FRAME_SIZE are defined. Just one must be set";
        throw NOMAD::InvalidParameter(__FILE__,__LINE__, s);
    }
    if ( iPS.isDefined() && ! iFS.isDefined())
        setAttributeValue("INITIAL_FRAME_SIZE", iPS );
 
    auto mPS = getAttributeValueProtected<NOMAD::ArrayOfDouble>("MIN_POLL_SIZE",false);
    auto mFS = getAttributeValueProtected<NOMAD::ArrayOfDouble>("MIN_FRAME_SIZE",false);
    if ( mFS.isDefined() && mPS.isDefined() )
    {
        std::string s = "Check: MIN_POLL_SIZE (deprecated) and MIN_FRAME_SIZE are defined. Just one must be set.";
        throw NOMAD::InvalidParameter(__FILE__,__LINE__, s);
    }
    if ( mPS.isDefined() && ! mFS.isDefined())
    {
        setAttributeValue("MIN_FRAME_SIZE", iPS );
    }
    
    
    
    /*--------------------*/
    /* Granular variables */
    /*--------------------*/
    setGranularityAndBBInputType();

    /*-----------------*/
    /* Fixed variables */
    /*-----------------*/
    setFixedVariables();
    
    /* ------------------------------------------------------------------*/
    /* Lower and upper bounds might have been adjusted when              */
    /* granularity, BBInputTypes, and fixed variables were set.          */
    /* Check X0 against lower and upper bounds now.                      */
    /* ------------------------------------------------------------------*/
    checkX0AgainstBounds();
    
    /*----------------------------*/
    /*       Poll and Mesh        */
    /*----------------------------*/
    
    setMinMeshParameters("MIN_MESH_SIZE");
    setMinMeshParameters("MIN_FRAME_SIZE");
    
    setInitialMeshParameters();
    
    // Verify X0 and MESH parameters are conform with granularity.
    checkX0ForGranularity();
    checkForGranularity("MIN_MESH_SIZE");
    checkForGranularity("MIN_FRAME_SIZE");
    checkForGranularity("INITIAL_MESH_SIZE");
    checkForGranularity("INITIAL_FRAME_SIZE");
    
    _toBeChecked = false;
    
}
// End checkAndComply()


// Set and adjust GRANULARITY and BB_INPUT_TYPE.
// Adjust LOWER_BOUND and UPPER_BOUND if needed.
void NOMAD::PbParameters::setGranularityAndBBInputType()
{
    const size_t n = getAttributeValueProtected<size_t>("DIMENSION",false);
    auto granularity = getAttributeValueProtected<NOMAD::ArrayOfDouble>("GRANULARITY",false);
    auto bbInputType = getAttributeValueProtected<NOMAD::BBInputTypeList>("BB_INPUT_TYPE",false);
    auto lb = getAttributeValueProtected<NOMAD::ArrayOfDouble>("LOWER_BOUND",false);
    auto ub = getAttributeValueProtected<NOMAD::ArrayOfDouble>("UPPER_BOUND",false);
    std::ostringstream oss;
    
    if (granularity.isDefined() && granularity.size() != n)
    {
        std::string err;
        err = "Warning: Parameter GRANULARITY resized from ";
        err += std::to_string(granularity.size()) + " to " + std::to_string(n);
        err += ". Values may be lost.";
        std::cerr << err << std::endl;

        granularity.resize(n);
    }
    
    if (!granularity.isDefined())
    {
        granularity = NOMAD::ArrayOfDouble(n, 0.0);
        setAttributeValue("GRANULARITY", granularity);
    }
    
    
    for (size_t i = 0; i < n; i++)
    {
        if (!granularity[i].isDefined())
        {
            granularity[i] = 0.0;
        }
        else if (granularity[i] < 0.0)
        {
            throw NOMAD::InvalidParameter(__FILE__,__LINE__, "Check: invalid granular variables (negative values)" );
        }
    }
    // Update attribute GRANULARITY
    setAttributeValue("GRANULARITY", granularity);
    // Now we will adjust GRANULARITY with BB_INPUT_TYPE.
    
    
    /*----------------------*/
    /*  Init BB_INPUT_TYPE  */
    /*----------------------*/
    if (bbInputType.empty())
    {
        bbInputType.resize(n);
        // By default, all CONTINUOUS.
        std::fill(bbInputType.begin(), bbInputType.end(), NOMAD::BBInputType::CONTINUOUS);
        
        setAttributeValue("BB_INPUT_TYPE", bbInputType);
    }
    
    // Fix BB_INPUT_TYPE in a list form using the 'all of the same type (*)' syntax (for example: *R)
    std::vector<NOMAD::BBInputType>::const_iterator it=bbInputType.begin();
    switch (*it)
    {
        case NOMAD::BBInputType::ALL_CONTINUOUS:
            bbInputType.resize(n);
            std::fill(bbInputType.begin(), bbInputType.end(), NOMAD::BBInputType::CONTINUOUS);
            setAttributeValue("BB_INPUT_TYPE", bbInputType);
            break;
        case NOMAD::BBInputType::ALL_INTEGER:
            bbInputType.resize(n);
            std::fill(bbInputType.begin(), bbInputType.end(), NOMAD::BBInputType::INTEGER);
            setAttributeValue("BB_INPUT_TYPE", bbInputType);
            break;
        case NOMAD::BBInputType::ALL_BINARY:
            bbInputType.resize(n);
            std::fill(bbInputType.begin(), bbInputType.end(), NOMAD::BBInputType::BINARY);
            setAttributeValue("BB_INPUT_TYPE", bbInputType);
            break;
        default:
            break;

    }

    if (!isSetByUser("BB_INPUT_TYPE"))
    {
        bbInputType.resize(n);
        setAttributeValue("BB_INPUT_TYPE", bbInputType);
    }

    // By this point, bbInputType must be of size n.
    if (n != bbInputType.size())
    {
        oss << "Error: BB_INPUT_TYPE " << bbInputType << " has dimension ";
        oss << bbInputType.size() << " which is different from ";
        oss << "problem dimension " << n;
        throw NOMAD::InvalidParameter(__FILE__,__LINE__, oss.str());
    }
    
    /*-----------------------------------------*/
    /*  Adjust granularity with BB_INPUT_TYPE  */
    /*-----------------------------------------*/
   
    size_t i = 0;
    for (it = bbInputType.begin(); it != bbInputType.end(); ++it, i++)
    {
        switch (*it)
        {
            case NOMAD::BBInputType::CONTINUOUS:
                // do nothing - keep granularity as is.
                break;
            case NOMAD::BBInputType::INTEGER:
                // Ensure we have an integer
                if (!granularity[i].isInteger())
                {
                    oss << "Check: Invalid granularity[";
                    oss << i;
                    oss << "] = " << granularity[i];
                    oss << " for BB_INPUT_TYPE INTEGER";
                    throw NOMAD::InvalidParameter(__FILE__, __LINE__, oss.str());
                }
                // ensure we have 1 (or could be more)
                granularity[i] = max(granularity[i], 1.0);
                break;
            case NOMAD::BBInputType::BINARY:
                // ensure we have exactly 1
                granularity[i] = 1.0;
                // set lower and upper bounds too
                lb[i] = 0.0;
                ub[i] = 1.0;
                break;
            default:
                break;
        }
    }
    // Update attributes GRANULARITY and bounds
    setAttributeValue("GRANULARITY", granularity);
    setAttributeValue("LOWER_BOUND", lb);
    setAttributeValue("UPPER_BOUND", ub);
    
}


// This -> should be set at read
// Set fixed variables to their provided values, or to their
// X0 values if no provided values
//
// This -> should be done here.
// Adjust lower and upper bounds to fixed variables.
void NOMAD::PbParameters::setFixedVariables()
{
    const size_t n = getAttributeValueProtected<size_t>("DIMENSION",false);
    // Get bounds
    auto lb = getAttributeValueProtected<NOMAD::ArrayOfDouble>("LOWER_BOUND",false);
    auto ub = getAttributeValueProtected<NOMAD::ArrayOfDouble>("UPPER_BOUND",false);
    // Get X0
    const auto x0s = getAttributeValueProtected<NOMAD::ArrayOfPoint>("X0",false);
    auto fixedVariable = getAttributeValueProtected<NOMAD::Point>("FIXED_VARIABLE",false);

    if (!fixedVariable.isDefined())
    {
        fixedVariable.resize(n);
    }
    else if (fixedVariable.size() != n)
    {
        std::string err = "Warning: Parameter FIXED_VARIABLE resized from ";
        err += std::to_string(fixedVariable.size()) + " to " + std::to_string(n);
        err += ". Values may be lost.";
        std::cerr << err << std::endl;

        fixedVariable.resize(n);
    }
    
    for (size_t x0index = 0; x0index < x0s.size(); x0index++)
    {
        auto x0 = x0s[x0index];
        for (size_t i = 0; i < n; i++)
        {
            if (fixedVariable[i].toBeDefined() && x0[i].isDefined())
            {
                // Set fixed variable value to X0.
                fixedVariable[i] = x0[i];
            }
            else if (fixedVariable[i].isDefined() && x0[i].isDefined())
            {
                // Verify the fixed variable value is the same as X0.
                if (fixedVariable[i] != x0[i])
                {
                    std::ostringstream oss;
                    oss << "Parameters check: fixed variable must be the same as x0 at index ";
                    oss << i << ":" << std::endl;
                    oss << "v["; oss << i << "] = " << fixedVariable[i];
                    oss << " != x0[" << i << "] = " << x0[i] << std::endl;
                    throw NOMAD::InvalidParameter(__FILE__, __LINE__, oss.str());
                }
            }

            if (fixedVariable[i].isDefined())
            {
                if (lb[i].isDefined() && fixedVariable[i] < lb[i])
                {
                    std::ostringstream oss;
                    oss << "Parameters check: fixed variable v under lower bound: v[";
                    oss << i << "] = " << fixedVariable[i];
                    oss << " < " << lb[i] << std::endl;
                    throw NOMAD::InvalidParameter(__FILE__, __LINE__, oss.str());
                }
                if (ub[i].isDefined() && fixedVariable[i] > ub[i])
                {
                    std::ostringstream oss;
                    oss << "Parameters check: fixed variable v over upper bound: v[";
                    oss << i << "] = " << fixedVariable[i];
                    oss << " > " << ub[i] << std::endl;
                    throw NOMAD::InvalidParameter(__FILE__, __LINE__, oss.str());
                }
            }
        }
    }

    // Update values
    setAttributeValue("FIXED_VARIABLE", fixedVariable);
    // Note: We could update upper bound and lower bound to fixed variable.
    // We choose not to do this. The subproblem will be on non-fixed variables
    // only. The fixed variables may be modified later to solve another
    // subproblem: the bounds will still be valid.
}


void NOMAD::PbParameters::checkX0AgainstBounds()
{
    const size_t n = getAttributeValueProtected<size_t>("DIMENSION",false);
    // Get bounds
    auto lb = getAttributeValueProtected<NOMAD::ArrayOfDouble>("LOWER_BOUND",false);
    auto ub = getAttributeValueProtected<NOMAD::ArrayOfDouble>("UPPER_BOUND",false);
    // Get X0
    const auto x0s = getAttributeValueProtected<NOMAD::ArrayOfPoint>("X0",false);
    
    for (size_t x0index = 0; x0index < x0s.size(); x0index++)
    {
        auto x0 = x0s[x0index];
        for ( size_t i = 0 ; i < n ; ++i )
        {
            // Check that x0s are within bounds when defined
            if (!x0[i].isDefined())
            {
                continue;
            }
            if (lb[i].isDefined())
            {
                if (x0[i] < lb[i])
                {
                    std::ostringstream oss;
                    oss << "Parameters check: x0 under lower bound: x0[";
                    oss << i << "] = " << x0[i];
                    oss << " < " << lb[i] << " " << x0.display() << std::endl;
                    throw NOMAD::InvalidParameter(__FILE__, __LINE__, oss.str());
                }
            }
            if (ub[i].isDefined())
            {
                if (x0[i] > ub[i])
                {
                    std::ostringstream oss;
                    oss << "Parameters check: x0 over upper bound: x0[";
                    oss << i << "] = " << x0[i];
                    oss << " > " << ub[i] << std::endl;
                    throw NOMAD::InvalidParameter(__FILE__, __LINE__, oss.str());
                }
            }
        }
    }
}


// Set and adjust MIN_MESH_SIZE and MIN_FRAME_SIZE.
void NOMAD::PbParameters::setMinMeshParameters(const std::string &paramName)
{
    const size_t n = getAttributeValueProtected<size_t>("DIMENSION",false);
    const auto granularity = getAttributeValueProtected<NOMAD::ArrayOfDouble>("GRANULARITY",false);
    
    // minArray = either min mesh size (\delta_min) or min frame size (\Delta min).
    auto minArray = getAttributeValueProtected<NOMAD::ArrayOfDouble>(paramName,false);
    
    if (!minArray.isDefined())
    {
        // Default values: granularity if it is > 0, epsilon otherwise.
        minArray = NOMAD::ArrayOfDouble(n, NOMAD::Double::getEpsilon());
        for (size_t i = 0 ; i < n ; ++i)
        {
            if (0.0 < granularity[i])
            {
                minArray[i] = granularity[i];
            }
        }
    }
    else
    {
        if (minArray.size() != n)
        {
            std::string err = "Check: dimension of parameter " + paramName;
            err += " resized from " + std::to_string(minArray.size()) + " to " + std::to_string(n);
            err += ". Values may be lost.";
            std::cerr << err << std::endl;

            minArray.resize(n);
        }
        
        for (size_t i = 0 ; i < n ; ++i)
        {
            if (minArray[i].isDefined() && minArray[i].todouble() <= 0.0)
            {
                std::string err = "Check: invalid value for parameter " + paramName;
                throw NOMAD::InvalidParameter(__FILE__, __LINE__, err);
            }
            else if (!minArray[i].isDefined()
                     || minArray[i].todouble() < NOMAD::Double::getEpsilon()
                     || (0.0 < granularity[i] && minArray[i].todouble() < granularity[i]) )
            {
                // Set default value
                minArray[i] = NOMAD::Double::getEpsilon();
                if (0.0 < granularity[i])
                {
                    minArray[i] = granularity[i];
                }
            }
        }
    }
    // Update attribute
    setAttributeValue(paramName, minArray);
}


void NOMAD::PbParameters::setInitialMeshParameters()
{
    // Set INITIAL_MESH_SIZE and INITIAL_FRAME_SIZE.
    // Take into account GRANULARITY, MIN_MESH_SIZE and MIN_FRAME_SIZE.
    const size_t n = getAttributeValueProtected<size_t>("DIMENSION",false);
    const auto minMeshSize = getAttributeValueProtected<NOMAD::ArrayOfDouble>("MIN_MESH_SIZE",false);
    const auto minFrameSize = getAttributeValueProtected<NOMAD::ArrayOfDouble>("MIN_FRAME_SIZE",false);
    const auto granularity = getAttributeValueProtected<NOMAD::ArrayOfDouble>("GRANULARITY",false);
    auto initialMeshSize = getAttributeValueProtected<NOMAD::ArrayOfDouble>("INITIAL_MESH_SIZE",false);
    auto initialFrameSize = getAttributeValueProtected<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE",false);
    // Remember the previous values to show a warning only if the values changed.
    auto initialMeshSize0 = initialMeshSize;
    const auto lb = getAttributeValueProtected<NOMAD::ArrayOfDouble>("LOWER_BOUND",false);
    const auto ub = getAttributeValueProtected<NOMAD::ArrayOfDouble>("UPPER_BOUND",false);
    const auto x0s = getAttributeValueProtected<NOMAD::ArrayOfPoint>("X0",false);
    bool warningInitialFrameSizeReset = true;
    
    // Basic checks
    if (initialMeshSize.isDefined() && initialMeshSize.size() != n)
    {
        std::string err = "Warning: Parameter INITIAL_MESH_SIZE resized from ";
        err += std::to_string(initialMeshSize.size()) + " to " + std::to_string(n);
        err += ". Values may be lost.";
        std::cerr << err << std::endl;
        initialMeshSize.resize(n);
    }
    
    if (initialFrameSize.isDefined() && initialFrameSize.size() != n)
    {
        std::string err = "Warning: Parameter INITIAL_MESH_SIZE resized from ";
        err += std::to_string(initialFrameSize.size()) + " to " + std::to_string(n);
        err += ". Values may be lost.";
        std::cerr << err << std::endl;
        initialFrameSize.resize(n);
    }
    
    if (initialMeshSize.isDefined() && initialFrameSize.isDefined())
    {
        //initialMeshSize will be redefined from initialFrameSize.
        initialMeshSize.reset(n);
    }
    
    if (!initialMeshSize.isDefined())
    {
        initialMeshSize = NOMAD::ArrayOfDouble(n , NOMAD::Double()) ;
        setAttributeValue("INITIAL_MESH_SIZE", initialMeshSize);
    }
    
    if (!initialFrameSize.isDefined())
    {
        initialFrameSize = NOMAD::ArrayOfDouble(n, NOMAD::Double()) ;
        setAttributeValue("INITIAL_FRAME_SIZE", initialFrameSize);
    }
    
    
    // initial mesh size or frame size:
    // --------------------------------
    for (size_t i = 0; i < n; ++i)
    {
        // Forces initial frame size from initial mesh size
        if (initialMeshSize[i].isDefined())
        {
            if (initialFrameSize[i].isDefined() && warningInitialFrameSizeReset)
            {
                warningInitialFrameSizeReset = false;
                std::cerr << "Warning: initial frame size reset from initial mesh." << std::endl;
            }
            initialFrameSize[i] = initialMeshSize[i] * pow(n, 0.5);
            // Adjust value with granularity
            initialFrameSize[i] = initialFrameSize[i].nextMult(granularity[i]);
            // Adjust value with minFrameSize
            if (initialFrameSize[i] < minFrameSize[i])
            {
                initialFrameSize[i] = minFrameSize[i];
            }
        }
        
        
        // Compute x0lb and x0ub for frame size initialization.
        NOMAD::Point x0lb(n);
        NOMAD::Point x0ub(n);
        for (size_t j = 0; j < n; j++)
        {
            for (size_t x0index = 0; x0index < x0s.size(); x0index++)
            {
                auto x0 = x0s[x0index];
                if (!x0lb[j].isDefined() || x0[j] < x0lb[j])
                {
                    x0lb[j] = x0[j];
                }
                if (!x0ub[j].isDefined() || x0[j] > x0ub[j])
                {
                    x0ub[j] = x0[j];
                }
            }
        }
            
        // default value for initial mesh/frame size
        if (!initialFrameSize[i].isDefined())
        {
            
            if (lb[i].isDefined() &&  ub[i].isDefined())
            {
                initialFrameSize.set(i, NOMAD::Double(0.1), true , lb[i],ub[i]);
                
            }
            else if (lb[i].isDefined() && x0lb[i].isDefined() && lb[i] != x0lb[i])
            {
                initialFrameSize[i] = (x0lb[i]-lb[i]) / 10.0;   // Case x0 < lb tested elsewhere
            }
            else if (ub[i].isDefined() && x0ub[i].isDefined() && ub[i] != x0ub[i])
            {
                initialFrameSize[i] = (ub[i]-x0ub[i]) / 10.0;   // Case x0 > ub tested elsewhere
            }
            else
            {
                if (x0lb[i].isDefined() && (x0lb[i].abs() > NOMAD::Double::getEpsilon() * 10.0))
                {
                    initialFrameSize[i] = x0lb[i].abs() / 10.0;
                }
                else
                {
                    initialFrameSize[i] = 1.0;
                }
            }
            // Adjust value with granularity
            initialFrameSize[i] = initialFrameSize[i].nextMult(granularity[i]);
            // Adjust value with minFrameSize
            if (initialFrameSize[i] < minFrameSize[i])
            {
                initialFrameSize[i] = minFrameSize[i];
            }
        }
        // Determine initial mesh size from initial frame size
        if (!initialMeshSize[i].isDefined())
        {
            initialMeshSize[i] = initialFrameSize[i] * pow(n, -0.5);
            // Adjust value with granularity
            initialMeshSize[i] = initialMeshSize[i].nextMult(granularity[i]);
            // Adjust value with minMeshSize
            if (initialMeshSize[i] < minMeshSize[i])
            {
                initialMeshSize[i] = minMeshSize[i];
            }
        }
    }

    if (_showWarningMeshSizeRedefined && initialMeshSize0.isDefined() && initialMeshSize0 != initialMeshSize)
    {
        std::string err = "Warning: initial mesh size reset from initial frame size.\n";
        err += "INITIAL_MESH_SIZE: " + initialMeshSize.display() + "\n";
        err += "INITIAL_FRAME_SIZE: " + initialFrameSize.display() + "\n";
        std::cerr << err;

        // Show the warning only once.
        _showWarningMeshSizeRedefined = false;
    }

    setAttributeValue("INITIAL_FRAME_SIZE", initialFrameSize);
    setAttributeValue("INITIAL_MESH_SIZE", initialMeshSize);
    
    if (!(minMeshSize <= initialMeshSize))
    {
        std::string err = "Check: initial mesh size is lower than min mesh size.\n";
        err += "INITIAL_MESH_SIZE " + initialMeshSize.display() + "\n";
        err += "MIN_MESH_SIZE " + minMeshSize.display();
        throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
    }
    if (!(minFrameSize <= initialFrameSize))
    {
        std::string err = "Check: initial frame size is lower than min frame size.\n";
        err += "INITIAL_FRAME_SIZE\t" + initialFrameSize.display() + "\n";
        err += "MIN_FRAME_SIZE\t\t" + minFrameSize.display();
        throw NOMAD::InvalidParameter(__FILE__,__LINE__, err);
    }
    
}


void NOMAD::PbParameters::checkX0ForGranularity() const
{
    auto x0s = getAttributeValueProtected<NOMAD::ArrayOfPoint>("X0", false);
    for (size_t x0index = 0; x0index < x0s.size(); x0index++)
    {
        auto x0 = x0s[x0index];
        if (!x0.toBeDefined())
        {
            checkForGranularity("X0", x0);
        }
    }
}


void NOMAD::PbParameters::checkForGranularity(const std::string &paramName) const
{
    // Assuming paramName is of type ArrayOfDouble.    
    NOMAD::ArrayOfDouble arrayToCheck = getAttributeValueProtected<NOMAD::ArrayOfDouble>(paramName,false);
    checkForGranularity(paramName, arrayToCheck);
}
    

void NOMAD::PbParameters::checkForGranularity(const std::string &paramName, const NOMAD::ArrayOfDouble &arrayToCheck) const
{
    const auto granularity = getAttributeValueProtected<NOMAD::ArrayOfDouble>("GRANULARITY",false);

    int index = -1;
    if (!arrayToCheck.isMultipleOf(granularity, index))
    {
        std::ostringstream oss;
        oss << std::setprecision(16);
        oss << "Check: Invalid granularity of parameter ";
        oss << paramName << " at index ";
        oss << index;
        oss << ": " << arrayToCheck[index];
        oss << " vs granularity value " << granularity[index];
        throw NOMAD::InvalidParameter(__FILE__, __LINE__, oss.str());
    }
}




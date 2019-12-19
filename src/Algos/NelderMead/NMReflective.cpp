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
    
#include "../../Algos/NelderMead/NMReflective.hpp"
#include "../../Algos/NelderMead/NMIteration.hpp"
#include "../../Algos/NelderMead/NMUpdate.hpp"


const NOMAD::Double deltaR = 1;

void NOMAD::NMReflective::init()
{
    _name = getAlgoName() + "Step";
    
    _currentStepType = _nextStepType = NMStepType::UNSET;
    
     _deltaE = _runParams->getAttributeValue<NOMAD::Double>("NM_DELTA_E");
     _deltaIC = _runParams->getAttributeValue<NOMAD::Double>("NM_DELTA_IC");
     _deltaOC = _runParams->getAttributeValue<NOMAD::Double>("NM_DELTA_OC");

    
    if ( _deltaE <= 1 )
        throw NOMAD::Exception(__FILE__,__LINE__,"Delta value deltaE not compatible with expansion");
    
    if ( _deltaOC < 0 || _deltaOC > 1 )
        throw NOMAD::Exception(__FILE__,__LINE__,"Delta value deltaOC not compatible with outside contraction");
    
    if ( _deltaIC > 0 )
        throw NOMAD::Exception(__FILE__,__LINE__,"Delta value deltaIC not compatible with inside contraction");
    
    auto nmOptimization = _runParams->getAttributeValue<bool>("NM_OPTIMIZATION");
    auto nmSearchRankEps = _runParams->getAttributeValue<NOMAD::Double>("NM_SEARCH_RANK_EPS");
    _rankEps = ( nmOptimization ) ? NOMAD::DEFAULT_EPSILON:nmSearchRankEps;
    
    verifyParentNotNull();
 

}

void NOMAD::NMReflective::setCurrentNMStepType ( NMStepType stepType )
{
    _currentStepType = stepType;
    
    switch ( _currentStepType ) {
        case NMStepType::REFLECT:
            _name = getAlgoName() + "Reflect";
            _delta = deltaR;
            break;
        case NMStepType::EXPAND:
            _name = getAlgoName() + "Expansion";
            _delta = _deltaE;
            break;
        case NMStepType::OUTSIDE_CONTRACTION:
            _name = getAlgoName() + "Outside Contraction";
            _delta = _deltaOC;
            break;
        case NMStepType::INSIDE_CONTRACTION:
            _name = getAlgoName() + "Inside Contraction";
            _delta = _deltaIC;
            break;
        default:
            throw NOMAD::Exception(__FILE__,__LINE__,"Only REFLECT, EXPAND, INSIDE_CONTRACTION and OUTSIDE_CONTRACTION are supported");
            break;
    }
}


void NOMAD::NMReflective::startImp()
{
    
    // Specific
    if ( _currentStepType == NMStepType::UNSET )
        throw NOMAD::Exception(__FILE__,__LINE__,"The NM step type must be set");
    
    // Create EvalPoints
    generateTrialPoints();
    
    verifyPointsAreOnMesh(getName());
    updatePointsWithFrameCenter();
    
}


bool NOMAD::NMReflective::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }
    
    // If trial point is in cache, it does not need eval. Stop NM step.
    if ( getNbEvalPointsThatNeededEval() == 0 )
        setStopReason( );
    else
    {
        
        // Update the eval points by step type
        switch ( _currentStepType )
        {
            case NMStepType::REFLECT:
                _xr = *(getTrialPoints().begin())->getX();
                break;
            case NMStepType::EXPAND:
                _xe = *(getTrialPoints().begin())->getX();;
                break;
            case NMStepType::OUTSIDE_CONTRACTION:
                _xoc = *(getTrialPoints().begin())->getX();;
                break;
            case NMStepType::INSIDE_CONTRACTION:
                _xic = *(getTrialPoints().begin())->getX();;
                break;
            default:
                throw NOMAD::Exception(__FILE__,__LINE__,"Current step must be REFLECT, EXPAND, OUTSIDE_CONTRACTION or INSIDE_CONTRACTION.");
                break;
        }
    }
    
    // Specific
    if ( ! _stopReasons->checkTerminate() )
        setNextNMStepType();
    
    // From IterationUtils
    postProcessing(getEvalType());
    
    return foundBetter;
}


void NOMAD::NMReflective::generateTrialPoints ()
{
    
    // The pb params handle only variables (fixed variables are not considered)
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    
    // The size of simplex must be large enough to perform reflexion
    size_t minYSize = n + 1;
    size_t YSize = _nmY->size();
    if ( YSize < minYSize )
    {
        AddOutputInfo("Not enough points in the simplex to generate points for " + _name +" (delta=" + _delta.tostring() +")");
        return;
    }
    
    AddOutputInfo("Generate point for " + _name +" (delta=" + _delta.tostring()+")");
    
    // Clear the previous trial points
    clearTrialPoints();

    // Determine the centroid of Y
    std::set<NOMAD::EvalPoint>::iterator it;
    std::set<NOMAD::EvalPoint>::iterator itBeforeEnd = _nmY->end();
    --itBeforeEnd;
    int i=0;
    NOMAD::Point yc(n,0);
    for (it = _nmY->begin() ; it != itBeforeEnd ; ++it, i++)
    {
        AddOutputInfo("y" + std::to_string(i) + ": " + (*it).display() );
        for (size_t k = 0 ; k < n ; ++k )
        {
            yc[k] += (*it)[k];
        }
    }
    yc *= 1.0/n;
    
    const NOMAD::Point & yn = *itBeforeEnd;
    AddOutputInfo("yn: " + yn.display() );
    AddOutputInfo("yc: " + yc.display() );
    
    NOMAD::Point d(n,0);
    for (size_t k = 0 ; k < n ; ++k )
    {
        d[k] = yc[k]-yn[k];
    }
    d *= _delta;

    // Creation of point
    NOMAD::Point xt(n);
    for (size_t k = 0 ; k < n ; ++k )
    {
        xt[k] = yc[k] + d[k];
    }

    auto lb = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto ub = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    
    if (snapPointToBoundsAndProjectOnMesh(xt, lb, ub, getIterationFrameCenter(), getIterationMesh()))
    {
        // New EvalPoint to be evaluated.
        // Add it to the list.
        bool inserted = insertTrialPoint(NOMAD::EvalPoint(xt));
        
        std::string s = "xr:";
        s += (inserted) ? " " : " not inserted: ";
        s += xt.display();
        AddOutputInfo(s);
    }
    
    if ( xt == yn )
    {
        clearTrialPoints();
        xt = NOMAD::Point();
        AddOutputInfo("No valid point generated: too close to yn.");
    }
    
    verifyPointsAreOnMesh(getName());
    updatePointsWithFrameCenter();
    
}

void NOMAD::NMReflective::setNextNMStepType()
{
    // Updated vector Y0 and Yn are needed to set the next step
    makeListY0();
    makeListYn();
    
    switch ( _currentStepType )
    {
        case NMStepType::REFLECT:
            setAfterReflect();
            break;
        case NMStepType::EXPAND:
            setAfterExpand();
            break;
        case NMStepType::OUTSIDE_CONTRACTION:
            setAfterOutsideContract();
            break;
        case NMStepType::INSIDE_CONTRACTION:
            setAfterInsideContract();
            break;
        default:
            throw NOMAD::Exception(__FILE__,__LINE__,"Current step must be REFLECT, EXPAND, OUTSIDE_CONTRACTION or INSIDE_CONTRACTION.");
            break;
    }

    // Unset the currectStepType.
    // The currentStepType must be set explicitely before each sequence Start, Run, End
    _currentStepType = NMStepType::UNSET;
}

void NOMAD::NMReflective::setAfterReflect( void )
{
    // Reflect is always the first step of NM iteration.
    // Depending on which zone x_r belongs, we perform an other step
    
    if ( _currentStepType != NMStepType::REFLECT )
        throw NOMAD::Exception(__FILE__,__LINE__,"The current step type should be REFLECT.");
    
    if ( ! _xr.isDefined() )
    {
        AddOutputDebug("The reflect point xr is not defined. Stop NM (no shrink).");
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    if ( pointDominatesY0( _xr ) )
    {
        // In NM-Mads paper: x_r belongs to the expansion zone. Next, expansion.
        AddOutputDebug("The reflect point xr: " + _xr.display() +" dominates Y0. Next perform expansion.");
        _nextStepType = NMStepType::EXPAND;
    }
    else if ( YnDominatesPoint( _xr ) )
    {
        // In NM-Mads paper: x_r belongs to the inside contraction zone. Next, perform inside contraction.
        AddOutputDebug("The reflect point xr: " + _xr.display() + " is dominated by Yn. Next, perform Inside Contraction.");
        _nextStepType = NMStepType::INSIDE_CONTRACTION;
    }
    else if ( pointDominatesPtsInY( _xr , 2 ) )
    {
        // In NM-Mads paper: x_r belongs to the reflection zone. Next, insert into simplex (3), update Y, Y0, Yn (1) and continue (2)
        AddOutputDebug("The reflect point xr: "  + _xr.display() + " dominates at least 2 points of Y.");
        
        _currentStepType = NMStepType::INSERT_IN_Y;
        if ( insertInY ( _xr ) )
        {
            AddOutputDebug("Insertion in Y is successful. NM iteration completed (no shrink).");
            _nextStepType = NMStepType::CONTINUE;
        }
        else
        {
            AddOutputDebug(" Cannot insert xr in Y. Perform shrink (if available).");
            _nextStepType = NMStepType::SHRINK;
        }
    }
    else if ( pointDominatesPtsInY( _xr , 1 ) || pointDominatesPtsInY( _xr , 0 ) )
    {
        // In NM-Mads paper: x_r belongs to the outside reflection zone. Next, perform outside contraction.
        AddOutputDebug("The reflect point xr:" + _xr.display() + " dominates 1 or 0 point of Y. Next, perform Outside Contraction.");
        _nextStepType = NMStepType::OUTSIDE_CONTRACTION;
    }
    else
    {
        setStopReason();
    }
}

void NOMAD::NMReflective::setAfterExpand( void )
{
    
    // Expand follows Reflect.
    // In NM-Mads paper: x_r belongs to the expansion zone. x_e has been evaluated. The best point between x_r and x_e is inserted in Y. NM iteration is completed.
    
    if ( _currentStepType != NMStepType::EXPAND )
        throw NOMAD::Exception(__FILE__,__LINE__,"The current step type should be EXPAND.");
    
    if ( ! _xe.isDefined() )
    {
        AddOutputDebug("The expansion point xe is not defined. Stop NM (no shrink).");
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    if ( ! _xr.isDefined() )
    {
        AddOutputDebug("The reflect point xr is not defined. Stop NM (no shrink).");
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    _currentStepType = NMStepType::INSERT_IN_Y;
    if ( insertInYBest( _xr , _xe ) )
    {
        AddOutputDebug("Insert in Y the best of xr and xe. NM iteration completed (no shrink). ");
        _nextStepType = NMStepType::CONTINUE;
    }
    else
    {
        // No point inserted.
        AddOutputDebug("The insertion in Y of the best of xr and xe dit not maintain a proper Y. Perform shrink (if available).");
        _nextStepType = NMStepType::SHRINK;
    }
}

void NOMAD::NMReflective::setAfterOutsideContract ( void )
{
    
    // Outside contraction follows Reflect.
    // In NM-Mads paper: x_r belongs to the outside contraction zone. x_oc has been evaluated. The best point between x_r and x_oc is inserted in Y. NM iteration is completed.
    if ( _currentStepType != NMStepType::OUTSIDE_CONTRACTION )
        throw NOMAD::Exception(__FILE__,__LINE__,"The current step type should be OUTSIDE_CONTRACTION.");
    
    if ( ! _xoc.isDefined() )
    {
        AddOutputDebug("The outside contraction point xoc is not defined. Stop NM (no shrink).");
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    if ( ! _xr.isDefined() )
    {
        AddOutputDebug("The reflect point xr is not defined. Stop NM (no shrink).");
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }

    _currentStepType = NMStepType::INSERT_IN_Y;
    if ( insertInYBest( _xr , _xoc ) )
    {
        AddOutputDebug("The insertion of the best of xr and xoc in Y is valid. NM iteration completed." );
        _nextStepType = NMStepType::CONTINUE;
        return;
    }
    else
    {
        AddOutputDebug( "The insertion of the best of xr and xoc in Y did not maintain a proper Y. Perform shrink (if available).");
        _nextStepType = NMStepType::SHRINK;
    }
    
}


void NOMAD::NMReflective::setAfterInsideContract ( void )
{
    
    // Inside contraction follows Reflect.
    // In NM-Mads paper: x_r belongs to the inside contraction zone. x_ic has been evaluated. If x_ic belongs to the inside contraction zone, stop NM. Otherwise insert x_ic in Y. NM iteration is completed.
    
    if ( _currentStepType != NMStepType::INSIDE_CONTRACTION )
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot set step after inside contraction if x_ic is not defined.");
    if ( ! _xic.isDefined() )
    {
        AddOutputDebug("The inside contraction point xic is not defined. Stop NM (no shrink).");
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    if ( YnDominatesPoint( _xic ) )
    {
        _nextStepType = NMStepType::SHRINK ;
        AddOutputDebug("Yn dominates xic: " + _xic.display() + " Next, perform Shrink.");
        return;
    }
    else
    {
        
        AddOutputDebug("The inside contraction point xic:" + _xic.display() + " is not dominated by Yn. Insert x_ic in Y." );
        
        _currentStepType = NMStepType::INSERT_IN_Y;
        if ( insertInY( _xic ) )
        {
            AddOutputDebug("Insertion in Y is successful. NM iteration completed (no shrink).");
            _nextStepType = NMStepType::CONTINUE;
        }
        else
        {
            // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
            AddOutputDebug("Cannot insert xic in Y. Perform Shrink if available." );
            _nextStepType = NMStepType::SHRINK;
        }
    }
}


bool NOMAD::NMReflective::insertInYBest( const NOMAD::Point & x1 , const NOMAD::Point & x2 )
{
    auto evalType = getEvalType();

    // Get the point from Cache
    NOMAD::EvalPoint X1;
    size_t findPt = NOMAD::CacheBase::getInstance()->find(x1,X1);
    
    if ( findPt == 0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"The trial point x1 is not in Cache" );
    
    bool x1EvalOK = ( X1.getEvalStatus(evalType) == NOMAD::EvalStatusType::EVAL_OK );
    
    NOMAD::EvalPoint X2;
    findPt = NOMAD::CacheBase::getInstance()->find(x2,X2);
    
    if ( findPt == 0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"The trial point x2 is not in Cache" );
    
    bool x2EvalOK = ( X2.getEvalStatus(evalType) == NOMAD::EvalStatusType::EVAL_OK );
    
    std::set<NOMAD::EvalPoint>::iterator itYn = _nmY->end();
    --itYn;
    
    AddOutputDebug("Insertion/Deletion of points in Y: ");
    
    if ( x1EvalOK )
    {
        std::pair<std::set<NOMAD::EvalPoint>::iterator,bool> ret_x1 = _nmY->insert ( X1 );
        
        // Insertion of x1 ok
        if ( ret_x1.second )
        {
            AddOutputDebug("Insertion in Y: ") ;
            AddOutputDebug(X1.display());
            
            std::pair<std::set<NOMAD::EvalPoint>::iterator,bool> ret_x2;
            
            bool same_x1_x2=false;
            
            //
            if ( ! x2EvalOK )
            {
                
                if ( std::distance(_nmY->begin(), ret_x1.first ) > std::distance(_nmY->begin(), itYn ) )
                {
                    _nmY->erase( ret_x1.first );
                    AddOutputDebug("Cannot perform insertion in Y. x1 is dominated by yn and x2 eval status is not ok. Delete x1.");
                    return false;
                }
                else
                    X2 = X1 ;  // Trick to continue only with X1 in the case X2 Eval is not ok. 
            }
            
            // Try to insert x2 only if x1!=x2
            if ( x2EvalOK && x1 != x2 )
                ret_x2 = _nmY->insert ( X2 );
            else
            {
                ret_x2.second = true;
                ret_x2.first = ret_x1.first;
                same_x1_x2 = true;
            }
            
            // Insertion of x2 ok
            if ( ret_x2.second )
            {
                AddOutputDebug("Insertion in Y: ") ;
                AddOutputDebug( X2.display());
                
                
                // if x1 and x2 are after yn --> remove both x1 and x2 ==> insertion failed
                if ( std::distance(_nmY->begin(), ret_x1.first ) > std::distance(_nmY->begin(), itYn ) && std::distance(_nmY->begin(), ret_x2.first ) > std::distance(_nmY->begin(), itYn ) )
                {
                    _nmY->erase( ret_x1.first );
                    if ( ! same_x1_x2 )
                        _nmY->erase( ret_x2.first );
                    
                    AddOutputDebug("Cannot perform insertion because both points are dominated by yn");
                    return false;
                }
                
                
                // if x1 after x2 --> keep x2 remove x1
                if ( std::distance(_nmY->begin(), ret_x1.first ) > std::distance(_nmY->begin(), ret_x2.first ) )
                {
                    AddOutputDebug("Delete from Y: ") ;
                    AddOutputDebug( X1.display() );
                    _nmY->erase( ret_x1.first );
                }
                else // if x2 after x1 --> keep x1 remove x2
                {
                    AddOutputDebug("Insertion in Y: ") ;
                    AddOutputDebug( X1.display());
                    if ( !same_x1_x2 )
                    {
                        AddOutputDebug("Delete from Y: ") ;
                        AddOutputDebug( X2.display() );
                        _nmY->erase( ret_x2.first );
                    }
                }
                
                // Remove yn
                std::set<NOMAD::EvalPoint>::iterator it = _nmY->end();
                --it;
                AddOutputDebug("Delete yn from Y: " + (*it).display());
                _nmY->erase ( it );
                
                // Update simplex characteristics (volumes and diameter)
                updateYCharacteristics();
                
                // Update vectors Y0 and Yn
                bool success = makeListY0();
                if ( !success )
                {
                    AddOutputDebug("Cannot make Y0. Let's continue." );
                    return false;
                }
                success = makeListYn();
                if ( !success )
                {
                    AddOutputDebug("Cannot make Yn. Let's continue" );
                    return false;
                }
                
                
                AddOutputDebug("After insertion best");
                displayYInfo();
                displayY0nInfo();
                
                if ( getRankDZ() != static_cast<int>(_nmY->size()) -1 )
                {
                    AddOutputInfo("Rank of DZ=[(y1-y0) (y2-y0) ... (yn-y0)] != n. Y is not a valid simplex. Let's continue. ");
                    return false;
                }
            }
            else
            {
                // Try to remove point from Y (in case it was there)
                _nmY->erase(ret_x2.first);
                
                // Update simplex characteristics (volumes and diameter)
                updateYCharacteristics();
                
                AddOutputDebug("After insertion cancelled");
                displayYInfo();
                displayY0nInfo();
                return false;
            }
        }
        else
        {
            // Try to remove point from Y (in case it was there)
            _nmY->erase(ret_x1.first);
            
            // Update simplex characteristics (volumes and diameter)
            updateYCharacteristics();
            
            AddOutputDebug("After insertion cancelled");
            displayYInfo();
            displayY0nInfo();
            return false;
        }
    }
    else
    {
        AddOutputDebug("Cannot insertion in Y because eval not OK: ") ;
        AddOutputDebug(X1.display());
    }
    return true;
}


bool NOMAD::NMReflective::insertInY( const NOMAD::Point & x )
{
    auto evalType = getEvalType();

    NOMAD::EvalPoint X;
    size_t findPt = NOMAD::CacheBase::getInstance()->find(x,X);
    
    if ( findPt == 0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"The trial point is not in Cache" );
    
    if ( X.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK )
    {
        AddOutputDebug("The trial point x:" + x.display() + " has not eval ok.");
        return false;
    }
    
    std::pair<std::set<NOMAD::EvalPoint>::iterator,bool> ret_x = _nmY->insert ( X );
    
    // Insertion of x ok
    if ( ret_x.second )
    {
        AddOutputDebug("Insertion in NM simplex: " + X.display() ) ;
        
        // Remove yn
        std::set<NOMAD::EvalPoint>::iterator it = _nmY->end();
        --it;
        AddOutputDebug("Delete from NM simplex: " + (*it).display() ) ;
        
        // The simplex is unchanged ==> insertion failed
        if ( it == ret_x.first )
        {
            AddOutputDebug("Inserted point is last ==> insertion not successful, simplex unchanged. Let's continue." );
            
            // Remove the inserted point
            _nmY->erase ( it );
            
            return false;
        }
        // Remove the last point
        _nmY->erase ( it );
        
        // Update the simplex characteristics (volumes and diameter)
        updateYCharacteristics();
        
        // Update vectors Y0 and Yn
        bool success = makeListY0();
        if ( ! success )
            return false;
        
        success = makeListYn();
        if ( ! success )
            return false;
        
        displayYInfo();
        displayY0nInfo();
        
        if ( getRankDZ() != static_cast<int>( _nmY->size()) -1 )
        {
            AddOutputDebug("Rank of DZ=[(y1-y0) (y2-y0) ... (yn-y0)] != n. Y is not a valid simplex. Let's continue. ");
            return false;
        }
    }
    else
    {
        // Try to remove point from Y (in case it was there)
        _nmY->erase(ret_x.first);
        
        // Update simplex characteristics (volumes and diameter)
        updateYCharacteristics();
        
        AddOutputDebug("Cannot insert point in Y. Point possibly already in Y.");
        displayYInfo();
        displayY0nInfo();
        return false;
    }
    return true;
}

bool NOMAD::NMReflective::pointDominatesY0( const NOMAD::Point & xt ) const
{
    auto evalType = getEvalType();
    
    if ( _nmY0.size()==0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"Y0 is empty" );
    
    // Get the point from Cache
    NOMAD::EvalPoint Xt;
    size_t findPt = NOMAD::CacheBase::getInstance()->find(xt,Xt);

    if ( findPt == 0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"The trial point is not in Cache" );

    if ( Xt.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK )
    {
        AddOutputDebug("The trial point xt:" + Xt.display() + " has not eval ok.");
        return false;
    }
    
    if ( Xt.getEval(evalType) == nullptr )
    {
        AddOutputDebug("The trial point xt:" + Xt.display() + " has no eval.");
        return false;
    }
    
    for (auto itY0 : _nmY0)
    {
        
        // Xt < y ?
        if (Xt.dominates(*itY0, evalType))
        {
            return true;
        }
        
    }
    return false;
}


bool NOMAD::NMReflective::YnDominatesPoint( const NOMAD::Point & xt ) const
{
    auto evalType = getEvalType();

    if ( _nmYn.size()==0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ," Yn is empty" );
    
    // Get the point from Cache
    NOMAD::EvalPoint Xt;
    size_t findPt = NOMAD::CacheBase::getInstance()->find(xt,Xt);
    
    if ( findPt == 0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"The trial point is not in Cache" );
    
    if ( Xt.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK )
    {
        AddOutputDebug("The trial point xt:" + Xt.display() + " has not eval ok.");
        return false;
    }
    
    if ( Xt.getEval(evalType) == nullptr )
    {
        AddOutputDebug("The trial point xt:" + Xt.display() + " has no eval.");
        return true;
    }
    
    // Without constraints, Yn contains a single point
    int flag = 0;
    for (auto itYn : _nmYn)
    {
        // y < xt ?
        if (itYn->dominates(Xt, evalType))
        {
            // One point of Yn dominates x
            flag = 1;
            break;
        }
    }
    
    if ( flag == 1 )
        return true;
    
    auto itYn = _nmYn[_nmYn.size()-1]; // Consider yn
    // no point of Yn dominates xt --> check if h(yn) < h(xt) --> Yn dominates
 
    // Case with EB constraints and a point from Yn is infeasible
    if ( ! itYn->getH(evalType).isDefined() )
        return false;
        
    // Test also case where xt has no value for h (case with EB constraint)
    if ( ! Xt.getH(evalType).isDefined() || itYn->getH(evalType) < Xt.getH(evalType) )
        return true;
    
    return false;
}


bool NOMAD::NMReflective::pointDominatesPtsInY( const NOMAD::Point & xt , size_t nbPointsToDominate ) const
{
    auto evalType = getEvalType();

    // Get the point from Cache
    NOMAD::EvalPoint Xt;
    size_t findPt = NOMAD::CacheBase::getInstance()->find(xt,Xt);
    
    if ( findPt == 0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"The trial point is not in Cache" );
    
    if ( Xt.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK )
    {
        AddOutputDebug("The trial point xt:" + Xt.display() + " has not eval ok.");
        return false;
    }
    
    if ( Xt.getEval(evalType) == nullptr )
    {
        AddOutputDebug("The trial point xt:" + Xt.display() + " has no eval.");
        return false;
    }
    
    
    
    size_t nDominates = 0;
    
    std::set<NOMAD::EvalPoint>::const_iterator itY = _nmY->begin();
    
    while ( itY != _nmY->end() && nDominates < nbPointsToDominate)
    {
        // xt < y ?
        if (Xt.dominates(*itY, evalType))
        {
            nDominates++;
        }
        
        ++itY;
        
    }
    return ( nDominates == nbPointsToDominate );
}


/*------------------------------------------------------------------*/
/*  Make vector of undominated points from Nelder-Mead simplex Y0   */
/*------------------------------------------------------------------*/
bool NOMAD::NMReflective::makeListY0 ()
{
    auto evalType = getEvalType();
    
    _nmY0.clear();
    
    std::set<NOMAD::EvalPoint>::const_iterator itx;
    std::set<NOMAD::EvalPoint>::const_iterator ity = _nmY->begin();
    
    size_t maxSizeY0 = _nmY->size();
    
    // NOTE: Yes, this uses the Copy Constructor.
    // Code before: _nmY0.push_back( &(*ity) );
    // For an unknown reason, using the raw pointers
    // eventually creates a seg fault on my side.
    // Keeping the shared_ptr until we can debug this.
    _nmY0.push_back(std::make_shared<NOMAD::EvalPoint>(*ity));
    ++ity;
    while( ity != _nmY->end() && _nmY0.size() < maxSizeY0 )
    {
        int flag =0;
        const NOMAD::EvalPoint & y = (*ity);
        itx = _nmY->begin();
        
        while ( itx != _nmY->end() )
        {
            const NOMAD::EvalPoint & x = (*itx);
            
            // Case x dominates y --> y cannot be included in Y0
            if (x.dominates(y, evalType))
            {
                flag = 1;
                break;
            }
            
            ++itx;
            
        }
        
        // No point dominates y --> y is included in Y0
        if (flag == 0)
        {
            _nmY0.push_back(std::make_shared<NOMAD::EvalPoint>(y));
        }
        
        ++ity;
    }
    if ( _nmY0.size() > 0 )
        return true;
    else
        return false;
}

/*----------------------------------------------------------------*/
/*  Make vector of dominated points from Nelder-Mead simplex Yn   */
/*----------------------------------------------------------------*/
bool NOMAD::NMReflective::makeListYn ()
{
    auto evalType = getEvalType();
    _nmYn.clear();
    
    std::set<NOMAD::EvalPoint>::const_iterator itx;
    std::set<NOMAD::EvalPoint>::const_iterator ity = _nmY->begin();
    
    while( ity != _nmY->end() )
    {
        int flag =0;
        const NOMAD::EvalPoint & y = (*ity);
        itx=_nmY->begin();
        while ( itx != _nmY->end() )
        {
            
            const NOMAD::EvalPoint & x = (*itx);
            
            // Case y dominates x --> y cannot be included in Yn
            if (y.dominates(x, evalType))
            {
                flag = 1;
                break;
            }
            ++itx;
            
        }
        
        // y dominates no points --> y is included in Y0
        if (0 == flag)
        {
            _nmYn.push_back(std::make_shared<NOMAD::EvalPoint>(y));
        }
        
        ++ity;
    }
    if ( _nmYn.size() > 0 )
        return true;
    else
        return false;
    
}


void NOMAD::NMReflective::displayY0nInfo() const
{
    
    AddOutputInfo("Number of points in Y0: " + std::to_string ( _nmY0.size()) );
    AddOutputInfo("Number of points in Yn: " + std::to_string ( _nmYn.size()) );
    
    NOMAD::OutputInfo dbgInfo("Display Y0 and Yn info", "The vector Y0 contains:", NOMAD::OutputLevel::LEVEL_DEBUG);
    
    for (auto y0 : _nmY0)
    {
        dbgInfo.addMsg(y0->display());
    }
    
    dbgInfo.addMsg ( "The vector Yn contains: ");
    for (auto yn : _nmYn)
    {
        dbgInfo.addMsg(yn->display());
    }

    NOMAD::OutputQueue::Add(std::move(dbgInfo));
    NOMAD::OutputQueue::Flush();
}


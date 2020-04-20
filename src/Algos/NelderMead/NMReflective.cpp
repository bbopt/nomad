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
    
#include "../../Algos/EvcInterface.hpp"
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
    
    // Update reference point according to the current step (if a trial point was evaluated).
    if (getTrialPoints().size() > 0)
    {
        // Update the eval points by step type
        switch ( _currentStepType )
        {
            case NMStepType::REFLECT:
                _xr = *(getTrialPoints().begin());
                break;
            case NMStepType::EXPAND:
                _xe = *(getTrialPoints().begin());
                break;
            case NMStepType::OUTSIDE_CONTRACTION:
                _xoc = *(getTrialPoints().begin());
                break;
            case NMStepType::INSIDE_CONTRACTION:
                _xic = *(getTrialPoints().begin());
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


void NOMAD::NMReflective::generateTrialPoints()
{
    // The pb params handle only variables (fixed variables are not considered)
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    // The size of simplex must be large enough to perform reflexion
    size_t minYSize = n + 1;
    size_t YSize = _nmY->size();
    if ( YSize < minYSize )
    {
        OUTPUT_INFO_START
        AddOutputInfo("Not enough points in the simplex to generate points for " + _name +" (delta=" + _delta.tostring() +")");
        OUTPUT_INFO_END
        return;
    }
    
    OUTPUT_INFO_START
    AddOutputInfo("Generate point for " + _name +" (delta=" + _delta.tostring()+")");
    OUTPUT_INFO_END
    
    // Clear the trial points.
    clearTrialPoints();

    // Determine the centroid of Y
    std::set<NOMAD::EvalPoint>::iterator it;
    std::set<NOMAD::EvalPoint>::iterator itBeforeEnd = _nmY->end();
    --itBeforeEnd;
    int i=0;
    NOMAD::Point yc(n,0);
    for (it = _nmY->begin() ; it != itBeforeEnd ; ++it, i++)
    {
        OUTPUT_INFO_START
        AddOutputInfo("y" + std::to_string(i) + ": " + (*it).display() );
        OUTPUT_INFO_END
        for (size_t k = 0 ; k < n ; ++k )
        {
            yc[k] += (*it)[k];
        }
    }
    yc *= 1.0/n;
    
    const NOMAD::Point & yn = *itBeforeEnd;
    OUTPUT_INFO_START
    AddOutputInfo("yn: " + yn.display() );
    AddOutputInfo("yc: " + yc.display() );
    OUTPUT_INFO_END
    
    NOMAD::Point d(n,0);
    for (size_t k = 0 ; k < n ; ++k )
    {
        d[k] = yc[k]-yn[k];
    }
    d *= _delta;

    // Creation of point
    NOMAD::EvalPoint xt(n);
    for (size_t k = 0 ; k < n ; ++k )
    {
        xt[k] = yc[k] + d[k];
    }

    auto lb = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto ub = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    
    if (snapPointToBoundsAndProjectOnMesh(xt, lb, ub))
    {
        bool inserted = insertTrialPoint(xt);
        
        OUTPUT_INFO_START
        std::string s = "xt:";
        s += (inserted) ? " " : " not inserted: ";
        s += xt.display();
        AddOutputInfo(s);
        OUTPUT_INFO_END
    }
    
    if (*xt.getX() == yn)
    {
        clearTrialPoints();
        OUTPUT_INFO_START
        AddOutputInfo("No valid point generated: too close to yn.");
        OUTPUT_INFO_END
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
    
    if ( getNbEvalPointsThatNeededEval() == 0 )
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("Cannot create a proper reflect point xr. Next, perform Inside Contraction.");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::INSIDE_CONTRACTION;
        return;
    }
    
    if (!_xr.NOMAD::Point::isDefined())
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The reflect point xr is not defined. Stop NM (no shrink).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    if (pointDominatesY0(_xr))
    {
        // In NM-Mads paper: x_r belongs to the expansion zone. Next, expansion.
        OUTPUT_DEBUG_START
        AddOutputDebug("The reflect point xr: " + _xr.display() +" dominates Y0. Next perform expansion.");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::EXPAND;
    }
    else if ( YnDominatesPoint( _xr ) )
    {
        // In NM-Mads paper: x_r belongs to the inside contraction zone. Next, perform inside contraction.
        OUTPUT_DEBUG_START
        AddOutputDebug("The reflect point xr: " + _xr.display() + " is dominated by Yn. Next, perform Inside Contraction.");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::INSIDE_CONTRACTION;
    }
    else if ( pointDominatesPtsInY( _xr , 2 ) )
    {
        // In NM-Mads paper: x_r belongs to the reflection zone. Next, insert into simplex (3), update Y, Y0, Yn (1) and continue (2)
        OUTPUT_DEBUG_START
        AddOutputDebug("The reflect point xr: "  + _xr.display() + " dominates at least 2 points of Y.");
        OUTPUT_DEBUG_END
        
        _currentStepType = NMStepType::INSERT_IN_Y;
        if ( insertInY ( _xr ) )
        {
            OUTPUT_DEBUG_START
            AddOutputDebug("Insertion in Y is successful. NM iteration completed (no shrink).");
            OUTPUT_DEBUG_END
            _nextStepType = NMStepType::CONTINUE;
        }
        else
        {
            OUTPUT_DEBUG_START
            AddOutputDebug(" Cannot insert xr in Y. Perform shrink (if available).");
            OUTPUT_DEBUG_END
            _nextStepType = NMStepType::SHRINK;
        }
    }
    else if ( pointDominatesPtsInY( _xr , 1 ) || pointDominatesPtsInY( _xr , 0 ) )
    {
        // In NM-Mads paper: x_r belongs to the outside reflection zone. Next, perform outside contraction.
        OUTPUT_DEBUG_START
        AddOutputDebug("The reflect point xr: " + _xr.display() + " dominates 1 or 0 point of Y. Next, perform Outside Contraction.");
        OUTPUT_DEBUG_END
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
    
    if (!_xe.NOMAD::Point::isDefined())
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The expansion point xe is not defined. Stop NM (no shrink).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    if (!_xr.NOMAD::Point::isDefined())
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The reflect point xr is not defined. Stop NM (no shrink).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    _currentStepType = NMStepType::INSERT_IN_Y;
    if ( insertInYBest( _xr , _xe ) )
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("Insert in Y the best of xr and xe. NM iteration completed (no shrink). ");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
    }
    else
    {
        // No point inserted.
        OUTPUT_DEBUG_START
        AddOutputDebug("The insertion in Y of the best of xr and xe dit not maintain a proper Y. Perform shrink (if available).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::SHRINK;
    }
}


void NOMAD::NMReflective::setAfterOutsideContract ( void )
{
    
    // Outside contraction follows Reflect.
    // In NM-Mads paper: x_r belongs to the outside contraction zone. x_oc has been evaluated. The best point between x_r and x_oc is inserted in Y. NM iteration is completed.
    if ( _currentStepType != NMStepType::OUTSIDE_CONTRACTION )
        throw NOMAD::Exception(__FILE__,__LINE__,"The current step type should be OUTSIDE_CONTRACTION.");
    
    if (!_xr.NOMAD::Point::isDefined())
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The reflect point xr is not defined. Stop NM (no shrink).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    if (getNbEvalPointsThatNeededEval() == 0)
    {
        if ( insertInY( _xr ) )
        {
            OUTPUT_DEBUG_START
            AddOutputDebug("Reflect point xr is successfully inserted in Y. Next perform Reflect.");
            OUTPUT_DEBUG_END
            
            _nextStepType = NMStepType::REFLECT;
        }
        else
            setStopReason();
        
        return;
    }
    
    
    if (!_xoc.NOMAD::Point::isDefined())
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The outside contraction point xoc is not defined. Stop NM (no shrink).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }

    _currentStepType = NMStepType::INSERT_IN_Y;
    if ( insertInYBest( _xr , _xoc ) )
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The insertion of the best of xr and xoc in Y is valid. NM iteration completed." );
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
        return;
    }
    else
    {
        OUTPUT_DEBUG_START
        AddOutputDebug( "The insertion of the best of xr and xoc in Y did not maintain a proper Y. Perform shrink (if available).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::SHRINK;
    }
    
}


void NOMAD::NMReflective::setAfterInsideContract ( void )
{
    
    // Inside contraction follows Reflect.
    // In NM-Mads paper: x_r belongs to the inside contraction zone. x_ic has been evaluated. If x_ic belongs to the inside contraction zone, stop NM. Otherwise insert x_ic in Y. NM iteration is completed.
    
    if ( _currentStepType != NMStepType::INSIDE_CONTRACTION )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot set step after inside contraction if x_ic is not defined.");
    }

    if (!_xic.NOMAD::Point::isDefined())
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The inside contraction point xic is not defined. Stop NM (no shrink).");
        OUTPUT_DEBUG_END
        _nextStepType = NMStepType::CONTINUE;
        setStopReason( );
        return;
    }
    
    if ( YnDominatesPoint( _xic ) )
    {
        _nextStepType = NMStepType::SHRINK;
        OUTPUT_DEBUG_START
        AddOutputDebug("Yn dominates xic: " + _xic.display() + " Next, perform Shrink.");
        OUTPUT_DEBUG_END
        return;
    }
    else
    {
        
        OUTPUT_DEBUG_START
        AddOutputDebug("The inside contraction point xic:" + _xic.display() + " is not dominated by Yn. Insert x_ic in Y." );
        OUTPUT_DEBUG_END
        
        _currentStepType = NMStepType::INSERT_IN_Y;
        if ( insertInY( _xic ) )
        {
            OUTPUT_DEBUG_START
            AddOutputDebug("Insertion in Y is successful. NM iteration completed (no shrink).");
            OUTPUT_DEBUG_END
            _nextStepType = NMStepType::CONTINUE;
        }
        else
        {
            // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
            OUTPUT_DEBUG_START
            AddOutputDebug("Cannot insert xic in Y. Perform Shrink if available." );
            OUTPUT_DEBUG_END
            _nextStepType = NMStepType::SHRINK;
        }
    }
}


bool NOMAD::NMReflective::insertInYBest(const NOMAD::EvalPoint& x1, const NOMAD::EvalPoint& x2)
{
    auto evalType = getEvalType();

    if (nullptr == x1.getEval(evalType))
    {
        throw NOMAD::Exception (__FILE__, __LINE__, "The trial point x1 does not have an evaluation: " + x1.display());
    }
    
    bool x1EvalOK = (x1.getEvalStatus(evalType) == NOMAD::EvalStatusType::EVAL_OK);
    
    if (nullptr == x2.getEval(evalType))
    {
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"The trial point x2 does not have an evaluation: " + x2.display());
    }
    
    bool x2EvalOK = (x2.getEvalStatus(evalType) == NOMAD::EvalStatusType::EVAL_OK);
    
    std::set<NOMAD::EvalPoint>::iterator itYn = _nmY->end();
    --itYn;
    
    OUTPUT_DEBUG_START
    AddOutputDebug("Insertion/Deletion of points in Y: ");
    OUTPUT_DEBUG_END
    
    if (x1EvalOK)
    {
        std::pair<std::set<NOMAD::EvalPoint>::iterator,bool> ret_x1 = _nmY->insert(x1);
        
        // Insertion of x1 ok
        if (ret_x1.second)
        {
            OUTPUT_DEBUG_START
            AddOutputDebug("Insertion in Y: ") ;
            AddOutputDebug(x1.display());
            OUTPUT_DEBUG_END
            
            std::pair<std::set<NOMAD::EvalPoint>::iterator,bool> ret_x2;
            
            bool same_x1_x2 = false;
            
            //
            if (!x2EvalOK)
            {
                if (std::distance(_nmY->begin(), ret_x1.first) > std::distance(_nmY->begin(), itYn))
                {
                    _nmY->erase( ret_x1.first );
                    OUTPUT_DEBUG_START
                    AddOutputDebug("Cannot perform insertion in Y. x1 is dominated by yn and x2 eval status is not ok. Delete x1.");
                    OUTPUT_DEBUG_END
                    return false;
                }
            }
            
            // Try to insert x2 only if x1 != x2
            if (x2EvalOK && x1 != x2)
            {
                ret_x2 = _nmY->insert(x2);
            }
            else
            {
                ret_x2.second = true;
                ret_x2.first = ret_x1.first;
                same_x1_x2 = true;
            }
            
            // Insertion of x2 ok
            if (ret_x2.second)
            {
                OUTPUT_DEBUG_START
                AddOutputDebug("Insertion in Y: ") ;
                AddOutputDebug(same_x1_x2 ? x1.display() : x2.display());
                OUTPUT_DEBUG_END
                
                
                // if x1 and x2 are after yn --> remove both x1 and x2 ==> insertion failed
                if (   std::distance(_nmY->begin(), ret_x1.first) > std::distance(_nmY->begin(), itYn)
                    && std::distance(_nmY->begin(), ret_x2.first) > std::distance(_nmY->begin(), itYn))
                {
                    _nmY->erase( ret_x1.first );
                    if (!same_x1_x2)
                    {
                        _nmY->erase( ret_x2.first );
                    }
                    
                    OUTPUT_DEBUG_START
                    AddOutputDebug("Cannot perform insertion because both points are dominated by yn");
                    OUTPUT_DEBUG_END
                    return false;
                }
                
                
                // if x1 after x2 --> keep x2 remove x1
                if (std::distance(_nmY->begin(), ret_x1.first) > std::distance(_nmY->begin(), ret_x2.first))
                {
                    OUTPUT_DEBUG_START
                    AddOutputDebug("Delete from Y: ") ;
                    AddOutputDebug(x1.display());
                    OUTPUT_DEBUG_END
                    _nmY->erase(ret_x1.first);
                }
                else // if x2 after x1 --> keep x1 remove x2
                {
                    if (!same_x1_x2)
                    {
                        OUTPUT_DEBUG_START
                        AddOutputDebug("Delete from Y: ") ;
                        AddOutputDebug(x2.display());
                        OUTPUT_DEBUG_END
                        _nmY->erase(ret_x2.first);
                    }
                }
                
                // Remove yn
                std::set<NOMAD::EvalPoint>::iterator it = _nmY->end();
                --it;
                OUTPUT_DEBUG_START
                AddOutputDebug("Delete yn from Y: " + (*it).display());
                OUTPUT_DEBUG_END
                _nmY->erase ( it );
                
                // Update simplex characteristics (volumes and diameter)
                updateYCharacteristics();
                
                // Update vectors Y0 and Yn
                bool success = makeListY0();
                if ( !success )
                {
                    OUTPUT_DEBUG_START
                    AddOutputDebug("Cannot make Y0. Let's continue." );
                    OUTPUT_DEBUG_END
                    return false;
                }
                success = makeListYn();
                if ( !success )
                {
                    OUTPUT_DEBUG_START
                    AddOutputDebug("Cannot make Yn. Let's continue" );
                    OUTPUT_DEBUG_END
                    return false;
                }
                
                
                OUTPUT_DEBUG_START
                AddOutputDebug("After insertion best");
                OUTPUT_DEBUG_END
                displayYInfo();
                displayY0nInfo();
                
                if ( getRankDZ() != static_cast<int>(_nmY->size()) -1 )
                {
                    OUTPUT_DEBUG_START
                    AddOutputInfo("Rank of DZ=[(y1-y0) (y2-y0) ... (yn-y0)] != n. Y is not a valid simplex. Let's continue. ");
                    OUTPUT_DEBUG_END
                    return false;
                }
            }
            else
            {
                // Try to remove point from Y (in case it was there)
                _nmY->erase(ret_x2.first);
                
                // Update simplex characteristics (volumes and diameter)
                updateYCharacteristics();
                
                OUTPUT_DEBUG_START
                AddOutputDebug("After insertion cancelled");
                OUTPUT_DEBUG_END
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
            
            OUTPUT_DEBUG_START
            AddOutputDebug("After insertion cancelled");
            OUTPUT_DEBUG_END
            displayYInfo();
            displayY0nInfo();
            return false;
        }
    }
    else
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("Cannot insertion in Y because eval not OK: ") ;
        AddOutputDebug(x1.display());
        OUTPUT_DEBUG_END
    }

    return true;
}


bool NOMAD::NMReflective::insertInY(const NOMAD::EvalPoint& x)
{
    auto evalType = getEvalType();

    if (x.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK)
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The trial point x:" + x.display() + " is not eval ok.");
        OUTPUT_DEBUG_END
        return false;
    }
    
    std::pair<std::set<NOMAD::EvalPoint>::iterator,bool> ret_x = _nmY->insert(x);
    
    // Insertion of x ok
    if ( ret_x.second )
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("Insertion in NM simplex: " + x.display()) ;
        OUTPUT_DEBUG_END
        
        // Remove yn
        std::set<NOMAD::EvalPoint>::iterator it = _nmY->end();
        --it;
        OUTPUT_DEBUG_START
        AddOutputDebug("Delete from NM simplex: " + (*it).display() ) ;
        OUTPUT_DEBUG_END
        
        // The simplex is unchanged ==> insertion failed
        if ( it == ret_x.first )
        {
            OUTPUT_DEBUG_START
            AddOutputDebug("Inserted point is last ==> insertion not successful, simplex unchanged. Let's continue." );
            OUTPUT_DEBUG_END
            
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
            OUTPUT_DEBUG_START
            AddOutputDebug("Rank of DZ=[(y1-y0) (y2-y0) ... (yn-y0)] != n. Y is not a valid simplex. Let's continue. ");
            OUTPUT_DEBUG_END
            return false;
        }
    }
    else
    {
        // Try to remove point from Y (in case it was there)
        _nmY->erase(ret_x.first);
        
        // Update simplex characteristics (volumes and diameter)
        updateYCharacteristics();
        
        OUTPUT_DEBUG_START
        AddOutputDebug("Cannot insert point in Y. Point possibly already in Y.");
        OUTPUT_DEBUG_END
        displayYInfo();
        displayY0nInfo();
        return false;
    }
    return true;
}


bool NOMAD::NMReflective::pointDominatesY0( const NOMAD::EvalPoint & xt ) const
{
    auto evalType = getEvalType();
    std::string s;
    
    if ( _nmY0.size()==0 )
        throw NOMAD::Exception ( __FILE__ , __LINE__ ,"Y0 is empty" );
    
    if (nullptr == xt.getEval(evalType))
    {
        s = "The evaluation for trial point xt = " + xt.display() + " was not found";
        throw NOMAD::Exception ( __FILE__ , __LINE__ , s);
    }

    if (xt.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK)
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The trial point xt: " + xt.display() + " is not eval ok.");
        OUTPUT_DEBUG_END
        return false;
    }
    
    for (auto evalPointY0 : _nmY0)
    {
        
        // xt < y ?
        if (xt.dominates(evalPointY0, evalType))
        {
            return true;
        }
        
    }
    return false;
}


bool NOMAD::NMReflective::YnDominatesPoint(const NOMAD::EvalPoint& xt) const
{
    auto evalType = getEvalType();

    if (_nmYn.size() == 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, " Yn is empty");
    }
    
    if (nullptr == xt.getEval(evalType))
    {
        throw NOMAD::Exception (__FILE__, __LINE__, "No evaluation for trial point " + xt.display());
    }
    
    if (xt.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK)
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The trial point xt: " + xt.display() + " is not eval ok.");
        OUTPUT_DEBUG_END
        return false;
    }
    
    // Without constraints, Yn contains a single point
    int flag = 0;
    for (auto evalPointYn : _nmYn)
    {
        // y < xt ?
        if (evalPointYn.dominates(xt, evalType))
        {
            // One point of Yn dominates x
            flag = 1;
            break;
        }
    }
    
    if (flag == 1)
    {
        return true;
    }
    
    auto evalPointYn = _nmYn[_nmYn.size()-1]; // Consider yn
    // no point of Yn dominates xt --> check if h(yn) < h(xt) --> Yn dominates
 
    // Case with EB constraints and a point from Yn is infeasible
    if (!evalPointYn.getH(evalType).isDefined())
    {
        return false;
    }
        
    // Test also case where xt has no value for h (case with EB constraint)
    if (!xt.getH(evalType).isDefined() || evalPointYn.getH(evalType) < xt.getH(evalType))
    {
        return true;
    }
    
    return false;
}


bool NOMAD::NMReflective::pointDominatesPtsInY(const NOMAD::EvalPoint& xt, size_t nbPointsToDominate) const
{
    auto evalType = getEvalType();

    if (nullptr == xt.getEval(evalType))
    {
        throw NOMAD::Exception (__FILE__, __LINE__, "No evaluation for trial point " + xt.display());
    }
    
    if (xt.getEvalStatus(evalType) != NOMAD::EvalStatusType::EVAL_OK)
    {
        OUTPUT_DEBUG_START
        AddOutputDebug("The trial point xt: " + xt.display() + " is not eval ok.");
        OUTPUT_DEBUG_END
        return false;
    }
    
    size_t nDominates = 0;
    
    std::set<NOMAD::EvalPoint>::const_iterator itY = _nmY->begin();
    
    while (itY != _nmY->end() && nDominates < nbPointsToDominate)
    {
        // xt < y ?
        if (xt.dominates(*itY, evalType))
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
    
    _nmY0.push_back(*ity);
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
            _nmY0.push_back(y);
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
    
    auto ity = _nmY->begin();
    
    while (ity != _nmY->end())
    {
        bool flag = false;
        const NOMAD::EvalPoint& y = (*ity);
        auto itx = _nmY->begin();
        while (itx != _nmY->end())
        {
            const NOMAD::EvalPoint & x = (*itx);
            
            // Case y dominates x --> y cannot be included in Yn
            if (y.dominates(x, evalType))
            {
                flag = true;
                break;
            }
            ++itx;
            
        }
        
        // y dominates no points --> y is included in Y0
        if (!flag)
        {
            _nmYn.push_back(y);
        }
        
        ++ity;
    }

    return (_nmYn.size() > 0);
}


void NOMAD::NMReflective::displayY0nInfo() const
{
    
    OUTPUT_INFO_START
    AddOutputInfo("Number of points in Y0: " + std::to_string ( _nmY0.size()) );
    AddOutputInfo("Number of points in Yn: " + std::to_string ( _nmYn.size()) );
    OUTPUT_INFO_END
    
    OUTPUT_DEBUG_START
    NOMAD::OutputInfo dbgInfo("Display Y0 and Yn info", "The vector Y0 contains:", NOMAD::OutputLevel::LEVEL_DEBUG);
    
    for (auto y0 : _nmY0)
    {
        dbgInfo.addMsg(y0.display());
    }
    
    dbgInfo.addMsg ( "The vector Yn contains: ");
    for (auto yn : _nmYn)
    {
        dbgInfo.addMsg(yn.display());
    }

    NOMAD::OutputQueue::Add(std::move(dbgInfo));
    NOMAD::OutputQueue::Flush();
    OUTPUT_DEBUG_END
}


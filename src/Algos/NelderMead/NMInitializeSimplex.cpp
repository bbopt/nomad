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

#include "../../Algos/NelderMead/NMInitializeSimplex.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"


void NOMAD::NMInitializeSimplex::init()
{
    _name = getAlgoName() + "Initialize Simplex";
    
    verifyParentNotNull();

}

bool NOMAD::NMInitializeSimplex::runImp()
{
    
    if ( _nmY == nullptr )
        NOMAD::Exception(__FILE__, __LINE__, "The simplex is not defined.");
    
    // create a simplex from EvalPoints in Cache
    if ( _nmY->size() == 0 )
        return createSimplexFromCache();
    else
    {
        AddOutputInfo("Simplex already initialized: " + std::to_string(_nmY->size()) + " points");
        return true;
    }
}


/*----------------------------------------------------------------------------------*/
/*  Create initial sets of points for Nelder-Mead within a radius of current best   */
/*----------------------------------------------------------------------------------*/
bool NOMAD::NMInitializeSimplex::createSimplexFromCache ( )
{
    auto evalType = getEvalType();

    auto iter = dynamic_cast<const NOMAD::NMIteration*>( NOMAD::Step::_parentStep );
    if ( nullptr == iter )
        NOMAD::Exception(__FILE__, __LINE__, "The simplex initialization must have a NMIteration Step has parent");
    
    const std::shared_ptr<NOMAD::EvalPoint> centerPt = iter->getFrameCenter();
    // Use center point of iteration, otherwise
    if ( centerPt == nullptr )
            NOMAD::Exception(__FILE__, __LINE__, "A center point must be defined.");
    
    // Clear the list of NM points
    _nmY->clear();
    
    // The pb params handle only variables (no need to consider fixed variables)
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    size_t minYSize = n + 1;
    
    // compute the include radius: points in Y must be at
    // a max distance of centerPt
    auto includeLength = _runParams->getAttributeValue<NOMAD::Double>("NM_SIMPLEX_INCLUDE_LENGTH");
    NOMAD::ArrayOfDouble includeRectangle(n, includeLength ) ;
    
    NOMAD::OutputInfo dbgInfo(_name,"Insertion of potential points to include in initial Y: ", NOMAD::OutputLevel::LEVEL_DEBUG);
    
    // If a mesh and include factor are supplied: the max distance is include factor times Delta
    // Else we use include length
    auto mesh = iter->getMesh(); // The mesh can be null
    auto includeFactor = _runParams->getAttributeValue<size_t>("NM_SIMPLEX_INCLUDE_FACTOR");
    if ( mesh != nullptr && includeFactor > 0 && includeFactor < NOMAD::INF_SIZE_T )
    {
        // The max distance is NM_search_include_factor times Delta
        
        includeRectangle = mesh->getDeltaFrameSize() ;
        
        if ( ! includeRectangle.isDefined() )
            NOMAD::Exception(__FILE__, __LINE__, "The frame size is not defined.");
        
        includeRectangle *= includeFactor ;
        
        dbgInfo.addMsg("The include rectangle: " + includeRectangle.display() );
    }
    if ( includeRectangle.max() == 0 )
        NOMAD::Exception(__FILE__, __LINE__, "The include rectangle has no volume");
    
    const NOMAD::ArrayOfDouble & bbo  = centerPt->getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble();
    size_t m = bbo.size();
    
    // The set of points initially included
    NOMAD::NMSimplexEvalPointSet T;
    
    // browse the cache:
    auto cache = NOMAD::CacheBase::getInstance().get();
    std::vector<NOMAD::EvalPoint> evalpointlist;
    
    cache->getAllPoints( evalpointlist );

    // variables used to limit display
    const size_t maxPointsToDisplay = 4;
    size_t nbPoints = 0;
    for ( const auto & cur : evalpointlist )
    {
        if ( cur.getEvalStatus(evalType) == NOMAD::EvalStatusType::EVAL_OK &&
            cur.getX()->size() == n             )
        {
            const NOMAD::ArrayOfDouble & bboCur  = cur.getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble();
            
            if ( bboCur.isComplete() && checkOutputs(bboCur, static_cast<int>(m) ) )
            {
                // the center point has been found and put in list
                if ( *centerPt == cur )
                {
                    NOMAD::EvalPoint Y ( cur );
                    T.insert ( Y );
                    if (nbPoints < maxPointsToDisplay)
                    {
                        dbgInfo.addMsg(Y.display());
                        nbPoints++;
                    }
                }
                // other points must be within the include region:
                else
                {
                    bool include = true;
                    for (size_t i = 0 ; i < n ; i++ )
                    {
                        NOMAD::Double val = cur[i] - (*centerPt)[i];
                        if ( val.abs() > includeRectangle[i]  )
                        {
                            include = false;
                            break;
                        }
                    }
                    if ( include )
                    {
                        
                        // Make sure to evaluate f or h for points in cache (important if cache is loaded from file)
                        // TODO
                        
                        NOMAD::EvalPoint Y ( cur );
                        std::pair<NMSimplexEvalPointSetIterator,bool> ret = T.insert ( Y );
                    
                        
                        if ( ! ret.second )
                        {
                            dbgInfo.addMsg("Cannot insert a point in Y: " + Y.display() );
                            break;
                        }
                        else
                        {
                            if (nbPoints <= maxPointsToDisplay)
                            {
                                dbgInfo.addMsg( ((nbPoints < maxPointsToDisplay) ? Y.display() : "...") );
                                nbPoints++;
                            }
                        }
                    }
                }
            }
        }
    }
    NOMAD::OutputQueue::Add(std::move(dbgInfo));

    AddOutputInfo("Number of potential points to include in initial Y: " + std::to_string(T.size()) );
    
    //
    auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( getAllStopReasons() );
    
    // not enough points (this is not counted as an error):
    if ( T.size() < minYSize )
    {
        nmStopReason->set ( NMStopType::INITIAL_FAILED );
        AddOutputInfo("Stop NM because not enough points in Y.");
        return false;
    }
    

    
    // Add points in simplex to obtain dim = n+1 and simplex affinely independant
    NOMAD::NMSimplexEvalPointSetIterator itT = T.begin();
    
    // For debugging
    NOMAD::OutputInfo dbgInfo2(_name,"Proceed to simplex creation", NOMAD::OutputLevel::LEVEL_DEBUG);
    
    // First point is always added
    _nmY->insert ( *itT );
    
    dbgInfo2.addMsg("k=0: Point z0:" + (*itT).display() ) ;
    dbgInfo2.addMsg(" ---> z0 KEPT in Y ");;
    
    int count_feasible = 0;
    
    if ( (*itT).getH(evalType).isDefined()
        && (*itT).isFeasible(evalType) )
        count_feasible = 1 ;
    
    itT++;
    int k=1;
    while ( itT != T.end() && ! nmStopReason->checkTerminate() )
    {
        
        if ( _nmY->size() == minYSize )
            break;
        
        dbgInfo2.addMsg("k=" + std::to_string(k) +": Point zk:" + (*itT).display()) ;
        
        std::pair<NMSimplexEvalPointSetIterator,bool> ret = _nmY->insert ( *itT );
        
        if ( ! ret.second )
        {
            nmStopReason->set ( NMStopType::INITIAL_FAILED );
            AddOutputInfo("Stop NM because cannot insert a point in Y.");
            break;
        }
        
        int rank = getRankDZ();
        
        if ( rank <= 0 )
        {
            nmStopReason->set ( NMStopType::INITIAL_FAILED );
            AddOutputInfo("Cannot get rank of DZ=[(y1-y0 (y2-y0) ...(yk-y0)].");
            break;
        }
        
        // Erase last point or not
        if ( rank != k )
        {
            dbgInfo2.addMsg(" ---> zk NOT KEPT in Y ");
            _nmY->erase( *itT );
            itT++;
        }
        else
        {
            dbgInfo2.addMsg( " ---> zk KEPT in Y " );
            if ( (*itT).isFeasible(evalType) )
                count_feasible++;
            k++;
            itT++;
        }
        
    }
    NOMAD::OutputQueue::Add(std::move(dbgInfo2));
    
    // not enough points or insufficient rank of simplex (this is not counted as an error):
    // Erase point
    if ( _nmY->size() < minYSize )
    {
        nmStopReason->set ( NMStopType::INITIAL_FAILED );
        AddOutputInfo( "Stop NM because not enough simplex points in Y." );
        return false;
    }
    if ( getRankDZ() != (int)n )
    {
        nmStopReason->set ( NMStopType::INITIAL_FAILED );
        AddOutputInfo( "Stop NM because rank of Y < n." );
        return false;
    }
    
    
    
    // Update simplex characteristics (volumes and diameter)
    updateYCharacteristics();
    
    displayYInfo();
    
    return true;
    
}


/*---------------------------------------------------------*/
/*  check evaluation point outputs before the integration  */
/*  into NM set (private)                                  */
/*---------------------------------------------------------*/
bool NOMAD::NMInitializeSimplex::checkOutputs ( const NOMAD::ArrayOfDouble & bbo , int m ) const
{
    
    if ( bbo.size() != (size_t)m )
        return false;
    
    for ( int i = 0 ; i < m ; ++i )
        if ( !bbo[i].isDefined() )
            return false;
    
    return true;
}

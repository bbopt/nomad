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

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/NelderMead/NMInitializeSimplex.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::NMInitializeSimplex::init()
{
    setStepType(NOMAD::StepType::NM_INITIALIZE_SIMPLEX);

    verifyParentNotNull();

}

bool NOMAD::NMInitializeSimplex::runImp()
{
    bool simplexCreated = false;
    if (nullptr == _nmY)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "The simplex is not defined.");
    }

    // Create a simplex from EvalPoints in Cache or in Barrier
    if (_nmY->empty())
    {
        simplexCreated = createSimplex();
    }
    else
    {
        OUTPUT_INFO_START
        AddOutputInfo("Simplex already initialized: " + std::to_string(_nmY->size()) + " points");
        OUTPUT_INFO_END
        simplexCreated = true;
    }

    return simplexCreated;
}


/*----------------------------------------------------------------------------------*/
/*  Create initial sets of points for Nelder-Mead within a radius of current best   */
/*----------------------------------------------------------------------------------*/
bool NOMAD::NMInitializeSimplex::createSimplex()
{
    auto evalType = NOMAD::EvcInterface::getEvaluatorControl()->getEvalType();
    auto computeType = NOMAD::EvcInterface::getEvaluatorControl()->getComputeType();

    auto iter = dynamic_cast<const NOMAD::NMIteration*>( NOMAD::Step::_parentStep );
    if (nullptr == iter)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "The simplex initialization must have a NMIteration Step as parent");
    }

    std::vector<NOMAD::EvalPoint> evalPointList;
    std::shared_ptr<NOMAD::EvalPoint> centerPt = nullptr;
    auto barrier = getMegaIterationBarrier();
    if (nullptr != barrier)
    {
        evalPointList = barrier->getAllPoints();
        centerPt = std::make_shared<NOMAD::EvalPoint>(evalPointList[0]);
    }
    if (nullptr == centerPt)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "A center point must be defined.");
    }
    if ( ! centerPt->getEval(evalType)->isBBOutputComplete() )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "A center point must have all output defined for functions (OBJ, PB, EB).");
    }

    // Clear the list of NM points
    _nmY->clear();

    // The pb params handle only variables (no need to consider fixed variables)
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    size_t minYSize = n + 1;

    // compute the include radius: points in Y must be at
    // a max distance of centerPt
    auto includeLength = _runParams->getAttributeValue<NOMAD::Double>("NM_SIMPLEX_INCLUDE_LENGTH");
    NOMAD::ArrayOfDouble includeRectangle(n, includeLength ) ;

    NOMAD::OutputInfo dbgInfo(getName(),"Insertion of potential points to include in initial Y: ", NOMAD::OutputLevel::LEVEL_DEBUG);

    // If a mesh and include factor are supplied: the max distance is include factor times Delta
    // Else we use include length
    auto mesh = iter->getMesh(); // The mesh can be null
    auto includeFactor = _runParams->getAttributeValue<size_t>("NM_SIMPLEX_INCLUDE_FACTOR");
    if ( mesh != nullptr && includeFactor > 0 && includeFactor < NOMAD::INF_SIZE_T )
    {
        // The max distance is NM_search_include_factor times Delta

        includeRectangle = mesh->getDeltaFrameSize() ;

        if ( ! includeRectangle.isDefined() )
            throw NOMAD::Exception(__FILE__, __LINE__, "The frame size is not defined.");

        includeRectangle *= (double)includeFactor ;

        OUTPUT_DEBUG_START
        dbgInfo.addMsg("The include rectangle: " + includeRectangle.display() );
        OUTPUT_DEBUG_END
    }
    if ( includeRectangle.max() == 0 )
        throw NOMAD::Exception(__FILE__, __LINE__, "The include rectangle has no volume");

    // The set of points initially included
    NOMAD::NMSimplexEvalPointSet T;

    if (NOMAD::EvcInterface::getEvaluatorControl()->getUseCache())
    {
        // browse the cache:
        NOMAD::CacheInterface cacheInterface(this);
        cacheInterface.getAllPoints(evalPointList);
    }
    else
    {
        if (nullptr != barrier)
        {
            evalPointList = barrier->getAllPoints();
        }
    }

    // variables used to limit display
    const size_t maxPointsToDisplay = 4;
    size_t nbPoints = 0;
    for ( const auto & cur : evalPointList )
    {
        if ( cur.getEvalStatus(evalType) == NOMAD::EvalStatusType::EVAL_OK &&
            cur.getX()->size() == n             )
        {
            if ( cur.getEval(evalType)->isBBOutputComplete() )
            {
                // the center point has been found and put in list
                if ( *centerPt == cur )
                {
                    NOMAD::EvalPoint Y ( cur );
                    T.insert ( Y );
                    if (nbPoints < maxPointsToDisplay)
                    {
                        OUTPUT_DEBUG_START
                        dbgInfo.addMsg(Y.display());
                        OUTPUT_DEBUG_END
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

                        // Issue #382: make sure to evaluate f or h for points in cache (important if cache is loaded from file) see
                        NOMAD::EvalPoint Y ( cur );
                        std::pair<NMSimplexEvalPointSetIterator,bool> ret = T.insert ( Y );


                        if ( ! ret.second )
                        {
                            OUTPUT_DEBUG_START
                            dbgInfo.addMsg("Cannot insert a point in Y (probably tied with another point): " + Y.display() );
                            OUTPUT_DEBUG_END
                        }
                        else
                        {
                            if (nbPoints <= maxPointsToDisplay)
                            {
                                OUTPUT_DEBUG_START
                                dbgInfo.addMsg( ((nbPoints < maxPointsToDisplay) ? Y.display() : "...") );
                                OUTPUT_DEBUG_END
                                nbPoints++;
                            }
                        }
                    }
                }
            }
        }
    }
    OUTPUT_DEBUG_START
    NOMAD::OutputQueue::Add(std::move(dbgInfo));
    OUTPUT_DEBUG_END

    OUTPUT_INFO_START
    AddOutputInfo("Number of potential points to include in initial Y: " + std::to_string(T.size()) );
    OUTPUT_INFO_END

    //
    auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( getAllStopReasons() );

    // not enough points (this is not counted as an error):
    if ( T.size() < minYSize )
    {
        nmStopReason->set ( NMStopType::INITIAL_FAILED );
        OUTPUT_INFO_START
        AddOutputInfo("Stop NM because not enough points in Y.");
        OUTPUT_INFO_END
        return false;
    }



    // Add points in simplex to obtain dim = n+1 and simplex affinely independant
    NOMAD::NMSimplexEvalPointSetIterator itT = T.begin();

    // For debugging
    NOMAD::OutputInfo dbgInfo2(getName(),"Proceed to simplex creation", NOMAD::OutputLevel::LEVEL_DEBUG);

    // First point is always added
    _nmY->insert ( *itT );

    OUTPUT_DEBUG_START
    dbgInfo2.addMsg("k=0: Point z0:" + (*itT).display() ) ;
    dbgInfo2.addMsg(" ---> z0 KEPT in Y ");;
    OUTPUT_DEBUG_END

    int count_feasible = 0;

    if ( (*itT).getH(evalType, computeType).isDefined()
        && (*itT).isFeasible(evalType, computeType) )
        count_feasible = 1 ;

    itT++;
    int k=1;
    while ( itT != T.end() && ! nmStopReason->checkTerminate() )
    {

        if ( _nmY->size() == minYSize )
            break;

        OUTPUT_DEBUG_START
        dbgInfo2.addMsg("k=" + std::to_string(k) +": Point zk:" + (*itT).display()) ;
        OUTPUT_DEBUG_END

        std::pair<NMSimplexEvalPointSetIterator,bool> ret = _nmY->insert ( *itT );

        if ( ! ret.second )
        {
            nmStopReason->set ( NMStopType::INITIAL_FAILED );
            OUTPUT_INFO_START
            AddOutputInfo("Stop NM because cannot insert a point in Y.");
            OUTPUT_INFO_END
            break;
        }

        int rank = getRankDZ();

        if ( rank <= 0 )
        {
            nmStopReason->set ( NMStopType::INITIAL_FAILED );
            OUTPUT_INFO_START
            AddOutputInfo("Cannot get rank of DZ=[(y1-y0 (y2-y0) ...(yk-y0)].");
            OUTPUT_INFO_END
            break;
        }

        // Erase last point or not
        if ( rank != k )
        {
            OUTPUT_DEBUG_START
            dbgInfo2.addMsg(" ---> zk NOT KEPT in Y ");
            OUTPUT_DEBUG_END
            _nmY->erase( ret.first );
            itT++;
        }
        else
        {
            OUTPUT_DEBUG_START
            dbgInfo2.addMsg( " ---> zk KEPT in Y " );
            OUTPUT_DEBUG_END
            if ( (*itT).isFeasible(evalType, computeType) )
                count_feasible++;
            k++;
            itT++;
        }

    }
    OUTPUT_DEBUG_START
    NOMAD::OutputQueue::Add(std::move(dbgInfo2));
    OUTPUT_DEBUG_END

    // not enough points or insufficient rank of simplex (this is not counted as an error):
    // Erase point
    if ( _nmY->size() < minYSize )
    {
        nmStopReason->set ( NMStopType::INITIAL_FAILED );
        OUTPUT_INFO_START
        AddOutputInfo( "Stop NM because not enough simplex points in Y." );
        OUTPUT_INFO_END
        return false;
    }
    if ( getRankDZ() != (int)n )
    {
        nmStopReason->set ( NMStopType::INITIAL_FAILED );
        OUTPUT_INFO_START
        AddOutputInfo( "Stop NM because rank of Y < n." );
        OUTPUT_INFO_END
        return false;
    }



    // Update simplex characteristics (volumes and diameter)
    updateYCharacteristics();

    displayYInfo();

    return true;

}

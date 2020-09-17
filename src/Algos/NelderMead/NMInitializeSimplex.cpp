
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/NelderMead/NMInitializeSimplex.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::NMInitializeSimplex::init()
{
    _name = getAlgoName() + "Initialize Simplex";

    verifyParentNotNull();

}

bool NOMAD::NMInitializeSimplex::runImp()
{
    bool simplexCreated = false;
    if (nullptr == _nmY)
    {
        NOMAD::Exception(__FILE__, __LINE__, "The simplex is not defined.");
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

    auto iter = dynamic_cast<const NOMAD::NMIteration*>( NOMAD::Step::_parentStep );
    if (nullptr == iter)
    {
        NOMAD::Exception(__FILE__, __LINE__, "The simplex initialization must have a NMIteration Step as parent");
    }

    const std::shared_ptr<NOMAD::EvalPoint> centerPt = iter->getFrameCenter();
    // Use center point of iteration, otherwise
    if (nullptr == centerPt)
    {
            NOMAD::Exception(__FILE__, __LINE__, "A center point must be defined.");
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

        OUTPUT_DEBUG_START
        dbgInfo.addMsg("The include rectangle: " + includeRectangle.display() );
        OUTPUT_DEBUG_END
    }
    if ( includeRectangle.max() == 0 )
        NOMAD::Exception(__FILE__, __LINE__, "The include rectangle has no volume");

    const NOMAD::ArrayOfDouble & bbo  = centerPt->getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble();
    size_t m = bbo.size();

    // The set of points initially included
    NOMAD::NMSimplexEvalPointSet T;

    std::vector<NOMAD::EvalPoint> evalpointlist;
    if (NOMAD::EvcInterface::getEvaluatorControl()->getUseCache())
    {
        // browse the cache:
        NOMAD::CacheInterface cacheInterface(this);
        cacheInterface.getAllPoints(evalpointlist);
    }
    else
    {
        auto barrier = getMegaIterationBarrier();
        if (nullptr != barrier)
        {
            evalpointlist = barrier->getAllPoints();
        }
    }

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

                        // Make sure to evaluate f or h for points in cache (important if cache is loaded from file)
                        // TODO

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
    NOMAD::OutputInfo dbgInfo2(_name,"Proceed to simplex creation", NOMAD::OutputLevel::LEVEL_DEBUG);

    // First point is always added
    _nmY->insert ( *itT );

    OUTPUT_DEBUG_START
    dbgInfo2.addMsg("k=0: Point z0:" + (*itT).display() ) ;
    dbgInfo2.addMsg(" ---> z0 KEPT in Y ");;
    OUTPUT_DEBUG_END

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
            if ( (*itT).isFeasible(evalType) )
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

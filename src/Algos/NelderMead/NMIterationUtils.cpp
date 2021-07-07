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
#include "../../Algos/NelderMead/NMIterationUtils.hpp"
#include "../../Math/MatrixUtils.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::NMIterationUtils::setStopReason ( ) const
{
    auto nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get ( _parentStep->getAllStopReasons() );

    if ( nmStopReason == nullptr )
        throw NOMAD::Exception(__FILE__,__LINE__,"NMReflect must have a NM stop reason.");


    switch ( _currentStepType )
    {
        case NOMAD::StepType::NM_REFLECT:
            nmStopReason->set( NOMAD::NMStopType::REFLECT_FAILED );
            break;
        case NOMAD::StepType::NM_EXPAND:
            nmStopReason->set( NOMAD::NMStopType::EXPANSION_FAILED );
            break;
        case NOMAD::StepType::NM_OUTSIDE_CONTRACTION:
            nmStopReason->set( NOMAD::NMStopType::OUTSIDE_CONTRACTION_FAILED);
            break;
        case NOMAD::StepType::NM_INSIDE_CONTRACTION:
            nmStopReason->set( NOMAD::NMStopType::INSIDE_CONTRACTION_FAILED );
            break;
        case NOMAD::StepType::NM_SHRINK:
            nmStopReason->set( NOMAD::NMStopType::SHRINK_FAILED );
            break;
        case NOMAD::StepType::NM_INSERT_IN_Y:
            nmStopReason->set( NOMAD::NMStopType::INSERTION_FAILED );
            break;
        default:
            nmStopReason->set( NOMAD::NMStopType::UNDEFINED_STEP );
            break;
    }
}


/*----------------------------------------------------------------------*/
/* Get the rank of the matrix DZ = [(y1-y0) (y2-y0) ... (yk-y0)]]       */
/*----------------------------------------------------------------------*/
int NOMAD::NMIterationUtils::getRankDZ( ) const
{
    if ( nullptr == _nmY )
        throw NOMAD::Exception(__FILE__, __LINE__, "The iteration utils must have a simplex to work with");

    // The dimension of DZ (k) is related to Y
    size_t k = _nmY->size() - 1 ;

    std::set<NOMAD::EvalPoint>::iterator itY = _nmY->begin();

    const NOMAD::Point & y0 = (*itY);
    const size_t dim = y0.size();

    // DZ : vector of yk-y0 (multidimensional array)
    double ** DZ = new double *[k];
    for (size_t i = 0 ; i < k ; ++i )
        DZ[i]=new double [dim];

    // For debugging
    std::ostringstream outDbg;
    OUTPUT_DEBUG_START
    outDbg << "The rank of DZ=[";
    OUTPUT_DEBUG_END


    // Create DZ
    ++itY;
    size_t j=0;
    while ( j < k )
    {
        OUTPUT_DEBUG_START
        outDbg << " (" ;
        OUTPUT_DEBUG_END
        for ( size_t i = 0; i < dim ; i++ )
        {
            // To get the rank we better use a scaled DZ
            // If Delta is not defined use 1 (no scaling)
            DZ[j][i] = ((*itY)[i].todouble() - y0[i].todouble()  );
            if (i < _Delta.size() && _Delta[i].isDefined())
            {
                DZ[j][i] /= _Delta[i].todouble();
            }
            OUTPUT_DEBUG_START
            outDbg << DZ[j][i] << " ";
            OUTPUT_DEBUG_END
        }
        j++;
        ++itY;
        OUTPUT_DEBUG_START
        outDbg << ")" ;
        OUTPUT_DEBUG_END

    }

    // Get the rank
    int rank= NOMAD::getRank(DZ , k , dim , _rankEps.todouble() );

    OUTPUT_DEBUG_START
    outDbg << " ] equals " << rank;

    NOMAD::OutputQueue::Add(outDbg.str(), NOMAD::OutputLevel::LEVEL_DEBUG);
    OUTPUT_DEBUG_END

    for (size_t i=0 ; i < k ; ++i)
        delete [] DZ[i];;
    delete [] DZ;

    return rank;
}

/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
void NOMAD::NMIterationUtils::updateYCharacteristics()
{
    if ( nullptr == _nmY )
        throw NOMAD::Exception(__FILE__, __LINE__, "The iteration utils must have a simplex to work with");

    // Update Y diameter
    // -----------------
    updateYDiameter();


    // Update Y volumes
    // ----------------
    _simplexVon = -1;
    _simplexVol = -1;

    std::set<NOMAD::EvalPoint>::iterator it1 = _nmY->begin();
    const size_t dim = (*it1).size();

    if ( _nmY->size() != dim + 1 )
        throw NOMAD::Exception(__FILE__, __LINE__, "Cannot get the volume of simplex Y when its dimension is not dimPb+1");

    const NOMAD::Point * y0 = (*it1).getX(); // y0: first element of Y

    if ( y0->size() != dim )
        throw NOMAD::Exception(__FILE__, __LINE__, "Cannot get the volume of simplex Y when dimension of an element is not dimPb");

    // Update volume
    //---------------

    // V : vector of yk-y0 ( for determinant, need square array (n x n) )
    double ** V = new double *[dim];
    for (size_t i = 0 ; i < dim ; i++ )
    {
        V[i] = new double [dim];
    }

    int j = 0;
    ++it1; // Go the second element of Y
    while ( it1 != _nmY->end() )
    {
        for ( size_t i = 0; i < dim ; i++ )
        {
            V[j][i] =  ( (*it1)[i].todouble() - (*y0)[i].todouble() ) ;
        }
        ++it1;
        j++;

    }

    // Get determinant of DZ
    double det;

    bool success = NOMAD::getDeterminant(V, det , dim);

    for ( size_t i = 0; i < dim ; i++ )
        delete [] V[i];
    delete [] V;

    if ( success )
    {
        OUTPUT_DEBUG_START
        NOMAD::OutputQueue::Add("The determinant of the matrix: det( [(y1-y0) (y2-y0) ... (ynf-y0)] ) = " + std::to_string(det), NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END

        double nfact = 1;

        // !n
        for ( size_t i=2 ; i < dim+1 ; i++)
        {
            nfact*=i;
        }

        _simplexVol = fabs(det) / nfact;  // Use fact(n) for volume

        if ( _simplexDiam > 0 )
            _simplexVon = _simplexVol / pow(_simplexDiam,dim) ;
        else
        {
            OUTPUT_DEBUG_START
            NOMAD::OutputQueue::Add("Cannot get the normalized volume of simplex Y because simplex diameter <=0. Let's continue. ", NOMAD::OutputLevel::LEVEL_DEBUG);
            OUTPUT_DEBUG_END
            return;
        }
    }
    else
    {
        OUTPUT_DEBUG_START
        NOMAD::OutputQueue::Add("Cannot get the volume of simplex Y because determinant failed. Continue", NOMAD::OutputLevel::LEVEL_DEBUG);
        OUTPUT_DEBUG_END
    }
    return ;
}


/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
void NOMAD::NMIterationUtils::updateYDiameter()
{
    std::set<NOMAD::EvalPoint>::const_iterator it1;
    std::set<NOMAD::EvalPoint>::iterator it2;

    _simplexDiam = 0;
    for (it1 = _nmY->begin(); it1 != _nmY->end(); ++it1)
    {
        it2 = it1;
        ++it2;
        while ( it2 != _nmY->end() )
        {
            const NOMAD::Direction DZ((*it1) - (*it2));
            const double lengthDZ = DZ.norm().todouble();

            if ( lengthDZ > _simplexDiam )
            {
                _simplexDiam = lengthDZ;
                _simplexDiamPt1 = &(*it1);
                _simplexDiamPt2 = &(*it2);
            }
            ++it2;
        }
    }

    return ;
}


/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
void NOMAD::NMIterationUtils::displayYInfo()const
{
    if ( nullptr == _nmY )
        throw NOMAD::Exception(__FILE__, __LINE__, "The iteration utils must have a simplex to work with");

    OUTPUT_DEBUG_START
    NOMAD::OutputInfo dbgInfo("NM iteration utils", "", NOMAD::OutputLevel::LEVEL_DEBUG );

    OUTPUT_INFO_START
    _parentStep->AddOutputInfo("Number of points in the simplex Y: " + std::to_string(_nmY->size()) );
    OUTPUT_INFO_END

    if ( _simplexVol > 0 )
        dbgInfo.addMsg("The volume of the simplex: " + std::to_string( _simplexVol ) );
    else
        dbgInfo.addMsg("WARNING: Evaluation of the simplex volume failed.");

    if ( _simplexDiam > 0 )
        dbgInfo.addMsg("The diameter of the simplex: " + std::to_string( _simplexDiam ) );
    else
        dbgInfo.addMsg( "WARNING: Evaluation of the simplex diameter failed.");

    if ( _simplexVon > 0 )
        dbgInfo.addMsg("The normalized volume of the simplex: " + std::to_string( _simplexVon ) );
    else
        dbgInfo.addMsg( "WARNING: Evaluation of the simplex normalized volume failed." );

    std::set<NOMAD::EvalPoint>::iterator itY;
    dbgInfo.addMsg("The simplex Y contains: ");
    for (itY =_nmY->begin(); itY != _nmY->end(); ++itY)
    {
        dbgInfo.addMsg( (*itY).display()) ;
    }
    getRankDZ();


    NOMAD::OutputQueue::Add(std::move(dbgInfo));
    NOMAD::OutputQueue::Flush();
    OUTPUT_DEBUG_END
}

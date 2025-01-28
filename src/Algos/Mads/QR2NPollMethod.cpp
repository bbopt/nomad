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

#include "../../Algos/Mads/QR2NPollMethod.hpp"
#include "../../Math/Direction.hpp"
#include "../../Math/MatrixUtils.hpp"

void NOMAD::QR2NPollMethod::init()
{
    setStepType(NOMAD::StepType::POLL_METHOD_QR_2N);
    verifyParentNotNull();

}

// Generate poll directions
void NOMAD::QR2NPollMethod::generateUnitPollDirections(std::list<NOMAD::Direction> &directions, const size_t n) const
{
    directions.clear();

    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);
    
    OUTPUT_DEBUG_START
    AddOutputDebug("Unit sphere direction: " + dirUnit.display());
    OUTPUT_DEBUG_END
    
    while (dirUnit[0] == 0)
    {
        NOMAD::Direction::computeDirOnUnitSphere(dirUnit);
    }
    
    // Matrix M
    auto M = new double*[n];
    for (size_t i = 0; i < n; ++i)
    {
        M[i] = new double [n];
        M[i][0] = dirUnit[i].todouble();
        for (size_t j = 1; j < n; ++j)
        {
            M[i][j] = (i == j)? 1.0:0.0;
        }
    }
    OUTPUT_DEBUG_START
    AddOutputDebug("M matrix for QR:");
    for (size_t i = 0; i < n; ++i)
    {
        NOMAD::ArrayOfDouble aod(n);
        for (size_t j = 0; j < n; ++j)
        {
            aod[j]=M[i][j];
        }
        AddOutputDebug(aod.display());
        
    }
    OUTPUT_DEBUG_END
    
    // Matrices Q and R
    auto Q = new double*[n];
    auto R = new double*[n];
    for (size_t i = 0; i < n; ++i)
    {
        Q[i] = new double [n];
        R[i] = new double [n];
    }
    
    std::string error_msg;
    bool success = NOMAD::qr_factorization (error_msg,M,Q,R,static_cast<int>(n),static_cast<int>(n));
    
    if ( !success || !error_msg.empty())
    {
        OUTPUT_INFO_START
        AddOutputInfo("QR decomposition for QR 2N poll method has failed");
        OUTPUT_INFO_END
    }
    
    OUTPUT_DEBUG_START
    AddOutputDebug("Direction after QR decomposition: ");
    OUTPUT_DEBUG_END
    
    // Ordering D_k alternates Qk and -Qk instead of [Q_k -Q_k]
    NOMAD::Direction dir(n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            dir[j] = Q[j][i];
        }
        OUTPUT_DEBUG_START
        AddOutputDebug(dir.display());
        OUTPUT_DEBUG_END
        
        directions.push_back(dir);
        directions.push_back(-dir);
    }
    
    OUTPUT_DEBUG_START
    NOMAD::OutputQueue::Flush();
    OUTPUT_DEBUG_END
    
    // Delete M, Q and R:
    for ( size_t i = 0 ; i < n ; ++i )
    {
        delete [] M[i];
        delete [] Q[i];
        delete [] R[i];
    }
    delete [] Q;
    delete [] R;
    delete [] M;
    
}

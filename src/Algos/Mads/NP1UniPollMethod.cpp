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

#include "../../Algos/Mads/NP1UniPollMethod.hpp"
#include "../../Math/Direction.hpp"

void NOMAD::NP1UniPollMethod::init()
{
    setStepType(NOMAD::StepType::POLL_METHOD_UNI_NPLUS1);
    verifyParentNotNull();
}


void NOMAD::NP1UniPollMethod::generateUnitPollDirections(std::list<Direction> &directions, const size_t n) const
{
    directions.clear();

    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);

    // Ortho MADS 2n
    // Householder Matrix
    // A vintage piece of code
    NOMAD::Direction** H = new NOMAD::Direction*[2*n];

    // Ordering D_k alternates Hk and -Hk instead of [H_k -H_k]
    std::list<Direction> vDirs;
    for (size_t i = 0; i < n; ++i)
    {
        vDirs.push_back(NOMAD::Direction(n, 0.0));
        H[i]   = &(vDirs.back());
        vDirs.push_back(NOMAD::Direction(n, 0.0));
        H[i+n] = &(vDirs.back());
    }
    // Householder transformations on the 2n directions on a unit n-sphere
    NOMAD::Direction::householder(dirUnit, true, H);

    // dir 0
    NOMAD::Direction dir0(*H[0]);
    for ( size_t i = 1 ; i < n ; ++i )
    {
        dir0=dir0+(*H[i]);
    }
    dir0*=-1.0/sqrt(double(n));
    directions.push_back(dir0);

    NOMAD::Double beta=(sqrt(double(n+1.0))-1.0)/sqrt(double(n));
    dir0*=beta;

    for ( size_t i = 0 ; i < n ; i++ )
    {
        NOMAD::Direction diri(*H[i]);
        diri*=sqrt(double(n+1));
        diri=diri+dir0;
        diri*=1.0/sqrt(double(n));

        directions.push_back(diri);
    }
    delete [] H;
}
 // end generateTrialPoints

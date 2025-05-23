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
#include "../Math/MathUtils.hpp"

#include <cmath>

bool NOMAD::roots_quadratic(const double q2, const double q1, const double q0,
                            double& r1, double& r2)
{
    if (q2 == 0)
    {
        if (q1 == 0)
        {
            r1 = r2 = 0;
            if (q0 != 0)
            {
                return false;
            }
        }
        else
        {
            r1 = r2 = -q0 / q1;
        }
        return true;
    }

    // q is quadratic
    const double tol = 1e-8; // =~ sqrt(eps(double))
    const double rhs = tol * q1 * q1;
    if (std::fabs(q0 * q2) > rhs)
    {
        const double rho = q1 * q1 - 4 * q2 * q0;
        if (rho < 0)
        {
            return false;
        }
        const double numd2 = -(q1 + std::copysign(sqrt(rho), q1)) / 2.0;
        r1 = numd2 / q2;
        r2 = q0 / numd2;
    }
    else
    {
        // Ill-conditioned quadratic
        r1 = -q1 / q2;
        r2 = 0;
    }

    // Improve accuracy with one newton iteration
    const size_t niter = 1;
    for (size_t i = 0; i < niter; ++i)
    {
        const double q = (q2 * r1 + q1) * r1 + q0;
        const double dq = 2 * q2 * r1 + q1;
        if (dq == 0)
            continue;

        r1 -= q / dq;
    }

    for (size_t i = 0; i < niter; ++i)
    {
        const double q = (q2 * r2 + q1) * r2 + q0;
        const double dq = 2 * q2 * r2 + q1;
        if (dq == 0)
            continue;

        r2 -= q / dq;
    }

    return true;
}

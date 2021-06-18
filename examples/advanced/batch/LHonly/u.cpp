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
#include <ctime>        // For time_t
#include <fstream>      // For ifstream
#include <iostream>
#include <cmath>        // For sqrt
#include <stdexcept>    // For logic_error
//#include <unistd.h>     // For usleep

const int n = 10;

int main (int argc, char **argv)
{
    bool eval_ok = false;

    // Remotely based on G2.
    double f = 1e+20, g1 = 1e+20, g2 = 1e+20;
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, prod1 = 1.0, prod2 = 1.0;
    double x[n];

    bool x0read = false;
    if (argc >= 2)
    {
        std::string x0file = argv[1];
        std::ifstream in (argv[1]);
        for (int i = 0; i < n; i++)
        {
            if (in.fail())
            {
                std::cerr << "Error reading file " << x0file << " for x0." << std::endl;
                x0read = false;
                break;
            }
            in >> x[i];
            x0read = true;
        }
        in.close();
    }

    if (x0read)
    {
        try
        {
            for (int i = 0; i < n ; i++)
            {
                //std::cerr << "sleep " << i << std::endl;
                //usleep(100000);
                sum1  += pow(cos(x[i]), 4);
                sum2  += x[i];
                sum3  += (i+1)*x[i]*x[i];
                prod1 *= pow(cos(x[i]), 2);
                if (prod2 != 0.0)
                {
                    if (x[i] == 0.0)
                    {
                        prod2 = 0.0;
                    }
                    else
                    {
                        prod2 *= x[i];
                    }
                }
            }

            g1 = -prod2 + 0.75;
            g2 = sum2 -7.5 * n;

            f = 10*g1 + 10*g2;
            if (0.0 != sum3)
            {
                f -= (sum1 -2*prod1) / std::abs(sqrt(sum3));
            }
            else
            {
                f = NAN;
            }
            // Scale
            if (!std::isnan(f))
            {
                f *= 1e-5;
            }

            eval_ok = !std::isnan(f);

        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }
    }

    std::cout << - f - 2000 << " " << f << std::endl;

    // Return 0 if eval_ok.
    // Note: returning a value other than 0 will make this evaluation
    // an EVAL_ERROR, meaning the point could be re-evaluated.
    // If eval_ok is true and f is NaN, then this evaluation is an
    // EVAL_FAILED, and the point will not be re-evaluated.
    return !eval_ok;
}

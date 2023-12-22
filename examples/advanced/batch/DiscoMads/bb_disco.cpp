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
#include <cmath>        // For sqrt
#include <ctime>        // For time_t
#include <fstream>      // For ifstream
#include <iostream>
#include <stdexcept>    // For logic_error
#include <string>

const int n = 2;

int main (int argc, char **argv)
{
    bool eval_ok = false;

    double f = 1e+20, c1 = 1e+20;
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
            double x0_f = 0;
            double y0_f = 10;
            double radius_f = 12;

            // Objective function
            if(pow(x[0]-x0_f,2)+pow(x[1]-y0_f,2)>pow(radius_f,2)){
                f=-0.025*x[1]+3;
            }
            else{
                f= 0.04*x[1];
            }
                
            // Constraint
            if (x[0]>0){
                c1=2;
            }
            else{
                c1=-0.1;
            }

            eval_ok=true;
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }
    }

    std::cout << f  << " " << c1 << std::endl;

    // Return 0 if eval_ok.
    // Note: returning a value other than 0 will make this evaluation
    // an EVAL_ERROR, meaning the point could be re-evaluated.
    // If eval_ok is true and f is NaN, then this evaluation is an
    // EVAL_FAILED, and the point will not be re-evaluated.
    return !eval_ok;
}

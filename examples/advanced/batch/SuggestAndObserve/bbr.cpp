/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0 has been created by                                        */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0 is owned by                               */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,            */
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
//Blackbox for evaluation

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#define DIMENSION 2

bool eval_xs(const std::vector<std::vector<double>> &xs, std::vector<double>& fxs)
{
    bool eval_ok = false;

    size_t nPoints = xs.size();

    double obj;

    try
    {
        for (size_t i = 0; i < nPoints ; i++)
        {
            // Rosenbrock 2D
            if (xs[i].size() != DIMENSION)
            {
                throw "Dimension should be " + DIMENSION;
            }
            obj = pow(10 * (xs[i][1] - pow(xs[i][0],2)), 2);
            obj += pow(1.0 - xs[i][0], 2);
            fxs.push_back(obj);
            std::cout << obj << std::endl;
        }

        eval_ok = true;
    }
    catch (std::exception &e)
    {
        std::string err("Exception: ");
        err += e.what();
        throw std::logic_error(err);
    }

    return eval_ok;
}


std::vector<std::vector<double>> readXFile(const std::string& xFileName)
{
    std::vector<std::vector<double>> xs;
    std::ifstream in(xFileName);

    // Get the input points
    while (!in.eof())
    {
        std::vector<double> x(DIMENSION);
        for (size_t i = 0; i < DIMENSION; i++)
        {
            in >> x[i];
        }
        xs.push_back(x);
    }
    in.close();
    xs.pop_back();

    return xs;
}


int main(int argc, char ** argv)
{
    try
    {
        if (argc < 2)
        {
            std::cerr << "Usage: " << argv[0] << " <xfile>" << std::endl;
            return 1;
        }

        auto xs = readXFile(argv[1]);

        std::vector<double> fxs;
        eval_xs(xs,fxs);
    }

    catch (std::exception &e)
    {
        std::cerr << "\n" << argv[0] << " has been interrupted (" << e.what() << ")\n\n";
    }

    return 0;
}

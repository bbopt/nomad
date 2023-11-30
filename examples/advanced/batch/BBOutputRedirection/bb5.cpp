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
//
//  BB_5
//
//  Created by CT
//
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;


// Blackbox manages sending the results in the nomad output file.
// The name of the output MUST be inputfilename.output.
// INPUT and OUTPUT files have a fixed name during optimization.
// Some verbose can be displayed on the standard output or standard error and will not be considered by Nomad
// The verbose is stored in a log file.
int main(int argc, const char ** argv) 
{
    
    double f = 1e20; 
    double x[2];
    if (argc >=2)
    {
		// Output filename MUST be input file name + ".output" 
        std::string outputfilename = std::string(argv[1]) + ".output";
        
        ifstream in (argv[1]); 
        for ( int i = 0 ; i < 2 ; i++ ) 
        {
          in >> x[i];
        }
        if ( in.fail() )
        {
            // std::cout << "Reading input file " << argv[1] << " fails. Nothing to output. Return execution status 1." << std::endl;
            return 1;
        }
        std::cout << "Input vector x ( " << x[0] << " " << x[1] << " )" << std::endl;
        
        f = pow (5 * x[0]-2 , 4) + pow (5 * x[0]-2, 2) * pow( x[1] , 2) +pow ( 3 * x[1] + 1 , 2);

        std::cout << "Objective computed successfully. F(x)=" << f << std::endl;
        
        ofstream outfile;
        outfile.open(outputfilename.c_str(), std::ofstream::trunc);
        if (outfile.fail())
        {
			std::cout << "Cannot open output file " << outputfilename.c_str() << std::endl;
        }
        else
        {
            std::cout << "Output objective value in file "<< outputfilename.c_str() << std::endl;
        }
        outfile << f << std::endl;
        outfile.close();
            
        return 0;
    }
    else
    {
        std::cout << "Input file name is not provided to the blackbox" << std::endl;
        return 1;
    }
       
  }
    

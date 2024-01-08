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
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
using namespace std;



int main ( int argc , char ** argv ) {
    
    double f = 1e20;
    
    int n = 12;
    const double pi = 3.1415927;
    f = 10.0*n;
    
    if ( argc >= 2 ) {

        double * x = new double[n];
        // read the point:
        ifstream in ( argv[1] );
        for (int i = 0 ; i < n ; ++i )
        {
            in >> x[i];
        }
        in.close();
        
        if ( !in.fail() )
            for ( int j = 0 ; j < n ; j++)
            {
                f += x[j]*x[j] - 10.0*cos(2.0*pi*x[j]);
            }
    }
    
    cout.setf(ios::fixed);
    cout.precision(15);
    cout << f << endl;
    return 0;
}
 



/*
int main ( int argc , char ** argv ) {

  double f = 1e20;

  if ( argc >= 2 ) {

    double   x , y;
    ifstream in ( argv[1] );

    in >> x >> y;
    
    in.close();

    if ( !in.fail() )
      f = 20 + x*x + y*y - 10*cos(6.283185307*x) - 10*cos(6.283185307*y);

  }

  cout.setf(ios::fixed);
  cout.precision(15);

  cout << f << endl;

  return 0;
}
 */

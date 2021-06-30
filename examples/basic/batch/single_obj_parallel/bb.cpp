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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

int main ( int argc , char ** argv )
{


	double f = 1e20, c1 = 1e20 , c2 = 1e20;
	std::vector< double * > X;

	if ( argc >= 2 )
	{

		ifstream in ( argv[1] );

		// Get the input vectors
		while( ! in.eof() )
		{
			double *x = new double [8];
			for ( int i = 0 ; i < 5 ; i++ ) 
				in >> x[i];

			X.push_back(x);  
		}

		in.close();
		X.pop_back();

		int nb_pts= (int)X.size();

		// Evaluate the points in parallel
#ifdef _OPENMP
#pragma omp parallel
#endif
		for (int i=0;i<nb_pts;++i)
		{
			// printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
			double *x_t=X[i];
			double c1 = 0.0 , c2 = 0.0;
			for ( int j = 0 ; j < 5 ; j++ ) 
			{
				//cout << x_t[j] << " ";
				c1 += pow(x_t[j]-1,2);
				c2 += pow(x_t[j]+1,2);
			}
			x_t[5] =  x_t[4];
			x_t[6] =  c1 - 25;
			x_t[7] =  25 - c2;
			// cout << x_t[6] << " " << x_t[7] << " " << endl;
		}

		// Output the points in the same order as in the input file
		for (int i=0 ; i < nb_pts ; ++i)
		{
			cout << X[i][5]  << " " << X[i][6]  << " " << X[i][7]  <<  endl;
			delete [](X[i]);
		}

	}
	return 0;

}	

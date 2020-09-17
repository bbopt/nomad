/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#define SGTELIB_PREDICT_VERSION "1.0  [Fev. 18th, 2019]"
// sgtelib_predict.exe X.txt Z.txt XX.txt model.txt


#include "sgtelib.hpp"
#include "Surrogate_Factory.hpp"
#include "Surrogate_Utils.hpp"

#include <stdio.h>
#include <string.h>

using namespace std;



char errstr[1024];

//Function Prototypes
void printSolverInfo();

// Main Entry Function
// -----------------------------------------------------------------
int main(int argc, char ** argv)
{
    
	if (argc < 5)
	{
		printInfo();
		return 0;
	}

    //Set Option Defaults
    int printLevel = 1;
    

    
    //Get the model definition in text format
	std::ifstream modelFile(argv[5]);
	std::string sgtelib_model= "";
	if (! modelFile.good())
	{
		std::cerr << "Error reading model file" << argv[5];
		modelFile.close();
		return false;
	}
    else
    {
		modelFile >> sgtelib_model;
		modelFile.close();
    }
       
	//Internal Vars
	size_t mX, nX, mZ, nZ, mXX, nXX;
    SGTELIB::Matrix X(argv[1]);
	X.set_name("X");
	X.display(std::cout);

    SGTELIB::Matrix Z(argv[2]);
	Z.set_name("Z");
	Z.display(std::cout);

    SGTELIB::Matrix XX(argv[3]);
	XX.set_name("XX");

	nXX = XX.get_nb_cols();
	nX = Z.get_nb_rows();
    SGTELIB::Matrix ZZ ("ZZ",mXX,nZ);

	//cout << "size of X: mRows=" << mX << " nCols=" << nX << endl;
	// cout << "size of Z: mRows=" << mZ << " nCols=" << nZ << endl;
	// cout << "size of XX: mRows=" << mXX << " nCols=" << nXX << endl;

	//Check Sizes
	if (mX == 0 || mZ == 0 || mXX == 0 || nX == 0 || nZ == 0 || nXX == 0)
		std::cerr << "error: size of X, Z and XX must be greater than 0" << std::endl;


	if (mX != mZ)
		std::cerr << "error: size of X and Z must be equal" << std::endl;

   
    //cout << " X nb_rows (nb points)= " << X.get_nb_rows() << std::endl;
    //cout << " X nb_cols (n) = " <<  X.get_nb_cols() << std::endl;
    //cout << " X nb_cols " << Z.get_nb_cols() << std::endl;
    
    try {
        SGTELIB::TrainingSet TS(X,Z); 
        SGTELIB::Surrogate * S = Surrogate_Factory(TS,sgtelib_model);
        S->build();
        
        TS.display(std::cout);

        S->predict(XX,&ZZ);
        ZZ.display(std::cout);
    } 
  catch(exception &e)
    {
        std::cerr << "SGTELIB Error:" << e.what();
    }
	
  return 0;
}

//Print Solver Information
void printInfo()
{
	std::cout
		<< "\n-----------------------------------------------------------" << std::endl
		<< " SGTELIB: 2.0.2 " << std::endl << "  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html" << std::endl
		<< "  - Source available from: https://github.com/bastientalgorn/sgtelib" << std::endl << std::end
		<< "Command line Interface C. Tribes 2018  " << std::endl
		<< "	usage: sgtelib_predict.exe X.txt Y.txt XX.txt model.txt" << std::endl
		<< "-----------------------------------------------------------" << std::endl;
}

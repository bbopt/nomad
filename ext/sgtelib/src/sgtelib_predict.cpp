
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

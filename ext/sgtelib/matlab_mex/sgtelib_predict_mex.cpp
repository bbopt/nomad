
#define SGTELIB_PREDICT_VERSION "1.0  [June. 12th, 2017]"
// 1) The default interface is:
//      [ZZ] = sgtelib_predict(X,Z,XX,model)


#include "mex.h"
#include "sgtelib.hpp"
#include "Surrogate_Factory.hpp"
#include "Surrogate_Utils.hpp"

#include <stdio.h>
#include <string.h>

using namespace std;


//Ctrl-C Detection
#ifdef __cplusplus
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);
#else
extern bool utIsInterruptPending();
extern void utSetInterruptPending(bool);
#endif


//PRHS Defines
#define pX      prhs[0]
#define pZ      prhs[1]
#define pXX     prhs[2]
#define pModel  prhs[3]


char errstr[1024];

//Function Prototypes
void printSolverInfo();
// int checkInputs(const mxArray *prhs[], int nrhs, mxArray *plhs[], int nlhs);
//void lower(char *str);
//double getStatus(int stat);

//cout Redirection
struct printfbuf : std::streambuf {
public:
    //Constructor
    printfbuf()
    {
        setp(m_buffer, m_buffer + s_size - 2);
    }
private:
    enum { s_size = 1024 }; //not sure on this size
    char m_buffer[s_size];
    int_type overflow(int_type c)
    {
        if (!traits_type::eq_int_type(c, traits_type::eof()))
        {
            *pptr() = traits_type::to_char_type(c);
            pbump(1);
        }
        return sync() != -1 ? traits_type::not_eof(c) : traits_type::eof();
    }
    
    int sync()
    {
        *pptr() = 0;
        mexPrintf(pbase());
        mexEvalString("drawnow;");
        setp(m_buffer, m_buffer + s_size - 2);
        return 0;
    }
};


static printfbuf buf;

// Main Entry Function
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    //Internal Vars
    size_t nX,nZ,nXX;
    
    //Redirect cout
    std::cout.rdbuf(&buf); //redirect buffer
    

    //Set Option Defaults
    int printLevel = 1;

    // Test input arguments
    if (nrhs < 4)
       mexErrMsgTxt("error: insufficient number of arguments");
    
    if( !mxIsDouble(pX) || mxIsComplex(pX) || mxIsEmpty(pX) )
        mexErrMsgTxt("X must be a real double column vector!");
    
    if( !mxIsDouble(pZ) || mxIsComplex(pZ) || mxIsEmpty(pZ) )
        mexErrMsgTxt("Z must be a real double column vector!");
    
    if( !mxIsDouble(pXX) || mxIsComplex(pXX) || mxIsEmpty(pXX) )
        mexErrMsgTxt("XX must be a real double column vector!");

    
    //Get Sizes
    nX = mxGetNumberOfElements(pX);
    
    // cout << nX << endl;
    
    nZ = mxGetNumberOfElements(pZ);
    
    //cout << nZ << endl;
    nXX =mxGetNumberOfElements(pXX);
    // cout << nXX << endl;
    
    //Check Sizes
    if ( nX == 0 || nZ == 0 || nXX == 0 )
        mexErrMsgTxt("error: size of X, Z and XX must be greater than 0");
    
    if ( nX!=nZ )
        mexErrMsgTxt("error: size of X, Z and XX must be equal");
    
    
    //Get Blackbox / Objective Function Handle
    const char *str=NULL;
    std::string model;
    if ( mxIsChar(pModel) )
    {
        int i=0;
        str = mxArrayToString(pModel);
        
        /* get the length of the input string */
        size_t buflen = (mxGetM(pModel) * mxGetN(pModel)) + 1;

        model = std::string(str);
        
        // cout << buflen <<endl;
    }
    else
    {
        mexErrMsgTxt("error reading model string");
    }
    
    
    
    if( str == NULL)
        mexErrMsgTxt("error reading model string (empty string)");

    
    
    if(printLevel)
    {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is SGTELIB v2.0.2 \n");
        mexPrintf(" Authors: B. Talgorn and S. Le Digabel\n");
        mexPrintf(" MEX Interface C. Tribes 2017 \n\n");
        mexPrintf("Variable properties:\n");
        mexPrintf(" # Decision Variables:               %4d\n",nX);
        mexPrintf(" # Number of dependant variable(s):  %4d\n",nZ);
        mexPrintf(" # Number of prediction variable(s):  %4d\n",nXX);
        cout << "Model name: " <<endl << "   " << model << endl;
        mexPrintf("------------------------------------------------------------------\n");
        mexEvalString("drawnow;"); //flush draw buffer
    }

    
    SGTELIB::Matrix X("X",nX,1);
    SGTELIB::Matrix Z("Z",nZ,1);
    SGTELIB::Matrix XX("XX",nXX,1);
    SGTELIB::Matrix ZZ ("ZZ",XX.get_nb_rows(),Z.get_nb_cols());
    
    double *x=mxGetPr(pX);
    double *z=mxGetPr(pZ);
    double *xx =mxGetPr(pXX);
    for (int i=0;i<nX;i++)
    {
        X.set_row(x[i],i);
        Z.set_row(z[i],i);
    }
    X.display(std::cout);
    Z.display(std::cout);
    
    for (int i=0;i<nXX;i++)
    {
        XX.set_row(xx[i],i);
    }
    XX.display(std::cout);
    
    try {
        SGTELIB::TrainingSet TS(X,Z);
        
        SGTELIB::Surrogate * S = Surrogate_Factory(TS,model);
        S->build();
        
        S->predict(XX,&ZZ);
        ZZ.display(std::cout);
    }
  catch(exception &e)
    {
        sprintf(errstr,"SGTELIB Error:\n\n%s",e.what());
        mexErrMsgTxt(errstr);
    }
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(nXX,1, mxREAL);
    
    double *zz = mxGetPr(plhs[0]);
    
    for(int j=0;j<nXX;j++)
        zz[j] = ZZ.get(j);

    
    //Return error control to default
    mexSetTrapFlag(0);
    
}

//Print Solver Information
void printSolverInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" SGTELIB: 2.0.2 \n");
    mexPrintf("  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html\n");
    mexPrintf("  - Source available from: https://github.com/bastientalgorn/sgtelib\n");
    mexPrintf("\n MEX Interface C. Tribes 2017  \n");
    mexPrintf("-----------------------------------------------------------\n");
}

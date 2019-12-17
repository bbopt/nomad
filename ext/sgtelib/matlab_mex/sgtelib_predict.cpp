/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
    size_t mX,nX,mZ,nZ,mXX,nXX;
    
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
    mX = mxGetM(pX);
    nX = mxGetN(pX);
    cout << "size of X: mRows=" << mX <<  " nCols=" << nX << endl;
    
    mZ = mxGetM(pZ);
    nZ = mxGetN(pZ);
    cout << "size of Z: mRows=" << mZ <<  " nCols=" << nZ << endl;

    mXX = mxGetM(pXX);
    nXX = mxGetN(pXX);
    cout << "size of XX: mRows=" << mXX <<  " nCols=" << nXX << endl;
    
    //Check Sizes
    if ( mX == 0 || mZ == 0 || mXX == 0 || nX == 0 || nZ == 0 || nXX == 0 )
        mexErrMsgTxt("error: size of X, Z and XX must be greater than 0");
    
    if ( nZ != 1 )
        mexErrMsgTxt("error: Z must be a scalar dependant variable  ");


    if ( mX!=mZ )
        mexErrMsgTxt("error: size of X and Z must be equal");
    
    
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
        mexPrintf(" MEX Interface C. Tribes 2018 \n\n");
        //mexPrintf("Variable properties:\n");
        //mexPrintf(" # Decision Variables:               %4d\n",nX);
        mexPrintf(" # Number of dependant variable(s):  %4d\n",nZ);
        mexPrintf(" # Number of prediction variable(s):  %4d\n",nXX);
        cout << "Model name: " <<endl << "   " << model << endl;
        mexPrintf("------------------------------------------------------------------\n");
        mexEvalString("drawnow;"); //flush draw buffer
    }

    
    SGTELIB::Matrix X("X",mX,nX);
    SGTELIB::Matrix Z("Z",mZ,nZ);
    SGTELIB::Matrix XX("XX",mXX,nXX);
    SGTELIB::Matrix ZZ ("ZZ",mXX,nZ);
    

    /*mxDouble *x  = mxGetDoubles(pX);
    mxDouble *z  = mxGetDoubles(pZ);
    mxDouble *xx = mxGetDoubles(pXX); */
        
    double *x = mxGetPr(pX); 
    int k=0;
    for (int j=0; j < nX ;j++)
    {
	for (int i=0; i < mX ; i++,k++)
	{	
           // std::cout << "i=" << i << " j="<< j << " X(i,j)=" << x[k] <<  std::endl;   
           X.set(i,j,(double)x[k]);
        }   
    }
 
   X.display(std::cout);
 

    double *z = mxGetPr(pZ);
   for (int i=0; i < mZ ; i++)
   {
      Z.set(i,0,(double)z[i]);
      //std::cout << "i=" << i << " Z(i,1)=" << z[i] <<  std::endl;
   }
 

   Z.display(std::cout);


    k=0;
    double *xx = mxGetPr(pXX);
    for (int j=0; j < nXX ;j++)
    {
        for (int i=0; i < mXX ; i++,k++)
        {	
           //std::cout << "i=" << i << " j="<< j << " XX(i,j)=" << xx[k] <<  std::endl;   
           XX.set(i,j,(double)xx[k]);
        } 
    }

    XX.display(std::cout);

    //cout << " X nb_rows (nb points)= " << X.get_nb_rows() << std::endl;
    //cout << " X nb_cols (n) = " <<  X.get_nb_cols() << std::endl;
    //cout << " X nb_cols " << Z.get_nb_cols() << std::endl;
    
    try {
        SGTELIB::TrainingSet TS(X,Z);
        
        SGTELIB::Surrogate * S = Surrogate_Factory(TS,model);
        S->build();
        
TS.display(std::cout);

        S->predict(XX,&ZZ);
        ZZ.display(std::cout);
    } 
  catch(exception &e)
    {
        sprintf(errstr,"SGTELIB Error:\n\n%s",e.what());
        mexErrMsgTxt(errstr);
    }



    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(mXX,1, mxREAL);
    
    double *zz = mxGetPr(plhs[0]);
    
    for(int i=0;i<mXX;i++)
        zz[i] = ZZ.get(i);

    
    //Return error control to default
    //mexSetTrapFlag(0);
  
}

//Print Solver Information
void printSolverInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" SGTELIB: 2.0.2 \n");
    mexPrintf("  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html\n");
    mexPrintf("  - Source available from: https://github.com/bastientalgorn/sgtelib\n");
    mexPrintf("\n MEX Interface C. Tribes 2018  \n");
    mexPrintf("-----------------------------------------------------------\n");
}

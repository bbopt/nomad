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

#define NOMADMEX_VERSION "1.4  [Feb 21th, 2022]"
//NOTE from Version 1.3 on this MEX file contains only the GERAD interface.
// The OPTI interface is no longer supported.
// Another change compared with previous versions: the parameters are passed as a list of strings (in the form struct('KEYWORD','value',...)). The keyword and value follow exactly the Nomad syntax.

// The GERAD interface is now:
//      [x,fval,exitflag,iter] = nomadOPt(fun,x0,lb,ub,params,callbackFun)


#include "mex.h"
#include "Nomad/nomad.hpp"
#include "Cache/CacheBase.hpp"
#include "Algos/EvcInterface.hpp"
#include <stdio.h>
#include <string.h>
#include <memory.h>


using namespace std;


//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 4 /* max nrhs for user function */

typedef struct {
    char f[FLEN];
    mxArray *plhs[2];
    mxArray *prhs[MAXRHS];
    int xrhs, nrhs;
    double *nlrhs;
} usrFcn;

// Evaluation callback structure
typedef struct {
	char f[FLEN];
	mxArray *plhs[1];
	mxArray *prhs[4];
	bool enabled;
} usrCallbackFcn;



//Ctrl-C Detection
#ifdef __cplusplus
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);
#else
extern bool utIsInterruptPending();
extern void utSetInterruptPending(bool);
#endif

//Argument Enums (in expected order of arguments)
enum {eFUN, eX0, eLB, eUB, ePARAM, eCB};

//PRHS Defines
#define pFUN    prhs[eFUN]
#define pX0     prhs[eX0]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pParam  prhs[ePARAM]
#define pCB     prhs[eCB]


//Function Prototypes
void printSolverInfo();
int checkInputs(const mxArray *prhs[], int nrhs, mxArray *plhs[], int nlhs);
void setNomadParams(const std::shared_ptr<NOMAD::AllParameters> p, const mxArray *params);
void lower(char *str);
double getStatus(int stat);

int counter_eval = 0;

//MATLAB Evaluator Class
class matlabEval : public NOMAD::Evaluator
{
private:
    usrFcn *fun;
    // TODO bool hasSur;
    // TODO bool hasAdditionalParam;
    size_t nobj;
    size_t ncon;

	usrCallbackFcn *callbackF;


public:
    //Constructor

    matlabEval(const std::shared_ptr<NOMAD::EvalParameters> p, usrFcn *_fun, size_t _nobj, size_t _ncon, usrCallbackFcn * _callbackF ) : NOMAD::Evaluator(p)
    {
        fun     = _fun;
        // hasSur  = false ; // TODO p.has_sgte();
        nobj    = _nobj;
        ncon    = _ncon;
		counter_eval = 0;
		callbackF = _callbackF;

    }
    //Destructor
    ~matlabEval(void) {}

// TODO
// std::vector<bool> eval_block(Block &block,
//                                          const Double &hMax,
//                                                                                   std::vector<bool> &countEval) const; 
//     bool eval_x(std::list<NOMAD::EvalPoint *> &x, const NOMAD::Double &h_max, std::list<bool> & list_count_eval )
//     {

//         char errstr[1024];
//         bool stop = false;
//         int i, j, m, n;
//         double  *fvals,*count;
//         mxLogical *sur;

//         m=static_cast<int>(x.size());
//         n=(*(x.begin()))->size();

//         if ( m !=list_count_eval.size())
//         {
//             mexPrintf("NomadMex Evaluator: inconsistent size of list" );
//             //Force exit
//             force_quit();
//             return false;
//         }


//         //Check for Ctrl-C
//         if ( utIsInterruptPending() )
//         {
//             utSetInterruptPending(false); /* clear Ctrl-C status */
//             mexPrintf("\nCtrl-C Detected. Exiting NOMAD...\n\n");
//             list_count_eval.assign(m,false);
//             force_quit();
//             return false;
//         }


//         fun->prhs[fun->xrhs] = mxCreateDoubleMatrix(m, n, mxREAL); //x
//         double *List_x = mxGetPr(fun->prhs[fun->xrhs]);
//         std::list<NOMAD::EvalPoint *>::iterator it_x=x.begin();
//         j=0;
//         for (it_x=x.begin();it_x!=x.end();++it_x,++j)
//             for(i=0;i<n;i++)
//                 List_x[i*m+j] = (*(*it_x))[i].value();


//         //Add Surrogate if present and requested
//         if( hasSur )
//         {
//             sur=mxGetLogicals(fun->prhs[fun->xrhs+1]);
//             ( x.front()->get_eval_type()==NOMAD::SGTE )? *sur=true:*sur=false;  // all evaluations in a list have the same eval_type
//         }


//         // Count eval for bbox
//         // The case where the evaluation is rejected by user (and should not be counted) is not managed in the matlab version
//         list_count_eval.assign(m,true);


//         //Call MATLAB Objective
//         try
//         {
//             // Use Trap to catch some errors on fun eval that are not properly catched as an exception


// #ifdef BLACKBOX_COUNT_EVAL
//             mxArray * except = mexCallMATLABWithTrap(2, fun->plhs, fun->nrhs, fun->prhs, fun->f);
// #else
//             mxArray * except = mexCallMATLABWithTrap(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
// #endif
//             if ( except != NULL )
//             {
//                 std::string error_message ( "\n +++++++++++++++++++++++++++++++++++++++++++++++ \n");
//                 error_message += " Error message captured from blackbox: \n";
//                 error_message += mxArrayToString(mxGetProperty(except,0,"message"));
//                 error_message += "....... \n Correct this error in the blackbox or handle exception.\n";
//                 error_message += " +++++++++++++++++++++++++++++++++++++++++++++++ \n";
//                 throw ( runtime_error(error_message.c_str()) );
//             }
//         }

//         //Note if these errors occur it is due to errors in MATLAB code, no way to recover?
//         catch(exception &e)
//         {
//             sprintf(errstr,"Unrecoverable Error from Objective / Blackbox Callback:\n%sExiting NOMAD...\n\n",e.what());
//             mexPrintf(errstr);
//             //Force exit
//             force_quit();
//             return false;
//         }
//         catch(...)
//         {
//             mexPrintf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
//             //Force exit
//             force_quit();
//             return false;
//         }

//         //Check we got the correct number of elements back
//         if(mxGetNumberOfElements(fun->plhs[0]) > (nobj+ncon)*m)
//             mexPrintf("Black box returns more elements than required. Please provide a BB_OUTPUT_TYPE consistent with your black box function or correct the black box function.");
//         else if(mxGetNumberOfElements(fun->plhs[0]) < (nobj+ncon)*m)
//         {
//             mexPrintf("Insufficient outputs provided by the black box function. Exiting NOMAD...\n\n");
//             force_quit();
//             return false;
//         }
//         else if(mxGetM(fun->plhs[0]) != m )
//         {
//             mexPrintf("Insufficient number of rows in the output of the black box function. The number of rows should be equal to the size of the block of evaluations. Exiting NOMAD...\n\n");
//             force_quit();
//             return false;
//         }
//         else if(mxGetN(fun->plhs[0]) != (nobj+ncon) )
//         {
//             mexPrintf("Insufficient number of columns in the output of the black box function. The number of columns should be the number of objectives plus the number of constraints. Exiting NOMAD...\n\n");
//             force_quit();
//             return false;
//         }

//         //Assign count_eval
// #ifdef BLACKBOX_COUNT_EVAL
//         count=mxGetPr(fun->plhs[1]);
//         if ( count != NULL )
// 	{
// 	     list_count_eval.clear();
// 	}
// #endif

//         //Assign bb output
//         fvals = mxGetPr(fun->plhs[0]);
//         j=0;
//         for (it_x=x.begin();it_x!=x.end();++it_x,++j)
// 	{
// #ifdef BLACKBOX_COUNT_EVAL
//             list_count_eval.push_back( ( ( count[j] == 1) ? true : false) ) ;
// #endif

//             for(i=0;i<(nobj+ncon);i++)
//                 (*it_x)->set_bb_output(i,fvals[m*i+j]);
// 	}

//         //Iteration Callback
//         if(iterF->enabled)
//         {

// 	    // 1 iter counter for block of evals
//             iterF->plhs[0] = NULL;
//             iterF->prhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); //x
//             double *iter = mxGetPr(iterF->prhs[1]);
//             iter[0] = double(citer);
 
//             // m evals of obj+con      
//             iterF->prhs[2] = mxCreateDoubleMatrix(m, nobj+ncon, mxREAL); //x
//             double *dup_fvals = mxGetPr(iterF->prhs[2]);
//             for (j=0;j<m;++j)
//             	for(i=0;i<(nobj+ncon);i++)
//                 	dup_fvals[i*m+j] = fvals[m*i+j];

//             // m eval points
//             iterF->prhs[3] = mxCreateDoubleMatrix(m, n, mxREAL); //x
//             double *List_x = mxGetPr(iterF->prhs[3]);
//             std::list<NOMAD::EvalPoint *>::iterator it_x=x.begin();
//             j=0;
//             for (it_x=x.begin();it_x!=x.end();++it_x,++j)
//                 for(i=0;i<n;i++)
//                     List_x[i*m+j] = (*(*it_x))[i].value();

//             try {

//             	mxArray * except = mexCallMATLABWithTrap(1, iterF->plhs, 4, iterF->prhs, iterF->f);

//             	if ( except != NULL )
//             	{
//                 	std::string error_message ( "\n +++++++++++++++++++++++++++++++++++++++++++++++ \n");
//                 	error_message += " Error message captured from callback >(iterfun): \n";
//                 	error_message += mxArrayToString(mxGetProperty(except,0,"message"));
//                	 	error_message += "....... \n Correct this error in the callback (iterfun) or handle exception.\n";
//                 	error_message += " +++++++++++++++++++++++++++++++++++++++++++++++ \n";
//              		   throw ( runtime_error(error_message.c_str()) );
//             	}
// 	    }
//             //Note if these errors occur it is due to errors in MATLAB code, no way to recover?
//             catch(exception &e)
//             {
//             	sprintf(errstr,"Unrecoverable Error from Objective / Blackbox Callback:\n%sExiting NOMAD...\n\n",e.what());
//             	mexPrintf(errstr);
//             	//Force exit
//             	force_quit();
//             	return false;
//             }
//             catch(...)
//             {
//             	mexPrintf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
//            	 //Force exit
//             	force_quit();
//             	return false;
//             }

//             //Collect return argument
//             stop = *(bool*)mxGetData(iterF->plhs[0]);

//         }

//         //Add Function Eval Counter
//         citer++;

//         // Clean up LHS Fun Ptr
//         mxDestroyArray(fun->plhs[0]);

// #ifdef BLACKBOX_COUNT_EVAL
// 	mxDestroyArray(fun->plhs[1]);
// #endif

//         //Check for iterfun stop
//         if(stop)
//         {
//             mexPrintf("\nIterFun Called Stop. Exiting NOMAD...\n\n");
//             force_quit();
//             return false;
//         }
//         else
//             return true;
//     }



    //Function + Constraint Evaluation
    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &h_max, bool &count_eval) const
    {
        char errstr[1024];
        bool stop = false;
        int i, n = static_cast<int>(x.size());
        double *xm, *fvals;
        // TODO sur mxLogical *sur;
        double *count=NULL;
        count_eval = true; 

		if (NOMAD::Step::getUserTerminate())
		{
			NOMAD::Step::setUserTerminate(); // Two calls to user terminate maybe necessary
			count_eval = false;
			return false;
		}

        //Check for Ctrl-C
        if ( utIsInterruptPending() )
        {
            utSetInterruptPending(false); /* clear Ctrl-C status */
            mexPrintf("\nCtrl-C Detected. Exiting NOMAD...\n\n");
            count_eval = false;
            NOMAD::Step::setUserTerminate();
            return false;
         }

         //Blackbox / Objective Evaluation
         xm = mxGetPr(fun->prhs[fun->xrhs]);
         for(i=0;i<n;i++)
             xm[i] = x[i].todouble();

// TODO
//         //Add Surrogate if present and requested
//         if( hasSur )
//         {
//             sur=mxGetLogicals(fun->prhs[fun->xrhs+1]);
//             (x.get_eval_type()==NOMAD::SGTE)? *sur=true:*sur=false;
//         }

         //Call MATLAB Objective
         try
         {
             // Use Trap to catch some errors on fun eval that are not properly catched as an exception
#ifdef BLACKBOX_COUNT_EVAL
            mxArray * except = mexCallMATLABWithTrap(2, fun->plhs, fun->nrhs, fun->prhs, fun->f);
#else
 			mxArray * except = mexCallMATLABWithTrap(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
#endif

            // Counting eval if Matlab bb has been called
            count_eval = true;
            if ( except != NULL )
            {
                std::string error_message ( "\n +++++++++++++++++++++++++++++++++++++++++++++++ \n");
                error_message += " Error message captured from blackbox: \n";
                error_message += mxArrayToString(mxGetProperty(except,0,"message"));
                error_message += "....... \n Correct this error in the blackbox or handle exception.\n";
                error_message += " +++++++++++++++++++++++++++++++++++++++++++++++ \n";
                throw ( runtime_error(error_message.c_str()) );
            }
         }
         //Note if these errors occur it is due to errors in MATLAB code, no way to recover?
         catch(exception &e)
         {
             sprintf(errstr,"Unrecoverable Error from Objective / Blackbox Callback:\n%sExiting NOMAD...\n\n",e.what());
             mexErrMsgTxt(errstr);
             //Force exit
             NOMAD::Step::setUserTerminate();
             return false;
         }
         catch(...)
         {
             mexPrintf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
             //Force exit
             NOMAD::Step::setUserTerminate();
             return false;
         }

         //Check we got the correct number of elements back

         if( mxGetNumberOfElements(fun->plhs[0]) > nobj+ncon )
             mexWarnMsgTxt("Black box returns more elements than required. Please provide a BB_OUTPUT_TYPE consistent with your black box function");
         else if( mxGetNumberOfElements(fun->plhs[0]) < nobj+ncon )
         {
             mexPrintf("Insufficient outputs provided by the black box function. Exiting NOMAD...\n\n");
             NOMAD::Step::setUserTerminate();
             return false;
         }
         //Assign bb output
         fvals = mxGetPr(fun->plhs[0]);

         //Assign count_eval
#ifdef BLACKBOX_COUNT_EVAL
         count=mxGetPr(fun->plhs[1]);
#endif

         if ( count != NULL )
 		{
             count_eval = (*count == 1) ? true : false ;
			 mxDestroyArray(fun->plhs[1]);
		 }

         std::string bboStr;
         for(i=0;i<(nobj+ncon);i++)
         {
             NOMAD::Double bbo = fvals[i];
             bboStr += bbo.tostring() ;   
         }
         x.setBBO(bboStr);
		 // mexPrintf("Function output: %s\n", bboStr);

		 //Iteration Callback
		 if (nullptr != callbackF && callbackF->enabled)
		 {
			 // mexPrintf("Callback function about to be called\n");

			 callbackF->plhs[0] = nullptr;

			 callbackF->prhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
			 callbackF->prhs[2] = mxCreateDoubleMatrix(1, (nobj + ncon), mxREAL);
			 callbackF->prhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);

			 memcpy(mxGetData(callbackF->prhs[1]), &counter_eval, sizeof(int));
			 memcpy(mxGetPr(callbackF->prhs[2]), fvals, (nobj + ncon) * sizeof(double));
			 memcpy(mxGetPr(callbackF->prhs[3]), xm, n * sizeof(double));
			 try {
				 mexCallMATLAB(1, callbackF->plhs, 4, callbackF->prhs, callbackF->f);
			 }
			 catch (...)
			 {
				 mexPrintf("Unrecoverable error from user eval callback. Callback function [stop] = cbFun(eval_counter,fevals,x) must be provided. Exiting NOMAD...\n\n");
				 //Force exit
				 NOMAD::Step::setUserTerminate();
				 return false;
			 }

			 //Collect return argument
			 stop = *(bool*)mxGetData(callbackF->plhs[0]);

			 //Clean up Ptr
			 mxDestroyArray(callbackF->plhs[0]);

		 }

		 //Add Function Eval Counter
		 counter_eval ++;

		 // Clean up LHS Fun Ptr
		 mxDestroyArray(fun->plhs[0]);

		 //Check for iterfun stop
		 if (stop)
		 {
			 mexPrintf("\nUser evaluation callback called for a stop. Exiting NOMAD with a CTRL-C signal ...\n\n");
			 NOMAD::Step::setUserTerminate();
		 }
		 
		 return true;

    }
};



class mxstreambuf : public std::streambuf {
public:
	mxstreambuf() {
		stdoutbuf = std::cout.rdbuf(this);
	}
	~mxstreambuf() {
		std::cout.rdbuf(stdoutbuf);
	}
protected:
	virtual std::streamsize xsputn(const char* s, std::streamsize n) override {
		mexPrintf("%.*s", n, s);
		return n;
	}
	virtual int overflow(int c = EOF) override {
		if (c != EOF) {
			mexPrintf("%.1s", &c);
		}
		return 1;
	}
private:
	std::streambuf *stdoutbuf;
};


// Main Entry Function
// -----------------------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// For cout redirection
	mxstreambuf mout;

	//Input Args
	usrFcn fun;
	double *x0, *lb = NULL, *ub = NULL;



	//Internal Vars
	size_t ndec;
	int i;
	size_t nobj = 1, ncon = 0;
	char errstr[1024]; //used for returning error info

	//Check user inputs
	if (!checkInputs(prhs, nrhs, plhs, nlhs))
		return;


	//Get Size
	ndec = mxGetNumberOfElements(pX0);

	//Get Blackbox / Objective Function Handle
	if (mxIsChar(pFUN))
	{
		if (mxGetString(pFUN, fun.f, FLEN) != 0)
			mexErrMsgTxt("\n Error reading objective name as string");
		fun.nrhs = 1;
		fun.xrhs = 0;
	}
	else
	{
		fun.prhs[0] = (mxArray*)pFUN;
		strcpy(fun.f, "feval");
		fun.nrhs = 2;
		fun.xrhs = 1;
	}
	fun.prhs[fun.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x


	usrCallbackFcn callbackFun;
	callbackFun.enabled = false;
	if (nrhs > eCB && !mxIsEmpty(pCB))
	{
		mexPrintf("Callback function is enabled");
		//Get Blackbox / Objective Function Handle
		if (mxIsChar(pCB))
		{
			if (mxGetString(pCB, callbackFun.f, FLEN) != 0)
				mexErrMsgTxt("\n Error reading eval callback function name as string");
		}
		else
		{
			callbackFun.prhs[0] = (mxArray*)pCB;
			strcpy(callbackFun.f, "feval");
		}
		callbackFun.enabled = true;

	}



    //Get x0
    x0 = mxGetPr(pX0);

    NOMAD::MainStep::resetComponentsBetweenOptimization();

    // Ready to set all parameters
    std::shared_ptr<NOMAD::AllParameters> p = std::make_shared<NOMAD::AllParameters>();

    //Setup ndec
    p->setAttributeValue("DIMENSION",ndec);

    //Warn if >1000
    if(ndec > 1000)
    {
        sprintf(errstr,"Warning: NOMAD is designed for problems with less than 1000 variables. Your model has %d.\nWhile unlikely, it is possible that NOMAD may not perform as intended on this problem.",static_cast<int>(ndec));
        mexWarnMsgTxt(errstr);
    }

    //Setup Lower Bounds
    if( nrhs > eLB && !mxIsEmpty(pLB) )
    {
        NOMAD::ArrayOfDouble aodLB(ndec);
        lb = mxGetPr(pLB);

        for(i=0;i<ndec;i++)
        {
            if(!mxIsInf(lb[i])) //if not initialized will not be used
            {
                aodLB[i] = lb[i];
            }
        }
        p->setAttributeValue("LOWER_BOUND",aodLB);
    }

    //Setup Upper Bounds
    if( nrhs > eUB && !mxIsEmpty(pUB) )
    {
        NOMAD::ArrayOfDouble aodUB(ndec);
        ub = mxGetPr(pUB);

        for(i=0;i<ndec;i++)
        {
            if(!mxIsInf(ub[i])) //if not initialized will not be used
            {
                aodUB[i] = ub[i];
            }
        }
        p->setAttributeValue("UPPER_BOUND",aodUB);
    }

    //Setup x0
    NOMAD::Point px0 (ndec);
    for(i=0;i<ndec;i++)
        px0[i] = x0[i];
    
    p->setAttributeValue("X0",px0);
    
    //Check NOMAD parameters
    try
    {
        //Set User Options
        if( nrhs > ePARAM && !mxIsEmpty(pParam) )
            setNomadParams(p,pParam);

        p->checkAndComply();
    }
    catch(exception &e)
    {
        sprintf(errstr,"\n NOMAD Parameter Error:\n%s",e.what());
        mexErrMsgTxt(errstr);
        return;
    }

    // Obtain number of objectives and constraints from parameters
    auto bboutputlist = p->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    nobj=NOMAD::getNbObj(bboutputlist);  // Only single objective is supported in Nomad 4
    ncon=bboutputlist.size()-nobj;
    // TODO
    //if ( p.has_sgte() )
    //{
    //     fun.prhs[fun.xrhs+1] = mxCreateLogicalMatrix(1,1); //extra logical indicating surrogate or not
    //     fun.nrhs++;
    //}

    // TODO
    // // Get additional bb param if specified
    // if ( nrhs > ePARAM && ! mxIsEmpty(pParam) )
    // {
    //     fun.nrhs++;

    //     if ( p.has_sgte() )
    //         fun.prhs[fun.xrhs+2]=mxDuplicateArray(pParam);
    //     else
    //         fun.prhs[fun.xrhs+1]=mxDuplicateArray(pParam);
    // }

    //Print Header
    if(p->getAttributeValue<int>("DISPLAY_DEGREE") > 1)
    {
         mexPrintf("\n------------------------------------------------------------------\n");
         mexPrintf(" This is NOMAD v%s\n",NOMAD_VERSION_NUMBER);
         mexPrintf(" Authors: C. Audet, S. Le Digabel, V. Rochon Montplaisir and C. Tribes\n");
         mexPrintf(" MEX Interface C. Tribes 2021 \n\n");
         mexPrintf(" Problem Properties:\n");
         mexPrintf(" # Decision Variables:               %4d\n",ndec);
         mexPrintf(" # Number of Objectives:             %4d\n",nobj);
         mexPrintf(" # Number of Nonlinear Constraints:  %4d\n",ncon);
         mexPrintf("------------------------------------------------------------------\n");
         mexEvalString("drawnow;"); //flush draw buffer
    }

    //Create evaluator and run mads based on number of objectives
    try
    {
        //NOMAD Vars
        NOMAD::MainStep mainStep;
        mainStep.setAllParameters(p); 

	    std::unique_ptr<matlabEval> ev(new matlabEval(p->getEvalParams(),&fun,nobj,ncon,&callbackFun));
        mainStep.setEvaluator(std::move(ev));

        mainStep.start();
        mainStep.run();
        mainStep.end();
    
    }
    catch(exception &e)
    {
        sprintf(errstr,"\n NOMAD run error:\n%s",e.what());
        mexErrMsgTxt(errstr);
    }

    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
    
	double *x= mxGetPr(plhs[0]);
	double *fval = mxGetPr(plhs[1]);
	double *hinf = mxGetPr(plhs[2]);
    double *exitflag = mxGetPr(plhs[3]);
    double *nfval = mxGetPr(plhs[4]);

    std::vector<NOMAD::EvalPoint> bf, bi ;
    NOMAD::EvalPoint sol;  

    size_t nbBestFeas = NOMAD::CacheBase::getInstance()->findBestFeas(bf, NOMAD::Point(ndec), NOMAD::EvalType::BB,NOMAD::ComputeType::STANDARD, nullptr);
    size_t nbBestInf  = NOMAD::CacheBase::getInstance()->findBestInf(bi, NOMAD::INF, NOMAD::Point(ndec), NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD,nullptr);

    if (nbBestFeas == 0)
    {
        if (nbBestInf == 0)
        { 
			mexPrintf("NO solution obtained\n");
            *exitflag = -1; //No solution
			*fval = 0;
			*hinf = 1;
			for (i = 0; i < ndec; i++)
				x[i] = x0[i];
            return;
        }
        else
        {
			mexPrintf("Infeasible solution obtained\n");
            sol = bi[0];
        }
    }
    else
    {
		mexPrintf("Feasible solution obtained\n");
        sol = bf[0];
    }

	// TODO Use exit flag to indicate the stopping criterion

    //Save x
    NOMAD::Point pt(*sol.getX());
    for(i=0;i<ndec;i++)
        x[i] = pt[i].todouble();
    
    //Save fval
    *fval = sol.getF().todouble();
    
    //Save hinf
    *hinf = sol.getH().todouble();

    //Save Status & Iterations
    *exitflag = 0;

	*nfval = (double) NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();

	// TODO return a structure with stats from run

    // Clean up of fun
    mxDestroyArray(fun.prhs[fun.xrhs]);
}



//Convert input string to lowercase
void lower(char *str)
{
    int i = 0;
    while(str[i])
    {
        str[i] = tolower(str[i]);
        i++;
    }
}


//User Input Checking + Version / Info / Help
int checkInputs(const mxArray *prhs[], int nrhs, mxArray *plhs[], int nlhs)
{
    size_t ndec;
    char *str = NULL;

    NOMAD::MainStep mainStep;

	//Redirect cout
	mxstreambuf mout;

    //MEX Display Version
    if (nrhs < 1)
    {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(NOMAD_VERSION_NUMBER);
        return 0;
    }

    

    //Check for display on options passed as structure
    if(nrhs == 1 && mxIsStruct(prhs[0]))
    {
        int i, no = mxGetNumberOfFields(prhs[0]);
        const char *field;
        //For all fields, display nomad help
        for(i=0;i<no;i++)
        {
            field = mxGetFieldNameByNumber(prhs[0],i);
            string st(field);
            mainStep.displayHelp(st);
        }
        return 0;
    }

    //Check for Version / Information / Help Request
    if(nrhs == 1 && mxIsChar(prhs[0]))
    {
        str = mxArrayToString(prhs[0]);
        //Check for Info Request
        if(!strcmp(str,"-I") || !strcmp(str,"-INFO") || !strcmp(str,"-i") || !strcmp(str,"-info"))
        {
            mainStep.displayInfo();
            mainStep.displayUsage("");

            return 0;
        }
        if(!strcmp(str,"-U") || !strcmp(str,"-USAGE") || !strcmp(str,"-u") || !strcmp(str,"-usage"))
        {
            mainStep.displayInfo();
            mainStep.displayUsage("");

            return 0;
        }
        //Check for Ver Request
        if(!strcmp(str,"-V") || !strcmp(str,"-v") || !strcmp(str,"-version"))
        {
            mainStep.displayVersion();
            mexPrintf("MEX Interface (GERAD) v%s\n",NOMADMEX_VERSION);

            return 0;
        }
        //Check for Help Request
        if (strcmp(str,"-H")<0 || strcmp(str,"-HELP")<0 || strcmp(str,"-h")<0 || strcmp(str,"-help")<0 )
        {
            const char * toks=" ";
            char *w = strtok(str,toks) ;
            for ( w = strtok(NULL,toks) ; w != NULL ; w = strtok(NULL,toks) )
            {
                mainStep.displayHelp(w);
            }

            return 0;
        }
    }

    //Otherwise assume we have a normal problem
    if(nrhs < 2)
        mexErrMsgTxt("You must supply at least 2 arguments to nomad!\n\nnomadOpt(fun,x0)\n");

    //Check Types
    if(!mxIsFunctionHandle(pFUN) && !mxIsChar(pFUN))
        mexErrMsgTxt("fun must be a function handle or function name!");

    if(!mxIsDouble(pX0) || mxIsComplex(pX0) || mxIsEmpty(pX0))
        mexErrMsgTxt("x0 must be a real double column vector!");

    //Get ndec
    ndec = mxGetNumberOfElements(prhs[1]);

    //Check Bounds
    if(nrhs > 2)
    {
        if(!mxIsDouble(pLB) || mxIsComplex(pLB))
            mexErrMsgTxt("lb must be a real double column vector!");
        if(nrhs > 3 && (!mxIsDouble(pUB) || mxIsComplex(pUB)))
            mexErrMsgTxt("ub must be a real double column vector!");
        //Check Sizes
        if(!mxIsEmpty(pLB) && (ndec != mxGetNumberOfElements(pLB)))
            mexErrMsgTxt("lb is not the same length as x0! Ensure they are both Column Vectors");
        if(nrhs > 3 && !mxIsEmpty(pUB) && (ndec != mxGetNumberOfElements(pUB)))
            mexErrMsgTxt("ub is not the same length as x0! Ensure they are both Column Vectors");
    }

    //Version check
	if (nrhs > eCB+1)
		mexErrMsgTxt("Calling NOMAD using the wrong syntax. Correct syntax: nomadOpt(bb,x0,lb,ub,params,bbCallBackFun). The last argument is optional.");

    //Return Continue
    return  1;

}

void setNomadParams(const std::shared_ptr<NOMAD::AllParameters> p, const mxArray *params)
{

	// For cout redirection
	mxstreambuf  mout;

    char strbuf[1024];
    int i, no = 0;
    mxArray *value;

	bool doAdd = false;

    if(!mxIsStruct(params))
    {
        mexErrMsgTxt("Params should be provided as Matlab structure: struct('KEY1','VAL1','KEY2','VAL2') "); 
    }

    if(params)
    {    
        no = mxGetNumberOfFields(params);
    }
    else
    {
        return;
    }
    
    //For each field, check if it's empty, if not, set it within NOMAD
    for(i=0;i<no;i++)
    {
       value = mxGetFieldByNumber(params,0,i);
       if(!mxIsChar(value))
       {
          mexErrMsgTxt("Params should be provided as Matlab structure with string value: struct('KEY1','VAL1','KEY2','VAL2') ");
       }
       const char * keyword = mxGetFieldNameByNumber(params,i);
       sprintf(strbuf,"%s %s",keyword,mxArrayToString(value));
        // std::cout<<strbuf << std::endl;
       p->readParamLine(strbuf);
    }       
}


double getStatus(int stat)
{

   // TODO in Nomad 4 we don't have a stop flag

    switch((int)stat)
    {
        case 5:         //mesh minimum
        case 8:         //min mesh size
        case 9:         //min poll size
        case 20:        //ftarget reached
        case 19:        //feas reached
            return 1;
            break;
        case 12:        //max time
        case 13:        //max bb eval
        case 14:        //max sur eval
        case 15:        //max evals
        case 16:        //max sim bb evals
        case 17:        //max iter
        case 23:        //max multi bb evals
        case 24:        //max mads runs
            return 0;
            break;
        case 10:        //max mesh index
        case 11:        //mesh index limits
        case 18:        //max consec fails
        case 25:        //stagnation
        case 26:        //no pareto
        case 27:        //max cache memory
            return -1;
        case 6:         //x0 fail
        case 7:         //p1_fail
            return -2;
        case 2:         //Unknown stop reason
            return -3;
        case 3:         //Ctrl-C
        case 4:         //User-stopped
            return -5;
        default:        //Not assigned flag
            return -3;
    }
}

//Print Solver Information
void printSolverInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" NOMAD: Nonlinear Optimization using the MADS Algorithm [v%s]\n",NOMAD_VERSION_NUMBER);
    mexPrintf("  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html\n");
    mexPrintf("  - Source available from: https://www.gerad.ca/nomad/\n");

    mexPrintf("\n MEX Interface C. Tribes 2021  \n");
    mexPrintf("-----------------------------------------------------------\n");
}

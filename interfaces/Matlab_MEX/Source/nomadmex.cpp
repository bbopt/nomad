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

#define NOMADMEX_VERSION "1.6  [Jan, 15th, 2025]"
//NOTE The OPTI interface is no longer supported.
// Another change compared with previous versions: the parameters are passed as a list of strings (in the form struct('KEYWORD','value',...)). The keyword and value follow exactly the Nomad syntax.

// The GERAD interface is now:
//      [x,fval,hinf,runflag,nfval]= nomadOPt(fun,x0,lb,ub,params,callbackFun)


#include "mex.h"
#include "Nomad/nomad.hpp"
#include "Cache/CacheBase.hpp"
#include "Algos/EvcInterface.hpp"
#include "Type/EvalType.hpp"
#include <stdio.h>
#include <string.h>
#include <memory.h>


using namespace std;


//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 5 /* max nrhs for user function */

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

int counter_eval = 0;

//MATLAB Evaluator Class (BB or SURROGATE)
class matlabEval : public NOMAD::Evaluator
{
private:
    usrFcn *fun;
    size_t nbbo;

    int indexBBOForCountEval;

    usrCallbackFcn *callbackF;


public:
    //Constructor

    matlabEval(const std::shared_ptr<NOMAD::EvalParameters> p, usrFcn *_fun, size_t _nbbo, NOMAD::EvalType evalType, usrCallbackFcn * _callbackF, int _indexBBOForCountEval ) : NOMAD::Evaluator(p, evalType)
    {
        fun     = _fun;
        nbbo    = _nbbo;
        counter_eval = 0;
        callbackF = _callbackF;
        indexBBOForCountEval = _indexBBOForCountEval ;

    }
    //Destructor
    ~matlabEval(void) {}

    //Function + Constraint Evaluation on a given block of eval points
    std::vector<bool> eval_block(NOMAD::Block &block,
                                                   const NOMAD::Double &hMax,
                                                   std::vector<bool> &countEval) const
    {
        
        size_t m = block.size();
        std::vector<bool> evalOk(m, false);
        countEval.resize(m, false);
        
        //mexPrintf("Optimization using block of points evaluation is not available\n");
        //NOMAD::Step::setUserTerminate();
        
        
        if (NOMAD::Step::getUserTerminate())
        {
            NOMAD::Step::setUserTerminate(); // Two calls to user terminate maybe necessary
            return evalOk;
        }
        
        //Check for Ctrl-C
        if ( utIsInterruptPending() )
        {
            utSetInterruptPending(false); /* clear Ctrl-C status */
            mexPrintf("\nCtrl-C Detected. Exiting NOMAD...\n\n");
            
            NOMAD::Step::setUserTerminate();
            return evalOk;
        }
        
        
        int n = static_cast<int>(block[0]->size());
        fun->prhs[fun->xrhs] = mxCreateDoubleMatrix(m, n, mxREAL); //x
        double *List_x = mxGetPr(fun->prhs[fun->xrhs]);
        
        mxLogical *sur;
        size_t j=0;
        for (size_t j = 0 ; j < block.size() ; j++)
        {
            for(size_t i=0;i<n;i++)
            {
                List_x[i*m+j] = (*block[j])[i].todouble();
            }
        }
        
        //Add a flag if surrogate eval
        if( NOMAD::EvalType::SURROGATE == _evalType )
        {
            sur=mxGetLogicals(fun->prhs[fun->xrhs+1]);
            *sur=true;
        }
        
        char errstr[1024];
        //Call MATLAB Objective
        try
        {
            // Count eval for bbox
            // An evaluation can be counted/unaccounted by setting BB_OUTPUT_TYPE CNT_EVAL and return the correct flag in the bb function (see example_count_eval).
            // Id exception is catched, we count the eval.
            countEval.assign(m,true);
            
            // Use call with Trap to catch some errors on fun eval that are not properly catched as an exception
            mxArray * except = mexCallMATLABWithTrap(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
            
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
            return evalOk;
        }
        catch(...)
        {
            mexPrintf("Unrecoverable Error from Objective / Blackbox Callback, Exiting NOMAD...\n\n");
            //Force exit
            NOMAD::Step::setUserTerminate();
            return evalOk;
        }
        
        
        //Check we got the correct number of elements back
        if(mxGetNumberOfElements(fun->plhs[0]) > (nbbo)*m)
            mexPrintf("Black box returns more elements than required. Please provide a BB_OUTPUT_TYPE consistent with your black box function or correct the black box function.");
        else if(mxGetNumberOfElements(fun->plhs[0]) < (nbbo)*m)
        {
            mexPrintf("Insufficient outputs provided by the black box function. Exiting NOMAD...\n\n");
            NOMAD::Step::setUserTerminate();
            return evalOk;
        }
        else if(mxGetM(fun->plhs[0]) != m )
        {
            mexPrintf("Insufficient number of rows in the output of the black box function. The number of rows should be equal to the size of the block of evaluations. Exiting NOMAD...\n\n");
            NOMAD::Step::setUserTerminate();
            return evalOk;
        }
        else if(mxGetN(fun->plhs[0]) != (nbbo) )
        {
            mexPrintf("Insufficient number of columns in the output of the black box function. The number of columns should be the number of objectives plus the number of constraints. Exiting NOMAD...\n\n");
            NOMAD::Step::setUserTerminate();
            return evalOk;
        }
        
        // Update the MatlabEvaluator counter
        counter_eval += m;
        
        //Assign bb output
        double *fvals = mxGetPr(fun->plhs[0]);
        j=0;
        for (size_t j =0 ; j < block.size() ; j++)
        {            
            std::string bboStr;
            for(size_t i=0;i<nbbo;i++)
            {
                NOMAD::Double bbo = fvals[m*i+j];
                bboStr += bbo.tostring() + " " ;
                // Modify eval count if flag (CNT_EVAL 0) is provided by the user.
                if (indexBBOForCountEval>=0 && i == indexBBOForCountEval)
                {
                    if (bbo == 0)
                    {
                        countEval[j]= false ;
                    }
                }
            }
            block[j]->setBBO(bboStr,_bbOutputTypeList, _evalType);
        }
        
        // We have passed the exception. The eval is ok.
        evalOk.assign(m,true);
        
        //Iteration Callback
        if (nullptr != callbackF && callbackF->enabled)
        {
            callbackF->plhs[0] = nullptr;
            
            callbackF->prhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            callbackF->prhs[2] = mxCreateDoubleMatrix(m, nbbo, mxREAL);
            callbackF->prhs[3] = mxCreateDoubleMatrix(m, n, mxREAL);
            
            memcpy(mxGetData(callbackF->prhs[1]), &counter_eval, sizeof(int));
            memcpy(mxGetPr(callbackF->prhs[2]), fvals, m*nbbo * sizeof(double));
            memcpy(mxGetPr(callbackF->prhs[3]), List_x, m*n * sizeof(double));
            try {
                mexCallMATLAB(1, callbackF->plhs, 4, callbackF->prhs, callbackF->f);
            }
            catch (...)
            {
                mexPrintf("Unrecoverable error from user eval callback. Callback function [stop] = cbFun(eval_counter,fevals,x) must be provided. Exiting NOMAD...\n\n");
                //Force exit
                NOMAD::Step::setUserTerminate();
                return evalOk;
            }
            
            //Collect return argument
            bool stop = *(bool*)mxGetData(callbackF->plhs[0]);
            //Check for iterfun stop
            if (stop)
            {
                mexPrintf("\nUser evaluation callback called for a stop. Exiting NOMAD with a CTRL-C signal ...\n\n");
                NOMAD::Step::setUserTerminate();
            }
        }
        
        // Clean up LHS Fun Ptr
        mxDestroyArray(fun->plhs[0]);
        
        return evalOk;
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
    size_t nbbo;
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
    nbbo=bboutputlist.size();

    int indexBBOForCountEval = -1; // Defautl is "no eval count eval"
    for (size_t i = 0; i < bboutputlist.size(); i++)
    {
        if (bboutputlist[i] == NOMAD::BBOutputType::Type::CNT_EVAL)
        {
            indexBBOForCountEval = i;
            break;
        }
    }

    size_t nobj=NOMAD::getNbObj(bboutputlist);  
    size_t ncon=NOMAD::getNbConstraints(bboutputlist);
    nbbo=bboutputlist.size();
    
    // Test for DiscoMads -> DiscoMads increases the number of constraints.
    if (p->getAttributeValue<bool>("DISCO_MADS_OPTIMIZATION"))
    {
        ncon--;
    }

    
    // Read fun output for surrogate use
    if ( p->mayUseSurrogate() )
    {
         fun.prhs[fun.xrhs+1] = mxCreateLogicalMatrix(1,1); //extra logical indicating surrogate or not
         fun.nrhs++;
    }

    //Print Header
    if(p->getAttributeValue<int>("DISPLAY_DEGREE") > 1)
    {
         mexPrintf("\n------------------------------------------------------------------\n");
         mexPrintf(" This is NOMAD v%s\n",NOMAD_VERSION_NUMBER);
         mexPrintf(" Authors: C. Audet, S. Le Digabel, V. Rochon Montplaisir and C. Tribes\n");
         mexPrintf(" MEX Interface C. Tribes 2023 \n\n");
         mexPrintf(" Problem Properties:\n");
         mexPrintf(" # Decision Variables:               %4d\n",ndec);
         mexPrintf(" # Number of Objectives:             %4d\n",nobj);
         mexPrintf(" # Number of Nonlinear Constraints:  %4d\n",ncon);
         mexPrintf("------------------------------------------------------------------\n");
         mexEvalString("drawnow;"); //flush draw buffer
    }

    int mainStepRunFlag = -3;
    
    //Create evaluator and run mads based on number of objectives
    try
    {
        //NOMAD Vars
        NOMAD::MainStep mainStep;
        mainStep.setAllParameters(p);

        auto evBB = std::make_unique<matlabEval>(p->getEvalParams(),&fun,nbbo,NOMAD::EvalType::BB, &callbackFun, indexBBOForCountEval);
        mainStep.addEvaluator(std::move(evBB));
        
        // The same matlab evaluator is used for surrgate and bb. A flag is passed to the matlab BB callback function
    
        if (p->mayUseSurrogate())
        {
            auto evSurrogate = std::make_unique<matlabEval>(p->getEvalParams(),&fun,nbbo,NOMAD::EvalType::SURROGATE, &callbackFun, indexBBOForCountEval);
            mainStep.addEvaluator(std::move(evSurrogate));
        }
        
        mainStep.start();
        mainStep.run();
        mainStep.end();
        
        mainStepRunFlag = mainStep.getRunFlag();
    
    }
    catch(exception &e)
    {
        sprintf(errstr,"\n NOMAD run error:\n%s",e.what());
        mexErrMsgTxt(errstr);
        return;
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
    double *runflag = mxGetPr(plhs[3]);
    double *nfval = mxGetPr(plhs[4]);

    std::vector<NOMAD::EvalPoint> bf, bi ;
    NOMAD::EvalPoint sol;

    size_t nbBestFeas = NOMAD::CacheBase::getInstance()->findBestFeas(bf, NOMAD::Point(ndec), NOMAD::defaultFHComputeType);
    size_t nbBestInf  = NOMAD::CacheBase::getInstance()->findBestInf(bi, NOMAD::INF, NOMAD::Point(ndec), NOMAD::defaultFHComputeType);

    // For now
    // If nbFeas > 0 we return a single best feasible point (no infeasible point)
    // Else (if nbFeas == 0) we return the least infeasible point with the smallest f (index 0, see findBestInf)
    // Maybe this could be generalized to show the best feasible point and all undominated infeasible points.
    // The same logic for Nomad C and Matlab interfaces and for PyNomad.
    if (nbBestFeas == 0)
    {
        if (nbBestInf == 0)
        {
            mexPrintf("NO solution obtained\n");
            *runflag = double(mainStepRunFlag) ;
            *fval = 0;
            *hinf = 1;
            for (i = 0; i < ndec; i++)
                x[i] = x0[i];
            return;
        }
        else
        {
            mexPrintf("Least infeasible solution obtained (with smallest f)\n");
            sol = bi[0];
        }
    }
    else
    {
        mexPrintf("Feasible solution obtained\n");
        sol = bf[0];
    }

    //Save x
    NOMAD::Point pt(*sol.getX());
    for(i=0;i<ndec;i++)
        x[i] = pt[i].todouble();
    
    //Save fval
    *fval = sol.getF(NOMAD::defaultFHComputeType).todouble();
    
    //Save hinf
    *hinf = sol.getH(NOMAD::defaultFHComputeType).todouble();

    //Run flag
    *runflag = double(mainStepRunFlag);

    *nfval = (double) NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();

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
            NOMAD::toupper(st);
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
                std::string st(w);
                NOMAD::toupper(st);
                mainStep.displayHelp(st);
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

//Print Solver Information
void printSolverInfo()
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" NOMAD: Nonlinear Optimization using the MADS Algorithm [v%s]\n",NOMAD_VERSION_NUMBER);
    mexPrintf("  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html\n");
    mexPrintf("  - Source available from: https://www.gerad.ca/nomad/\n");

    mexPrintf("\n MEX Interface C. Tribes updated in 2025  \n");
    mexPrintf("-----------------------------------------------------------\n");
}

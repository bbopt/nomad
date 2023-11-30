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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverOptimize.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Eval/ComputeSuccessType.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"
#include "../../Type/EvalSortType.hpp"
#include "../../Math/MatrixUtils.hpp"

#include "../../../ext/sgtelib/src/Surrogate_PRS.hpp"
#include "QPSolverOptimize.hpp"

NOMAD::EvalPointPtr NOMAD::QPSolverOptimize::_prevFeasRefCenter = nullptr;
NOMAD::EvalPointPtr NOMAD::QPSolverOptimize::_prevInfeasRefCenter = nullptr;
NOMAD::Point NOMAD::QPSolverOptimize::_prevFeasXopt = NOMAD::Point();
NOMAD::Point NOMAD::QPSolverOptimize::_prevInfeasXopt = NOMAD::Point();

void NOMAD::QPSolverOptimize::init()
{
    setStepType(NOMAD::StepType::MODEL_OPTIMIZE);
    verifyParentNotNull();
    
    // Set the bounds and fixed variables from the model
    setModelBoundsAndFixedVar();
    
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    _bbot = evc->getCurrentEvalParams()->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    _m = static_cast<int>(_bbot.size());
    _nbCons = static_cast<int>(getNbConstraints(_bbot));
    
    _quadModelMaxEval = evc->getEvaluatorControlGlobalParams()->getAttributeValue<size_t>("QUAD_MODEL_MAX_EVAL");

    if ( _modelFixedVar.nbDefined() == _modelFixedVar.size() )
    {
        OUTPUT_INFO_START
        std::ostringstream oss;
        oss << "Effective dimension is null. No QuadModelOptimize" << std::endl;
        AddOutputInfo(oss.str(), NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
        OUTPUT_INFO_END

        return;
    }
    
    
    /// Access to model coefficients for Surrogate_PRS
    if (nullptr == _model)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "QPSolverOptimize: a model is required (nullptr)");
    }
    _model->check_ready(__FILE__,__FUNCTION__,__LINE__);
    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    if ( nullptr != surrogate_prs)
    {
        SGTELIB::Matrix model_coef = surrogate_prs->get_alpha();
        SGTELIB::Matrix  model_monomes = surrogate_prs->get_PRS_monomes(static_cast<int>(_n), 2);
        
        if (model_coef.get_nb_cols() != _m)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Number of cols in polynom coefficients do not match number of models required (bbo)");
        }
        
        OUTPUT_INFO_START
        std::ostringstream os;
        model_coef.display(os);
        NOMAD::OutputQueue::Add(os.str(), _displayLevel);
        OUTPUT_INFO_END
    }
    
    _verbose = _runParams->getAttributeValue<bool>("QP_verbose");
    _verboseFull = _runParams->getAttributeValue<bool>("QP_verboseFull");
    
}


void NOMAD::QPSolverOptimize::startImp()
{
    auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");
    _displayLevel = (std::string::npos != modelDisplay.find("O"))
        ? NOMAD::OutputLevel::LEVEL_INFO
        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

    
    OUTPUT_INFO_START
    std::string s;
    s = "QUAD_MODEL_MAX_EVAL: " + std::to_string(_quadModelMaxEval);
    AddOutputInfo(s, _displayLevel);
    s = "BBOT: " + NOMAD::BBOutputTypeListToString(NOMAD::Algorithm::getBbOutputType());
    AddOutputInfo(s, _displayLevel);
    OUTPUT_INFO_END

    generateTrialPoints();
}

bool NOMAD::QPSolverOptimize::runImp()
{
    std::string s;
    bool foundBetter = false;
    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
        
        if (_modelFixedVar.nbDefined() > 0)
        {
            NOMAD::EvalPointSet evalPointSet;
            for (auto trialPoint : _trialPoints)
            {
                evalPointSet.insert(trialPoint.makeFullSpacePointFromFixed(_modelFixedVar));
            }
            _trialPoints.clear();
            _trialPoints = evalPointSet;
        }
        // Update barrier
        postProcessing();
        
        // If the oracle point cannot be evaluated, the optimization has failed.
        if (_success==NOMAD::SuccessType::NO_TRIALS)
        {
            auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );
            qmsStopReason->set( NOMAD::ModelStopType::NO_NEW_POINTS_FOUND);
        }
        
    }
    
    
    return foundBetter;
}



void NOMAD::QPSolverOptimize::generateTrialPointsImp()
{

    bool runOk = false;

    // We need a base point.
    // QPSolverOptimize can be called by QPSolverSearchMethod OR as a standalone optimization
    NOMAD::Point X_k;
    if (nullptr == _iterAncestor )
    {
        throw NOMAD::Exception(__FILE__,__LINE__,getName() + " must have an Iteration ancestor.");
    }

    auto refCenter = getParentOfType<QuadModelIteration*>()->getRefCenter();
    X_k = *(refCenter->getX());
    
    if (_verboseFull)
    {
        SGTELIB::Matrix out = getModelOut(X_k);
        std::cout <<"Model output for RefCenter: " << std::endl;
        out.display(std::cout);
        auto bbo = refCenter->getEval(NOMAD::EvalType::BB)->getBBOutput();
        std::cout<< "Blackbox output of RefCenter: " << bbo.getBBO()<<std::endl;
    }
    bool feasRefCenter = false;
    
    //_model->display_trainingset();
    
    if (refCenter->isFeasible(NOMAD::EvalType::BB))
    {
        if ( nullptr != _prevFeasRefCenter &&
            _prevFeasRefCenter->NOMAD::ArrayOfDouble::isDefined() &&
            X_k == *(_prevFeasRefCenter->getX()) &&
            _prevFeasXopt.isComplete()  )
        {
            X_k = _prevFeasXopt;
        }
        _prevFeasRefCenter = refCenter;
        feasRefCenter = true;
    }
    else if (! refCenter->isFeasible(NOMAD::EvalType::BB))
    {
        if (nullptr != _prevInfeasRefCenter &&
            _prevInfeasRefCenter->NOMAD::ArrayOfDouble::isDefined() &&
            X_k == *(_prevInfeasRefCenter->getX()) &&
            _prevInfeasXopt.isComplete() )
        {
            X_k = _prevInfeasXopt;
        }
        _prevInfeasRefCenter = refCenter;
        feasRefCenter = false;
    }

    
    OUTPUT_INFO_START
    std::string s = "Model base X:" +X_k.display() ;
    
    AddOutputInfo(s, _displayLevel);
    OUTPUT_INFO_END
    if (_m == 1 && _bbot[0].isObjective())
    {
        runOk = solveBCQP(X_k);
    }
    else
    {
        
        double MeshSize = 0.0;

        auto quadModelIter = getParentOfType<NOMAD::QuadModelIteration*>(false);
        if (nullptr == quadModelIter)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"QPSolverOptimize must have a quadModelIteration as parent");
        }
        
        auto mesh = quadModelIter->getMesh();
        if ( nullptr != mesh)
        {
            // auto FrameSize = mesh->getDeltaFrameSize().max().todouble();
            MeshSize = mesh->getdeltaMeshSize().max().todouble();
            auto meshIndex = mesh->getMeshIndex();
            _verbose && std::cout << " meshIndex=" << meshIndex << std::endl;
            // std::cout << " frame size=" << FrameSie.max() << " mesh size=" << MeshSie.max() << std::endl;
        }


        SGTELIB::Matrix Gk("Gk", static_cast<int>(_n), 1);
        getModelGrad(&Gk, X_k);
        double ng = Gk.norm();

        SGTELIB::Matrix Y("Y", static_cast<int>(_nbCons), 1);
        Y.fill(1.0);
        SGTELIB::Matrix Hlag = getModelLagHessian(X_k, Y);
        SGTELIB::Matrix svdHLag = Hlag.get_singular_values();
        double sing_val_min = svdHLag.min();
        double condHessian;
        if (sing_val_min > 0)
        {
            condHessian = svdHLag.max() / sing_val_min;
        }
        else
        {
            condHessian = INF;
        }
        

        auto maxIter = static_cast<int>(_runParams->getAttributeValue<size_t>("QP_maxIter"));
        auto tolDistDX = _runParams->getAttributeValue<Double>("QP_tolDistDX").todouble();
        auto atol =_runParams->getAttributeValue<Double>("QP_absoluteTol").todouble();
        auto rtol = _runParams->getAttributeValue<Double>("QP_relativeTol").todouble();
        auto tolMesh = _runParams->getAttributeValue<Double>("QP_tolMesh").todouble();
        auto tolCond = _runParams->getAttributeValue<Double>("QP_tolCond").todouble();
        if (MeshSize > 0)
        {
            atol = std::min(atol, MeshSize * tolMesh);
            rtol = std::min(rtol, MeshSize * tolMesh);
        }
        if (sing_val_min > 0)
        {
            atol = std::max(atol, condHessian * tolCond);
            rtol = std::max(rtol, condHessian * tolCond);
        }
        double ng0 = ng;
        double tol = atol + ng0 * rtol;

        
        
        // Specific to augmented Lagrangian
        //

        auto mu0 = _runParams->getAttributeValue<Double>("QP_AugLag_mu0").todouble(); // TODO: * std::min(std::sqrt(tol), 1.0);
        auto muDecrease = _runParams->getAttributeValue<Double>("QP_AugLag_muDecrease").todouble();

        auto eta0 = _runParams->getAttributeValue<Double>("QP_AugLag_eta0").todouble();
        auto omega0 = _runParams->getAttributeValue<Double>("QP_AugLag_omega0").todouble(); // TODO: * std::min(std::sqrt(tol), 1.0);

        auto successRatio = _runParams->getAttributeValue<Double>("QP_AugLag_successRatio").todouble();
        auto maxIterInner = _runParams->getAttributeValue<size_t>("QP_AugLag_maxIterInner");
        auto tolDistDXInner = _runParams->getAttributeValue<Double>("QP_AugLag_tolDistDXInner").todouble();
        auto maxSuccessivFail = _runParams->getAttributeValue<size_t>("QP_AugLag_maxSuccessivFail");

        auto SelectAlgo = _runParams->getAttributeValue<size_t>("QP_SelectAlgo");
        if (SelectAlgo == 0)
        {
            _verbose && std::cout << "Run solveAugLag (n=" << _n << ", m=" << _nbCons << ")" << std::endl;
            _verbose && std::cout << "atol=" << atol << " rtol=" << rtol << " tol=" << tol << " cond(H)=" << condHessian << " mesh=" << MeshSize << std::endl;
            runOk = solveAugLag(X_k, maxIter, tolDistDX, atol, rtol, mu0, muDecrease, eta0, omega0, successRatio, maxIterInner, tolDistDXInner, maxSuccessivFail);
        }
        else if (SelectAlgo == 1)
        {
            _verbose && std::cout << "Run solveTRIPM (n=" << _n << ", m=" << _nbCons << ")" << std::endl;
            _verbose && std::cout << "atol=" << atol << " rtol=" << rtol << " tol=" << tol << " cond(H)=" << condHessian << std::endl;
            runOk = solveTRIPM(X_k, maxIter, tolDistDX, atol, rtol, mu0, muDecrease, maxIterInner, _verbose, _verboseFull);
        }
        else if (SelectAlgo == 2)
        {
            runOk = solveL1AugLag(X_k);
        }
        else if (SelectAlgo == 3)
        { // feasibility only
            _verbose && std::cout << "Run feasibility check (n=" << _n << ", m=" << _nbCons << ")" << std::endl;
            _verbose && std::cout << "atol=" << atol << " rtol=" << rtol << " tol=" << tol << " cond(H)=" << condHessian << std::endl;

            SGTELIB::Matrix cons("cons", _nbCons, 1);
            getModelCons(&cons, X_k);

            SGTELIB::Matrix lvar("lvar", _n + _nbCons, 1);
            SGTELIB::Matrix uvar("uvar", _n + _nbCons, 1);
            SGTELIB::Matrix XS("XS", _n + _nbCons, 1);
            SGTELIB::Matrix p("p", _n + _nbCons, 1);

            bool strict = getStrictFeasiblePoint(X_k, XS, lvar, uvar, cons);
            _verbose && std::cout << " strict feasibility found? " << strict << std::endl;
            if (strict)
            {
                solveLM(X_k, XS, lvar, uvar, cons, mu0, tol, maxIterInner, tolDistDXInner, false, _verbose);
            }
        }
        else
        {
            // Not implemented
        }
    }

    if (!runOk)
    {
        OUTPUT_INFO_START
        std::string s = "Solver run NOT OK";
        AddOutputInfo(s);
        OUTPUT_INFO_END
        
        auto modelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        modelStopReasons->set(NOMAD::ModelStopType::MODEL_OPTIMIZATION_FAIL);
    }
    else
    {
        if (X_k.isComplete())
        {
            // Keep the solution returned by solver
            if (feasRefCenter)
            {
                _prevFeasXopt = X_k;
            }
            else
            {
                _prevInfeasXopt = X_k;
            }
            
            // New EvalPoint to be evaluated.
            // Add it to the list (local or in Search method).
            bool inserted = insertTrialPoint(NOMAD::EvalPoint(X_k));

            OUTPUT_INFO_START
            std::string s = "xt:";
            s += (inserted) ? " " : " not inserted: ";
            s += X_k.display() + " \n";
            if (inserted)
            {
                SGTELIB::Matrix out = getModelOut(X_k);
                std::ostringstream os;
                out.display(os);
                s += "Output(xt) = " + os.str();
            }
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
    }

}


// Set the bounds and the extra fixed variables (when lb=ub) of the model.
void NOMAD::QPSolverOptimize::setModelBoundsAndFixedVar()
{
    // When optWithScaleBounds is true, the training set is scaled with some generating directions. Warning: points are not necessarilly in [0,1]^n
    const SGTELIB::Matrix & X = _trainingSet->get_matrix_X();
    
    _n = static_cast<int>(_pbParams->getAttributeValue<size_t>("DIMENSION"));
    
    if (_n != X.get_nb_cols())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "QPSolverOptimize::setModelBounds() dimensions do not match");
    }
    
    int nbDim = X.get_nb_cols();
    int nbPoints = X.get_nb_rows();
    
    // Build model bounds & detect fixed variables
    NOMAD::Double lb;
    NOMAD::Double ub;
    bool isFixed = false;
    for (int j = 0; j < nbDim; j++)
    {
        lb = _modelLowerBound[j];
        ub = _modelUpperBound[j];
        for (int p = 0; p < nbPoints; p++)
        {
            // Use a comparison on regular double for better precision
            
            auto xpj = NOMAD::Double(X.get(p,j));
            if (lb.isDefined())
            {
                if (xpj.todouble() < lb.todouble())
                    lb = xpj;
            }
            else
            {
                lb = xpj;
            }
            if (ub.isDefined())
            {
                if (xpj.todouble() > ub.todouble())
                    ub = xpj;
            }
            else
            {
                ub = xpj;
            }
        }
        isFixed = false;
        // Comparison of Double at epsilon
        if (lb == ub)
        {
            _modelFixedVar[j] = ub;
            lb = NOMAD::Double(); // undefined
            ub = NOMAD::Double();
            isFixed = true;
        }
        
        
        if (!_optWithScaledBounds)
        {
            _modelLowerBound[j] = lb;
            _modelUpperBound[j] = ub;
        }
        else
        {
            // When scaling we force the bounds to be [0,1] whatever the training set except if we have a fixed variable.
            // If optWithScaledBounds is true and we detect a fixed variable, the modelFixedVar is not necessarily within [0,1] but the model center and the fixed variable must be equal.
            if (isFixed)
            {
                _modelLowerBound[j] = _modelUpperBound[j] = NOMAD::Double();
                _modelCenter[j] = _modelFixedVar[j];
            }
            // If not fixed, the model bounds and model center are pretermined.
            else
            {
                _modelLowerBound[j] = 0;
                _modelUpperBound[j] = 1;
                _modelCenter[j] = 0.5;
            }
        }
    }
    
    if (!_optWithScaledBounds)
    {
        // Detect the model center of the bounds
        // Scale the bounds around the model center
        auto reduction_factor = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR");
        for (int j = 0; j < nbDim; j++)
        {
            lb = _modelLowerBound[j];
            ub = _modelUpperBound[j];
            if (lb.isDefined() && ub.isDefined())
            {
                // The model center is the bounds middle point
                _modelCenter[j] = (lb + ub)/2.0;
                
                // Scale the bounds with respect to the bounds
                lb = _modelCenter[j] + (lb-_modelCenter[j])/reduction_factor;
                ub = _modelCenter[j] + (ub-_modelCenter[j])/reduction_factor;
                
                // Comparison of Double at epsilon
                if (lb == ub)
                {
                    _modelFixedVar[j] = ub;
                    _modelCenter[j] = ub;
                    lb = NOMAD::Double(); // undefined
                    ub = NOMAD::Double();
                    
                }
            }
            else
            {
                _modelCenter[j] = _modelFixedVar[j];
            }
            
            _modelLowerBound[j] = lb;
            _modelUpperBound[j] = ub;
        }
    }
    
    OUTPUT_INFO_START
    std::string s = "model lower bound: " + _modelLowerBound.display();
    AddOutputInfo(s);
    s = "model upper bound: " + _modelUpperBound.display();
    AddOutputInfo(s);
    s = "model center: " + _modelCenter.display();
    AddOutputInfo(s);
    OUTPUT_INFO_END
    
} // end setModelBounds

bool NOMAD::QPSolverOptimize::update(NOMAD::Point & X, const SGTELIB::Matrix & Y, const double a)
{
    for (int i = 0 ; i < Y.get_nb_rows() ; i++)
    {
        X[i] += a * Y.get(i, 0);
    }

    return true;
}

bool NOMAD::QPSolverOptimize::update(NOMAD::Point & Xup, const NOMAD::Point & X, const SGTELIB::Matrix & Y, const double a)
{
    for (int i = 0 ; i < Y.get_nb_rows() ; i++)
    {
        Xup[i] = X[i] + a * Y.get(i, 0);
    }

    return true;
}

void NOMAD::QPSolverOptimize::solve_TR_constrained_QP(
    SGTELIB::Matrix * d,
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g,
    SGTELIB::Matrix & grad,
    const double Delta )
{
    const int n = X.get_nb_rows(); // length of X and active

    lencheck(n, g);
    sizecheck(n, n, H);

    if (Delta < 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Invalid Delta value <0.");
    }

    bool * active = new bool [n];
    for (int i=0; i < n; ++i)
    {
        active[i] = 0;
    }

    solve_TR_constrained_QP(d, X, H, g, grad, active, Delta);

    delete [] active;
}

void NOMAD::QPSolverOptimize::solve_TR_constrained_QP(
    SGTELIB::Matrix * d,
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g,
    SGTELIB::Matrix & grad,
    const bool * active,
    const double Delta,
    bool verbose )
{

    int n = X.get_nb_rows(); // length of X and active
    int nfree = n - sum(active, n);

    lencheck(n, g);
    sizecheck(n, n, H);

    // Pre-allocations
    double nd, a_unb, slope;
    d->fill(0);

    _verbose && std::cout << "Starting solve_TR_constrained_QP with delta=" << Delta << " nfree=" << nfree << std::endl;

    // q(x) = 0.5 * x' * H * x + g' * x
    getModelGrad(&grad, X, H, g); // Hx + g

    SGTELIB::Matrix HW = matrix_subset(H, active); // HW("HW", nfree, nfree); // fixed
    HW.set_name("HW");
    SGTELIB::Matrix gW = vector_subset(grad, active); // SGTELIB::Matrix gW("gW", nfree, 1); // depend on X
    gW.set_name("gW");

    ////////// LDLt
    // init matrices for LDLt
    double ** M = new double *[nfree];
    double ** L = new double *[nfree];
    double ** D = new double *[nfree];
    for (int i = 0 ; i < nfree ; ++i )
    {
        M[i] = new double[nfree];
        L[i] = new double[nfree];
        D[i] = new double[nfree];
        for (int j = 0 ; j < nfree ; ++j )
        {
            M[i][j] = HW.get(i, j);
            L[i][j] = 0;
            D[i][j] = 0;
        }
    }
    int * pp = new int [nfree];
    for (int i = 0 ; i < nfree ; ++i )
    {
        pp[i] = 0;
    }
    std::string error_msg;

    bool success = NOMAD::LDLt_decomposition(error_msg, M, L, D, pp, nfree);
    if (!success)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt decomposition");
    }

    double eigmin = NOMAD::FindSmallestEigenvalue(D, nfree);
    _verbose && std::cout << " smallest eigenvalue= " << eigmin << std::endl;

    if (eigmin > 0) // positive definite
    {
        bool newton_success = Convex_TR_QP(d, grad, gW, H, HW, pp, D, L, active, Delta, _verbose);

        if (!newton_success)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Convex trust-region solve failure");
        }
    }
    else
    {
        _verbose && std::cout << "Not positive definite. Delta= " << Delta << " l=" << eigmin << std::endl;
        
        SGTELIB::Matrix bk("bk", nfree, 1); // depend on X
        bool successInverseIteration = true;

        successInverseIteration = InverseIteration(&bk, HW, eigmin, nfree, 1E-12);

        if (!successInverseIteration)
        { // TODO: Need a back-up plan
            std::cout << "Error InverseIteration";
            d->fill(0.0);
            // throw NOMAD::Exception(__FILE__, __LINE__, "Error InverseIteration");
        }
        else 
        {
            nd = bk.norm();
            slope = SGTELIB::Matrix::dot(gW, bk);
            if ((Delta < 1E15))
            {
                a_unb = Delta / nd;
            }
            else
            {
                a_unb = 1000 * abs(slope) / nd;
            }

            int ki = 0;
            d->fill(0);
            for (int i = 0 ; i < n ; ++i )
            {
                    if (!active[i])
                    {
                        d->set(i, 0, a_unb * bk.get(ki, 0));
                        ki ++;
                    }
            }

            if (ki != nfree)
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "Error dimension");
            }
        }
    }

    for (int i=0 ; i<nfree ; i++)
    {
        delete [] M[i];
    }
    delete [] M;
    for (int j=0 ; j<nfree ; j++)
    {
        delete [] L[j];
    }
    delete [] L;
    for (int j=0 ; j<nfree ; j++)
    {
        delete [] D[j];
    }
    delete [] D;
    delete [] pp;
}

bool NOMAD::QPSolverOptimize::solveBCQP(
    NOMAD::Point & X,
    const int max_iter,
    const double atol,
    const double rtol,
    bool verbose )
{

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    
    SGTELIB::Matrix Y("x0", _n, 1);
    Y.fill(0);
    double g0 = surrogate_prs->getModelObj(Y);
    SGTELIB::Matrix g("g", _n, 1); // rhs of the quadratic
    g = surrogate_prs->getModelGrad(Y);
    SGTELIB::Matrix H("H", _n, _n);
    H = surrogate_prs->getModelHessian(Y);

    SGTELIB::Matrix Xm("X", _n, 1);
    SGTELIB::Matrix lvar("lvar", _n, 1);
    SGTELIB::Matrix uvar("uvar", _n, 1);
    for (int i=0; i < _n; ++i)
    {
        double lb = ((_modelLowerBound[i].isDefined())? _modelLowerBound[i].todouble(): NOMAD::M_INF);
        double ub = ((_modelUpperBound[i].isDefined())? _modelUpperBound[i].todouble(): NOMAD::INF);
        Xm.set(i, 0, X[i].todouble());
        lvar.set(i, 0, lb);
        uvar.set(i, 0, ub);
    }

    bool success = solveBCQP(Xm, H, g, g0, lvar, uvar, max_iter, atol, rtol, verbose);

    for (int i=0; i < _n; ++i)
    {
        X[i] = Xm.get(i, 0);
    }

    return success;
}

void NOMAD::QPSolverOptimize::projectedGradient(
    SGTELIB::Matrix & X,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g,
    const double g0,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    bool * active_l,
    bool * active_u,
    bool * working,
    SGTELIB::Matrix & d_k,
    SGTELIB::Matrix & Temp,
    const double tol,
    const int max_iter,
    bool verbose )
{
    const int n = X.get_nb_rows();

    // Pre-allocation
    SGTELIB::Matrix armijo_Xp("armijo_Xp", n, 1);

    bool OKbis = false;
    // bool ProjectedStop = false;
    double slope;
    int nfree = n - sum(active_l, n) - sum(active_u, n);
    int nfreep;
    double qm = getModelObj(X, H, g, g0);
    getModelGrad(&d_k, X, H, g);
    double qmp;

    double a_max=1.0, a_k = 1.0;

    int k = 0;

    while (!OKbis && (k < max_iter))
    {
        d_k.multiply(-1.0);
        project_bounds(d_k, working);
        slope = -d_k.normsquare();

        a_k = 1; // min(a_k, 1) is the inital value for Armijo
        a_k = projected_armijo(X, H, g, g0, lvar, uvar, d_k, qm, slope, armijo_Xp, Temp, a_k); // a_k start from previous value
        d_k.multiply(a_k);
        X.add(d_k);
        
        // X is not necessarily feasible here.
        snapToBounds(X, lvar, uvar);

        qmp = getModelObj(X, H, g, g0);
        getModelGrad(&d_k, X, H, g);

        active_bounds(X, lvar, uvar, active_l, active_u);
        nfreep = n - sum(active_l, n) - sum(active_u, n);

        verbose && std::cout << "  Projected-gradient k=" << k << " f(x)=" << qmp << " |A|=" << nfreep << " " << nfree;
        verbose && std::cout << " |d|=" << d_k.norm() << " amax=" << a_max << " ak=" << a_k << std::endl;
        k ++;

        OKbis = (qmp >= qm) || (nfree == nfreep);
        nfree = nfreep;
        qm = qmp;
    }
    // return X updated
}

bool NOMAD::QPSolverOptimize::solveBCQP(
    SGTELIB::Matrix & X,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g,
    const double g0,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    const int max_iter,
    const double atol,
    const double rtol,
    bool verbose )
{
    const int n = X.get_nb_rows();

    lencheck(n, X);
    lencheck(n, g);
    lencheck(n, lvar);
    lencheck(n, uvar);
    sizecheck(n, n, H);

    // Check that X is feasible
    bool feasible = true;
    for (int i=0; i < n; ++i)
    {
        feasible = feasible && (X.get(i, 0) >= lvar.get(i, 0)) && (X.get(i, 0) <= uvar.get(i, 0));
    }
    if (!feasible)
    {
        verbose && std::cout << "solveBCQP assertion error: X not feasible. Compute projection." << std::endl;
        snapToBounds(X, lvar, uvar);
    }
    bool bound_compat = true;
    feasible = true;
    for (int i=0; i < n; ++i)
    {
        feasible = feasible && (X.get(i, 0) >= lvar.get(i, 0)) && (X.get(i, 0) <= uvar.get(i, 0));
        bound_compat = bound_compat && (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!feasible || !bound_compat)
        {
            verbose && std::cout << lvar.get(i, 0) << " " << X.get(i, 0) << " " << uvar.get(i, 0);
            throw NOMAD::Exception(__FILE__, __LINE__, "solveBCQP assertion error: Error compatibility lower and upper bound");
        }
    }

    // Pre-allocation
    SGTELIB::Matrix d_k("dk", n, 1);
    SGTELIB::Matrix Xp("Xp", n, 1);
    SGTELIB::Matrix Grad("Grad", n, 1);
    SGTELIB::Matrix armijo_Xp("armijo_Xp", n, 1);
    SGTELIB::Matrix Temp("Temp-vector", n, 1);

    bool * active_l = new bool [n];
    bool * active_u = new bool [n];
    bool * binding = new bool [n];
    bool * working = new bool [n]; // subset of active_l OR active_u

    // Initialization
    double nd = d_k.norm(); // \|dk\|
    double f0 = getModelObj(X, H, g, g0);
    double fk, fktemp, fkptemp ;
    getModelGrad(&Grad, X, H, g);
    double ng0 = Grad.norm();

    active_bounds(X, lvar, uvar, active_l, active_u);
    binding_bounds(Grad, active_l, active_u, binding);

    for (int i=0; i < n; ++i)
    {
        working[i] = false;
    }

    double a_max, a_k, qm;

    // Start iterating
    int k = 0;
    verbose && std::cout << "  k=" << k << " f(x0)=" << f0;
    verbose && std::cout << " |W|" << sum(working, n) << " |L|" << sum(active_l, n) << " |U|" << sum(active_u, n) << std::endl;
    bool OK = false;
    while (!OK && (k < max_iter))
    {
        // First, we compute a sequence of projected gradient steps
        // Update active_l and active_u
        projectedGradient(X, H, g, g0, lvar, uvar, active_l, active_u, working, Grad, Temp, atol + ng0 * rtol, n);
        qm = getModelObj(X, H, g, g0);
        getModelGrad(&Grad, X, H, g); // TODO: projectedGradient return -Grad.
        
        // Optimize within the working set
        solve_unconstrained_QP(&d_k, X, H, g, Grad, working);

        if (d_k.has_nan()) {
            throw NOMAD::Exception(__FILE__, __LINE__, "d_k contains NaN");
        }

        // check that it worked
        nd = d_k.norm();
        Xp = X; Xp.add(d_k); // Xp = X + d_k
        bool unconstrainedSucess = (getModelObj(X, H, g) > getModelObj(Xp, H, g)) || (nd <= atol);
        if (!unconstrainedSucess)
        {
            d_k.fill(0);
            // verbose && std::cout << getModelObj(X, H, g) << " " << getModelObj(SGTELIB::Matrix::add(X, d_k), H, g) << " " << getModelObj(SGTELIB::Matrix::add(X, SGTELIB::Matrix::product(d_k, 1000.0)), H, g);
            // verbose && std::cout  << " |d|=" << nd << std::endl;
            // throw NOMAD::Exception(__FILE__, __LINE__, "solve_unconstrained_QP: unexpected failure.");
        }

        a_max = max_step_bounds(X, lvar, uvar, d_k);

        double slope = SGTELIB::Matrix::dot(d_k, Grad);
        if (slope >= 0)
        {
            // verbose && std::cout << "  Use negative curvature direction" << std::endl;
            Xp = d_k; Xp.multiply(a_max); Xp.add(X); // Xp = X + d_k a_max
            fkptemp = getModelObj(Xp, H, g);
            fktemp = getModelObj(X, H, g);
            if (fkptemp <= fktemp)
            {
                a_k = a_max;
            }
            else
            {
                a_k = 0;
                // a_k = projected_armijo(X, H, g, g0, lvar, uvar, d_k, qm, slope); // check a_max automatically
                // a_k = std::min(a_k, a_max);
            }
        }
        else
        { // If no negative curvature: we do a line search!
            if (a_max > 1E-15) // 1E-15 is the smallest value in projected_armijo
            {
                a_k = projected_armijo(X, H, g, g0, lvar, uvar, d_k, qm, slope, armijo_Xp, Temp, a_max); // check a_max automatically
            }
            else
            {
                a_k = 0;
            }
            // a_k = std::min(a_k, a_max);
        }

        // Update iterate
        d_k.multiply(a_k);
        X.add(d_k);

        getModelGrad(&Grad, X, H, g);
        fk = getModelObj(X, H, g, g0);

        if (fk > qm)
        {
            //throw NOMAD::Exception(__FILE__, __LINE__, "non-minimizing");
            verbose && std::cout << " This is weird: " << a_k << " " << a_max << " " << fk << " < " << qm << " = " << f0 << std::endl;
            // throw NOMAD::Exception(__FILE__, __LINE__, "non-minimizing");
        }

        bool update = false; // switch to `true` if working has been updated
        if (a_max <= 1) // need to update active sets and the working set with a new constraint
        {
            for (int i = 0; i < n; ++i)
            { // We will correct if needed
                if (std::fabs(X.get(i, 0) - lvar.get(i, 0)) <= 1E-15)
                {
                    X.set(i, 0, lvar[i]);
                }
                if (std::fabs(X.get(i, 0) - uvar.get(i, 0)) <= 1E-15)
                {
                    X.set(i, 0, uvar[i]);
                }
                active_l[i] = (X.get(i, 0) == lvar.get(i, 0));
                active_u[i] = (X.get(i, 0) == uvar.get(i, 0));
                if (!update && !working[i] && (active_l[i] || active_u[i]))
                {
                    working[i] = true;
                    update = true;
                }
            }
        }

        snapToBounds(X, lvar, uvar);
        feasible = true;
        for (int i=0; i < n; ++i)
        {
            feasible = feasible && (X.get(i, 0) >= lvar.get(i, 0)) && (X.get(i, 0) <= uvar.get(i, 0));
            if (!feasible)
            {
                std::cout << lvar.get(i, 0) << " " << X.get(i, 0) << " " << uvar.get(i, 0) << " a=" << a_k << " amax=" << a_max << std::endl;
                std::cout << active_l[i] << active_u[i] << " " << nd << std::endl;
                std::cout << "Error compatibility lower and upper bound" << std::endl;
                // throw NOMAD::Exception(__FILE__, __LINE__, "solveBCQP assertion error (4): Error compatibility lower and upper bound");
            }
        }

        binding_bounds(Grad, active_l, active_u, binding);

        if ((nd < atol + ng0 * rtol) && !update) // we should not remove a constraint we just added
        {
            // If OK = false, then working has been updated
            OK = check_subset_binding_update(working, binding, _n);
        }
        else if ((!unconstrainedSucess || a_k == 0) && !update)
        {
            OK = true; // no progress
        }
        k++;
        verbose && std::cout << "  k=" << k << " (" << unconstrainedSucess << ") " << " f(xk)=" << fk << " |d|=" << nd << " a(amax)= " << a_k << " ( " << a_max << " )";
        verbose && std::cout << " |W|" << sum(working, n) << " |L|" << sum(active_l, n) << " |U|" << sum(active_u, n) << " OK? " << OK << std::endl;
    }

    bool success = getModelObj(X, H, g, g0) < f0;

    delete [] active_l;
    delete [] active_u;
    delete [] binding;
    delete [] working;
    
    return success;
}

bool NOMAD::QPSolverOptimize::check_subset_binding_update(
    bool * working,
    const bool * binding,
    size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        if (working[i] && !binding[i])
        {
            working[i] = false;
            return false;
        }
    }  
    return true;
}

double NOMAD::QPSolverOptimize::projected_armijo(
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g,
    const double g0,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    const SGTELIB::Matrix & d,
    const double fk,
    const double slope,
    SGTELIB::Matrix & Xp,
    SGTELIB::Matrix & gradientF_kp,
    const double t_max )
{
    int n = X.get_nb_rows();

    // Linesearch parameters
    double armijo_tol = 1E-4; // > 0
    double t_small = 1E-15; // small
    double t_decrease = 2.5; // > 1
    double wolfe_tol = 0.9999; // < 1
    int bW_max = 5; // nonnegative integer
    double t_increase = 5; // > 1

    bool good_grad = false; // true if gradient has been updated. TODO: use it to save a gradient call.

    // Pre-allocation
    lencheck(n, Xp);
    lencheck(n, gradientF_kp);

    // Initialize
    double tk = std::min(1.0, t_max);
    
    Xp = d; Xp.multiply(tk); Xp.add(X); // Xp = X + t_k d
    snapToBounds(Xp, lvar, uvar);
    double fkp = getModelObj(Xp, H, g, g0);
    getModelGrad(&gradientF_kp, Xp, H, g);
    double slope_t = SGTELIB::Matrix::dot(d, gradientF_kp);

    int nbW = 0;
    while ((slope_t < wolfe_tol * slope) // Wolfe condition
           && (fkp <= fk - armijo_tol * tk * fabs(slope)) // Armijo condition
           && (nbW < bW_max)
           && (tk <= t_max))
    {
        tk *= t_increase;
        Xp = d; Xp.multiply(tk); Xp.add(X); // Xp = X + t_k d
        snapToBounds(Xp, lvar, uvar);
        getModelGrad(&gradientF_kp, Xp, H, g);
        fkp = getModelObj(Xp, H, g, g0);
        slope_t = SGTELIB::Matrix::dot(d, gradientF_kp);
        nbW ++;
        good_grad = true;
    }
    nbW=0;

    bool armijo = fkp <= fk - armijo_tol * tk * fabs(slope);
    // Enrich Armijo's condition with Hager & Zhang numerical trick
    // armijo = armijo || ((fkp <= fk + 4E-10 * abs(fk)) && (slope_t <= fact * slope));
    while (!armijo && tk > t_small)
    {
        tk /= t_decrease;
        Xp = d; Xp.multiply(tk); Xp.add(X); // Xp = X + t_k d
        snapToBounds(Xp, lvar, uvar);
        fkp = getModelObj(Xp, H, g, g0);
        nbW ++;

        armijo = fkp <= fk - armijo_tol * tk * fabs(slope);
        //good_grad = false;
        //if (!armijo && (fkp <= fk + 4E-10 * abs(fk)))
        //{
        //    getModelGrad(&gradientF_kp, Xp, H, g);
        //    slope_t = SGTELIB::Matrix::dot(d, gradientF_kp);
        //    armijo = slope_t <= fact * slope;
        //    good_grad = true;
        //}
    }

    if (!armijo)
    {
        return 0.0;
    }

    return tk;
}

// Solve method with outter loop and inner loop
bool NOMAD::QPSolverOptimize::solveL1AugLag(
    NOMAD::Point & X_k,
    const int max_iter,
    const double atol,
    const double rtol)
{


    double fk = getModelObj(X_k);
    // double fk_old;
    
    SGTELIB::Matrix gradientLag_k("gradientLag_k", _n, 1);
    getModelGrad(&gradientLag_k, X_k);
    double ng = gradientLag_k.norm();
    double ng0 = ng;
    double tol = atol + ng0 * rtol;
    std::cout << "Start solveL1AugLag with tol=" << tol << std::endl;

    SGTELIB::Matrix hessianLag_k("hessianLag_k", _n, _n);
    SGTELIB::Matrix invHessianLag("invHessianLag_k", _n, _n);
    
    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X_k);

    bool* active = new bool [_nbCons]; // true for indices of active constraints
    bool* feasible = new bool [_nbCons]; // true for indices of strictly feasible constraints
    bool* infeasible = new bool [_nbCons]; // true for indices of ineasible constraints

    getModelActiveCons(cons, tol, active);
    int nbActive = sum(active, _nbCons);
    getModelFeasibleCons(cons, tol, feasible);
    getModelInfeasibleCons(cons, tol, infeasible);

    // Lagrange multipliers
    SGTELIB::Matrix multiplier_k("multiplier_k", _nbCons, 1);
    multiplier_k.fill(1.0);
    SGTELIB::Matrix active_multiplier_k("active_multiplier_k", nbActive, 1); // TODO: bad...
    multiplier_k.fill(0.0);

    SGTELIB::Matrix Jacobian_k("Jacobian_k", _nbCons, _n);
    // std::cout << _nbCons << " " << _n << " " << _m << std::endl;
    Jacobian_k = getModelJacobian(X_k);
    Jacobian_k.display(std::cout);

    SGTELIB::Matrix activeJacobian_k = getModelActiveJacobian(Jacobian_k, active);

    // Compute pseudo-gradient:
    SGTELIB::Matrix pseudoGradient_k("pseudoGradient_k", _n, 1);
    for (int w=0 ; w < _nbCons ; w++ )
    {
        multiplier_k.set(w, 0, infeasible[w]); // 1 for constraints active
    }
    pseudoGradient_k = getModelLagGradient(X_k, multiplier_k);

    active_multiplier_k = SGTELIB::Surrogate_PRS::compute_multiplier(pseudoGradient_k, activeJacobian_k);
    double dual_norm = compute_dual_residual(pseudoGradient_k, activeJacobian_k, active_multiplier_k);

    bool innerSuccess;
    // bool innerFailure;
    innerSuccess = (dual_norm < tol) && isFeasible(cons, tol);

    // Outer loop parameters
    SGTELIB::Matrix lambda_l("lambda_l", _nbCons, 1);
    lambda_l.fill(0.0);
    double mu_l = 1.0;
    double eta_l = 1.0;
    double omega_l = 1.0;

    bool unbounded_subpb = false; // unbounded subproblem

    // Inner loop pre-allocation
    NOMAD::Point X_km1(X_k), X_om1(X_k), X_can(X_k);
    SGTELIB::Matrix h_k("h_k", _n, 1); // horizontal step
    SGTELIB::Matrix v_k("v_k", _n, 1); // vertical step


    size_t iterOutterLoop = 0;
    double distXOutterLoop = INF;
    double F_k =INF, F_km1 = INF;
    const size_t maxIterInnerLoop = 20, maxIterOutterLoop = 10;
    
    const double toleranceF = 1E-12, tolerance_distDX = 1E-12 ;
    size_t quadModelNbEval = 0;

    F_k = getPenalizedL1AugLagModelObj(X_k, cons, lambda_l, mu_l);
    gradientLag_k = getModelLagGradient(X_k, lambda_l);
    double ngproj = check_optimality_bounds(X_k, gradientLag_k);

    double outerSuccess;
    outerSuccess = (ngproj <= tol) && isFeasible(cons, tol);

    std::cout << " |grad|= " << ng << " |Proj(x - grad) - x|= " << ngproj << " P=" << F_k << " f=" << fk << " |c|=" << cons.norm() << std::endl;

    bool outerFailure = (distXOutterLoop <= tolerance_distDX) || (iterOutterLoop >= maxIterOutterLoop) || (quadModelNbEval >= _quadModelMaxEval);

    // outer iteration
    while (!outerFailure && !outerSuccess)
    {
        
        X_om1 = X_k;
        
        // Initialization inner loop
        size_t iterInnerLoop = 0;
        double innerPrecision = 1.0;
        getModelCons(&cons, X_k);
        std::cout << _nbCons << " " << _n << " " << _m << std::endl;
        Jacobian_k = getModelJacobian(X_k);
        getModelActiveCons(cons, innerPrecision, active);
        getModelFeasibleCons(cons, innerPrecision, feasible);
        getModelInfeasibleCons(cons, innerPrecision, infeasible);
        double deltaF = INF , distXkXkm1 = INF;
        bool innerSuccess = false;
        bool innerFailure = false;
        bool horizontalSuccess;
        bool verticalSuccess;

        while (!innerFailure && !innerSuccess)
        {
            F_km1 = F_k;
            X_km1 = X_k;
            X_can = X_k;

            horizontalSuccess = compute_horizontal_step(X_k, h_k, Jacobian_k, active, feasible, infeasible, lambda_l, mu_l);
            update(X_can, h_k);
            // Check optimality
            getModelCons(&cons, X_can);
            F_k = getPenalizedL1AugLagModelObj(X_can, cons, lambda_l, mu_l);
            std::cout << _nbCons << " " << _n << " " << _m << std::endl;
            Jacobian_k = getModelJacobian(X_can);
            // multiplier_k = SGTELIB::Surrogate_PRS::compute_multiplier(pseudoGradient_k, activeJacobian_k); // TODO: recompute multiplier
            dual_norm = check_inner_success(X_can, Jacobian_k, multiplier_k, lambda_l, mu_l, active, infeasible);
            std::cout << " H_can: (l=" << iterInnerLoop << ") dualnorm= " << dual_norm << " |h|= " << h_k.norm() << " Pk=" << F_k << " |c|= " << cons.norm() << std::endl;
            bool innerFirstOrder = (dual_norm <= tol) || (h_k.norm() <= tol);
            if (innerFirstOrder)
            { // Increase feasibility
//                std::cout << "Compute vertical step" << std::endl;
                verticalSuccess = compute_vertical_step(X_k, v_k, activeJacobian_k, cons, active);
                update(X_can, v_k); // X_k + h_k + v_k
                getModelCons(&cons, X_can);
                F_k = getPenalizedL1AugLagModelObj(X_can, cons, lambda_l, mu_l);
                std::cout << " V: (l=" << iterInnerLoop << ") |v|= " << v_k.norm() << " Pk=" << F_k << " ||c||= " << cons.norm() << std::endl;
                bool decreasePenalty = F_k < F_km1; // =...
                if (decreasePenalty)
                {
                    X_k = X_can;
                }
                else
                {
                    innerPrecision = innerPrecision / 2;
                    getModelCons(&cons, X_k);
                    getModelActiveCons(cons, innerPrecision, active);
                    getModelFeasibleCons(cons, innerPrecision, feasible);
                    getModelInfeasibleCons(cons, innerPrecision, infeasible);
                    horizontalSuccess = compute_horizontal_step(X_k, h_k, Jacobian_k, active, feasible, infeasible, lambda_l, mu_l);
                    double beta_k = piecewise_line_search(X_k, h_k, active , feasible, infeasible, lambda_l, mu_l);
                    update(X_k, h_k, beta_k);
                }

            }
            else
            {
                double beta_k = piecewise_line_search(X_k, h_k, active , feasible, infeasible, lambda_l, mu_l);
                update(X_k, h_k, beta_k);
            }
            // TODO
            
            // X_k.snapToBounds(_modelLowerBound, _modelUpperBound);

            // Check for success
            getModelCons(&cons, X_k);
            F_k = getPenalizedL1AugLagModelObj(X_k, cons, lambda_l, mu_l);
            std::cout << _nbCons << " " << _n << " " << _m << std::endl;
            Jacobian_k = getModelJacobian(X_k);
            dual_norm = check_inner_success(X_k, Jacobian_k, multiplier_k, lambda_l, mu_l, active, infeasible);
            innerSuccess = (dual_norm < tol) && isFeasible(cons, tol);
                   
            // Check for failure:
            deltaF = fabs(F_k-F_km1);
            distXkXkm1 = NOMAD::Point::dist(X_k,X_km1).todouble();
            // normGradLag_k = gradLag_k.norm();
            std::cout << " Inner: (l=" << iterInnerLoop << ") Pl= " << F_k << " |c|= " << cons.norm() << " |L|=" << dual_norm << " |dF|=" << deltaF << std::endl;

            iterInnerLoop++;
            quadModelNbEval++;
            // innerFailure = verticalSuccess && horizontalSuccess;
            innerFailure = (deltaF <= toleranceF) || (distXkXkm1 <= tolerance_distDX) || (iterInnerLoop >= maxIterInnerLoop) || (quadModelNbEval >= _quadModelMaxEval);
        
        } // end of inner loop - X_k has been updated


        getModelCons(&cons, X_k);

        F_k = getPenalizedL1AugLagModelObj(X_k, cons, lambda_l, mu_l);
        unbounded_subpb = (F_k < - 1/tol); 

        if (isFeasible(cons, eta_l))
        { // no update of mu_l
            for (int i = 0; i < _nbCons ; i++)
            {
                double li = lambda_l.get(i, 0);
                if (fabs(cons.get(i, 0)) <= tol)
                {
                    active[i] = true;
                    feasible[i] = false;
                    infeasible[i] = false;
                    lambda_l.set(i, 0, li + min(li, 0).todouble());
                }
                else if (cons.get(i, 0) < -tol)
                {
                    active[i] = false;
                    feasible[i] = true;
                    infeasible[i] = false;
                    lambda_l.set(i, 0, 0); // no update of lambda_l -- should be 0 I think. TODO: Make sure?
                }
                else
                {
                    active[i] = false;
                    feasible[i] = false;
                    infeasible[i] = true;
                    lambda_l.set(i, 0, li - 1 / mu_l);                    
                }
            }
            eta_l = eta_l * pow( mu_l, 0.9);
            omega_l = std::max(omega_l * mu_l, 1E-15);
        }
        else
        { // no update of lambda_l
            mu_l = mu_l / 10;
            eta_l = pow( mu_l, 0.1) / 10;
            omega_l = mu_l;
            if (unbounded_subpb)
            {
                X_k = X_km1;
            }
        }

        X_k.snapToBounds(_modelLowerBound, _modelUpperBound);

        // Check optimality
        gradientLag_k = getModelLagGradient(X_k, lambda_l);
        ngproj = check_optimality_bounds(X_k, gradientLag_k); // ||Proj(x - ∇L) - x||
        outerSuccess = (ngproj < tol) && isFeasible(cons, tol);
        
        iterOutterLoop++; // l = l + 1;
        fk = getModelObj(X_k);
        std::cout << "k= " << iterOutterLoop <<  " |Proj(x - grad) - x|= " << ngproj  << " f(x)= " << fk << " |c(x)|=" << cons.norm() << std::endl;

        // Check failure
        distXOutterLoop = NOMAD::Point::dist(X_k, X_om1).todouble();
        outerFailure = !unbounded_subpb && (distXOutterLoop <= tolerance_distDX);
        // outerFailure = outerFailure || innerFailure;
        outerFailure = outerFailure || (iterOutterLoop >= maxIterOutterLoop); //|| (quadModelNbEval >= _quadModelMaxEval);
        if (outerFailure)
        {
            std::cout << "Early stop: |d|=" << distXOutterLoop << " inner? " << innerFailure << " unbounded? " << unbounded_subpb << std::endl;
        }
        
    }
    
    delete [] active;
    delete [] feasible;
    delete [] infeasible;

    return outerSuccess;
}

// f(X) - lambda_i cx_i + max(cx_i, 0)/mu
double NOMAD::QPSolverOptimize::getPenalizedL1AugLagModelObj(
    const Point & X,
    const SGTELIB::Matrix & cons,
    const SGTELIB::Matrix & lambda,
    double mu) const
{

    lencheck(_nbCons, cons);
    lencheck(_nbCons, lambda);

    double lag;
    lag = getModelLag(X, lambda);

    for (int i = 0; i < _nbCons ; i++)
    {
        if (cons.get(i, 0) > 0)
        {
            lag += cons.get(i, 0) / mu;
        }
    }

    return lag;
}

double NOMAD::QPSolverOptimize::getPenalizedL1AugLagModelObj(
    const Point & X,
    const SGTELIB::Matrix & lambda,
    double mu) const
{
    SGTELIB::Matrix cons("cons", static_cast<int>(_nbCons), 1);
    getModelCons(&cons, X);
    return getPenalizedL1AugLagModelObj(X, cons, lambda, mu);
}

double NOMAD::QPSolverOptimize::check_inner_success(
    NOMAD::Point & X,
    const SGTELIB::Matrix & Jacobian_k,
    SGTELIB::Matrix & multiplier_k, // Modified
    const SGTELIB::Matrix & lambda,
    double mu,
    bool * active,
    bool * infeasible) const
{
    const int nbVar = _n;
    const int nbActive = sum(active, _nbCons);

    SGTELIB::Matrix activeJacobian (  "activeJacobian", nbActive, nbVar);
    SGTELIB::Matrix temp_multiplier(  "temp_multiplier", nbActive, 1);
    SGTELIB::Matrix pseudoGradient (  "pseudoGradient", nbVar, 1);

    for (int w=0 ; w < _nbCons ; w++ )
    {
        multiplier_k.set(w, 0, - int(infeasible[w]) / mu + lambda.get(w, 0)); // 1 for constraints active
    }
    pseudoGradient = getModelLagGradient(X, multiplier_k);
    activeJacobian = getModelActiveJacobian(Jacobian_k, active);
    temp_multiplier = SGTELIB::Surrogate_PRS::compute_multiplier(pseudoGradient, activeJacobian);

    return compute_dual_residual(pseudoGradient, activeJacobian, temp_multiplier);
}

double NOMAD::QPSolverOptimize::piecewise_line_search(
    const Point &X,
    const SGTELIB::Matrix &d,
    const bool * active,
    const bool * feasible,
    const bool * infeasible,
    const SGTELIB::Matrix & lambda,
    double mu,
    double small_gamma, // = 1E-20
    double gamma_update, // = 1.5
    double delta /* = 1E-4 // Pk < (P0 - delta) */ ) const
{
    
    // Allocations
    // const int nbActive = sum(active, ncon);
    // double slope;
    double ak;
    bool* Ik = new bool [_nbCons];
    NOMAD::Point X_k(X);
    X_k = X;

    SGTELIB::Matrix multiplier_k = get_pseudo_multiplier(active, feasible, infeasible, lambda, mu);
    SGTELIB::Matrix pseudoGradient = getModelLagGradient(X, multiplier_k);

    SGTELIB::Matrix gamma (  "gamma", _nbCons, 1);
    gamma.fill(0.0);
    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X);
    SGTELIB::Matrix Jacobian = getModelJacobian(X);
    SGTELIB::Matrix jprod = SGTELIB::Matrix::product(Jacobian, d);

    // Step 1:
    ak = SGTELIB::Matrix::dot(d, pseudoGradient); // < 0
    if (ak >= 0)
    {
        std::cout << "piecewise_line_search: error slope should be negative." << std::endl;
        return 0.0;
    }
    for (int i = 0; i < _nbCons; ++i)
    {
        if (!active[i])
        {
            gamma[i] = -cons[i] / jprod[i];
        }
        Ik[i] = (gamma[i] > 0) && (!active[i]);
    }

    bool unbounded = (sum(Ik, _nbCons) == 0); // Step 2
    if (unbounded)
    {
        return 0.0;
    }

    int lk;
    double gamma_lk = 1;
    int k=0;
    bool OK = false; // (ak >= 0)
    while (!OK && !unbounded)
    { // Step 3
        gamma_lk = INF;
        lk = -1;
        for (int i = 0; i < _nbCons; ++i)
        {
            if (Ik[i])
            {
                if (gamma[i] < gamma_lk)
                {
                    lk = i;
                    gamma_lk = gamma[i];
                }
            }
        }
        if (lk == -1)
        {
            unbounded = sum(Ik, _nbCons) == 0;
            OK = true; // TODO: add error
            std::cout << "piecewise_line_search: step 3 failure." << std::endl;
        }
        else
        {
            ak += fabs(jprod[lk]);
            OK = (ak >= 0);
            Ik[lk] = false;
            unbounded = (!OK) && (sum(Ik, _nbCons) == 0);
            k++;
        }
    }
    for (int i = 0 ; i < d.get_nb_rows() ; i++)
    {
        X_k[i] = X[i] +  gamma_lk * d.get(i, 0);
    }
    
    double P0 = getPenalizedL1AugLagModelObj(X, cons, lambda, mu);
    getModelCons(&cons, X_k);
    double Pk = getPenalizedL1AugLagModelObj(X_k, cons, lambda, mu);
    OK = Pk < (P0 - delta);
//    std::cout << "L-S P0=" << P0 << " Pk=" << Pk << " g=" << gamma_lk << std::endl;
    while (!OK)
    {
        // cubic interpolation...
        gamma_lk /= gamma_update; // Armijo instead
        for (int i = 0 ; i < d.get_nb_rows() ; i++)
        {
            X_k[i] = X[i] +  gamma_lk * d.get(i, 0);
        }
        getModelCons(&cons, X_k);
        Pk = getPenalizedL1AugLagModelObj(X_k, cons, lambda, mu);
        OK = (Pk < (P0 - delta)) || (gamma_lk <= small_gamma);
//        std::cout << "L-S P0=" << P0 << " Pk=" << Pk << " g=" << gamma_lk << std::endl;
    }

    if (gamma_lk <= small_gamma)
    {
        std::cout << "piecewise_line_search: no sufficient decrease found." << std::endl;
    }
    
    delete [] Ik;

    return gamma_lk;
}

bool NOMAD::QPSolverOptimize::compute_vertical_step(
    const Point &X,
    SGTELIB::Matrix &v_k,
    const SGTELIB::Matrix &activeJacobian_k,
    const SGTELIB::Matrix &cons,
    const bool * active) const
{
    const int ncon = static_cast<int>(_nbCons);
    const int nbActive = activeJacobian_k.get_nb_rows();
    SGTELIB::Matrix activeCons("activeCons", nbActive, 1);
    int k=0;
    for (int i=0; i < ncon; i++)
    {
        if (active[i])
        {
            activeCons.set(i, 0, -cons.get(k,0));
            k ++;
        }
    }
    if (nbActive != k)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Active jacobian number of rows do not match active indices.");
    }
    v_k = SGTELIB::Matrix::solve_least_squares_SVD(activeJacobian_k, activeCons);
    return true;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::get_pseudo_multiplier(
    const bool * active,
    const bool * feasible,
    const bool * infeasible,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{
    const int ncon = static_cast<int>(_nbCons);
    SGTELIB::Matrix multiplier_k("multiplier_k", ncon, 1);
    for (int j=0; j < ncon; ++j)
    {
        if (infeasible[j])
        {
            multiplier_k.set(j, 0, lambda.get(j, 0) - 1 / mu);
        }
        else
        {
            multiplier_k.set(j, 0, lambda.get(j, 0));
        }
    }
    return multiplier_k;
}

bool NOMAD::QPSolverOptimize::compute_horizontal_step(
    const Point &X,
    SGTELIB::Matrix &h_k,
    const SGTELIB::Matrix &Jacobian_k,
    const bool * active,
    const bool * feasible,
    const bool * infeasible,
    const SGTELIB::Matrix &lambda_l,
    const double mu_l) const
{
    SGTELIB::Matrix activeJacobian_k = getModelActiveJacobian(Jacobian_k, active);
    SGTELIB::Matrix Z = activeJacobian_k.null_space();
    SGTELIB::Matrix multiplier_k = get_pseudo_multiplier(active, feasible, infeasible, lambda_l, mu_l);

    SGTELIB::Matrix Hlag_k = getModelLagHessian(X, multiplier_k);
    SGTELIB::Matrix ZLZ = SGTELIB::Matrix::product(Z.transpose(), Hlag_k, Z);
    SGTELIB::Matrix invZLZ = ZLZ.SVD_inverse();

    SGTELIB::Matrix Glag_k = getModelLagGradient(X, multiplier_k);
    SGTELIB::Matrix ZL = SGTELIB::Matrix::product(Z.transpose(), Glag_k);
    ZL.multiply(-1);
    SGTELIB::Matrix delta_k = SGTELIB::Matrix::product(invZLZ, ZL);
    h_k = SGTELIB::Matrix::product(Z, delta_k);

    double slope = SGTELIB::Matrix::dot(h_k, Glag_k);
    if (slope >= 0)
    {
        SGTELIB::Matrix::inplace_product(h_k, Z, ZL);
    }

    return true;
}

//*****************************************************************************
//
// Augmented Lagrangian algorithm and functions
//
//*****************************************************************************

bool NOMAD::QPSolverOptimize::solveAugLag(
    NOMAD::Point & X,
    const int max_iter,
    const double tolDistDX, // 1E-15 ;
    const double atol,
    const double rtol,
    const double mu0,
    const double muDecrease,
    const double eta0,
    const double omega0,
    const double successRatio, // 0.99
    const size_t maxIterInner, // = 50;
    const double tolDistDXInner, // = 1E-15;
    const size_t maxSuccessivFail // = 3;
   )
{

    const int nbVar = _n + _nbCons; // we use slack variables here

    double fk = getModelObj(X);
    double f0 = fk;
    SGTELIB::Matrix Gk("Gk", nbVar, 1);
    getModelGrad(&Gk, X);
    double ng = Gk.norm();
    double ng0 = ng;


    double tol = atol + ng0 * rtol;

    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X);

    SGTELIB::Matrix lvar("lvar", nbVar, 1);
    lvar.fill(0.0);
    SGTELIB::Matrix uvar("uvar", nbVar, 1);
    uvar.fill(INF);
    SGTELIB::Matrix XS("XS", nbVar, 1);
    for (int i=0; i < _n; ++i)
    {
        double lb = ((_modelLowerBound[i].isDefined())? _modelLowerBound[i].todouble(): NOMAD::M_INF);
        double ub = ((_modelUpperBound[i].isDefined())? _modelUpperBound[i].todouble(): NOMAD::INF);
        lvar.set(i, 0, lb);
        uvar.set(i, 0, ub);
        XS.set(i, 0, X[i].todouble());
    }

    for (int j=0; j < _nbCons; ++j)
    {
        XS.set(j + _n, 0, -cons.get(j, 0)); // S = -cons
    }


    snapToBounds(XS, lvar, uvar);
    bool bound_compat = true;
    bool feasible = true;
    for (int i=0; i < nbVar; ++i)
    {
        bound_compat = bound_compat && (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!bound_compat)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solveAugLag assertion error: Error compatibility lower and upper bound");
        }

        feasible = feasible && (XS.get(i, 0) >= lvar.get(i, 0));
        feasible = feasible && (XS.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            std::cout << XS.get(i, 0) - lvar.get(i, 0) << " " << uvar.get(i, 0) - XS.get(i, 0) << std::endl;
            std::cout << "Error compatibility lower and upper bound" << std::endl;
            throw NOMAD::Exception(__FILE__, __LINE__, "solveAugLag assertion error: Error XS is not feasible");
        }
    }

    double mu_l = mu0;
    double mu_decrease = muDecrease;
    double eta_l = eta0;
    double omega_l = omega0;

    solveLM(X, XS, lvar, uvar, cons, mu_l, omega_l, 30, 1E-15, false, _verbose);

    SGTELIB::Matrix cslack("c+s", _nbCons, 1);
    for (int j=0; j < _nbCons; ++j)
    {
        cslack.set(j, 0, XS.get(j + _n, 0) + cons.get(j, 0));
    }
    double cx = cslack.norm();
    double cxp = cx;

    fk = getModelObj(X);
    getModelGrad(&Gk, X);

    // Next iterate
    SGTELIB::Matrix XSp("XSp", nbVar, 1);
    Point Xp(X);

    // Outer loop parameters
    SGTELIB::Matrix lambda_l("lambda_l", _nbCons, 1);
    lambda_l.fill(0.0);

    double Pk = getAugLagModelObj(XS, cons, fk, lambda_l, mu_l);
    double Pkp;
    SGTELIB::Matrix GradPk("GradAugLag", nbVar, 1);
    SGTELIB::Matrix HessPk("HessAugLag", nbVar, nbVar);
    getAugLagModelGrad(&GradPk, XS, lambda_l, mu_l);

    int iterOuterLoop = 0;

    double distXOutterLoop = INF;

    SGTELIB::Matrix DualFeas("DualFeas", nbVar, 1);
    double ngproj = check_optimality_bounds(XS, GradPk, lvar, uvar, DualFeas);
    bool outerSuccess = ngproj < tol; // nonlinear equality constraints feasible
    bool outerFailure = (distXOutterLoop <= tolDistDX) || (iterOuterLoop >= max_iter);

    _verbose && std::cout << "Outer ("<< iterOuterLoop <<")  : |G|= " << GradPk.norm();
    _verbose && std::cout << " |Proj(x - grad) - x|= " << ngproj << " P=" << Pk;
    _verbose && std::cout << " f=" << fk << " |c+s|=" << cx;
    _verbose && std::cout << " mu=" << mu_l << " omega=" << omega_l << " eta=" << eta_l << std::endl;

    int innerResult;
    bool innerSuccess;
    size_t successivFailure = 0;
    size_t successivAcceptable = 0;
    size_t successivBeforeUpdate = 2;

    // outer iteration
    while (!outerFailure && !outerSuccess)
    {
        // Solve non-linear bound-constrained subproblem
        innerResult = solve_bound_AugLag(XSp, XS, lvar, uvar, lambda_l, omega_l, mu_l, GradPk, HessPk, maxIterInner, tolDistDXInner, _verboseFull);
        innerSuccess = innerResult > 0;

        for (int i=0; i < _n; ++i)
        {
            Xp[i] = XSp.get(i, 0);
        }
        getModelCons(&cons, Xp);
        for (int j=0; j < _nbCons; ++j)
        {
            cslack.set(j, 0, XSp.get(j + _n, 0) + cons.get(j, 0));
        }
        cxp = cslack.norm();
        fk = getModelObj(Xp);
        Pkp = getAugLagModelObj(XSp, cons, fk, lambda_l, mu_l);
        
        // Update parameters
        if ((cxp <= eta_l) && innerSuccess)
        {
            cslack.multiply(- 1 / mu_l); // lambda_l.add(SGTELIB::Matrix::product(cslack, -1 / mu_l));
            lambda_l.add(cslack);

            eta_l = std::max(eta_l * std::pow(mu_l, 0.9), atol);
            if (innerResult == 1)
            {
                successivAcceptable = 0;
                omega_l *= mu_l;
            }
            else if (successivAcceptable >= successivBeforeUpdate)
            {
                successivAcceptable += 1;
                omega_l *= std::sqrt(mu_l);
            }
            else
            {
                successivAcceptable += 1; // Try re-run from new point
            }
            omega_l = std::max(omega_l, std::max(atol * atol, 1E-15));
        }
        else if (innerSuccess) // not feasible
        {
            if (innerResult == 1)
            {
                successivAcceptable = 0;
                mu_l /= mu_decrease;
            }
            else if (successivAcceptable >= successivBeforeUpdate)
            {
                successivAcceptable += 1;
                mu_l /= sqrt(mu_decrease);
            }
            else
            {
                successivAcceptable += 1; // Try re-run from new point
            }
            eta_l = std::max(std::pow(mu_l, 0.1), atol);
            omega_l = mu_l;
        }
        else
        { // no success and not feasible
            successivAcceptable = 0;
            mu_l /= mu_decrease;
            eta_l = std::max(std::pow(mu_l, 0.1), atol);
            omega_l = mu_l;
        }

        // Check optimality and accept new iterate
        getAugLagModelGrad(&GradPk, XSp, lambda_l, mu_l);
        ngproj = check_optimality_bounds(XSp, GradPk, lvar, uvar, DualFeas);
        outerSuccess = (ngproj < tol) && (cxp < tol);

        if (innerSuccess || (Pkp < successRatio * Pk)) // we accept this step
        {
            distXOutterLoop = NOMAD::Point::dist(X, Xp).todouble(); // diff XSp, XS
            XS = XSp;
            X = Xp;
            Pk = Pkp;
            cx = cxp;
            successivFailure = 0;
        }
        else
        { // for verbose only
            size_t feasibilityPhase = solveLM(Xp, XSp, lvar, uvar, cons, mu_l, omega_l, 30, 1E-15, false, _verbose);
            if (feasibilityPhase > 0 )
            {
                XS = XSp;
                X = Xp;
                for (int j=0; j < _nbCons; ++j)
                {
                    cslack.set(j, 0, XS.get(j + _n, 0) + cons.get(j, 0));
                }
                cxp = cslack.norm();
                cx = cxp;
                fk = getModelObj(X);
                successivFailure = 0;
            }
            else
            {
                successivFailure += 1;
            }
            getAugLagModelGrad(&GradPk, XSp, lambda_l, mu_l);
            ngproj = check_optimality_bounds(XSp, GradPk, lvar, uvar, DualFeas);
        }

        iterOuterLoop += 1;
        outerFailure = (iterOuterLoop >= max_iter) || (distXOutterLoop <= tolDistDX);
        outerFailure = outerFailure || (mu_l <= atol / mu_decrease);
        outerFailure = outerFailure || (successivFailure >= maxSuccessivFail);
    
        _verbose && std::cout << "Outer ("<< iterOuterLoop <<") " << innerResult;
        _verbose && std::cout << ": |Proj(x - grad) - x|= " << ngproj << " P=" << Pk << " f=" << fk ;
        _verbose && std::cout << " |c+s|=" << cx << " mu=" << mu_l << " omega=" << omega_l << " eta=" << eta_l << std::endl;

    }

    // Ending output
    _verbose && std::cout << "End of solveAugLag" << std::endl;
    _verbose && std::cout << "f(x0)=" << f0 << " f(x*)=" << fk <<std::endl;
    _verbose && std::cout << "|c(x*) + s*|=" << cxp << " tol=" << tol << std::endl;
    if (outerSuccess)
    {
        _verbose && std::cout << "success" << " : " << ngproj << " < " << tol << std::endl;
    }
    else if (outerFailure)
    {
        _verbose && std::cout << " failure" << " iter=" << iterOuterLoop << std::endl;
        _verbose && std::cout << " dist=" << distXOutterLoop << " <=" << tolDistDX << std::endl;
        _verbose && std::cout << " too small parameters? " << mu_l * mu_l << "<=" << atol << std::endl;
        _verbose && std::cout << " successivFailure=" << successivFailure << std::endl;
    }
    else
    {
        _verbose && std::cout << "unknown stopping";
    }
    return true;
}

int NOMAD::QPSolverOptimize::solve_bound_AugLag(
    SGTELIB::Matrix & XSp,
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    const SGTELIB::Matrix & lambda,
    double omega,
    double mu,
    SGTELIB::Matrix & GradPk,
    SGTELIB::Matrix & HessPk,
    const size_t maxIterInnerLoop, // = 50
    const double tolerance_distDX, // = 1E-15
    bool verbose )
{
    const int nvar = XS.get_nb_rows();
    const int ncon = static_cast<int>(_nbCons);

    lencheck(nvar, XS);
    lencheck(nvar, XSp);
    lencheck(nvar, lvar);
    lencheck(nvar, uvar);
    lencheck(ncon, lambda);

    if (mu <= 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error mu must be positive");
    }
    if (omega < 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error omega must be non-negative");
    }

    // initialization:
    XSp = XS;
    snapToBounds(XSp, lvar, uvar);
    bool bound_compat = true;
    bool feasible = true, dfeasible;
    for (int i=0; i < nvar; ++i)
    {
        bound_compat = bound_compat && (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!bound_compat)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solve_bound_AugLag assertion error: Error compatibility lower and upper bound");
        }

        feasible = feasible && (XSp.get(i, 0) >= lvar.get(i, 0));
        feasible = feasible && (XSp.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solve_bound_AugLag assertion error: Error XS is not feasible");
        }
    }

    // Inner loop parameters
    double ared, pred, nd;
    const double epsilon_1 = 0.05;
    const double epsilon_2 = 0.9;
    const double gamma_1 = 0.5; // 0.01;
    const double gamma_2 = 2; // 100;

    double delta = 1.0; // trust region initial radius
    const double smallestDelta = 1E-15;
    double largestDelta = 1E15;

    const size_t limitUnsuccessful = 40;
    const size_t maxIterBCQP = 50;

    // Pre-allocations
    SGTELIB::Matrix Xcan("Xcan", nvar, 1);
    SGTELIB::Matrix d("d", nvar, 1);
    SGTELIB::Matrix DualFeas("DualFeas", nvar, 1);
    SGTELIB::Matrix dlvar("dlvar", nvar, 1);
    SGTELIB::Matrix duvar("duvar", nvar, 1);
    
    double Pk = getAugLagModelObj(XSp, lambda, mu);
    getAugLagModelHess(&HessPk, XSp, lambda, mu);
    
    size_t iterInnerLoop = 0;
    size_t successivUnsuccessful = 0;
    double distXInnerLoop = INF;

    double ngproj = check_optimality_bounds(XS, GradPk, lvar, uvar, DualFeas);

    double tolSuccess;
    if (omega >= 1)
    {
        tolSuccess = omega;
    }
    else
    {
        tolSuccess = omega * (1 + ngproj); // if omega >=1, this is always satisfied.
    }
    bool innerSuccess = ngproj < tolSuccess; // nonlinear equality constraints feasible
    bool subpbSuccess = false;
    bool innerFailure = (distXInnerLoop <= tolerance_distDX) || (iterInnerLoop >= maxIterInnerLoop);

    double atol_BCQP, rtol_BCQP; // TODO: fine-tune these 2 in function of omega, 1E-7, 1E-15 and norm of Hessian.
    atol_BCQP = std::min(1E-8, omega); // by default: omega
    rtol_BCQP = 1E-15 * HessPk.norm() / nvar; // 1E-7 by default

    verbose && std::cout << "Inner 0 ("<< iterInnerLoop <<")  : |Proj(x - grad) - x|= " << ngproj << " P=" << Pk << " theo. tol=" << rtol_BCQP << std::endl;
    
    bool success = innerSuccess; // overall success
    
    innerSuccess = false; // ChT To enter the loop once

    while (!innerFailure && !innerSuccess)
    {
        // compute d;
        if (successivUnsuccessful == 0)
        {
            getAugLagModelGrad(&GradPk, XSp, lambda, mu);
            getAugLagModelHess(&HessPk, XSp, lambda, mu);
        }

        dfeasible = true;
        for (int j=0; j < nvar; ++j)
        {
            dlvar.set(j, 0, std::max(lvar.get(j, 0) - XSp.get(j, 0), -delta));
            duvar.set(j, 0, std::min(uvar.get(j, 0) - XSp.get(j, 0), delta));
            bound_compat = bound_compat && (dlvar.get(j, 0) <= duvar.get(j, 0));

            feasible = feasible && (dlvar.get(j, 0) <= 0) && (duvar.get(j, 0) >= 0);
            dfeasible = dfeasible && (dlvar.get(j, 0) <= d.get(j, 0)) && (duvar.get(j, 0) >= d.get(j, 0));

            // Init scaling
            //scaling.set(j, j, 1 / sqrt(abs(HessPk.get(j, j)) + 1));
            //unscaling.set(j, j, sqrt(abs(HessPk.get(j, j)) + 1));
        }
        
        if (!bound_compat)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solve_bound_AugLag assertion error (d): Error compatibility lower and upper bound");
        }

        if ((successivUnsuccessful == 0) || !dfeasible)
        { // we updated x OR the last direction is still feasible
            d.fill(0);
            if (!feasible)
            {
                // throw NOMAD::Exception(__FILE__, __LINE__, "solve_bound_AugLag assertion error: Error d is not feasible");
                std::cerr << "solve_bound_AugLag assertion error: Error d is not feasible " << delta << std::endl;
                //dlvar.display(std::cout);
                //duvar.display(std::cout);
                //XSp.display(std::cout);
                return false;
            }

            /* Try scaling
            HessPk.display(std::cout); GradPk.display(std::cout);
            GradPk = SGTELIB::Matrix::product(scaling, GradPk);
            HessPk = SGTELIB::Matrix::product(scaling, HessPk, scaling);
            HessPk.display(std::cout); GradPk.display(std::cout);
            dlvar = SGTELIB::Matrix::product(scaling, dlvar);
            duvar = SGTELIB::Matrix::product(scaling, duvar);
            */
            
            subpbSuccess = solveBCQP(d, HessPk, GradPk, 0.0, dlvar, duvar, maxIterBCQP, atol_BCQP, rtol_BCQP); // rtol_BCQP, specify verbose

            /* Try scaling
                    d.display(std::cout);
                    d = SGTELIB::Matrix::product(unscaling, d);
                    d.display(std::cout);
            */
        }
        else if (!dfeasible)
        {
            // throw NOMAD::Exception(__FILE__, __LINE__, "solve_bound_AugLag assertion error: Error d is not feasible");
            std::cerr << "solve_bound_AugLag assertion error: Error d is not feasible" << std::endl;
            return false;
        }

        pred = getModelObj(d, HessPk, GradPk);
        if (pred > 0)
        {
            std::cerr << "Assertion error: prediction " << pred << " > 0" <<std::endl;
            return false;
        }

        // Trust-region update
        Xcan = XSp; Xcan.add(d);
        nd = d.norm();

        ared = compute_AugLag_TR_ared(XSp, Xcan, lambda, mu);

        bool pointAccepted = false;
        // Ratio: ared / pred;
        if (subpbSuccess && (ared <= pred * epsilon_1)) // r >= epsilon_1
        {
            XSp = Xcan; // accept the point
            pointAccepted = true;
            if (ared <= pred * epsilon_2) // r >= epsilon_2
            {
                delta = std::min(gamma_2 * delta, std::max(1 / omega, largestDelta)); // std::max(d.norm(), delta);
            }
            success = true;
            successivUnsuccessful = 0;
        }
        else
        {
            delta = std::max(gamma_1 * std::min(delta, nd), smallestDelta);
            successivUnsuccessful += 1;
        }

        // Check optimality
        getAugLagModelGrad(&GradPk, XSp, lambda, mu);
        ngproj = check_optimality_bounds(XSp, GradPk, lvar, uvar, DualFeas);
        innerSuccess = ngproj < tolSuccess;

        iterInnerLoop += 1;
        distXInnerLoop = nd;
        innerFailure = !subpbSuccess || (iterInnerLoop >= maxIterInnerLoop) || (distXInnerLoop <= tolerance_distDX);
        innerFailure = innerFailure || (successivUnsuccessful > limitUnsuccessful);
    
        Pk = getAugLagModelObj(XSp, lambda, mu);
        verbose && std::cout << "Inner ("<< iterInnerLoop <<") " << subpbSuccess;
        verbose && std::cout << " " << pointAccepted;
        verbose && std::cout << " |Proj(x - grad) - x|= " << ngproj;
        verbose && std::cout << " P=" << Pk << " |d|=" << distXInnerLoop;
        verbose && std::cout  << " delta=" << delta << " r=";
        if (pred == 0)
        {
            verbose && std::cout << ared << " / " << pred << std::endl;
        } else
        {
            verbose && std::cout << ared / pred << std::endl;
        }
    }

    int result;
    // Success if at least one successful iteration
    // bool success = (iterInnerLoop >= 2) || !innerFailure;
    if (innerSuccess)
    {
        result = 1;
    }
    else if (!success)
    {
        verbose && std::cout << "Trust-region failure: return XS" << std::endl; 
        XSp = XS;
        result = 0;
    }
    else
    {
        result = 2; // acceptable
    }
    return result;
}

double NOMAD::QPSolverOptimize::compute_AugLag_TR_ared(
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & XSp,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{

    const int nbVar = _n + _nbCons;

    lencheck(nbVar, XS);
    lencheck(nbVar, XSp);
    lencheck(_nbCons, lambda);

    double ared;
    ared = getAugLagModelObj(XSp, lambda, mu) - getAugLagModelObj(XS, lambda, mu);

    return ared;
}

//*****************************************************************************
//
// Augmented Lagrangian with slack variables model getters
//
//*****************************************************************************

double NOMAD::QPSolverOptimize::getAugLagModelObj(
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{
    SGTELIB::Matrix X("X", _n, 1);
    for (int i=0; i < _n; ++i)
    {
        X.set(i, 0, XS.get(i, 0));
    }

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    SGTELIB::Matrix cons = surrogate_prs->getModelCons(X.transpose());
    double lag = surrogate_prs->getModelObj(X.transpose());

    return getAugLagModelObj(XS, cons, lag, lambda, mu);
}

double NOMAD::QPSolverOptimize::getAugLagModelObj(
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & cons,
    double fx,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{
    const int nbVar = _n + _nbCons;
    
    lencheck(nbVar, XS);
    lencheck(_nbCons, lambda);
    lencheck(_nbCons, cons);

    double lag = fx;

    double cpsi; // cons + S
    for (int i=0; i < _nbCons; i++)
    {
        cpsi = XS.get(i + _n, 0) + cons.get(i, 0);
        lag -= lambda[i] * cpsi;
        lag += 1/(2 * mu) * std::pow(cpsi, 2);
    }

    return lag;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getAugLagModelGrad(
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{
    const int nbVar = _n + _nbCons;
    SGTELIB::Matrix lagGrad("lagGradient", nbVar, 1);
    getAugLagModelGrad(&lagGrad, XS, lambda, mu);
    return lagGrad;
}

void NOMAD::QPSolverOptimize::getAugLagModelGrad(
    SGTELIB::Matrix * lagGrad,
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{

    const int nbVar = _n + _nbCons;

    lencheck(nbVar, XS);

    SGTELIB::Matrix X("X", _n, 1);
    SGTELIB::Matrix S("S", _nbCons, 1);

    for (int i=0; i < _n; ++i)
    {
        X.set(i, 0, XS.get(i, 0));
    }
    for (int j=0; j < _nbCons; ++j)
    {
        S.set(j, 0, XS.get(j + _n, 0));
    }

    lencheck(_nbCons, lambda);

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);

    lagGrad->fill(0.0);

    SGTELIB::Matrix lagGradS("temp", _nbCons, 1);
    surrogate_prs->getModelCons(&lagGradS, X.transpose());
    lencheck(_nbCons, lagGradS);
    lagGradS.add(S);
    lagGradS.multiply(-1 / mu);
    lagGradS.add(lambda);

    // Derivative w.r.t. X
    SGTELIB::Matrix lagGradX("tempX", _n, 1);
    SGTELIB::Matrix Mpredict_grad ("grad_predict", _nbCons + 1, _n);
    SGTELIB::Matrix Jx ("Jx", _nbCons, _n);
    surrogate_prs->getModelLagGrad(&lagGradX, &Mpredict_grad, &Jx, X.transpose(), lagGradS);

    // Derivative w.r.t. S
    for (int i=0; i < _n; i++)
    {
        lagGrad->set(i, 0, lagGradX.get(i, 0));
    }
    for (int i=0; i < _nbCons; i++)
    {
        lagGrad->set(i + _n, 0, -lagGradS.get(i, 0));
    }
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getAugLagModelHess(
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{
    int nbVar = _n + _nbCons;
    SGTELIB::Matrix lagHess("lagHess", nbVar, nbVar);
    getAugLagModelHess(&lagHess, XS, lambda, mu);
    return lagHess;
}

void NOMAD::QPSolverOptimize::getAugLagModelHess(
    SGTELIB::Matrix * lagHess,
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lambda,
    double mu ) const
{
    const int nbVar = _n + _nbCons;

    lencheck(nbVar, XS);

    SGTELIB::Matrix X("X", _n, 1);
    SGTELIB::Matrix S("S", _nbCons, 1);

    for (int i=0; i < _n; ++i)
    {
        X.set(i, 0, XS.get(i, 0));
    }
    for (int j=0; j < _nbCons; ++j)
    {
        S.set(j, 0, XS.get(j + _n, 0));
    }

    lencheck(_nbCons, lambda);

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);

    lagHess->fill(0.0);
    SGTELIB::Matrix lagGradS("temp", _nbCons, 1);
    surrogate_prs->getModelCons(&lagGradS, X.transpose());
    lagGradS.add(S);
    lagGradS.multiply(-1 / mu);
    lagGradS.add(lambda);

    // Derivative w.r.t. X
    SGTELIB::Matrix lagHessX = surrogate_prs->getModelLagHessian(X.transpose(), lagGradS);

    SGTELIB::Matrix Jac = surrogate_prs->getModelJacobian(X.transpose());
    sizecheck(_nbCons, _n, Jac);

    SGTELIB::Matrix JtJ = SGTELIB::Matrix::product(Jac.transpose(), Jac);
    JtJ.multiply(1 / mu);
    lagHessX.add(JtJ);

    for (int i=0; i < _n; i++)
    {
        for (int j=0; j < _n; j++)
        {
            lagHess->set(i, j, lagHessX.get(i, j));
        }
    }
    for (int i=0; i < _nbCons; i++)
    {
        for (int j=0; j < _n; j++)
        {
            lagHess->set(i + _n, j, Jac.get(i, j) / mu);
            lagHess->set(j, i + _n, Jac.get(i, j) / mu);
        }
    }
    for (int i=0; i < _nbCons; i++)
    {
        lagHess->set(i + _n, i + _n, 1 / mu);
    }
}

bool NOMAD::QPSolverOptimize::getStrictFeasiblePoint(
    NOMAD::Point &X,
    SGTELIB::Matrix & XS,
    SGTELIB::Matrix & lvar,
    SGTELIB::Matrix & uvar,
    const SGTELIB::Matrix & cX )
{
    double lb, ub, mid, xi;
    for (int i=0; i < _n; ++i)
    {
        lb = ((_modelLowerBound[i].isDefined())? _modelLowerBound[i].todouble(): -NOMAD::INF);
        ub = ((_modelUpperBound[i].isDefined())? _modelUpperBound[i].todouble(): NOMAD::INF);
        lvar.set(i, 0, lb);
        uvar.set(i, 0, ub);
        xi = X[i].todouble();
        if ((xi <= lb) || (xi >= ub))
        {
            if ((_modelLowerBound[i].isDefined()) && !(_modelUpperBound[i].isDefined()))
            {
                xi = lb + 0.5;
            }
            else if(!(_modelLowerBound[i].isDefined()) && (_modelUpperBound[i].isDefined()))
            {
                xi = ub - 0.5;
            }
            else if((_modelLowerBound[i].isDefined()) && (_modelUpperBound[i].isDefined()))
            {
                mid = uvar.get(i, 0) - lvar.get(i, 0);
                xi = lb + mid / 2;
            }
            else
            {
                xi = 0.0;
            }
        }
        XS.set(i, 0, xi);
    }

    for (int j=0; j < _nbCons; ++j)
    {
        lvar.set(j + _n, 0, 0.0);
        uvar.set(j + _n, 0, INF);
        XS.set(j + _n, 0, std::max(-cX.get(j, 0), 0.5)); // S = -cons
    }

    return check_strict_feasible(XS, lvar, uvar);
}

bool NOMAD::QPSolverOptimize::solveTRIPM(
    NOMAD::Point &X,
    const int max_iter, // 30
    const double tolDistDX, // -1
    const double atol, // 1E-7
    const double rtol, // 1E-7
    const double mu0, // 0.5
    const double muDecrease, // 2
    const size_t maxIterInner, // 40
    bool verbose, // true
    bool verbose_SolverBarrier) // true
{

    const int nbVar = _n + _nbCons; // we use slack variables here

    double mu = mu0;
    if (muDecrease <= 1)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, " muDecrease must be > 1");
    }
    double mu_decrease = muDecrease;
    double tol_mu = mu;
    double smallest_tol_mu = atol / 100;

    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X);

    SGTELIB::Matrix lvar("lvar", nbVar, 1);
    SGTELIB::Matrix uvar("uvar", nbVar, 1);
    SGTELIB::Matrix XS("XS", nbVar, 1);
    SGTELIB::Matrix p("p", nbVar, 1);

    getStrictFeasiblePoint(X, XS, lvar, uvar, cons);
    solveLM(X, XS, lvar, uvar, cons, mu, tol_mu, 30, 1E-15, true, verbose);

    double fk = getModelObj(X);
    double f0 = fk;
    SGTELIB::Matrix Gk("Gk", nbVar, 1);
    getModelGrad(&Gk, X);
    double ng = Gk.norm();
    double ng0 = ng;
    SGTELIB::Matrix Jx = getModelJacobian(X);

    bool bound_compat = true;
    bool feasible = true;
    for (int i=0; i < nbVar; ++i)
    {
        bound_compat = bound_compat && (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!bound_compat)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solveAugLag assertion error: Error compatibility lower and upper bound");
        }

        feasible = feasible && (XS.get(i, 0) >= lvar.get(i, 0));
        feasible = feasible && (XS.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            std::cout << XS.get(i, 0) - lvar.get(i, 0) << " " << uvar.get(i, 0) - XS.get(i, 0) << std::endl;
            std::cout << "Error compatibility lower and upper bound" << std::endl;
            throw NOMAD::Exception(__FILE__, __LINE__, "solveAugLag assertion error: Error XS is not feasible");
        }
    }

    SGTELIB::Matrix cslack("c+s", _nbCons, 1);
    for (int j=0; j < _nbCons; ++j)
    {
        cslack.set(j, 0, XS.get(j + _n, 0) + cons.get(j, 0));
    }
    double cx = cslack.norm();
    double cxp = cx;

    // Next iterate
    SGTELIB::Matrix XSp("XSp", nbVar, 1);
    Point Xp(X);

    int innerResult;
    bool innerSuccess;

    double tolDistDXInner = 1E-15;
    size_t maxSuccessivFail = 3;

    // Lagrange multiplier estimate
    SGTELIB::Matrix lambda("lambda", _nbCons, 1);
    compute_slack_multiplier(lambda, XS, Jx, Gk, 0);

    int iterOuterLoop = 0;
    double distXOutterLoop = INF;
    size_t successivFailure = 0;
    size_t successivAcceptable = 0;
    size_t successivBeforeUpdate = 2;

    bool outerSuccess, outerFailure;
    double res = errorTRIPM(XS, lvar, uvar, lambda, cslack, 0);
    double tol = atol + std::max(ng0, res) * rtol;

    outerSuccess = (res <= tol);
    outerFailure = (distXOutterLoop <= tolDistDX) || (iterOuterLoop >= max_iter);
 
    verbose && std::cout << "Outer ("<< iterOuterLoop <<"): ";
    verbose && std::cout << " |E(x,s,y;0)|= " << res;
    verbose && std::cout << " f=" << fk << " |c+s|=" << cx;
    verbose && std::cout << " mu=" << mu << " e_mu=" << tol_mu << std::endl;

    // outer iteration
    while (!outerFailure && !outerSuccess)
    {
        innerResult = solver_barrier(Xp, XSp, p, cslack, XS, lvar, uvar, lambda, Gk, cons, Jx, mu, tol_mu, maxIterInner, tolDistDXInner, verbose_SolverBarrier);
        innerSuccess = innerResult > 0;

        for (int i=0; i < _n; ++i)
        {
            Xp[i] = XSp.get(i, 0);
        }
        getModelCons(&cons, Xp);
        for (int j=0; j < _nbCons; ++j)
        {
            cslack.set(j, 0, XSp.get(j + _n, 0) + cons.get(j, 0));
        }
        cxp = cslack.norm();
        fk = getModelObj(Xp);

        if (innerSuccess)
        {
            distXOutterLoop = NOMAD::Point::dist(X, Xp).todouble(); // diff XSp, XS
            XS = XSp;
            X = Xp;
            cx = cxp;
            successivFailure = 0;
            if (innerResult == 1)
            {
                successivAcceptable = 0;
                mu /= mu_decrease;
                tol_mu /= mu_decrease;
            }
            else if ((successivFailure > 0) || (successivAcceptable >= successivBeforeUpdate))
            { // slower decrease
                successivAcceptable = 0;
                mu /= sqrt(mu_decrease);
                tol_mu /= sqrt(mu_decrease);
            }
            else
            {
                // Try re-run from new point.
                successivAcceptable += 1;
            }
            res = errorTRIPM(XS, lvar, uvar, lambda, cslack, 0);
            outerSuccess = (res <= tol);
        }
        else
        {
            successivAcceptable = 0;
            successivFailure += 1;

            if (cxp > tol_mu)
            {
                size_t feasibilityPhase = solveLM(Xp, XSp, lvar, uvar, cons, mu, tol_mu, 30, 1E-15, true, verbose);
                if (feasibilityPhase > 0 )
                {
                    XS = XSp;
                    X = Xp;
                    for (int j=0; j < _nbCons; ++j)
                    {
                        cslack.set(j, 0, XS.get(j + _n, 0) + cons.get(j, 0));
                    }
                    cxp = cslack.norm();
                    cx = cxp;
                    fk = getModelObj(X);
                    res = errorTRIPM(XS, lvar, uvar, lambda, cslack, 0);
                    outerSuccess = (res <= tol);
                    successivFailure = 0;
                }
                else
                {
                    mu /= mu_decrease;
                    tol_mu /= mu_decrease;
                }
            }
            else
            {
                mu /= mu_decrease;
                tol_mu /= mu_decrease;
            }        
        }
        tol_mu = std::max(smallest_tol_mu, tol_mu);

        iterOuterLoop += 1;
        outerFailure = (distXOutterLoop <= tolDistDX) || (iterOuterLoop >= max_iter);
        outerFailure = outerFailure || (mu <= atol / mu_decrease);
        outerFailure = outerFailure || (successivFailure >= maxSuccessivFail);

        verbose && std::cout << "Outer ("<< iterOuterLoop <<") " << innerResult;
        verbose && std::cout << ": |E(x,s,y;0)|= " << res << " f=" << fk ;
        verbose && std::cout << " |c+s|=" << cx << " mu=" << mu << " e_mu=" << tol_mu << std::endl;
    
    }

    // Ending output
    verbose && std::cout << "End of solveTR-IPM" << std::endl;
    verbose && std::cout << "f(x0)=" << f0 << " f(x*)=" << fk <<std::endl;
    verbose && std::cout << "|c(x*) + s*|=" << cxp << " tol=" << tol << std::endl;
    if (outerSuccess)
    {
        verbose && std::cout << "success" << " : " << res << " < " << tol << std::endl;
    }
    else if (outerFailure)
    {
        verbose && std::cout << " failure" << " iter=" << iterOuterLoop << std::endl;
        verbose && std::cout << " dist=" << distXOutterLoop << " <=" << tolDistDX << std::endl;
        verbose && std::cout << " too small parameters? " << mu << "<=" << atol << std::endl;
        verbose && std::cout << " successivFailure=" << successivFailure << std::endl;   
    }
    else
    {
        verbose && std::cout << "unknown stopping";
    }

    bool success = (fk < f0);
    return success;
}

void NOMAD::QPSolverOptimize::compute_slack_multiplier(
    SGTELIB::Matrix & y,
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & Jx,
    const SGTELIB::Matrix & Gx,
    const double mu )
{

    const int nbVar = _n + _nbCons;

    lencheck(_nbCons, y);
    lencheck(nbVar, XS);
    sizecheck(_nbCons, _n, Jx);

    SGTELIB::Matrix W("W", nbVar, _nbCons);
    SGTELIB::Matrix bls("bls", nbVar, 1);
    for (int i=0; i < _n; i++)
    {
        for (int j=0; j < _nbCons; j++)
        {
            W.set(i, j, Jx.get(j, i));
        }
        bls.set(i, 0, Gx.get(i, 0));
    }
    for (int i=0; i < _nbCons; i++)
    {
        for (int j=0; j < _nbCons; j++)
        {
            if (i == j)
            {
                W.set(_n + i, j, XS.get(_n + i));
            }
            else
            {
                W.set(_n + i, j, 0);
            }
        }
        bls.set(_n + i, 0, -mu);
    }
    // LS solve of ||Wy - bls||
    y = SGTELIB::Matrix::solve_least_squares_SVD(W, bls);

    // enforce sign of y
    double si;
    for (int i=0; i < _nbCons; i++)
    {
        if (y.get(i, 0) >= 0)
        {
            si = XS.get(_n + i, 0);
            y.set(i, 0, std::min(-1E-3, -mu / si)); // TODO: or max 
        }
    }

}

double NOMAD::QPSolverOptimize::errorTRIPM(
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    const SGTELIB::Matrix & lambda,
    const SGTELIB::Matrix & cslack,
    const double mu )
{

    const int nbVar = _n + _nbCons;

    lencheck(nbVar, XS);
    lencheck(_nbCons, lambda);
    lencheck(_nbCons, cslack);

    SGTELIB::Matrix X("X", _n, 1);

    for (int i=0; i < _n; ++i)
    {
        X.set(i, 0, XS.get(i, 0));
    }

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);

    // Derivative w.r.t. X
    SGTELIB::Matrix lagGradX("tempX", _n, 1);
    SGTELIB::Matrix Mpredict_grad ("grad_predict", _nbCons + 1, _n);
    SGTELIB::Matrix Jx ("Jx", _nbCons, _n);
    surrogate_prs->getModelLagGrad(&lagGradX, &Mpredict_grad, &Jx, X.transpose(), lambda);

    SGTELIB::Matrix dual_feas = SGTELIB::Matrix("dual_feas", _n, 1);

    for (int i = 0 ; i < _n ; i++)
    {
        // lagGradX.set(i, 0, lagGradX.get(i, 0) - mu / (X.get(i, 0) - lvar.get(i, 0)) + mu / (uvar.get(i, 0) - X.get(i, 0)));
        dual_feas.set(i, 0, X.get(i, 0) - lagGradX.get(i, 0));
        if (dual_feas.get(i, 0) < lvar.get(i, 0))
        {
            dual_feas.set(i, 0, lvar.get(i, 0));
        }
        else if (uvar.get(i, 0) < dual_feas.get(i, 0))
        {
            dual_feas.set(i, 0, uvar.get(i, 0));
        }
        else
        {
            // no snap
        }
        dual_feas.set(i, 0, dual_feas.get(i, 0) - X.get(i, 0));
    }

    for (int i = 0 ; i < _n ; i++)
    {
        lagGradX.set(i, 0, lagGradX.get(i, 0) + mu / (X.get(i, 0) - lvar.get(i, 0)) - mu / (uvar.get(i, 0) - X.get(i, 0)));
    }

    // |-Sy - mu|
    double lagGradS = 0;
    for (int i=0; i < _nbCons; i++)
    {
        lagGradS += std::pow(-XS.get(i + _n, 0) * lambda.get(i, 0) - mu, 2);
    }
    lagGradS = std::sqrt(lagGradS);

/*
    std::cout << " |G_s|=" << lagGradS;
    std::cout << " |c + s|=" << cslack.norm();
    std::cout << " |G_x|=" << dual_feas.norm() << " |newG_x|=" << lagGradX.norm() << std::endl;
    // XS.display(std::cout); lvar.display(std::cout); uvar.display(std::cout); lambda.display(std::cout);
*/    
    return std::max(lagGradS, std::max(cslack.norm(), dual_feas.norm()));
}

size_t NOMAD::QPSolverOptimize::solveLM(
    NOMAD::Point & X,
    SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    SGTELIB::Matrix & cX,
    const double Fx,
    const double tol,
    const size_t maxIterInner,
    const double tolDistDXInner,
    const bool checkStrict, // if true, it checks wether the point is strictly feasible.
    // const double Delta0, // Initial trust-region radius
    const bool verbose)
{

    const int nbVar = _n + _nbCons; // we use slack variables here

    // Initialize X
    Point Xp(X);
    SGTELIB::Matrix XSp = XS;
    checkStrict && check_strict_feasible(XSp, lvar, uvar);
    SGTELIB::Matrix Xcan("Xcan", nbVar, 1);

    // Check optimality:
    SGTELIB::Matrix cons("cx", _nbCons, 1);
    getModelCons(&cons, X);
    SGTELIB::Matrix cslack("cx+s", _nbCons, 1);
    for (int j=0; j < _nbCons; ++j)
    {
        cslack.set(j, 0, XSp.get(j + _n, 0) + cons.get(j, 0));
    }
    SGTELIB::Matrix Jx = getModelJacobian(X);

    // Pre-allocation for the normal step
    SGTELIB::Matrix W("W", _nbCons, _nbCons + _n);
    SGTELIB::Matrix WtW("W", _nbCons + _n, _nbCons + _n);
    SGTELIB::Matrix WtWr("WtW*r", _nbCons + _n, 1);
    SGTELIB::Matrix wq("vxs", _nbCons + _n, 1);
    SGTELIB::Matrix vxs("vxs", _nbCons + _n, 1);
    SGTELIB::Matrix zer("zer", _nbCons + _n, 1);
    zer.fill(0.0);
    SGTELIB::Matrix temp("-temp-", _nbCons + _n, 1);
    SGTELIB::Matrix r("r", _nbCons, 1);
    SGTELIB::Matrix cxp("cxp", _nbCons, 1);
    SGTELIB::Matrix checkslack("cxp+sp", _nbCons, 1);
    double f_normal_model;
    double backtrack_length;
    double vsi, vxi, xi, li, ui;

    // Trust-region should handle it, but...
    double normal_step_regularization = tol; // regularizer of the normal equation in normal step

    // W is the Jacobian of the residual
    for (int i = 0; i < _nbCons; i++)
    {
        for (int j = 0; j < _n; j++)
        {
            W.set(i, j, Jx.get(i, j));
        }
        for (int j = 0; j < _nbCons; j++)
        {
            if (i == j)
            {
                W.set(i, j + _n, 1.0);
            }
            else
            {
                W.set(i, j + _n, 0.0);
            }
        }
    }
    SGTELIB::Matrix::inplace_product(WtW, W.transpose(), W);
    for (int i = 0; i < nbVar; i++)
    {
        WtW.set(i, i, WtW.get(i, i) + normal_step_regularization); // regularization term
    }
    SGTELIB::Matrix::inplace_product(wq, W.transpose(), cslack);

    // Trust-region parameters:
    double ared, pred;
    const double epsilon_1 = 1E-8; // trust-region successful ratio
    const double epsilon_2 = 0.9; // trust-region very successful ratio
    const double gamma_1 = 0.5; // trust-region decrease factor
    const double gamma_2 = 2; // trust-region increase factor

    double Delta = 10000; // trust-region initial radius
    double smallestDelta = 1E-15;
    double largestDelta = 1E15;

    double tau = 0.5;

    double res = cslack.norm();
    double res0 = res;
    bool success = (res <= Fx) || (wq.norm() <= tol);

    size_t iterInnerLoop = 0;
    size_t successivUnsuccessful = 0;
    double distXInnerLoop = INF;
    bool failure = false;
    failure = failure || (distXInnerLoop <= tolDistDXInner);
    failure = failure || (iterInnerLoop >= maxIterInner);

    _verbose && std::cout << " Feas. " << iterInnerLoop << " |c+s|=" << res;
    _verbose && std::cout << " ? " << tol << " D=" << Delta << std::endl;

    while (!failure && !success)
    {

        // Compute the normal step (vx, vs):
        solve_TR_constrained_QP(&vxs, zer, WtW, wq, temp, Delta);
        f_normal_model = getModelObj(vxs, WtW, wq, 0.5 * std::pow(cslack.norm(), 2));

        // backtrack to satisfy vs >= -tau
        backtrack_length = 1;
        for (int i = 0; i < _nbCons; i++)
        {
            vsi = vxs.get(_n + i, 0);
            if (vsi < -tau)
            {
                backtrack_length = std::min(backtrack_length, - tau / vsi);
            }
        }
        // backtrack to satisfy vx + tau * x >= tau * l and tau * u >= vx + tau * x
        for (int i = 0; i < _n; i++)
        {
            vxi = vxs.get(i, 0);
            xi = XSp.get(i, 0);
            li = lvar.get(i, 0);
            ui = uvar.get(i, 0);
            if (vxi == 0)
            {
                if (xi < li)
                {
                    std::cout << " lvar issue: " << std::fabs(li - xi) << std::endl;
                }
                if (xi > ui)
                {
                    std::cout << " uvar issue: " << std::fabs(ui - xi) << std::endl;
                }
            }
            else
            {
                if (vxi < tau * (li - xi))
                {
                    backtrack_length = std::min(backtrack_length, tau * (li - xi) / vxi);
                }
                if (vxi > tau * (ui - xi))
                {
                    backtrack_length = std::min(backtrack_length, tau * (ui - xi) / vxi);
                }
            }
        }
        vxs.multiply(backtrack_length);

        // Update the candidate
        for (int i = 0; i < _n; i++)
        {
            Xcan.set(i, 0, XSp.get(i, 0) + vxs.get(i, 0));
            Xp[i] = X[i] + vxs.get(i, 0);
        }
        for (int i = 0; i < _nbCons; i++)
        {
            Xcan.set(i + _n, 0, XSp.get(i + _n, 0) + vxs.get(i + _n, 0));
        }

        // Make sure the bounds are satisfied strictly
        checkStrict && check_strict_feasible(Xcan, lvar, uvar);

        if (f_normal_model < 0)
        {
            _verbose && std::cout << " solver normal step: |v|=" << vxs.norm() << " f(v)=" << f_normal_model << " |c(xv) + sv|=" << checkslack.norm() << " |c(x) + s|=" << cslack.norm() << std::endl;
            // throw NOMAD::Exception(__FILE__, __LINE__, " NOT POSSIBLE");
        }

        // Compute the residual: r = Jx * vx + vs + (cx + s)
        SGTELIB::Matrix::inplace_product(r, W, vxs);
        r.add(cslack);
        SGTELIB::Matrix::inplace_product(WtWr, W.transpose(), r);

        getModelCons(&cxp, Xp);
        for (int j=0; j < _nbCons; ++j)
        {
            checkslack.set(j, 0, Xcan.get(j + _n, 0) + cxp.get(j, 0));
        }

        ared = cslack.norm() - checkslack.norm(); // > 0 if that works
        pred = cslack.norm() - r.norm();

        if ((ared >= pred * epsilon_1) && (pred > 0))
        { // r >= epsilon_1
            // Accept x and s steps
            XSp = Xcan;
            distXInnerLoop = NOMAD::Point::dist(X, Xp).todouble();
            X = Xp;

            // Increase Delta
            if (ared >= pred * epsilon_2) // r >= epsilon_2
            {
                Delta = std::min(gamma_2 * Delta, largestDelta); // std::max(d.norm(), delta);
            }

            successivUnsuccessful = 0;

            // Check optimality:
            cons = cxp;
            cslack = checkslack;
            Jx = getModelJacobian(X);
            for (int i = 0; i < _nbCons; i++)
            {
                for (int j = 0; j < _n; j++)
                {
                    W.set(i, j, Jx.get(i, j));
                }
            }
            SGTELIB::Matrix::inplace_product(WtW, W.transpose(), W);
            for (int i = 0; i < nbVar; i++)
            {
                WtW.set(i, i, WtW.get(i, i) + normal_step_regularization); // regularization term
            }
            SGTELIB::Matrix::inplace_product(wq, W.transpose(), cslack);
        }
        else
        {
            if (pred < 0)
            {
                // Oh...
            }
            // Decrease Delta
            Delta = std::max(gamma_1 * std::min(Delta, vxs.norm()), smallestDelta);
            successivUnsuccessful += 1;
        }

        // Check optimality
        res = cslack.norm();
        success = (res <= Fx) || (WtWr.norm() <= tol);

        iterInnerLoop += 1;
        failure = failure || (distXInnerLoop <= tolDistDXInner);
        failure = failure || (iterInnerLoop >= maxIterInner);
        failure = failure || (backtrack_length == 0);

        _verbose && std::cout << " Feas. " << iterInnerLoop << " |c+s|=" << res;
        _verbose && std::cout << " ? " << tol << " D=" << Delta << " |v|=" << vxs.norm();
        _verbose && std::cout << " ared=" << ared << " pred=" << pred << " r=" << ared/pred;
        _verbose && std::cout << " back=" << backtrack_length << std::endl;
    }

    size_t resultFeasibility;
    if (success)
    {
        resultFeasibility = 1; // success
        XS = XSp;
        cX = cons;
    }
    else if (res < res0)
    {
        resultFeasibility = 2; // improved
        XS = XSp;
        cX = cons;
    }
    else
    {
        resultFeasibility = 0; // failed
        for (int i = 0; i< _n; i++)
        {
            X[i] = XS.get(i, 0);
        }
    }

    return resultFeasibility;
}

int NOMAD::QPSolverOptimize::solver_barrier(
    NOMAD::Point & X,
    SGTELIB::Matrix & XSp,
    SGTELIB::Matrix & p,
    SGTELIB::Matrix & cslack,
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    SGTELIB::Matrix & lambda,
    SGTELIB::Matrix & Gx,
    SGTELIB::Matrix & cons,
    SGTELIB::Matrix & Jx,
    const double mu,
    const double tol_mu,
    const size_t maxIterInner,
    const double tolDistDXInner,
    const bool verbose,
    const bool verbosePCG)
{

    const int nbVar = _n + _nbCons; // we use slack variables here

    // Initialize X
    XSp = XS;
    check_strict_feasible(XSp, lvar, uvar);

    SGTELIB::Matrix Xcan("Xcan", nbVar, 1);
    Point Xp(X);
    SGTELIB::Matrix r("r", _nbCons, 1);

    // For the normal step:
    SGTELIB::Matrix W("W", _nbCons, _nbCons + _n);
    SGTELIB::Matrix Wscal("W", _nbCons, _nbCons + _n);
    SGTELIB::Matrix WtW("W", _nbCons + _n, _nbCons + _n);
    SGTELIB::Matrix wq("wq", _nbCons + _n, 1);
    SGTELIB::Matrix vxs("vxs", _nbCons + _n, 1);
    SGTELIB::Matrix zer("zer", _nbCons + _n, 1);
    SGTELIB::Matrix x0PCG("x0PCG", _nbCons + _n, 1);
    SGTELIB::Matrix r0("r0", _nbCons, 1);
    r0.fill(0.0);
    zer.fill(0.0);
    SGTELIB::Matrix temp("-temp-", _nbCons + _n, 1);
    SGTELIB::Matrix vs("vs", _nbCons, 1);
    SGTELIB::Matrix vx("vx", _n, 1);
    SGTELIB::Matrix Jxvx("Jxvx", _nbCons, 1);

    // Compute p
    double np = p.norm();
    SGTELIB::Matrix Q("Q", _nbCons + _n, _nbCons + _n);
    Q.fill(0);
    SGTELIB::Matrix qc("qc", _nbCons + _n, 1);
    SGTELIB::Matrix ps("ps", _nbCons, 1);
    SGTELIB::Matrix px("px", _n, 1);
    double psi, si, pxi, vxi, xi, li, ui, vsi;

    ////////////////////////////////////////////////////////////
    // Barrier solver parameters:
    double ared, pred;
    const double epsilon_1 = 1E-8; // trust-region successful ratio
    const double epsilon_2 = 0.9; // trust-region very successful ratio
    const double gamma_1 = 0.5; // trust-region decrease factor
    const double gamma_2 = 2; // trust-region increase factor

    double Delta = 1; // trust-region initial radius
    double smallestDelta = 1E-15;
    double largestDelta = 1E15;

    const double tau = 0.995;
    double rho = 0.5;

    double normal_step_regularization = 1E-7; // regularizer of the normal equation in normal step
    double Delta_normal_step_factor = 0.8; // Factor of Delta used in normal step
    double min_tol_PCG = 1E-10;
    double tol_PCG;
    bool success_PCG;
    double den, nu;
    double small_p = 1E-15; // below this value `p` is considered 0 (declare success).
    ////////////////////////////////////////////////////////////

    size_t iterInnerLoop = 0;
    double distXInnerLoop = INF;
    size_t successivUnsuccessful = 0;

    bool innerSuccess, innerFailure;
    compute_slack_multiplier(lambda, XSp, Jx, Gx, mu);
    double res = errorTRIPM(XSp, lvar, uvar, lambda, cslack, mu);
    innerSuccess = (res <= tol_mu);
    bool success = innerSuccess;
    innerFailure = false;
    innerFailure = innerFailure || (distXInnerLoop <= tolDistDXInner);
    innerFailure = innerFailure || (iterInnerLoop >= maxIterInner);

    verbose && std::cout << "Inner 0 (-): ";
    verbose && std::cout << " |E(x,s,y; mu)|= " << res << " mu=" << mu;
    verbose && std::cout << " |c+s|=" << cslack.norm();
    verbose && std::cout << " mu=" << mu << " e_mu=" << tol_mu << std::endl;

    while (!innerFailure && !innerSuccess)
    {
        // Compute the normal step (vx, vs):
        for (int i = 0; i < _nbCons; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                W.set(i, j, Jx.get(i, j));
            }
            for (int j = 0; j < _nbCons; j++)
            {
                if (i == j)
                {
                    W.set(i, j + _n, XSp.get(_n + i, 0));
                }
                else
                {
                    W.set(i, j + _n, 0.0);
                }
            }
        }
        SGTELIB::Matrix::inplace_product(WtW, W.transpose(), W);
        for (int i = 0; i < nbVar; i++)
        {
            WtW.set(i, i, WtW.get(i, i) + normal_step_regularization); // regularization term
        }
        SGTELIB::Matrix::inplace_product(wq, W.transpose(), cslack);
        vxs.fill(0.0);
        solve_TR_constrained_QP(&vxs, zer, WtW, wq, temp, Delta_normal_step_factor * Delta);
        double f_normal_model = getModelObj(vxs, WtW, wq, 0.5 * std::pow(cslack.norm(), 2));

        ///////////////////////////////////////////////////////////////        
        for (int i = 0; i< _n; i++)
        {
            Xcan.set(i, 0, XSp.get(i, 0) + vxs.get(i, 0));
            Xp[i] = X[i] + vxs.get(i, 0);
        }
        for (int i = 0; i < _nbCons; i++)
        {
            Xcan.set(i + _n, 0, XSp.get(i + _n, 0) + vxs.get(i + _n, 0));
        }
        SGTELIB::Matrix checkcons = getModelCons(Xp);
        SGTELIB::Matrix checkslack("checkslack", _nbCons, 1);
        for (int j=0; j < _nbCons; ++j)
        {
            checkslack.set(j, 0, XSp.get(j + _n, 0) + checkcons.get(j, 0));
        }

        if (f_normal_model < 0)
        {
                    std::cout << " solver normal step: |v|=" << vxs.norm() << " f(v)=" << f_normal_model << " |c(xv) + sv|=" << checkslack.norm() << " |c(x) + s|=" << cslack.norm() << std::endl;
        //            WtW.display(std::cout);
        //            wq.display(std::cout);
        //            cslack.display(std::cout);
        //            vxs.display(std::cout);
                    // throw NOMAD::Exception(__FILE__, __LINE__, " NOT POSSIBLE");
        }
        ///////////////////////////////////////////////////////////////
        

        // backtrack to satisfy vs >= -tau / 2
        double backtrack_length = 1;
        for (int i = 0; i < _nbCons; i++)
        {
            vsi = vxs.get(_n + i, 0);
            if (vsi < -tau / 2)
            {
                backtrack_length = std::min(backtrack_length, - tau / (2 * vsi));
            }
        }
        // backtrack to satisfy vx + tau/2 * x >= tau/2 * l and tau/2 * u >= vx + tau/2 * x
        for (int i = 0; i < _n; i++)
        {
            vxi = vxs.get(i, 0);
            xi = XSp.get(i, 0);
            li = lvar.get(i, 0);
            ui = uvar.get(i, 0);
            if (vxi == 0)
            {
                if (xi < li)
                {
                    std::cout << " lvar issue: " << std::fabs(li - xi) << std::endl;
                }
                if (xi > ui)
                {
                    std::cout << " uvar issue: " << std::fabs(ui - xi) << std::endl;
                }
            }
            else
            {
                if (vxi < tau * (li - xi))
                {
                    backtrack_length = std::min(backtrack_length, tau * (li - xi) / (2 * vxi));
                }
                if (vxi > tau * (ui - xi))
                {
                    backtrack_length = std::min(backtrack_length, tau * (ui - xi) / (2 * vxi));
                }
            }
        }
        vxs.multiply(backtrack_length);

        // Compute(px, ptildes) with projected CG:
        SGTELIB::Matrix HLag = getModelLagHessian(X, lambda);
        for (int i = 0; i < _n; i++)
        {
            xi = XSp.get(i, 0);
            li = lvar.get(i, 0);
            ui = uvar.get(i, 0);
            for (int j = 0; j < _n; j++)
            {
                if (i == j)
                {
                    Q.set(i, j, HLag.get(i, j) + mu *(1 / std::pow(xi - li, 2) + 1 / std::pow(ui - xi, 2)));
                }
                else
                {
                    Q.set(i, j, HLag.get(i, j));
                }
            }
            qc.set(i, 0, Gx.get(i, 0) - mu / (xi - li) + mu / (ui - xi));
        }
        for (int i = 0; i < _nbCons; i++)
        {
            double si = XSp.get(_n + i, 0);
            Q.set(i + _n, i + _n, lambda.get(i, 0) * si);
            qc.set(i + _n, 0, -mu);
        }
        
        for (int i = 0; i < _nbCons; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                Wscal.set(i, j, Jx.get(i, j));
            }
            for (int j = 0; j < _nbCons; j++)
            {
                if (i == j)
                {
                    Wscal.set(i, j + _n, 1);
                }
                else
                {
                    Wscal.set(i, j + _n, 0.0);
                }
            }
        }

        SGTELIB::Matrix::inplace_product(r, W, vxs);
        x0PCG.fill(0.0);
        tol_PCG = std::min(tol_mu, min_tol_PCG);
        success_PCG = projected_conjugate_gradient(p, W, r, Q, qc, Delta, x0PCG, tol_PCG, verbosePCG);
        
        if (!success_PCG)
        {
            // std::cout << "Check PCG |Wp - r|=" << SGTELIB::Matrix::sub(SGTELIB::Matrix::product(W, p), r).norm() << std::endl;
            // std::cout << "Check PCG |Qp + q|=" << SGTELIB::Matrix::add(SGTELIB::Matrix::product(Q, p), qc).norm() << std::endl;
            // p.display(std::cout);
            // SGTELIB::Matrix::add(SGTELIB::Matrix::product(Q, p), qc).display(std::cout);
            // W.display(std::cout);
            std::cout << " PCG failed, we should call off. |x0|=" << x0PCG.norm() << " d=" << Delta << std::endl;
        }
        
        // backtrack to satisfy ptildes >= -tau
        if (p.norm() > Delta)
        {
            backtrack_length = Delta / p.norm();
        }
        else
        {
            backtrack_length = 1;
        }
        for (int i = 0; i < _nbCons; i++)
        {
            psi = p.get(_n + i, 0);
            if (psi < -tau)
            {
                backtrack_length = std::min(backtrack_length, - tau / psi);
            }
        }
        // backtrack to satisfy px + tau * x >= tau * l and tau * u >= px + tau * x
        for (int i = 0; i < _n; i++)
        {
            pxi = p.get(i, 0);
            xi = XSp.get(i, 0);
            li = lvar.get(i, 0);
            ui = uvar.get(i, 0);
            // std::cout << " (n=" << n << ") i=" << i << " li=" << li << " ui=" << ui << " xi=" << xi << " pxi + xi=" << pxi + xi;
            if (pxi == 0)
            {
                if (xi < li)
                {
                    std::cout << " lvar issue: " << std::fabs(li - xi) << std::endl;
                }
                if (xi > ui)
                {
                    std::cout << " uvar issue: " << std::fabs(ui - xi) << std::endl;
                }
            }
            else
            {
                if (pxi < tau * (li - xi))
                {
                    backtrack_length = std::min(backtrack_length, tau * (li - xi) / pxi);
                }
                if (pxi > tau * (ui - xi))
                {
                    backtrack_length = std::min(backtrack_length, tau * (ui - xi) / pxi);
                }
            }
            //std::cout << " b=" << backtrack_length << std::endl;
        }
        if (backtrack_length <= 0)
        {
            std::cout << " backtrack length=" << backtrack_length << " tau=" << tau << std::endl;
        }
        p.multiply(backtrack_length);

        // ps = S ptildes:
        for (int i = 0; i< _nbCons; i++)
        {
            si = XSp.get(_n + i, 0);
            p.set(i + _n, 0, p.get(i + _n, 0) * si);
        }
        np = p.norm();
        // std::cout << "Wscalp - Wv =" << SGTELIB::Matrix::sub(r, SGTELIB::Matrix::product(W, vxs)).norm() << std::endl; // Should be 0
    
        // Compute the residual: r = Jx * px + S * ps + (cx + s)
        SGTELIB::Matrix::inplace_product(r, W, p);
        r.add(cslack);
        double nr = r.norm();

        // Update nu_l:
        den = cslack.norm() - nr + tol_mu;
        if (den == 0)
        {
            nu = 0;
        }
        else if (den < 0)
        {
            // We try a feasibility step only, in this case
            verbose && std::cout << " p do not increase feasibility" << den <<", we try using v: ";
            nr = checkslack.norm();
            den = cslack.norm() - checkslack.norm();
            verbose && std::cout << den << std::endl;
            if (den <= 0)
            {
                verbose && std::cout << " normal step: |v|=" << vxs.norm() << " f(v)=" << f_normal_model << " |c(xv) + sv|=" << checkslack.norm() << " |c(x) + s|=" << cslack.norm();
                verbose && std::cout << "Wscalp - Wv=" << SGTELIB::Matrix::sub(SGTELIB::Matrix::product(Wscal, p), SGTELIB::Matrix::product(W, vxs)).norm();
                verbose && std::cout << " Wp= " << SGTELIB::Matrix::product(Wscal, p).norm() << std::endl;
                nu = INF;
            }
            else
            { // den > 0
                p = vxs;
                nu = 1E15;
            }
        }
        else
        {
            nu = getModelObj(p, Q, qc) / ((1 - rho) * den);
            // nu = std::max(nu, -nu); // std::max(nu, smallest_nu);
        }

        // Trust-region update
        for (int i = 0; i< _n; i++)
        {
            Xcan.set(i, 0, XSp.get(i, 0) + p.get(i, 0));
            Xp[i] = X[i] + p.get(i, 0);
        }
        for (int i = 0; i < _nbCons; i++)
        {
            Xcan.set(i + _n, 0, XSp.get(i + _n, 0) + p.get(i + _n, 0));
        }
        
        // Compute ared (> 0):
        ared = merit_function_barrier(X, XSp, lvar, uvar, mu, nu);
        ared -= merit_function_barrier(Xp, Xcan, lvar, uvar, mu, nu);

        // Compute pred (> 0):
        pred = nu * den;
        pred -= getModelObj(p, Q, qc);

        if ((ared >= pred * epsilon_1) && (pred > 0))
        { // r >= epsilon_1
            // Accept x and s steps
            XSp = Xcan;
            distXInnerLoop = NOMAD::Point::dist(X, Xp).todouble();
            X = Xp;

            // Increase Delta
            if (ared >= pred * epsilon_2) // r >= epsilon_2
            {
                Delta = std::min(gamma_2 * Delta, std::max(1 / mu, largestDelta)); // std::max(d.norm(), delta);
            }

            success = true;
            successivUnsuccessful = 0;

            // Check optimality:
            getModelCons(&cons, X);
            for (int j=0; j < _nbCons; ++j)
            {
                cslack.set(j, 0, XSp.get(j + _n, 0) + cons.get(j, 0));
            }
            getModelGrad(&Gx, X);
            Jx = getModelJacobian(X);
            compute_slack_multiplier(lambda, XSp, Jx, Gx, mu);
            res = errorTRIPM(XSp, lvar, uvar, lambda, cslack, mu);
        }
        else
        {
            if (pred < 0)
            {
                verbose && std::cout << "quad=" << getModelObj(p, Q, qc) << " <= " << getModelObj(zer, Q, qc);
                verbose && std::cout << " |Wscalp|=" << SGTELIB::Matrix::product(Wscal, p).norm() << " |Wp|=" << SGTELIB::Matrix::product(W, p).norm() << " |Wvxs|=" << SGTELIB::Matrix::product(W, vxs).norm();
                verbose && std::cout << " |r|=" << nr << " |rv|=" << SGTELIB::Matrix::add(SGTELIB::Matrix::product(W, vxs), cslack).norm();
                verbose && std::cout << " |c+s|=" << cslack.norm() << " nu=" << nu << std::endl;
                verbose && std::cout << "pred=" << nu * cslack.norm();
                verbose && std::cout << " + " << - getModelObj(p, Q, qc);
                verbose && std::cout << " + " << - nu * nr << std::endl;
            }
            // Decrease Delta
            Delta = std::max(gamma_1 * std::min(Delta, np), smallestDelta);
            successivUnsuccessful += 1;
        }
        iterInnerLoop += 1;

        verbose && std::cout << "Inner " << iterInnerLoop << " (" << successivUnsuccessful << "):";
        verbose && std::cout << " |E(x,s,y;mu)|=" << res << " |r|=" << nr << " |c+s|=" << cslack.norm();
        verbose && std::cout << " nu=" << nu;
        if (pred != 0)
        {
            verbose && std::cout << " ared/pred=" << ared / pred;
        } else
        {
            verbose && std::cout << " ared/pred=" << ared << "/" << pred;
        }
        verbose && std::cout << " Delta=" << Delta << " |p|=" << np << " |v|=" << vxs.norm() << std::endl;

        innerSuccess = (res <= tol_mu) || (np <= small_p);
        innerFailure = innerFailure || (distXInnerLoop <= tolDistDXInner);
        innerFailure = innerFailure || (iterInnerLoop >= maxIterInner);
    }

    int resultSolverBarrier;
    if (res <= tol_mu)
    {
        resultSolverBarrier = 1;
    }
    else if (np <= small_p)
    {
        resultSolverBarrier = 2;
    }
    else if (success)
    {
        // at least one step has been made
        resultSolverBarrier = 3;
    }
    else
    {
        resultSolverBarrier = 0;
    }
    return resultSolverBarrier;
}

double NOMAD::QPSolverOptimize::merit_function_barrier(
    NOMAD::Point & X,
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    const double mu,
    const double nu )
{

    double fx = getModelObj(X);

    check_strict_feasible(XS, lvar, uvar);

    double xi, si, li, ui;

    double res = 0;
    for (int i = 0; i < _nbCons; i++)
    {
        si = XS.get(i + _n, 0);
        res -= mu * std::log(si);
    }
    
    for (int i = 0; i < _n; i++)
    {
        xi = XS.get(i, 0);
        ui = uvar.get(i, 0);
        li = lvar.get(i, 0);
        res -= mu * std::log(xi - li);
        res -= mu * std::log(ui - xi);
    }

    SGTELIB::Matrix cons("cons", _nbCons, 1);
    SGTELIB::Matrix cslack("cslack", _nbCons, 1);

    getModelCons(&cons, X);
    for (int j=0; j < _nbCons; ++j)
    {
        cslack.set(j, 0, XS.get(j + _n, 0) + cons.get(j, 0));
    }

    return fx + res + nu * cslack.norm();
}

bool NOMAD::QPSolverOptimize::check_strict_feasible(
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar
)
{
    double xi, li, ui;
    bool strict_feasible = true;
    for (int i = 0; i < _n; i++)
    {
        xi = X.get(i, 0);
        ui = uvar.get(i, 0);
        li = lvar.get(i, 0);
        strict_feasible = strict_feasible && (xi > li) && (xi < ui);
    }

    if (!strict_feasible)
    {
        X.display(std::cout);
        lvar.display(std::cout);
        uvar.display(std::cout);
        throw NOMAD::Exception(__FILE__, __LINE__, X.get_name() + " is not strictly feasible.");
    }

    return strict_feasible;
}

//*****************************************************************************
//
// Optimality functions
//
//*****************************************************************************

// Check that ||Proj(x - ∇f) - x|| <= tol
double NOMAD::QPSolverOptimize::check_optimality_bounds(
    NOMAD::Point & X_k,
    const SGTELIB::Matrix & gradientF_k)
{

    NOMAD::Point temp(X_k);
    for (int i = 0 ; i < _n ; i++)
    {
        temp[i] = X_k[i] - gradientF_k.get(i, 0);
    }
    temp.snapToBounds(_modelLowerBound, _modelUpperBound);
    temp = temp - X_k;

    Double norm_temp = 0;
    for (int i = 0 ; i < _n ; i++)
    {
        norm_temp += temp[i].pow2();
    }
    
    return norm_temp.sqrt().todouble();
}

double NOMAD::QPSolverOptimize::check_optimality_bounds(
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & gradient,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar)
{

    int n = X.get_nb_rows();
    SGTELIB::Matrix dual_feas = SGTELIB::Matrix("dual_feas", n, 1);
    return check_optimality_bounds(X, gradient, lvar, uvar, dual_feas);
}

double NOMAD::QPSolverOptimize::check_optimality_bounds(
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & gradient,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    SGTELIB::Matrix & dual_feas)
{

    int n = X.get_nb_rows();
    if (lvar.get_nb_rows() != n || uvar.get_nb_rows() != n || gradient.get_nb_rows() != n)
    {
        std::string err = "check_optimality_bounds: ";
        err += "Inconsistent dimension for bounds. Expecting ";
        err += std::to_string(n);
        err += " but sizes are " + std::to_string(lvar.get_nb_rows());
        err += " and " + std::to_string(uvar.get_nb_rows()) + ".";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (int i = 0 ; i < n ; i++)
    {
        dual_feas.set(i, 0, X.get(i, 0) - gradient.get(i, 0));
    }
    snapToBounds(dual_feas, lvar, uvar);
    dual_feas.sub(X);
    
    return dual_feas.norm();
}

double NOMAD::QPSolverOptimize::compute_dual_residual(
    const SGTELIB::Matrix & Grad_k,
    const SGTELIB::Matrix & Jacobian_k,
    const SGTELIB::Matrix & multiplier_k) const
{
    const int nvar = Jacobian_k.get_nb_cols();
    const int ncon = Jacobian_k.get_nb_rows();

    lencheck(nvar, Grad_k);
    lencheck(nvar, multiplier_k);
    sizecheck(ncon, nvar, Jacobian_k);

    if (Jacobian_k.has_nan()) {
        throw NOMAD::Exception(__FILE__, __LINE__, "Jacobian_k contains NaN");
    }

    // SGTELIB::Matrix residual("residual", static_cast<int>(nvar), 1);
    SGTELIB::Matrix residual = SGTELIB::Matrix::product(Jacobian_k.transpose(), multiplier_k);
    residual.sub(Grad_k);

    return residual.norm();
}

//*****************************************************************************
//
// Basic optimization functions
//
//*****************************************************************************

bool NOMAD::QPSolverOptimize::Convex_TR_QP(
    SGTELIB::Matrix * d,
    const SGTELIB::Matrix & g, // gradient
    const SGTELIB::Matrix & gW, // gradient active
    const SGTELIB::Matrix & H, // hessian
    const SGTELIB::Matrix & HW, // hessian active
    int * pp,
    double ** D,
    double ** L,
    const bool * active, // length n
    const double Delta,
    const bool verbose
)
{
    bool success = true;

    int n = g.get_nb_rows();
    int nfree = n - sum(active, n);

    lencheck(n, *d);
    sizecheck(n, n, H);
    sizecheck(nfree, nfree, HW);
    lencheck(n, g);
    lencheck(nfree, gW);

    double * sol = new double [nfree];

    bool solve_success = ComputeNewtonDirection(gW, pp, D, L, sol, nfree);
    if (!solve_success)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt solve");
    }

    d->fill(0);
    int ki = 0;
    for (int i = 0 ; i < n ; ++i )
    {
        if (!active[i])
        {
            d->set(i, 0, sol[ki]);
            ki += 1;
        }
    }
    if (ki != nfree)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error dimension");
    }

    double slope = SGTELIB::Matrix::dot(g, *d);
    if (slope > 0)
    {
        verbose && std::cout << "Numerical issue Newton direction is not positive definite, slope= " << slope << std::endl;
    }

    double nd = d->norm();

    // Check if we solve a "fake" trust-region problem.
    if ((Delta < 1E15) && (nd > Delta))
    {
        verbose && std::cout << " Newton direction is not inside the trust-region: " << nd << " >= " << Delta << std::endl;
        d->multiply(Delta / nd);
    }
    verbose && std::cout << "|d|= " << nd << " slope= " << slope;

    delete [] sol;

    return success;
}

bool NOMAD::QPSolverOptimize::ComputeNewtonDirection(
    const SGTELIB::Matrix g, 
    int * pp,
    double ** D,
    double ** L, // const SGTELIB::Matrix H,
    double * sol,
    int n
)
{
    lencheck(n, g);

    bool success = true;
    string error_msg;

    double * rhs = new double [n];

    for (int i = 0 ; i < n ; ++i )
    {
        rhs[i] = - g.get(i, 0);
        sol[i] = 0;
    }

    success = NOMAD::ldl_solve(error_msg, D, L, rhs, sol, pp, n);

    delete [] rhs;

    return success;
}

bool NOMAD::QPSolverOptimize::InverseIteration(
    SGTELIB::Matrix * sol,
    const SGTELIB::Matrix & HW,
    const double eigmin,
    const int nfree,
    const double tol, // > 0
    const bool verbose
)
{

    lencheck(nfree, *sol);
    sizecheck(nfree, nfree, HW);

    // We compute the eigenvector corresponding to it
    SGTELIB::Matrix bk("bk", nfree, 1); // depend on X
    bk.fill(1.0 / nfree);

    SGTELIB::Matrix bkp("bkp", nfree, 1);
    bkp.fill(0);

    double Ck = 1.0, Ckp;
    double fix_point;
    SGTELIB::Matrix HWp = HW;
    for (int i = 0 ; i < nfree ; ++i )
    {
        HWp.set(i, i, HWp.get(i, i) - eigmin + 1E-7);
    }

    SGTELIB::Matrix invHWp = HWp.SVD_inverse();
    if (invHWp.has_nan())
    {
        return false;
    }

    SGTELIB::Matrix invHWpbk = SGTELIB::Matrix::product(invHWp, bk);

    bool OK = false;
    size_t count = 0;
    const size_t max_count = 1000;
    while (!OK && count < max_count )
    {
        // b_{k+1} = H^{-1} * b_k * C_k
        SGTELIB::Matrix::inplace_product(bkp, invHWpbk, Ck);
        if (bkp.has_nan())
        {
            return false;
        }

        fix_point = SGTELIB::Matrix::sub(bkp, bk).norm();

        bk = bkp;
        SGTELIB::Matrix::inplace_product(invHWpbk, invHWp, bk);

        verbose && std::cout << fix_point << " Ck=" << Ck;
        verbose && std::cout << " |bk|=" << bk.norm() << " |bkp|=" << bkp.norm() << std::endl;

        if (invHWpbk.norm() <= 0.0) {
            return false;
        }

        Ckp = 1 / invHWpbk.norm();

        // OK = (fix_point <= tol) || (fabs(Ck - Ckp) <= tol);
        OK = (fix_point <= 1E-7) || (fabs(Ck - Ckp) <= tol);
        Ck = Ckp;

        if (bk.has_nan() || (invHWpbk.norm() == 0))
        {
            return false;
        }
        count++;
    }

    // sol = bk;
    for (int k = 0; k < nfree; ++k)
    {
        sol->set(k, 0, bk.get(k, 0));
    }

    return true;
}

bool NOMAD::QPSolverOptimize::projected_conjugate_gradient (
    SGTELIB::Matrix & x ,
    const SGTELIB::Matrix & A ,
    const SGTELIB::Matrix & b ,
    const SGTELIB::Matrix & G ,
    const SGTELIB::Matrix & c ,
    const double delta ,
    const SGTELIB::Matrix & x0 ,
    const double tol,
    const bool verbose)
{
    const int n = x0.get_nb_rows();
    const int m = b.get_nb_rows();
    // m should be less than n
    const int npm = n + m;

    sizecheck(m, n, A);
    lencheck(m, b);
    sizecheck(n, n, G);
    lencheck(n, c);

    // other choices are possible for H
    SGTELIB::Matrix H = SGTELIB::Matrix::identity(n);

    ////////// LDLt
    // init matrices for LDLt
    double ** M = new double *[npm];
    double ** L = new double *[npm];
    double ** D = new double *[npm];
    for (int i = 0 ; i < npm ; ++i )
    {
        M[i] = new double[npm];
        L[i] = new double[npm];
        D[i] = new double[npm];
        for (int j = 0 ; j < npm ; ++j )
        {
            if ((i < n) & (j < n)) {
                M[i][j] = H.get(i, j);
            } else if ((i < n)) { // j >= n
                M[i][j] = A.get(j - n, i);
            } else if (j < n) { // i >= n
                M[i][j] = A.get(i - n, j);
            } else { // i >= n, j >= n
                M[i][j] = 0;
            }
            L[i][j] = 0;
            D[i][j] = 0;
        }
    }
    int * pp = new int [npm];
    for (int i = 0 ; i < npm ; ++i )
    {
        pp[i] = 0;
    }
    std::string error_msg;

    bool success = NOMAD::LDLt_decomposition(error_msg, M, L, D, pp, npm);
    if (!success)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt decomposition");
    }

    double * rhs = new double [npm];
    double * sol = new double [npm];

    // Initial guess
    x = x0;
    
    // Just for information:
    SGTELIB::Matrix Ax = SGTELIB::Matrix::product(A, x);
    SGTELIB::Matrix f = b;
    f.sub(Ax);
    double feas = f.norm();

    // If x0 is not feasible:
    if ((feas > 1E-15) | (x0.norm() > delta)) {
        for (int i = 0 ; i < n ; ++i )
        {
            rhs[i] = 0;
            sol[i] = 0;
        }
        for (int i = n ; i < npm ; ++i )
        {
            rhs[i] = b.get(i - n, 0);
            sol[i] = 0;
        }
        success = NOMAD::ldl_solve(error_msg, D, L, rhs, sol, pp, npm);
        if (!success) {
            throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt solve");
        }
        for (int i = 0 ; i < n ; ++i )
        {
            x.set(i, 0, sol[i]);
        }
    }

    // r = Gx + c
    SGTELIB::Matrix r = SGTELIB::Matrix::product(G, x);
    r.add(c);
    // Init r+
    SGTELIB::Matrix rp = r;

    // g = Pr
    SGTELIB::Matrix g("g", n, 1);
    for (int i = 0 ; i < n ; ++i )
    {
        rhs[i] = r.get(i, 0);
        sol[i] = 0;
    }
    for (int i = n ; i < npm ; ++i )
    {
        rhs[i] = 0;
        sol[i] = 0;
    }
    success = NOMAD::ldl_solve(error_msg, D, L, rhs, sol, pp, npm);
    if (!success) {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt solve");
    }
    for (int i = 0 ; i < n ; ++i )
    {
        g.set(i, 0, sol[i]);
    }
    
    // Init g+
    SGTELIB::Matrix gp = g;

    // d = -g
    SGTELIB::Matrix d = g;
    d.multiply(-1);

    // Gd = G * d
    SGTELIB::Matrix Gd = SGTELIB::Matrix::product(G, d);

    // Loop:
    double alpha, beta, tau;
    double xtd, nd2, xd2;

    // Stopping criteria
    const size_t max_iter = npm * 2;
    size_t iter = 0;
    double rg = SGTELIB::Matrix::dot(r, g);
    double rgp;
    double dtGd = SGTELIB::Matrix::dot(d, Gd);
    bool bounds = x.norm() < delta; // false if we do not satisfy bounds

    bool OK = (rg < tol) | (dtGd <= 0) | !bounds;

    verbose && std::cout << "PCG-It " << iter << " :";
    verbose && std::cout << " rg=" << rg << " tol=" << tol;
    verbose && std::cout << " dtGd=" << dtGd;
    verbose && std::cout << " nx=" << x.norm() << " <= delta=" << delta;
    verbose && std::cout << " feas=" << feas << std::endl;

    while (!OK & (iter < max_iter)) {
        alpha = rg / dtGd;

        // tau is the largest step to the boundary
        xtd = SGTELIB::Matrix::dot(x, d);
        nd2 = d.normsquare();
        xd2 = x.normsquare();
        tau = (-xtd + sqrt(xtd + nd2 * (std::pow(delta, 2) - xd2))) / nd2;
        alpha = std::min(alpha, tau);

        // x = x + alpha * d
        x.add(SGTELIB::Matrix::product(d, alpha));

        // r+ = r + alpha * G * d
        rp = r;
        Gd.multiply(alpha);
        rp.add(Gd);

        // g+ = Pr+
        for (int i = 0 ; i < n ; ++i )
        {
            rhs[i] = rp.get(i, 0);
            sol[i] = 0;
        }
        for (int i = n ; i < npm ; ++i )
        {
            rhs[i] = 0;
            sol[i] = 0;
        }
        success = NOMAD::ldl_solve(error_msg, D, L, rhs, sol, pp, npm);
        if (!success) {
            throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt solve");
        }
        for (int i = 0 ; i < n ; ++i )
        {
            gp.set(i, 0, sol[i]);
        }

        // beta = r+' g+ / r'g
        rgp = SGTELIB::Matrix::dot(rp, gp);
        beta = rgp / rg;

        // d = -(g+) + beta * d
        d.multiply(beta);
        d.sub(gp);

        g = gp;
        r = rp;
        rg = rgp;
        SGTELIB::Matrix::inplace_product(Gd, G, d);
        dtGd = SGTELIB::Matrix::dot(d, Gd);

        // Just for information:
        SGTELIB::Matrix::inplace_product(Ax, A, x);
        f = b; f.sub(Ax); feas = f.norm();

        bounds = x.norm() < delta;
        OK = (rg < tol) | (dtGd <= 0) | !bounds;
        iter += 1;

        verbose && std::cout << "PCG-It " << iter << " :";
        verbose && std::cout << " rg=" << rg << " tol=" << tol;
        verbose && std::cout << " dtGd=" << dtGd;
        verbose && std::cout << " nx=" << x.norm() << " <= delta=" << delta;
        verbose && std::cout << " feas=" << feas << std::endl;

    }

    for (int i=0 ; i<npm ; i++)
    {
        delete [] M[i];
    }
    delete [] M;
    for (int j=0 ; j<npm ; j++)
    {
        delete [] L[j];
    }
    delete [] L;
    for (int j=0 ; j<npm ; j++)
    {
        delete [] D[j];
    }
    delete [] D;
    delete [] pp;
    delete [] rhs;
    delete [] sol;

    if ((iter == 0) && !bounds) {
        // TODO: do something here
        verbose && std::cout << " Initial guess does not satisfy trust-region constraint" << std::endl;
        // return false;
        x.add(d);
        x.multiply(delta / d.norm());
    }

    if ((iter == 0) && (dtGd <= 0)) {
        verbose && std::cout << " PCG stopped at iteration 0 with negative curvature, do a projected gradient step:" << std::endl;
        x.add(d);
    }

    return true;
}

//*****************************************************************************
//
// Model getters with NOMAD points
//
//*****************************************************************************

// Get all model outputs on point x
SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelOut(const NOMAD::Point & x) const
{

    
    
    // Verify there is at least one point to evaluate
    if (!x.isComplete())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Evaluator: eval_x called with undefined eval point");
    }
    
    
    
    // Init the matrices for prediction
    // Creation of matrix for input / output of SGTELIB model
    SGTELIB::Matrix Mpredict (  "M_predict", 1, static_cast<int>(_m));
    SGTELIB::Matrix Xpredict("X_predict", 1, static_cast<int>(_n));
    
    std::string s = "X =" + x.display();
    NOMAD::OutputQueue::Add(s, _displayLevel);
    
    // Set the input matrix
    for (int i = 0; i < _n; i++)
    {
        Xpredict.set(0, static_cast<int>(i), x[i].todouble());
    }
     
    // ------------------------- //
    //   Output Prediction    //
    // ------------------------- //
    NOMAD::OutputQueue::Add("Predict with quadratic formulation... ", _displayLevel);
    
    
    // Unfortunately, Sgtelib is not thread-safe.
    // For this reason we have to set part of the eval_x code to critical.
#ifdef _OPENMP
#pragma omp critical(SgtelibEvalBlock)
#endif // _OPENMP
    {
        _model->check_ready(__FILE__,__FUNCTION__,__LINE__);
        
        _model->predict(Xpredict, &Mpredict);
        NOMAD::OutputQueue::Add("ok", _displayLevel);
    }
    
    return Mpredict;
}

double NOMAD::QPSolverOptimize::getModelObj(const NOMAD::Point &X) const
{
    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, static_cast<int>(i), X[i].todouble());
    }
    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    double fx = surrogate_prs->getModelObj(XX);

    return fx;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelGrad(const NOMAD::Point &X) const
{
    SGTELIB::Matrix Gx("Gx", static_cast<int>(_n), 1);
    getModelGrad(&Gx, X);
    return Gx;
}

void NOMAD::QPSolverOptimize::getModelGrad(SGTELIB::Matrix * Gx, const NOMAD::Point &X) const
{
    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, i, X[i].todouble());
    }

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    SGTELIB::Matrix Mpredict_grad (  "grad_predict", static_cast<int>(_m), static_cast<int>(_n));
    surrogate_prs->getModelGrad(Gx, &Mpredict_grad, XX);
    lencheck(_n, *Gx);
}

// Get hessian of model for output j on point x
SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelHessian(const NOMAD::Point & x, int j) const
{
        
    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, i, x[i].todouble());
    }

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    SGTELIB::Matrix Hx = surrogate_prs->getModelHessian(XX, j);
    sizecheck(_n, _n, Hx);
    
    return Hx;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelHessian(const NOMAD::Point &X) const
{   
        
    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, i, X[i].todouble());
    }

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    SGTELIB::Matrix Hx = surrogate_prs->getModelHessian(XX);
    sizecheck(_n, _n, Hx);

    return Hx;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelCons(const NOMAD::Point & X_k) const
{
    SGTELIB::Matrix cons("cons", static_cast<int>(_nbCons), 1);
    getModelCons(&cons, X_k);
    return cons;
}

void NOMAD::QPSolverOptimize::getModelCons(SGTELIB::Matrix * cons, const NOMAD::Point & X_k) const
{

    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, i, X_k[i].todouble());
    }
    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    surrogate_prs->getModelCons(cons, XX);
    lencheck(_nbCons, *cons);
}

// Get gradient of all model outputs on point x
SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelJacobian(const NOMAD::Point & X) const
{

    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < static_cast<int>(_n); i++)
    {
        XX.set(0, i, X[i].todouble());
    }
    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    SGTELIB::Matrix Jx = surrogate_prs->getModelJacobian(XX);
    sizecheck(_nbCons, _n, Jx);
    
    return Jx;
}

double NOMAD::QPSolverOptimize::getModelLag(const NOMAD::Point &x,
                                            const SGTELIB::Matrix &multiplier,
                                            const double &sigma) const
{
    lencheck(_nbCons, multiplier);

    double lag = sigma * getModelObj(x);

    if (_nbCons > 0) {
        SGTELIB::Matrix cx = getModelCons(x);
        lencheck(_nbCons, cx);

        lag -= SGTELIB::Matrix::dot(cx, multiplier);

    }

    return lag;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelLagGradient(
    const NOMAD::Point &x,
    const SGTELIB::Matrix &multiplier,
    const double &sigma ) const
{
    const int nbVar = _n;

    SGTELIB::Matrix lagGradient("lagGradient", nbVar, 1);
    lagGradient.fill(0.0);
    
    SGTELIB::Matrix outGradient("tmp", nbVar, 1);
    SGTELIB::Matrix modelJacobian = getModelJacobian(x);

    lencheck(_nbCons, multiplier);
    sizecheck(_nbCons, nbVar, modelJacobian);
    
    getModelGrad(&outGradient, x);
    outGradient.multiply(sigma) ;
    lagGradient.add(outGradient);

    SGTELIB::Matrix::inplace_product(outGradient, modelJacobian.transpose(), multiplier);
    outGradient.multiply(-1) ;
    lagGradient.add(outGradient);

    return lagGradient;
}


SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelLagHessian(
    const NOMAD::Point &X,
    const SGTELIB::Matrix &multiplier,
    const double &sigma) const
{

    SGTELIB::Matrix XX("X_k", 1, _n);
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, i, X[i].todouble());
    }
    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    SGTELIB::Matrix lagHx = surrogate_prs->getModelLagHessian(XX, multiplier);
    sizecheck(_n, _n, lagHx);

    return lagHx;
}

//*****************************************************************************
//
// Model getters for quadratic functionals
//
//*****************************************************************************
double NOMAD::QPSolverOptimize::getModelObj(
    const SGTELIB::Matrix & x,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g,
    const double g0) const
{
    int n = x.get_nb_rows();

    lencheck(n, x);
    lencheck(n, g);
    sizecheck(n, n, H);

    // Do not use  Sgtelib matrix operations for faster computation
    double sum = g0;
    for (int i = 0 ; i < n ; i++)
    {
        sum += g.get(i,0)*x.get(i,0);
        double phxi = 0;
        for (int j = 0 ; j < n ; j++)
        {
            phxi += H.get(i,j)*x.get(j,0);
        }
        sum += 0.5*x.get(i,0)*phxi;
    }
    return  sum;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelGrad(
    const SGTELIB::Matrix & x,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g
) const
{
    const int n = x.get_nb_rows();

    lencheck(n, x);
    lencheck(n, g);
    sizecheck(n, n, H);

    SGTELIB::Matrix Gx("Gx", n, 1);
    getModelGrad(&Gx, x, H, g);
    lencheck(n, Gx);
    return Gx;
}

void NOMAD::QPSolverOptimize::getModelGrad(
    SGTELIB::Matrix * Gx,
    const SGTELIB::Matrix & x,
    const SGTELIB::Matrix & H,
    const SGTELIB::Matrix & g
) const
{
    SGTELIB::Matrix::inplace_product(*Gx, H, x);
    Gx->add(g);
}

//*****************************************************************************
//
// Utils
//
//*****************************************************************************

int NOMAD::QPSolverOptimize::sum(
    const bool * indices,
    int len) const
{
    if (len < 0 )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: len should be > 0");
    }
    int sum = 0;
    for (int i = 0; i < len; i++)
    {
        if (indices[i])
        {
            sum += 1;
        }
    }
    return sum;
}

void NOMAD::QPSolverOptimize::lencheck(const int n, const SGTELIB::Matrix & X) const
{
    if (X.get_nb_rows() != n || X.get_nb_cols() != 1 )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, X.get_name() + " has wrong dimensions!");
    }
}

void NOMAD::QPSolverOptimize::sizecheck(const int m, const int n, const SGTELIB::Matrix & X) const
{
    if (X.get_nb_rows() != m || X.get_nb_cols() != n )
    {   
        std::cout << X.get_nb_rows() << " != " << m << " and " << X.get_nb_cols() << " != " << n << std::endl;
        throw NOMAD::Exception(__FILE__, __LINE__, X.get_name() + " has wrong dimensions!");
    }
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::vector_subset(
    const SGTELIB::Matrix & X,
    const bool * active)
{
    const int n = X.get_nb_rows();
    const int nfree = n - sum(active, n);
    SGTELIB::Matrix Xsub("Xsub", nfree, 1);
    Xsub.fill(0);

    int ki = 0;
    for (int i=0; i < n; ++i)
    {
        if (!active[i])
        {
            Xsub.set(ki, 0, X.get(i, 0));
            ki += 1;
        }
    }
    if (ki != nfree)
    {
        throw NOMAD::
        Exception(__FILE__, __LINE__, "Error dimension");
    }

    return Xsub;
}

void NOMAD::QPSolverOptimize::vector_broadcast(
    SGTELIB::Matrix * X,
    const SGTELIB::Matrix & Xsub,
    const bool * active)
{
    const int nsub = Xsub.get_nb_rows();
    int ki = 0;

    for (int i = 0 ; i < nsub ; ++i )
    {
        if (!active[i])
        {
            X->set(i, 0, Xsub.get(ki, 0));
            ki += 1;
        }
    }
    if (ki != nsub)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error dimension");
    }

}

SGTELIB::Matrix NOMAD::QPSolverOptimize::matrix_subset(
    const SGTELIB::Matrix & X,
    const bool * active)
{
    int n = X.get_nb_rows();
    int nfree = n - sum(active, n);
    SGTELIB::Matrix Xsub("Xsub", nfree, nfree);

    int ki = 0;
    for (int i=0; i < n; ++i)
    {
        int kj = 0;
        if (!active[i])
        {
            for (int j=0; j < n; ++j)
            {
                if (!active[j])
                {
                    Xsub.set(ki, kj, X.get(i, j));
                    kj += 1;
                }
            }
            ki += 1;
        }
    }
    if (ki != nfree)
    {
        throw NOMAD::
        Exception(__FILE__, __LINE__, "Error dimension");
    }

    return Xsub;
}

double NOMAD::QPSolverOptimize::max_step_bounds(const NOMAD::Point & X, const SGTELIB::Matrix & d)
{
    double t_max = INF;
    double gamma = INF;
    if (d.get_nb_rows() != _n)
    {
        
        // should return an error
    }
    // X should be feasible
    for (int i = 0; i < _n ; ++i)
    {
        double di = d.get(i, 0);
        if (di > 0 && (_modelUpperBound[i] - X[i]) > 0)
        {
            gamma = (_modelUpperBound[i] - X[i]).todouble()/std::fabs(di);
        }
        else if (di < 0 && (X[i] - _modelLowerBound[i]) > 0)
        {
            gamma = (X[i] - _modelLowerBound[i]).todouble()/std::fabs(di);
        }
        else
        {
            // nothing
        }
        if (gamma < t_max)
        {
            t_max = gamma;
        }
    }
    return t_max;
}

double NOMAD::QPSolverOptimize::max_step_bounds(
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    const SGTELIB::Matrix & d)
{
    int n = X.get_nb_rows();

    lencheck(n, X);
    lencheck(n, lvar);
    lencheck(n, uvar);
    lencheck(n, d);

    bool bound_compat = true;
    bool feasible = true;
    for (int i=0; i < n; ++i)
    {
        feasible = true;
        bound_compat = bound_compat && (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!bound_compat)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: Error compatibility lower and upper bound");
        }

        feasible = feasible && (X.get(i, 0) >= lvar.get(i, 0));
        feasible = feasible && (X.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            // Temp for debugging
            // throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: Error X is not feasible");
            std::cout << lvar.get(i, 0) << " " << X.get(i, 0) << " " << uvar.get(i, 0) << std::endl;
            std::cout << "Assertion error: Error X is not feasible" << std::endl;
        }
    }

    double t_max = INF;
    double gamma = INF;
    double di;

    for (int i = 0; i < n ; ++i)
    {
        di = d.get(i, 0);
        if (di > 0)
        {
            gamma = (uvar.get(i, 0) - X.get(i, 0))/std::fabs(di);
        }
        else if (di < 0)
        {
            gamma = (X.get(i, 0) - lvar.get(i, 0))/std::fabs(di);
        }
        else
        {
            // nothing
        }
        if (gamma < t_max)
        {
            t_max = gamma;
        }
    }

    return t_max;
}

void NOMAD::QPSolverOptimize::project_bounds(NOMAD::Point & X_k, SGTELIB::Matrix & d_k) 
{
    for (int i = 0 ; i < _n ; ++i )
    {
        if ( _modelLowerBound[i].isDefined() && X_k[i] == _modelLowerBound[i] && d_k.get(i, 0) < 0 )
        {
            d_k.set(i, 0, 0);
        }
        else if (_modelLowerBound[i].isDefined() && X_k[i] == _modelUpperBound[i] && d_k.get(i, 0) > 0 )
        {
            d_k.set(i, 0, 0);
        }
        else {
            // no change
        }
    }

}

void NOMAD::QPSolverOptimize::project_bounds(SGTELIB::Matrix & d_k, bool * active) 
{
    for (int i = 0 ; i < _n ; ++i )
    {
        if ( active[i] )
        {
            d_k.set(i, 0, 0);
        }
        else {
            // no change
        }
    }

}

/*
Update the direction d to project over the bounds
*/
void NOMAD::QPSolverOptimize::project_bounds(
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    SGTELIB::Matrix & d_k)
{
    int n = d_k.get_nb_rows();
    lencheck(n, X);
    lencheck(n, lvar);
    lencheck(n, uvar);
    lencheck(n, d_k);

    bool bound_compat = true;
    bool feasible = true;
    for (int i=0; i < n; ++i)
    {
        bound_compat = bound_compat && (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!bound_compat)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: Error compatibility lower and upper bound");
        }

        feasible = feasible && (X.get(i, 0) >= lvar.get(i, 0));
        feasible = feasible && (X.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: Error X is not feasible");
        }
    }

    for (int i = 0 ; i < n ; ++i )
    {
        if ( X.get(i, 0) == lvar.get(i, 0) && d_k.get(i, 0) < 0 )
        {
            d_k.set(i, 0, 0);
        }
        else if (X.get(i, 0) == uvar.get(i, 0) && d_k.get(i, 0) > 0 )
        {
            d_k.set(i, 0, 0);
        }
        else {
            // no change
        }
    }
}

void NOMAD::QPSolverOptimize::active_bounds(
    const SGTELIB::Matrix & X,
    bool * active_l,
    bool * active_u)
{

    lencheck(_n, X);

    for (int i = 0; i < _n; ++i)
    {
        active_l[i] = (X.get(i, 0) == _modelLowerBound[i]);
        active_u[i] = (X.get(i, 0) == _modelUpperBound[i]);
    }    
}

void NOMAD::QPSolverOptimize::active_bounds(
    const SGTELIB::Matrix & X,
    const SGTELIB::Matrix & lvar,
    const SGTELIB::Matrix & uvar,
    bool * active_l, // length n
    bool * active_u, // length n
    double tol)
{
    const int n = X.get_nb_rows();
    lencheck(n, X);
    lencheck(n, lvar);
    lencheck(n, uvar);

    if ( tol < 0 )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter error tol should be nonnegative");
    }

    for (int i = 0; i < n; ++i)
    {
        active_l[i] = (fabs(X.get(i, 0) - lvar.get(i, 0)) < tol);
        active_u[i] = (fabs(X.get(i, 0) - uvar.get(i, 0)) < tol);
    }    
}

/*
Return true for indices where one of the two cases hold:
- the lower bound is active and grad >= 0;
- the upper bound is active and grad <= 0.
*/
void NOMAD::QPSolverOptimize::binding_bounds(
    SGTELIB::Matrix & Grad,
    const bool * active_l,
    const bool * active_u,
    bool * binding)
{
    const int n = Grad.get_nb_rows();
    for (int i = 0; i < n; ++i)
    {
        binding[i] = (active_l[i] && (Grad[i] >= 0)) || (active_u[i] && (Grad[i] <= 0));
    }    
}

void NOMAD::QPSolverOptimize::getModelActiveCons(
    const SGTELIB::Matrix &cons,
    const double tol,
    bool * indices // length _nbCons
) const
{

    lencheck(_nbCons, cons);

    if (tol < 0 )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: tol should be > 0");
    }

    for (int i = 0; i < static_cast<int>(_nbCons) ; i++)
    {
        indices[i] = (fabs(cons.get(i, 0)) <= tol);
    }
}

void NOMAD::QPSolverOptimize::getModelFeasibleCons(
    const SGTELIB::Matrix &cons,
    const double tol,
    bool * indices // length _nbCons
) const
{

    lencheck(_nbCons, cons);

    if (tol < 0 )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: tol should be > 0");
    }

    for (int i = 0; i < static_cast<int>(_nbCons) ; i++)
    {
        indices[i] = (cons.get(i, 0) < -tol);
    }
}

void NOMAD::QPSolverOptimize::getModelInfeasibleCons(
    const SGTELIB::Matrix &cons,
    const double tol,
    bool * indices // length _nbCons
) const
{

    lencheck(_nbCons, cons);

    if (tol < 0 )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: tol should be > 0");
    }

    for (int i = 0; i < static_cast<int>(_nbCons) ; i++)
    {
        indices[i] = (cons.get(i, 0) > tol);
    }
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::feasibility(const SGTELIB::Matrix &cons) const
{
    SGTELIB::Matrix feas("feas", static_cast<int>(_nbCons ), 1);
    feasibility(&feas, cons);
    return feas;
}

void NOMAD::QPSolverOptimize::feasibility(SGTELIB::Matrix * feas, const SGTELIB::Matrix &cons) const
{

    lencheck(_nbCons, cons);
    
    for (int i = 0; i < static_cast<int>(_nbCons) ; i++)
    {
        feas[i] = max(cons.get(i, 0), 0).todouble();
    }
}

bool NOMAD::QPSolverOptimize::isFeasible(
    const SGTELIB::Matrix & cons,
    const double tol
) const 
{

    lencheck(_nbCons, cons);

    if (tol < 0 )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: tol should be > 0");
    }
    
    for (int i = 0; i < static_cast<int>(_nbCons) ; i++)
    {
        if (cons.get(i, 0) > tol)
        {
            return false;
        }
    }
    return true;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelActiveJacobian(
    const SGTELIB::Matrix & Jacobian,
    const bool * indices // length _nbCons
) const
{
    const int n = static_cast<int>(_n);
    const int ncon =  static_cast<int>(_nbCons);
    const int nbActive = sum(indices, ncon);
    SGTELIB::Matrix activeJacobian (  "activeJacobian", nbActive, n);
    getModelActiveJacobian(&activeJacobian, Jacobian, indices);
    return activeJacobian;
}

void NOMAD::QPSolverOptimize::getModelActiveJacobian(
    SGTELIB::Matrix * activeJacobian,
    const SGTELIB::Matrix & Jacobian,
    const bool * indices // length _nbCons
) const
{
    const int n = static_cast<int>(_n);
    const int ncon =  static_cast<int>(_nbCons);

    sizecheck(ncon, n, Jacobian);

    int nbActive = sum(indices, ncon);

    int k = 0;
    for (int i = 0; i < ncon; i++)
    {
        if (indices[i])
        {
            for (int j = 0; j < n; j++)
            {
                activeJacobian->set(k, j, Jacobian.get(i, j));
            }
            k++;
        }
    }

    if (k != nbActive)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error");
    }
}


void NOMAD::QPSolverOptimize::snapToBounds(SGTELIB::Matrix & X, const SGTELIB::Matrix & lowerBound, const SGTELIB::Matrix & upperBound) const
{
    
    const int n = X.get_nb_rows();
    
    if (lowerBound.get_nb_rows() != n || upperBound.get_nb_rows() != n)
    {
        std::string err = "snapToBounds: ";
        err += "Inconsistent dimension for bounds or for X. Expecting a vector";
        err += std::to_string(n);
        err += " but sizes are " + std::to_string(lowerBound.get_nb_rows());
        err += " and " + std::to_string(upperBound.get_nb_rows()) + ".";
        throw Exception(__FILE__, __LINE__, err);
    }

    for (int i = 0; i < n; ++i)
    {
        if (X.get(i, 0) < lowerBound.get(i, 0))
        {
            X.set(i, 0, lowerBound.get(i, 0));
        }
        else if (upperBound.get(i, 0) < X.get(i, 0))
        {
            X.set(i, 0, upperBound.get(i, 0));
        }
        else
        {
            // no snap
        }
    }

    
}

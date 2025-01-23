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
#include "../../Algos/QPSolverAlgo/DoglegTRSolver.hpp"
#include "../../Algos/QPSolverAlgo/L1AugLagSolver.hpp"
#include "../../Algos/QPSolverAlgo/ProjectedConjugateGradientSolver.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverOptimize.hpp"
#include "../../Algos/QPSolverAlgo/QPSolverAlgo.hpp"
#include "../../Algos/QPSolverAlgo/TRIPMSolver.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Math/MathUtils.hpp"
#include "../../Math/MatrixUtils.hpp"

#include "../../../ext/sgtelib/src/Surrogate_PRS.hpp"

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
    
    // Get the model box size limit.
    // Compare with the actual model bounds to enable or not the search at start
    _modelBoxSizeLimit = _runParams->getAttributeValue<NOMAD::Double>("QP_SEARCH_MODEL_BOX_SIZE_LIMIT");
    
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

    // Test size of model box bounds
    bool modelBoxOk = false;
    for (int i = 0; i < _n; i++)
    {
        if ( _modelLowerBound[i].isDefined() && _modelUpperBound[i].isDefined() && _modelUpperBound[i] - _modelLowerBound[i] > _modelBoxSizeLimit )
        {
            modelBoxOk = true;
            break;
        }
    }
    // Model box is too small, no point generation.
    if (!modelBoxOk)
    {
        OUTPUT_INFO_START
        std::string s;
        s = "Bounds are too tight. Do not perform QP Solver point generation.";
        AddOutputInfo(s, _displayLevel);
        OUTPUT_INFO_END
        return;
    }

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
    bool foundBetter = false;
    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
        
        if (_modelFixedVar.nbDefined() > 0)
        {
            NOMAD::EvalPointSet evalPointSet;
            for (const auto& trialPoint : _trialPoints)
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
    // QPSolverOptimize can be called by QPSolverSearchMethod OR as a standalone optimization
    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,getName() + " must have an Iteration ancestor.");
    }

    // We need a base point.
    const auto refCenter = getParentOfType<QuadModelIteration*>()->getRefCenter();
    NOMAD::Point X_k = *(refCenter->getX());
    if (_verboseFull)
    {
        SGTELIB::Matrix out = getModelOut(X_k);
        std::cout <<"Model output for RefCenter: " << std::endl;
        out.display(std::cout);
        auto bbo = refCenter->getEval(NOMAD::EvalType::BB)->getBBOutput();
        std::cout<< "Blackbox output of RefCenter: " << bbo.getBBO() << std::endl;
    }

    auto isInBounds = [](const NOMAD::ArrayOfDouble& x,
                         const NOMAD::ArrayOfDouble& lb,
                         const NOMAD::ArrayOfDouble& ub) -> bool
    {
        for (size_t i = 0; i < x.size(); ++i)
        {
            if (((lb[i].isDefined()) && (x[i] < lb[i])) ||
                ((ub[i].isDefined()) && (x[i] > ub[i])))
            {
                return false;
            }
        }
        return true;
    };

    bool feasRefCenter = false;
    const auto computeType = getMegaIterationBarrier()->getFHComputeType();
    if (refCenter->isFeasible(computeType))
    {
        if ( nullptr != _prevFeasRefCenter &&
            _prevFeasRefCenter->NOMAD::ArrayOfDouble::isDefined() &&
            X_k == *(_prevFeasRefCenter->getX()) &&
            _prevFeasXopt.isComplete() &&
            isInBounds(_prevFeasXopt, _modelLowerBound, _modelUpperBound))
        {
            X_k = _prevFeasXopt;
        }
        _prevFeasRefCenter = refCenter;
        feasRefCenter = true;
    }
    else if (! refCenter->isFeasible(computeType))
    {
        if (nullptr != _prevInfeasRefCenter &&
            _prevInfeasRefCenter->NOMAD::ArrayOfDouble::isDefined() &&
            X_k == *(_prevInfeasRefCenter->getX()) &&
            _prevInfeasXopt.isComplete() &&
            isInBounds(_prevInfeasXopt, _modelLowerBound, _modelUpperBound))
        {
            X_k = _prevInfeasXopt;
        }
        _prevInfeasRefCenter = refCenter;
        feasRefCenter = false;
    }

    
    OUTPUT_INFO_START
    std::string s = "Model base X:" + X_k.display();
    AddOutputInfo(s, _displayLevel);
    OUTPUT_INFO_END
    bool runOk = false;
    if (_m == 1 && _bbot[0].isObjective())
    {
        runOk = solveBCQP(X_k);
    }
    else
    {
        auto quadModelIter = getParentOfType<NOMAD::QuadModelIteration*>(false);
        if (nullptr == quadModelIter)
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"QPSolverOptimize must have a quadModelIteration as parent");
        }
        
        const auto mesh = quadModelIter->getMesh();
        double MeshSize = 0.0;
        if (nullptr != mesh)
        {
            // auto FrameSize = mesh->getDeltaFrameSize().max().todouble();
            MeshSize = mesh->getdeltaMeshSize().max().todouble();
            const auto meshIndex = mesh->getMeshIndex();
            _verbose && std::cout << " meshIndex=" << meshIndex << std::endl;
        }

        SGTELIB::Matrix Gk("Gk", static_cast<int>(_n), 1);
        getModelGrad(&Gk, X_k);
        const double ng = Gk.norm();

        SGTELIB::Matrix Y("Y", static_cast<int>(_nbCons), 1);
        Y.fill(1.0);
        SGTELIB::Matrix HLag = getModelLagHessian(X_k, Y);
        SGTELIB::Matrix svdHLag = HLag.get_singular_values();
        const double sing_val_min = svdHLag.min();
        const double condHessian = (sing_val_min > 0) ? svdHLag.max() / sing_val_min : NOMAD::INF;

        // Fix tolerances
        const auto tolMesh = _runParams->getAttributeValue<Double>("QP_tolMesh").todouble();
        const auto tolCond = _runParams->getAttributeValue<Double>("QP_tolCond").todouble();
        auto atol =_runParams->getAttributeValue<Double>("QP_absoluteTol").todouble();
        auto rtol = _runParams->getAttributeValue<Double>("QP_relativeTol").todouble();
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
        const double ng0 = ng;
        const double tol = atol + ng0 * rtol;

        const auto maxIter = static_cast<int>(_runParams->getAttributeValue<size_t>("QP_maxIter"));
        const auto tolDistDX = _runParams->getAttributeValue<Double>("QP_tolDistDX").todouble();

        // Parameters specific to augmented Lagrangian
        const auto mu0 = _runParams->getAttributeValue<Double>("QP_AugLag_mu0").todouble();
        const auto muDecrease = _runParams->getAttributeValue<Double>("QP_AugLag_muDecrease").todouble();

        const auto eta0 = _runParams->getAttributeValue<Double>("QP_AugLag_eta0").todouble();
        const auto omega0 = _runParams->getAttributeValue<Double>("QP_AugLag_omega0").todouble();

        const auto successRatio = _runParams->getAttributeValue<Double>("QP_AugLag_successRatio").todouble();
        const auto maxIterInner = _runParams->getAttributeValue<size_t>("QP_AugLag_maxIterInner");
        const auto tolDistDXInner = _runParams->getAttributeValue<Double>("QP_AugLag_tolDistDXInner").todouble();
        const auto maxSuccessiveFail = _runParams->getAttributeValue<size_t>("QP_AugLag_maxSuccessivFail");

        const auto SelectAlgo = _runParams->getAttributeValue<size_t>("QP_SelectAlgo");
        if (SelectAlgo == 0)
        {
            _verbose && std::cout << "Run solveAugLag (n=" << _n << ", m=" << _nbCons << ")" << std::endl;
            _verbose && std::cout << "atol=" << atol << " rtol=" << rtol << " tol=" << tol << " cond(H)=" << condHessian << " mesh=" << MeshSize << std::endl;
            runOk = solveAugLag(X_k, maxIter, tolDistDX, atol, rtol, mu0, muDecrease, eta0, omega0, successRatio, maxIterInner, tolDistDXInner, maxSuccessiveFail);
        }
        else if (SelectAlgo == 1 || SelectAlgo == 2)
        {
            // Extract matrix and lower and upper bounds
            SGTELIB::Matrix QPModel = computeQPModelMatrix();
            const int nvar = _trainingSet->get_nvar();
            SGTELIB::Matrix lb("lb", nvar, 1);
            SGTELIB::Matrix ub("ub", nvar, 1);
            SGTELIB::Matrix x("x", nvar, 1);
            int k = 0;
            for (int i = 0; i < _n; ++i)
            {
                if (_trainingSet->get_X_nbdiff(i) <= 1)
                    continue;

                const double lbi = _modelLowerBound[i].isDefined() ? _modelLowerBound[i].todouble() : NOMAD::M_INF;
                const double ubi = _modelUpperBound[i].isDefined() ? _modelUpperBound[i].todouble() : NOMAD::INF;
                lb.set(k, 0, lbi);
                ub.set(k, 0, ubi);
                x.set(k, 0, X_k[i].todouble());
                k++;
            }
            if (SelectAlgo == 1)
            {
                _verbose && std::cout << "Run solveTRIPM (n=" << _n << ", m=" << _nbCons << ")" << std::endl;
                _verbose && std::cout << "atol=" << atol << " rtol=" << rtol << " tol=" << tol << " cond(H)=" << condHessian << std::endl;

                TRIPMSolver tripm_solver{mu0, muDecrease, 1e-12,
                                         80, 90, 0};

                auto status = tripm_solver.solve(x, QPModel, lb, ub);
                runOk = status != TRIPMSolverStatus::PARAM_ERROR &&
                        status != TRIPMSolverStatus::MATRIX_DIMENSIONS_FAILURE &&
                        status != TRIPMSolverStatus::NUM_ERROR &&
                        status != TRIPMSolverStatus::BOUNDS_ERROR;
            }
            else if (SelectAlgo == 2)
            {
                L1AugLagSolver l1_auglag_solver {1e-12, 50, 25, 0};
                auto status = l1_auglag_solver.solve(x, QPModel, lb, ub);
                runOk = status != L1AugLagSolverStatus::PARAM_ERROR &&
                        status != L1AugLagSolverStatus::MATRIX_DIMENSIONS_FAILURE &&
                        status != L1AugLagSolverStatus::NUM_ERROR &&
                        status != L1AugLagSolverStatus::BOUNDS_ERROR;
            }

            if (runOk)
            {
                k = 0;
                for (int i = 0; i < _n; ++i)
                {
                    if (_trainingSet->get_X_nbdiff(i) <= 1)
                        continue;

                    X_k[i] = x.get(k, 0);
                    k++;
                }
            }
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

            const bool strict = getStrictFeasiblePoint(X_k, XS, lvar, uvar, cons);
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

    if (!runOk) {
        OUTPUT_INFO_START
            std::string s = "Solver run NOT OK";
            AddOutputInfo(s);
        OUTPUT_INFO_END

        auto modelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_stopReasons);
        modelStopReasons->set(NOMAD::ModelStopType::MODEL_OPTIMIZATION_FAIL);
        return;
    }

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
        const bool inserted = insertTrialPoint(NOMAD::EvalPoint(X_k));

        OUTPUT_INFO_START
        std::string s = "xt:";
        s += (inserted) ? " " : " not inserted: ";
        s += X_k.display() + "\n";
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


// Set the bounds and the extra fixed variables (when lb=ub) of the model.
void NOMAD::QPSolverOptimize::setModelBoundsAndFixedVar()
{
    // When optWithScaleBounds is true, the training set is scaled with some generating directions. Warning: points are not necessarily in [0,1]^n
    const SGTELIB::Matrix & X = _trainingSet->get_matrix_X();
    
    _n = static_cast<int>(_pbParams->getAttributeValue<size_t>("DIMENSION"));
    
    if (_n != X.get_nb_cols())
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "QPSolverOptimize::setModelBounds() dimensions do not match");
    }
    
    const int nbDim = X.get_nb_cols();
    const int nbPoints = X.get_nb_rows();
    
    // Build model bounds and detect fixed variables
    //NOMAD::Double lb;
    //NOMAD::Double ub;
    bool isFixed = false;
    for (int j = 0; j < nbDim; j++)
    {
        NOMAD::Double lb = _modelLowerBound[j];
        NOMAD::Double ub = _modelUpperBound[j];
        for (int p = 0; p < nbPoints; p++)
        {
            // Use a comparison on regular double for better precision
            const auto xpj = NOMAD::Double(X.get(p,j));
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
            // When scaling we force the bounds to be in [0,1] whatever the training set except if we have a fixed variable.
            // If optWithScaledBounds is true, and we detect a fixed variable, the modelFixedVar is not necessarily within [0,1]
            // but the model center and the fixed variable must be equal.
            if (isFixed)
            {
                _modelLowerBound[j] = _modelUpperBound[j] = NOMAD::Double();
                _modelCenter[j] = _modelFixedVar[j];
            }
            // If not fixed, the model bounds and model center are predetermined.
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
        const auto reduction_factor = _runParams->getAttributeValue<NOMAD::Double>("QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR");
        for (int j = 0; j < nbDim; j++)
        {
            NOMAD::Double lb = _modelLowerBound[j];
            NOMAD::Double ub = _modelUpperBound[j];
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
        active[i] = false;
    }

    solve_TR_constrained_QP(d, X, H, g, grad, active, Delta);

    delete [] active;
}

void NOMAD::QPSolverOptimize::solve_TR_constrained_QP(
    SGTELIB::Matrix* d,
    const SGTELIB::Matrix& X,
    const SGTELIB::Matrix& H,
    const SGTELIB::Matrix& g,
    SGTELIB::Matrix& grad,
    const bool* active,
    const double Delta,
    const bool verbose)
{
    const int n = X.get_nb_rows();
    const int nfree = n - sum(active, n);

    // Check dimension compatibility
    lencheck(n, g);
    sizecheck(n, n, H);

    _verbose && std::cout << "Starting solve_TR_constrained_QP with delta=" << Delta << " nfree=" << nfree << std::endl;

    // Pre-allocation of active sub-matrices of the Hessian and the gradient
    getModelGrad(&grad, X, H, g); // Hx + g
    SGTELIB::Matrix gW = vector_subset(grad, active); // NB: depend on X
    gW.set_name("gW");
    SGTELIB::Matrix HW = matrix_subset(H, active);
    HW.set_name("HW");

    ////////// LDLt
    // Initialize matrices for LDLt decomposition. The LDLt decomposition is more generic than
    // the Cholesky decomposition, as for the former, the matrix needs to be positive-definite.
    // At the end of the initialization, M = HW.
    auto M = new double*[nfree];
    auto L = new double*[nfree];
    auto D = new double*[nfree];
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
    auto pp = new int[nfree];
    for (int i = 0 ; i < nfree ; ++i )
    {
        pp[i] = 0;
    }

    std::string error_msg;
    const bool success = NOMAD::LDLt_decomposition(error_msg, M, L, D, pp, nfree);
    if (!success)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt decomposition");
    }

    // NB: contrary to the title of this function, eigmin is not the minimum eigenvalue
    // of the hessian matrix, but it enables to indicate if the matrix is positive definite.
    const double eigmin = NOMAD::FindSmallestEigenvalue(D, nfree);
    _verbose && std::cout << " smallest eigenvalue= " << eigmin << std::endl;

    if (eigmin > 0) // positive definite
    {
        d->fill(0);
        const bool newton_success = Convex_TR_QP(d, grad, gW, H, HW, pp, D, L, active, Delta, _verbose);
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

        // See https://en.wikipedia.org/wiki/Inverse_iteration for an explanation
        // Compute an eigenvector associated to the eigmin value.
        // One could also use the Rayleigh quotient iteration algorithm, which
        // is more efficient and robust that the inverse iteration algorithm...
        successInverseIteration = InverseIteration(&bk, HW, eigmin, nfree, 1E-12);

        if (!successInverseIteration)
        {
            std::cerr << "Error InverseIteration" << std::endl;
            d->fill(0.0);
        }
        else 
        {
            // Move in the direction of the eigenvector associated to the minimum
            // eigenvalue.
            const double nd = bk.norm();
            const double slope = SGTELIB::Matrix::dot(gW, bk);
            const double a_unb = (Delta < 1e15) ? Delta / nd : 1000 * std::abs(slope) / nd;

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

    // Free memory
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

    const bool success = solveBCQP(Xm, H, g, g0, lvar, uvar, max_iter, atol, rtol, verbose);

    for (int i=0; i < _n; ++i)
    {
        X[i] = Xm.get(i, 0);
    }

    return success;
}


void NOMAD::QPSolverOptimize::projectedGradient(
    SGTELIB::Matrix& X,
    const SGTELIB::Matrix& H,
    const SGTELIB::Matrix& g,
    const double g0,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    bool* active_l,
    bool* active_u,
    SGTELIB::Matrix& d_k,
    SGTELIB::Matrix& Temp,
    const double kappa,
    const double tol,
    const int max_iter,
    const bool verbose)
{
    const int n = X.get_nb_rows();

    // Pre-allocation
    SGTELIB::Matrix armijo_Xp("armijo_Xp", n, 1);

    // active constraints
    auto active = new bool[n];
    for (int i = 0; i < n; ++i)
    {
        active[i] = active_l[i] || active_u[i];
    }
    auto activep = new bool[n];

    double qm = getModelObj(X, H, g, g0);
    getModelGrad(&d_k, X, H, g);

    const double a_max = 1.0;
    double a_k = 1.0;

    // Initialize the value between qm and qmp
    double diffqmqmp = 0;

    bool runOk = false;
    int k = 0;
    while (!runOk && (k < max_iter))
    {
        // Make dk a descent direction.
        d_k.multiply(-1.0);

        // Apply an Armijo projection linesearch. NB: at the end of the linesearch, a_k is >= 0
        const double slope = -d_k.normsquare();
        a_k = 1; // min(a_k, 1) is the initial value for Armijo linesearch.
        a_k = projected_armijo(X, H, g, g0, lvar, uvar, d_k, qm, slope, armijo_Xp, Temp, a_k); // a_k start from previous value

        // X := P_Omega(X + a_k d_k), where Omega = {y : lvar <= y <= uvar}
        d_k.multiply(a_k);
        X.add(d_k);
        snapToBounds(X, lvar, uvar);

        // Compute active sets
        active_bounds(X, lvar, uvar, active_l, active_u);
        int nfree = n;
        for (int i = 0; i < n; ++i)
        {
            activep[i] = active_l[i] || active_u[i];
            if (activep[i])
            {
                nfree --;
            }
        }

        // Update current iterates.
        const double qmp = getModelObj(X, H, g, g0);
        getModelGrad(&d_k, X, H, g);
        verbose && std::cout << "  Projected-gradient k =" << k << " f(x) = " << qmp << " |A| = " << nfree << " " << nfree;
        verbose && std::cout << " |d| = " << d_k.norm() << " amax = " << a_max << " ak = " << a_k << std::endl;

        // Compare active sets and update them
        bool checkActive = true;
        for (int i = 0; i < n; ++i)
        {
            if (active[i] != activep[i])
            {
                checkActive = false;
            }
            active[i] = activep[i];
        }

        // Compute stopping criteria.
        runOk = ((qm - qmp) <= kappa * diffqmqmp) || checkActive;
        diffqmqmp = std::max(qm - qmp, diffqmqmp);
        qm = qmp;
        k ++;
    }
    delete [] active;
    delete [] activep;
    // return X updated
}

bool NOMAD::QPSolverOptimize::conjugateGradient(
    SGTELIB::Matrix& X,
    const SGTELIB::Matrix& H,
    const SGTELIB::Matrix& g,
    const int maxIter,
    const double xi,
    const double atol,
    const double rtol,
    const bool verbose) const
{
    const int n = X.get_nb_rows();

    // Check dimension compatibility
    lencheck(n, g);
    sizecheck(n, n, H);

    // Pre-allocation of initial vectors
    SGTELIB::Matrix v("v", n, 1);
    for (int i = 0; i < n; i++)
    {
        v[i] = X[i];
    }

    // Pre-allocation of matrices used in the conjugate gradient algorithm
    SGTELIB::Matrix Ask("Ask", n, 1); // The matrix vector product A sk
    SGTELIB::Matrix rk("rk", n, 1);   // The current residual
    SGTELIB::Matrix sk("sk", n, 1);   // The conjugate gradient direction
    SGTELIB::Matrix Avk("Avk", n, 1); // The matrix vector product A vk
    SGTELIB::Matrix::inplace_product(Avk, H, v);

    // Initialization
    // Compute initial residual vector r0 := gW + HW v0
    rk = g;
    rk.add(Avk);

    // Initialize the conjugate direction s0 := -r0
    sk = rk;
    sk.multiply(-1.0);

    // Initialize gamma := rk^t rk
    double gamma =  SGTELIB::Matrix::dot(rk, rk);
    double sNormSquare = gamma;

    // Define tolerance
    double rNorm = std::sqrt(gamma);
    const double tol = atol + rtol * rNorm;

    // Declare some variables for verbose information
    double sAs =0 ;

    // Define stopping criteria
    bool solved = rNorm <= tol;
    bool zeroCurvature = false;

    double diffqmqmp = 0;
    double qmp = 0.5 * SGTELIB::Matrix::dot(Avk, v) + SGTELIB::Matrix::dot(v, g);
    double qm = qmp;

    // Main iteration
    int iter = 0;
    while (!((iter >= maxIter) || solved || zeroCurvature))
    {
        // Compute sk^t HW sk
        SGTELIB::Matrix::inplace_product(Ask, H, sk);
        sAs = SGTELIB::Matrix::dot(Ask, sk);

        // Negative curvature detection
        if (sAs <= atol * atol * sNormSquare)
        {
            if (std::abs(sAs) <= atol * sNormSquare)
            {
                zeroCurvature = true;
            }

            // When at iteration 0, set v := gW (the negative gradient)
            // otherwise, the iterate at the previous iteration will be
            // returned
            if (iter == 0)
            {
                v = g;
                v.multiply(-1.0);
            }

            solved = true;
        }

        if (zeroCurvature || solved)
            continue;

        // Compute alpha := rk^t rk / sk^t A sk, where A = HW
        double alpha = gamma / sAs;

        // vk := vk + alpha * sk
        for (int i = 0; i < n; ++i)
        {
           v[i] += alpha * sk[i];
        }

        // rk := rk + alpha * A sk, where A = HW
        for (int i = 0; i < n; ++i)
        {
            rk[i] += alpha * Ask[i];
        }

        double gammap = SGTELIB::Matrix::dot(rk, rk);
        rNorm = std::sqrt(gammap);

        // Compute stopping criterion
        solved = rNorm <= tol;
        if (solved)
        {
            continue;
        }

        // Update parameters
        // beta := rk+1^t rk+1 / rk^t rk
        const double beta = gammap / gamma;

        // sk+1 := -rk+1 + beta * sk
        sk.multiply(beta);
        sk.sub(rk);

        gamma = gammap;

        SGTELIB::Matrix::inplace_product(Avk, H, v);

        // Compute special criterion for sufficient conjugate gradient decrease
        qmp = 0.5 * SGTELIB::Matrix::dot(Avk, v) + SGTELIB::Matrix::dot(v, g);
        solved = solved || ((qm - qmp <= xi * diffqmqmp));
        diffqmqmp = std::max(qm - qmp, diffqmqmp);
        qm = qmp;
        iter += 1;
    }

    verbose && std::cout << "CG tol: " << tol;
    verbose && std::cout << " CG total niter: " << iter;
    verbose && std::cout << " CG residual norm:" << rNorm;
    verbose && solved && sAs <= 0 && std::cout << " Non positive curvature detected";
    verbose && std::cout << std::endl;

    // Compute solution in full dimension
    for (int i = 0; i < n; i++)
    {
            X[i] = v[i];
    }
    return (solved && sAs <= 0);
}

bool NOMAD::QPSolverOptimize::solveBCQP(
    SGTELIB::Matrix& X,
    const SGTELIB::Matrix& H,
    const SGTELIB::Matrix& g,
    const double g0,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    const int max_iter,
    const double atol,
    const double rtol,
    const bool verbose)
{
    const int n = X.get_nb_rows();

    // Check dimension compatibility
    lencheck(n, g);
    lencheck(n, lvar);
    lencheck(n, uvar);
    sizecheck(n, n, H);

    // Check that X is feasible
    bool feasible = true;
    for (int i=0; i < n; ++i)
    {
        if (X.get(i, 0) < lvar.get(i, 0) || X.get(i, 0) > uvar.get(i, 0))
        {
            feasible = false;
            break;
        }
    }
    if (!feasible)
    {
        verbose && std::cout << "solveBCQP assertion error: X not feasible. Compute projection." << std::endl;
        snapToBounds(X, lvar, uvar);
    }

    feasible = true;
    for (int i=0; i < n; ++i)
    {
        feasible = feasible && (X.get(i, 0) >= lvar.get(i, 0)) && (X.get(i, 0) <= uvar.get(i, 0));
        const bool areBoundsCompatible = lvar.get(i, 0) <= uvar.get(i, 0);
        if (!feasible || !areBoundsCompatible)
        {
            verbose && std::cout << lvar.get(i, 0) << " " << X.get(i, 0) << " " << uvar.get(i, 0);
            throw NOMAD::Exception(__FILE__, __LINE__, "solveBCQP assertion error: Error compatibility lower and upper bound");
        }
    }

    // Pre-allocation
    SGTELIB::Matrix d_k("dk", n, 1);
    SGTELIB::Matrix Xm("Xm", n, 1);
    SGTELIB::Matrix Xp("Xp", n, 1);
    SGTELIB::Matrix Grad("Grad", n, 1);
    SGTELIB::Matrix XmPXGrad("XmPXGrad", n, 1);
    SGTELIB::Matrix armijo_Xp("armijo_Xp", n, 1);
    SGTELIB::Matrix Temp("Temp-vector", n, 1);

    std::vector<bool> currentFree(n);
    std::vector<bool> currentActiveL(n);
    std::vector<bool> currentActiveU(n);

    // Initialisation of active, working and binding sets of constraints.
    auto active_l = new bool[n];
    auto active_u = new bool[n];
    auto active = new bool[n]; // subset of active_l OR active_u
    active_bounds(X, lvar, uvar, active_l, active_u);
    for (int i = 0; i < n; ++i)
    {
        active[i] = active_l[i] || active_u[i];
    }

    auto binding = new bool[n];
    binding_bounds(Grad, active_l, active_u, binding);

    // Initialization
    const double f0 = getModelObj(X, H, g, g0);
    getModelGrad(&Grad, X, H, g);
    const double ng0 = Grad.norm_inf();

    XmPXGrad = X;
    XmPXGrad.sub(Grad);
    snapToBounds(XmPXGrad, lvar, uvar);
    XmPXGrad.multiply(-1.0);
    XmPXGrad.add(X);

    double a_k;
    const double kappa = 0.1;

    // Start iterating
    int k = 0;
    verbose && std::cout << "  k = " << k << " f(x0) = " << f0;
    verbose && std::cout << " |W| = " << n - sum(active, n) << " |L| = " << sum(active_l, n) << " |U| = " << sum(active_u, n) << std::endl;
    bool OK = false;
    while (!OK && (k < max_iter))
    {
        // Check stopping criterion
        const double nXmPXGrad = XmPXGrad.norm_inf();
        OK = nXmPXGrad <= 1e-7 * ng0;
        if (OK)
        {
            continue;
        }

        // Save the candidate X to compute stopping criterion
        Xm = X;
        double qCurrent = getModelObj(X, H, g, g0);

        // Save the active sets at the beginning of iteration to compute stopping criterion
        for (int i = 0; i < n; ++i)
        {
            currentActiveL[i] = active_l[i];
            currentActiveU[i] = active_u[i];
        }

        // Generate a sequence of projected gradient steps. Update active_l and active_u at the same time.
        projectedGradient(X, H, g, g0, lvar, uvar, active_l, active_u, Grad, Temp, kappa, atol + ng0 * rtol, n);

        // Compute active set
        for (int i = 0; i < n; ++i)
        {
            active[i] = active_l[i] || active_u[i];
        }

        // When the set of working variables is empty, stop the algorithm.
        if (sum(active, n) == n)
        {
            OK = true;
            continue;
        }

        getModelGrad(&Grad, X, H, g);

        // Allocate "working" submatrices.
        for (int i = 0; i < n; ++i)
        {
            currentFree[i] = !active[i];
        }
        SGTELIB::Matrix dz("dz", n - sum(active, n), 1);
        dz.fill(0);
        SGTELIB::Matrix gz = vector_subset(Grad, active);
        SGTELIB::Matrix ZHZ = matrix_subset(H, active);

        // Generate a sequence of conjugate gradient steps within the working set.
        double xi = 1e-3;
        bool hasNegativeCurvature = conjugateGradient(dz, ZHZ, gz, 120, xi, 1e-7, 1e-7);

        // Compute direction
        d_k.fill(0);
        int kfree = 0;
        for (int i = 0; i < n; ++i)
        {
            if (currentFree[i])
            {
                d_k[i] = dz[kfree];
                kfree ++;
            }
        }

        if (d_k.has_nan()) {
            throw NOMAD::Exception(__FILE__, __LINE__, "d_k contains NaN");
        }

        const double a_max = max_step_bounds(X, lvar, uvar, d_k);
        if (hasNegativeCurvature)
        {
            // d_k is a direction of negative curvature: we find the smallest step gamma such that
            // X + gamma d_k is at the boundary of Omega^{k,j} = {y : lvar <= y <= uvar}.
            // We need the objective function to decrease (to compute the predicted reduction in
            // the trust-region test), we then add an additionnal test.

            // Xp = X + a_max d_k
            Xp = d_k; Xp.multiply(a_max); Xp.add(X);
            if (getModelObj(Xp, H, g) <= getModelObj(X,H,g))
            {
                // We accept the step.
                a_k = a_max;
            }
            else
            {
                // We perform a backtracking strategy, starting from a_max and decreasing along the line.
                const double qm = getModelObj(X, H, g, g0);
                double step = a_max;
                const int nbMaxTrials = 10;
                int niter = 0;
                bool finished = false;
                while (!finished)
                {
                    Xp = d_k; Xp.multiply(step); Xp.add(X);
                    snapToBounds(Xp, lvar, uvar);
                    const double qXp = getModelObj(Xp, H, g);
                    Xp.sub(X);
                    const double slope = SGTELIB::Matrix::dot(Grad, Xp);
                    if (qXp <= qm + 1e-4 * slope)
                    {
                        finished = true;
                        continue;
                    }

                    step /= 3.0;
                    niter += 1;
                    finished = niter >= nbMaxTrials;
                }
                a_k = (getModelObj(Xp, H, g) <= getModelObj(X,H,g)) ? step : 0;
            }
        }
        else
        {
            // d_k is not a direction of negative curvature: we execute a projected Armijo linesearch.
            // 1e-15 is the smallest value in projected_armijo
            const double slope = SGTELIB::Matrix::dot(d_k, Grad);
            const double qm = getModelObj(X, H, g, g0);
            a_k = (a_max > 1e-15) ? projected_armijo(X, H, g, g0, lvar, uvar, d_k, qm, slope, armijo_Xp, Temp, a_max)
                                  : 0;
            // a_k = std::min(a_k, a_max);
        }

        // Update the incumbent. Normally, it should remain in Omega^{k,j} = {y: lvar <= y <= uvar}
        if (a_k > 0)
        {
            for (int i = 0; i < n; ++i)
            {
                X.set(i, 0, X[i] + a_k * d_k[i]);
            }

            snapToBounds(X, lvar, uvar);

            // Update the set of active constraints
            for (int i = 0; i < n; ++i)
            { // We will correct if needed (see for example SolverTools.jl)
                if (std::abs(X.get(i, 0) - lvar.get(i, 0)) <= 1E-15)
                {
                    X.set(i, 0, lvar[i]);
                }
                if (std::abs(X.get(i, 0) - uvar.get(i, 0)) <= 1E-15)
                {
                    X.set(i, 0, uvar[i]);
                }
                active_l[i] = (X.get(i, 0) == lvar.get(i, 0));
                active_u[i] = (X.get(i, 0) == uvar.get(i, 0));
            }
            for (int i = 0; i < n; ++i)
            {
                active[i] = active_l[i] || active_u[i];
            }
            getModelGrad(&Grad, X, H, g);
        }

        // Compute the set of binding constraints
        binding_bounds(Grad, active_l, active_u, binding);

        bool areBindingActiveIdentical = true;
        for (int i = 0; i < n; ++i)
        {
            if (active[i] != binding[i])
            {
                areBindingActiveIdentical = false;
                break;
            }
        }

        // Reapply the conjugate gradient method with warm-start from dz with tightening
        // stopping criterion
        if (areBindingActiveIdentical)
        {
            xi = 1e-5;
            hasNegativeCurvature = conjugateGradient(dz, ZHZ, gz, 120, xi);
            d_k.fill(0);
            kfree = 0;
            for (int i = 0; i < n; ++i)
            {
                if (currentFree[i])
                {
                    d_k[i] = dz[kfree];
                    kfree++;
                }
            }

            if (d_k.has_nan()) {
                throw NOMAD::Exception(__FILE__, __LINE__, "d_k contains NaN");
            }

            const double a_max = max_step_bounds(X, lvar, uvar, d_k);
            if (hasNegativeCurvature)
            {
                // d_k is a direction of negative curvature: we find the smallest step gamma size such that
                // X + gamma d_k is at the boundary of Omega^{k,j} = {y : lvar <= y <= uvar}.
                // We need the objective function to decrease (to compute the predicted reduction in
                // the trust-region test), we then add an additionnal test.

                // Xp = X + a_max d_k
                Xp = d_k; Xp.multiply(a_max); Xp.add(X);
                if (getModelObj(Xp, H, g) <= getModelObj(X,H,g))
                {
                    // We accept the step.
                    a_k = a_max;
                }
                else
                {
                    // We perform a backtracking strategy, starting from a_max and decreasing along the line.
                    const double qm = getModelObj(X, H, g, g0);
                    double step = a_max;
                    const int nbMaxTrials = 10;
                    int niter = 0;
                    bool finished = false;
                    while (!finished)
                    {
                        Xp = d_k; Xp.multiply(step); Xp.add(X);
                        snapToBounds(Xp, lvar, uvar);
                        const double qXp = getModelObj(Xp, H, g);
                        Xp.sub(X);
                        const double slope = SGTELIB::Matrix::dot(Grad, Xp);
                        if (qXp <= qm + 1e-4 * slope)
                        {
                            finished = true;
                            continue;
                        }

                    step /= 3.0;
                    niter += 1;
                    finished = niter >= nbMaxTrials;
                    }
                    a_k = (getModelObj(Xp, H, g) <= getModelObj(X,H,g)) ? step : 0;
                }
            }
            else
            {
                // d_k is not a direction of negative curvature: we execute a projected Armijo linesearch.
                // 1e-15 is the smallest value in projected_armijo
                const double slope = SGTELIB::Matrix::dot(d_k, Grad);
                const double qm = getModelObj(X, H, g, g0);
                a_k = (a_max > 1e-15) ? projected_armijo(X, H, g, g0, lvar, uvar, d_k, qm, slope, armijo_Xp, Temp, a_max)
                    : 0;
                // a_k = std::min(a_k, a_max);
            }

            if (a_k > 0)
            {
                // Update the incumbent. Normally, it should remain in Omega^{k,j} = {y: lvar <= y <= uvar}
                for (int i = 0; i < n; ++i)
                {
                    X.set(i, 0, X[i] + a_k * d_k[i]);
                }
                snapToBounds(X, lvar, uvar);

                // Update the set of active constraints
                for (int i = 0; i < n; ++i)
                { // We will correct if needed (see for example SolverTools.jl)
                    if (std::abs(X.get(i, 0) - lvar.get(i, 0)) <= 1E-15)
                    {
                        X.set(i, 0, lvar[i]);
                    }
                    if (std::abs(X.get(i, 0) - uvar.get(i, 0)) <= 1E-15)
                    {
                        X.set(i, 0, uvar[i]);
                    }
                    active_l[i] = (X.get(i, 0) == lvar.get(i, 0));
                    active_u[i] = (X.get(i, 0) == uvar.get(i, 0));
                }
                for (int i = 0; i < n; ++i)
                {
                    active[i] = active_l[i] || active_u[i];
                }
                getModelGrad(&Grad, X, H, g);
            }
        }

        // When the set of working variables is empty, stop the algorithm.
        if (sum(active, n) == n)
        {
            OK = true;
            continue;
        }

        XmPXGrad = X;
        XmPXGrad.sub(Grad);
        snapToBounds(XmPXGrad, lvar, uvar);
        XmPXGrad.multiply(-1.0);
        XmPXGrad.add(X);

        const double fk = getModelObj(X, H, g, g0);

        // We consider there is no progress when the distance between the previous point
        // and the new iterate has not changed or the decreasing in the objective function
        // has stalled and the active set of constraints has not changed.
        Xm.sub(X);
        bool areActiveSetsChanged = false;
        for (int i = 0; i < n; ++i)
        {
            if ((active_l[i] != currentActiveL[i]) || (active_u[i] != currentActiveU[i]))
            {
                areActiveSetsChanged = true;
                break;
            }
        }
        OK = !areActiveSetsChanged && ((Xm.norm() <= 1e-9) || (std::abs(fk - qCurrent) <= 1e-9));

        k++;
        verbose && std::cout << "  k = " << k << " f(xk) = " << fk << " |d| = " << d_k.norm() << " a(amax) = " << a_k << " ( " << a_max << " )";
        verbose && std::cout << " |W| = " << n - sum(active, n) << " |L| = " << sum(active_l, n) << " |U| = " << sum(active_u, n) << " OK? " << OK << std::endl;
    }

    const bool success = getModelObj(X, H, g, g0) < f0;

    delete [] active_l;
    delete [] active_u;
    delete [] active;
    delete [] binding;

    return success;
}



bool NOMAD::QPSolverOptimize::check_subset_binding_update(
    bool* working,
    const bool* binding,
    const size_t n)
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
    const SGTELIB::Matrix& X,
    const SGTELIB::Matrix& H,
    const SGTELIB::Matrix& g,
    const double g0,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    const SGTELIB::Matrix& d,
    const double fk,
    const double slope,
    SGTELIB::Matrix& Xp,
    SGTELIB::Matrix& gradientF_kp,
    const double t_max)
{
    const int n = X.get_nb_rows();

    // A classical algorithm to compute a step-length satisfying Armijo/Wolfe conditions can be found at
    //
    // "Line Search Algorithms with Guaranteed Sufficient Decrease",
    // by J.J. More, D.J. Thuente, ACM Transactions on Mathematical Software, 20 (1994), Issue 3, pp 286–307
    //
    // Here, we use a procedure described in SolverTools.jl, which is an improved Armijo linesearch.

    // Linesearch parameters
    const double armijo_tol = 1E-4; // > 0
    const double t_small = 1E-15; // small
    const double t_decrease = 2.5; // > 1
    const double wolfe_tol = 0.9999; // < 1
    const int bW_max = 5; // non-negative integer
    // const int bA_max = 10; // non-negative integer
    const double t_increase = 5; // > 1

    bool good_grad = false; // true if gradient has been updated.

    // Check compatibility dimension
    lencheck(n, Xp);
    lencheck(n, gradientF_kp);

    // Initialization: starting steplength
    double tk = std::min(1.0, t_max);

    // Xp = P_Omega(X + tk d) where Omega = {y , lvar <= y <= uvar}
    Xp = d; Xp.multiply(tk); Xp.add(X);
    snapToBounds(Xp, lvar, uvar);

    double fkp = getModelObj(Xp, H, g, g0);
    getModelGrad(&gradientF_kp, Xp, H, g);
    double slope_t = SGTELIB::Matrix::dot(d, gradientF_kp);

    // First try to increase tk to satisfy Wolfe conditions.
    // NB: we do not enter often in this loop.
    int nbW = 0;
    bool wolfeCond = slope_t < wolfe_tol * slope; // NB: we use the strong Wolfe condition.
    bool armijoCond = fkp <= fk - armijo_tol * tk * std::abs(slope);
    while (wolfeCond && armijoCond && (nbW < bW_max) && (tk <= t_max))
    {
        tk *= t_increase;

        // Xp = P_Omega(X + tk d) where Omega = {y , lvar <= y <= uvar}
        Xp = d; Xp.multiply(tk); Xp.add(X);
        snapToBounds(Xp, lvar, uvar);

        // Recompute Wolfe and Armijo conditions.
        getModelGrad(&gradientF_kp, Xp, H, g);
        fkp = getModelObj(Xp, H, g, g0);
        slope_t = SGTELIB::Matrix::dot(d, gradientF_kp);
        wolfeCond = slope_t < wolfe_tol * slope;
        armijoCond = fkp <= fk - armijo_tol * tk * std::abs(slope);

        nbW ++;
        good_grad = true;
    }

    // Then try to satisfy Armijo's conditions.
    int nbA = 0;
    armijoCond = fkp <= fk - armijo_tol * tk * std::abs(slope);
    while (!armijoCond && tk > t_small)
    {
        tk /= t_decrease;
        Xp = d; Xp.multiply(tk); Xp.add(X); // Xp = X + t_k d
        snapToBounds(Xp, lvar, uvar);
        fkp = getModelObj(Xp, H, g, g0);

        armijoCond = fkp <= fk - armijo_tol * tk * fabs(slope);
        nbA ++;
    }

    if (!armijoCond)
    {
        return 0.0;
    }

    return tk;
}

// Solve method with outer loop and inner loop
bool NOMAD::QPSolverOptimize::solveL1AugLag(
    NOMAD::Point& X_k,
    const int max_iter,
    const double atol,
    const double rtol,
    const bool verbose)
{
    
    // Compute stopping tolerance
    SGTELIB::Matrix gradientLag_k("gradientLag_k", _n, 1);
    getModelGrad(&gradientLag_k, X_k);
    double ng = gradientLag_k.norm();
    const double ng0 = ng;
    const double tol = atol + ng0 * rtol;
    verbose && std::cout << "Start solveL1AugLag with tol = " << tol << std::endl;

    // Compute initial starting point and bounds
    SGTELIB::Matrix lvar("lvar", _n, 1);
    lvar.fill(0.0);
    SGTELIB::Matrix uvar("uvar", _n, 1);
    uvar.fill(INF);
    SGTELIB::Matrix XS("XS", _n, 1);
    for (int i=0; i < _n; ++i)
    {
        const double lb = ((_modelLowerBound[i].isDefined())? _modelLowerBound[i].todouble(): NOMAD::M_INF);
        const double ub = ((_modelUpperBound[i].isDefined())? _modelUpperBound[i].todouble(): NOMAD::INF);
        lvar.set(i, 0, lb);
        uvar.set(i, 0, ub);
        XS.set(i, 0, X_k[i].todouble());
    }

    // Check bound compatibility and feasibility
    snapToBounds(XS, lvar, uvar);
    for (int i = 0; i < _n; ++i)
    {
        const bool areBoundsCompatible = (lvar.get(i, 0) < uvar.get(i, 0));
        if (!areBoundsCompatible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solveL1AugLag assertion error: Error compatibility lower and upper bound");
        }

        const bool feasible = (XS.get(i, 0) >= lvar.get(i, 0)) && (XS.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            std::cout << XS.get(i, 0) - lvar.get(i, 0) << " " << uvar.get(i, 0) - XS.get(i, 0) << std::endl;
            std::cout << "Error compatibility lower and upper bound" << std::endl;
            throw NOMAD::Exception(__FILE__, __LINE__, "solveL1AugLag assertion error: Error XS is not feasible");
        }
    }

    // Particular case : when the gap between lvar and uvar is too small, we immediately leave the procedure.
    for (int i = 0; i < _n; ++i)
    {
        if (std::abs(lvar.get(i, 0) - uvar.get(i, 0)) <= 1e-8)
        {
            return false;
        }
    }

    // Allocate memory for the set of active, strictly feasible constraints
    // and infeasible constraints.
    auto active = new bool [_nbCons]; // true for indices of active constraints
    auto feasible = new bool [_nbCons]; // true for indices of strictly feasible constraints
    auto infeasible = new bool [_nbCons]; // true for indices of infeasible constraints
    auto active_lb = new bool [_n]; // true when X_k is close to its lower bounds
    auto active_ub = new bool [_n]; // true when X_k is close to its upper bounds

    // Allocate memory for inner loop iterates
    NOMAD::Point X_km1(X_k), X_om1(X_k), X_can(X_k);
    SGTELIB::Matrix Jacobian_k("Jacobian_k", _nbCons, _n);
    SGTELIB::Matrix h_k("h_k", _n, 1); // horizontal step
    SGTELIB::Matrix v_k("v_k", _n, 1); // vertical step
    SGTELIB::Matrix pseudoGradient_k("pseudoGradient_k", _n, 1); // pseudo-gradient
    SGTELIB::Matrix multiplier_k("multiplier_k", _nbCons, 1); // multipliers to compute pseudo-gradient
    SGTELIB::Matrix ZtransposePseudoGrad_k;
    SGTELIB::Matrix activeMultiplier_k;
    SGTELIB::Matrix innerTolBounds("innerTolBounds", _n, 1); // inner tolerance for bound constraints

    // Outer loop augmented Lagrangian multipliers.
    SGTELIB::Matrix lambda_l("lambda_l", _nbCons, 1);
    lambda_l.fill(0.0);
    double mu_l = 1.0;
    double eta_l = 1.0;
    double omega_l = 1.0;

    // Initialize outer loop parameters
    size_t iterOuterLoop = 0;
    double distXOuterLoop = NOMAD::INF;
    const size_t maxIterInnerLoop = 20, maxIterOuterLoop = 10;
    // const double toleranceF = 1e-12;
    const double tolerance_distDX = 1e-12;
    size_t quadModelNbEval = 0;

    // Initialize outerSuccess and outerFailure conditions
    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X_k);
    gradientLag_k = getModelLagGradient(X_k, lambda_l);
    double ngproj = check_optimality_bounds(X_k, gradientLag_k);
    bool outerSuccess = (ngproj <= tol) && isFeasible(cons, tol);
    bool outerFailure = (distXOuterLoop <= tolerance_distDX) || (iterOuterLoop >= maxIterOuterLoop) || (quadModelNbEval >= _quadModelMaxEval);

    double fk = getModelObj(X_k);
    double F_k = getPenalizedL1AugLagModelObj(X_k, cons, lambda_l, mu_l);
    double F_km1 = F_k;
    verbose && std::cout << " |grad| = " << ng << " |Proj(x - grad) - x| = " << ngproj << " P = " << F_k << " f = " << fk << " |c| = " << cons.norm() << std::endl;

    // Outer iteration
    while (!outerFailure && !outerSuccess)
    {
        X_om1 = X_k;

        // Initialize inner loop parameters
        size_t iterInnerLoop = 0;
        double innerPrecision = 1.0;
        double innerTolerance = 1.0;
        for (int i = 0; i < _n; ++i)
        {
            if ((lvar.get(i, 0) == NOMAD::M_INF) || (uvar.get(i, 0) == NOMAD::INF))
            {
                innerTolBounds.set(i, 0, 0.1); // 10 times lower than innerPrecision
            }
            else
            {
                innerTolBounds.set(i, 0, std::max(0.1 * (uvar.get(i, 0) - lvar.get(i, 0)), 1e-8));
            }
        }

        // Adjust tolerance for bound constraints
        bool areActiveBoundConstraintsAboveTol = check_active_bound_constraints(X_k, active_lb, active_ub,
                                                                                innerTolBounds, lvar, uvar);
        // No need to continue, the algorithm is stuck in a ``corner'' of the feasible space
        // or has reached the given tolerance
        if (!areActiveBoundConstraintsAboveTol)
        {
            break;
        }

        // Compute sets of constraints for the inner loop
        getModelCons(&cons, X_k);

        // Adjust tolerance for constraints
        int nbActiveBounds = 0;
        for (int i = 0; i < _n; ++i)
        {
            if (active_lb[i] || active_ub[i])
            {
                nbActiveBounds += 1;
            }
        }
        bool areConstraintsAboveTol = check_active_constraints(cons, nbActiveBounds,
                                                               active, feasible, infeasible, innerPrecision);
        // No need to continue, too many constraints are active
        if (!areConstraintsAboveTol)
        {
            break;
        }
        innerTolerance = innerPrecision;

        // Stopping criteria conditions for the inner loop

        // Inner success
        // 1- Compute the pseudo-gradient of the penalized augmented lagrangian
        multiplier_k = get_pseudo_multiplier(active, feasible, infeasible, lambda_l, mu_l);
        pseudoGradient_k = getModelLagGradient(X_k, multiplier_k);

        // 2- Compute a nullspace matrix for active Jacobian constraints. It is the identity matrix if
        // the set of active constraints is empty.
        int nbActive = sum(active, _nbCons);
        bool areActiveConstraints = nbActive + nbActiveBounds > 0;
        Jacobian_k = getModelJacobian(X_k);
        SGTELIB::Matrix activeJacobian_k = getModelActiveJacobian(Jacobian_k, active);
        // Complete with bound constraints
        for (int i = 0; i < _n; ++i)
        {
            SGTELIB::Matrix gradActiveBoundCst("grad", 1, _n);
            gradActiveBoundCst.fill(0.0);
            if (active_lb[i])
            {
                gradActiveBoundCst.set(0, i, -1.0);
                activeJacobian_k.add_rows(gradActiveBoundCst);
            }
            if (active_ub[i])
            {
                gradActiveBoundCst.set(0, i, 1.0);
                activeJacobian_k.add_rows(gradActiveBoundCst);
            }
        }
        SGTELIB::Matrix Zk = activeJacobian_k.null_space();
        if (areActiveConstraints)
        {
            ZtransposePseudoGrad_k = SGTELIB::Matrix::product(Zk.transpose(), pseudoGradient_k);
        }
        else
        {
            ZtransposePseudoGrad_k = pseudoGradient_k;
        }

        // Inner iteration
        bool innerSuccess = ZtransposePseudoGrad_k.norm() <= tol && isFeasible(cons, tol);
        bool innerFailure = false;
        while (!innerFailure && !innerSuccess)
        {
            F_km1 = F_k;
            X_km1 = X_k;
            X_can = X_k;

            // Compute a horizontal step (minimize the change in the penalty function subject to active
            // constraints)
            compute_horizontal_step(X_k, h_k, Jacobian_k, active, feasible, infeasible,
                                                                   active_lb, active_ub, lambda_l, mu_l);

            bool computeStrengthenedStep = true;
            bool computeActiveSets = true;
            if (!areActiveConstraints || ZtransposePseudoGrad_k.norm() > innerTolerance)
            {
                // Apply a linesearch along h_k
                const double beta_k = piecewise_line_search(X_k, h_k, active , feasible, infeasible,
                                                            active_lb, active_ub, lvar, uvar, lambda_l, mu_l);
                update(X_k, h_k, beta_k);
                computeStrengthenedStep = false;
            }
            else
            {
                // Compute active multipliers
                activeMultiplier_k = SGTELIB::Surrogate_PRS::compute_multiplier(pseudoGradient_k, activeJacobian_k);

                // Compute a drop constraint step if possible
                SGTELIB::Matrix d;
                const bool dropConstraintSuccess = compute_drop_constraint_step(d, activeJacobian_k, activeMultiplier_k, pseudoGradient_k, mu_l);

                if (dropConstraintSuccess)
                {
                    // Apply a linesearch along direction d
                    const double beta_k = piecewise_line_search(X_k, d, active , feasible, infeasible,
                                                                active_lb, active_ub, lvar, uvar, lambda_l, mu_l);
                    update(X_k, h_k, beta_k);
                    computeStrengthenedStep = false;
                }
                else
                {
                    update(X_can, h_k); // X_can := X_k + h_k

                    // Increase feasibility by computing a vertical step. It is computed using active constraint
                    // gradients at X_k (and NOT X_k + h_k).
                    compute_vertical_step(X_k, v_k, activeJacobian_k, cons, active,
                                                                       active_lb, active_ub, lvar, uvar);
                    update(X_can, v_k); // X_can := X_k + h_k + v_k

                    getModelCons(&cons, X_can);
                    F_k = getPenalizedL1AugLagModelObj(X_can, cons, lambda_l, mu_l);
                    verbose && std::cout << " V: (l=" << iterInnerLoop << ") |v| = " << v_k.norm() << " Pk = " << F_k << " ||c|| = " << cons.norm() << std::endl;
                    // const bool decreasePenalty = F_k < F_km1; // This criterion is the simpler
                    const bool decreasePenalty = F_k - F_km1 <= - 0.01 * (activeJacobian_k.norm() + std::pow(ZtransposePseudoGrad_k.norm(), 2));
                    if (decreasePenalty)
                    {
                        // Accept the candidate
                        X_k = X_can;
                        computeStrengthenedStep = false;
                        computeActiveSets = false;
                    }
                }
            }

            if (computeStrengthenedStep)
            {
                // Recompute the set of active, feasible and infeasible constraints with a lower threshold.
                innerPrecision = innerPrecision / 2;
                innerTolerance = innerTolerance / 2;
                getModelCons(&cons, X_k);
                getModelActiveCons(cons, innerPrecision, active);
                getModelFeasibleCons(cons, innerPrecision, feasible);
                getModelInfeasibleCons(cons, innerPrecision, infeasible);

                // Recompute the set of active bound constraints
                for (int i = 0; i < _n; ++i)
                {
                    innerTolBounds[i] /= 2.0;
                    innerTolBounds[i] = std::max(innerTolBounds[i], 1e-8);
                }
                for (int i = 0; i < _n; ++i)
                {
                    active_lb[i] = std::abs(X_k[i].todouble() - lvar[i]) <= innerTolBounds[i];
                    active_ub[i] = std::abs(X_k[i].todouble() - uvar[i]) <= innerTolBounds[i];
                }

                // Recompute a horizontal step
                compute_horizontal_step(X_k, h_k, Jacobian_k, active, feasible, infeasible,
                                        active_lb, active_ub, lambda_l, mu_l);
                const double beta_k = piecewise_line_search(X_k, h_k, active , feasible, infeasible,
                                                            active_lb, active_ub, lvar, uvar, lambda_l, mu_l);
                update(X_k, h_k, beta_k);
            }

            X_k.snapToBounds(_modelLowerBound, _modelUpperBound);

            // Check for success
            getModelCons(&cons, X_k);
            if (computeActiveSets)
            {
                // 1- Recompute the different set of constraints used in the inner loop.
                getModelActiveCons(cons, innerPrecision, active);
                getModelFeasibleCons(cons, innerPrecision, feasible);
                getModelInfeasibleCons(cons, innerPrecision, infeasible);

                areActiveBoundConstraintsAboveTol = check_active_bound_constraints(X_k, active_lb, active_ub,
                                                                                   innerTolBounds, lvar, uvar);
                // No need to continue, the algorithm is stuck in a ``corner'' of the feasible space
                // or has reached the given tolerance
                if (!areActiveBoundConstraintsAboveTol)
                {
                    break;
                }

                nbActiveBounds = 0;
                for (int i = 0; i < _n; ++i)
                {
                    if (active_lb[i] || active_ub[i])
                    {
                        nbActiveBounds += 1;
                    }
                }

                areConstraintsAboveTol = check_active_constraints(cons, nbActiveBounds, active, feasible, infeasible,
                                                                  innerPrecision);
                // No need to continue, all constraints are below tolerance
                if (!areConstraintsAboveTol)
                {
                    break;
                }
                innerTolerance = innerPrecision;
            }

            // 2- Compute stopping criterion based on pseudo-gradient
            multiplier_k = get_pseudo_multiplier(active, feasible, infeasible, lambda_l, mu_l);
            pseudoGradient_k = getModelLagGradient(X_k, multiplier_k);

            nbActive = sum(active, _nbCons);
            nbActiveBounds = 0;
            for (int i = 0; i < _n; ++i)
            {
                if (active_lb[i] || active_ub[i])
                {
                    nbActiveBounds += 1;
                }
            }
            areActiveConstraints = nbActive + nbActiveBounds > 0;
            Jacobian_k = getModelJacobian(X_k);
            activeJacobian_k = getModelActiveJacobian(Jacobian_k, active);
            // Complete with bound constraints
            for (int i = 0; i < _n; ++i)
            {
                SGTELIB::Matrix gradActiveBoundCst("grad", 1, _n);
                gradActiveBoundCst.fill(0.0);
                if (active_lb[i])
                {
                    gradActiveBoundCst.set(0, i, -1.0);
                    activeJacobian_k.add_rows(gradActiveBoundCst);
                }
                if (active_ub[i])
                {
                    gradActiveBoundCst.set(0, i, 1.0);
                    activeJacobian_k.add_rows(gradActiveBoundCst);
                }
            }
            Zk = activeJacobian_k.null_space();
            if (areActiveConstraints)
            {
                ZtransposePseudoGrad_k = SGTELIB::Matrix::product(Zk.transpose(), pseudoGradient_k);
            }
            else
            {
                ZtransposePseudoGrad_k = pseudoGradient_k;
            }

            innerSuccess = ZtransposePseudoGrad_k.norm() <= tol && isFeasible(cons, tol);
            if (areActiveConstraints && innerSuccess)
            {
                activeMultiplier_k = SGTELIB::Surrogate_PRS::compute_multiplier(pseudoGradient_k, activeJacobian_k);

                for (int i = 0; i < activeMultiplier_k.get_nb_cols(); ++i)
                {
                    if (activeMultiplier_k[i] < 0 || activeMultiplier_k[i] > (1.0 / mu_l))
                    {
                        innerSuccess = false;
                        break;
                    }
                }
            }

            // Check for failure
            const double distXkXkm1 = NOMAD::Point::dist(X_k,X_km1).todouble();
            // normGradLag_k = gradLag_k.norm();
            F_k = getPenalizedL1AugLagModelObj(X_k, cons, lambda_l, mu_l);
            const double deltaF = std::abs(F_k - F_km1);
            // double dual_norm = 0;
            verbose && std::cout << " Inner: (l=" << iterInnerLoop << ") Pl = " << F_k << " |c| = " << cons.norm() << " |L| = " << ZtransposePseudoGrad_k.norm() << " |dF| = " << deltaF;
            verbose && std::cout << " inner precision = " << innerPrecision << " bounds precision (inf) = " << innerTolBounds.norm_inf() << std::endl;
            verbose && std::cout << " |cons|_epsilon = " << nbActive << " nb bds active = " << nbActiveBounds << std::endl;

            iterInnerLoop++;
            quadModelNbEval++;
            // innerFailure = verticalSuccess && horizontalSuccess;
            innerFailure = (deltaF <= 1e-6) || (distXkXkm1 <= tolerance_distDX) || (iterInnerLoop >= maxIterInnerLoop) || (quadModelNbEval >= _quadModelMaxEval);

        } // end of inner loop: Xk has been updated

        getModelCons(&cons, X_k);
        getModelActiveCons(cons, tol, active);
        getModelFeasibleCons(cons, tol, feasible);
        getModelInfeasibleCons(cons, tol, infeasible);

        multiplier_k = get_pseudo_multiplier(active, feasible, infeasible, lambda_l, mu_l);
        pseudoGradient_k = getModelLagGradient(X_k, multiplier_k);

        F_k = getPenalizedL1AugLagModelObj(X_k, cons, lambda_l, mu_l);
        bool unbounded_subpb = (F_k < - 1/tol);

        // Update Lagrange multipliers
        if (isFeasible(cons, eta_l))
        { // no update of mu_l

            // Compute active multipliers
            Jacobian_k = getModelJacobian(X_k);
            activeJacobian_k = getModelActiveJacobian(Jacobian_k, active);
            activeMultiplier_k = SGTELIB::Surrogate_PRS::compute_multiplier(pseudoGradient_k, activeJacobian_k);

            int activeCstId = 0;
            for (int i = 0; i < _nbCons ; i++)
            {
                const double li = lambda_l.get(i, 0);
                if (std::abs(cons.get(i, 0)) <= tol)
                {
                    active[i] = true;
                    feasible[i] = false;
                    infeasible[i] = false;
                    lambda_l.set(i, 0, li + activeMultiplier_k[activeCstId]);
                    activeCstId++;
                }
                else if (cons.get(i, 0) < -tol)
                {
                    active[i] = false;
                    feasible[i] = true;
                    infeasible[i] = false;
                    lambda_l.set(i, 0, li); // no update of lambda_l
                }
                else
                {
                    active[i] = false;
                    feasible[i] = false;
                    infeasible[i] = true;
                    lambda_l.set(i, 0, li - 1 / mu_l);
                }
            }
            eta_l = eta_l * pow(mu_l, 0.9);
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

        iterOuterLoop++; // l = l + 1;
        fk = getModelObj(X_k);
        verbose && std::cout << "k = " << iterOuterLoop <<  " |Proj(x - grad) - x| = " << ngproj  << " f(x) = " << fk << " |c(x)| = " << cons.norm() << std::endl;

        // Check failure
        for (int i = 0; i < _n; ++i)
        {
            active_lb[i] = std::abs(X_k[i].todouble() - lvar[i]) <= 1e-8;
            active_ub[i] = std::abs(X_k[i].todouble() - uvar[i]) <= 1e-8;
        }
        int nbActiveLb = sum(active_lb, _n);
        int nbActiveUb = sum(active_ub, _n);

        distXOuterLoop = NOMAD::Point::dist(X_k, X_om1).todouble();
        outerFailure = (!unbounded_subpb && (distXOuterLoop <= tolerance_distDX)) || (nbActiveLb + nbActiveUb) >= _n;
        outerFailure = outerFailure || (iterOuterLoop >= maxIterOuterLoop);
        if (outerFailure)
        {
            verbose && std::cout << "Early stop: |d| = " << distXOuterLoop << " inner? " << innerFailure << " unbounded? " << unbounded_subpb << std::endl;
        }

    }

    delete [] active;
    delete [] feasible;
    delete [] infeasible;
    delete [] active_lb;
    delete [] active_ub;

    return outerSuccess;
}

bool NOMAD::QPSolverOptimize::check_active_bound_constraints(const NOMAD::Point &X,
                                                             bool *active_lb, bool *active_ub,
                                                             SGTELIB::Matrix &tolBounds,
                                                             const SGTELIB::Matrix &lvar,
                                                             const SGTELIB::Matrix &uvar) const
{
    // Before starting, detect if upper and lower bound constraints are active for the same variable
    // with the lowest tolerance
    for (int i = 0; i < _n; ++i)
    {
        if (std::abs(X[i].todouble() - lvar[i]) <= 1e-8 &&
            std::abs(X[i].todouble() - uvar[i]) <= 1e-8)
        {
            return false;
        }
    }

    // Update active bound sets
    for (int i = 0; i < _n; ++i)
    {
        active_lb[i] = std::abs(X[i].todouble() - lvar[i]) <= tolBounds[i];
        active_ub[i] = std::abs(X[i].todouble() - uvar[i]) <= tolBounds[i];
    }

    // Special check: further reduce the tolerance of active bound constraints to prevent
    // the nullspace matrix derived from the active Jacobian constraints from being empty.
    int nbActiveLb = sum(active_lb, _n);
    int nbActiveUb = sum(active_ub, _n);
    while (nbActiveLb + nbActiveUb >= _n)
    {
        // 1- Detect if all bound constraints are active (which means either active_lb is true
        // or active_ub is true for the smallest tolerance)
        bool areActiveBoundConstraintsBelowTol = true;
        for (int i = 0; i < _n; ++i)
        {
            if (std::abs(X[i].todouble() - lvar[i]) > 1e-8 &&
                std::abs(X[i].todouble() - uvar[i]) > 1e-8)
            {
                areActiveBoundConstraintsBelowTol = false;
                break;
            }
        }

        // No need to continue, the algorithm is stuck in a ``corner'' of the feasible space
        // or has reached the given tolerance
        if (areActiveBoundConstraintsBelowTol)
        {
            break;
        }

        // 2- Decrease the inner tolerance for all bound constraints
        for (int i = 0; i < _n; ++i)
        {
            tolBounds[i] /= 2.0;
            tolBounds[i] = std::max(tolBounds[i], 1e-8);
        }

        // 3- Recompute the set of active bound constraints
        for (int i = 0; i < _n; ++i)
        {
            active_lb[i] = std::abs(X[i].todouble() - lvar[i]) <= tolBounds[i];
            active_ub[i] = std::abs(X[i].todouble() - uvar[i]) <= tolBounds[i];
        }
        nbActiveLb = sum(active_lb, _n);
        nbActiveUb = sum(active_ub, _n);
    }
    return nbActiveUb + nbActiveLb < _n;
}

bool NOMAD::QPSolverOptimize::check_active_constraints(const SGTELIB::Matrix& cons,
                                                       const int nbActiveBounds,
                                                       bool* active,
                                                       bool* feasible,
                                                       bool* infeasible,
                                                       double& innerPrecision) const
{
    getModelActiveCons(cons, innerPrecision, active);
    getModelFeasibleCons(cons, innerPrecision, feasible);
    getModelInfeasibleCons(cons, innerPrecision, infeasible);
    int nbActive = sum(active, _nbCons);
    while (nbActiveBounds + nbActive >= _n)
    {
        // Detect if all active constraints are below tolerance
        bool areConstraintsBelowTol = true;
        for (int i = 0; i < _nbCons; ++i)
        {
            if (active[i] && std::abs(cons.get(i, 0)) > 1e-5)
            {
                areConstraintsBelowTol = false;
                break;
            }
        }
        if (areConstraintsBelowTol)
        {
            return false;
        }

        // Reduce innerPrecision and recompute set of constraints
        innerPrecision /= 2.0;
        innerPrecision = std::max(innerPrecision, 1e-5);
        getModelActiveCons(cons, innerPrecision, active);
        getModelFeasibleCons(cons, innerPrecision, feasible);
        getModelInfeasibleCons(cons, innerPrecision, infeasible);

        nbActive = sum(active, _nbCons);
    }
    return true;
}

bool NOMAD::QPSolverOptimize::compute_drop_constraint_step(
    SGTELIB::Matrix& d,
    const SGTELIB::Matrix& activeJacobian,
    const SGTELIB::Matrix& activeMultipliers,
    const SGTELIB::Matrix& pseudoGradient,
    const double mu) const
{
    // 1- Check if there exists potential constraints to drop
    int cstToDropId = -1;
    for (int i = 0; i < activeMultipliers.get_nb_cols(); ++i)
    {
        if (activeMultipliers[i] < 0 || activeMultipliers[i] > (1.0 / mu))
        {
            cstToDropId = i;
            break;
        }
    }
    if (cstToDropId < 0)
        return false;

    // There exists a potential constraint to drop: compute the corresponding new direction.
    // 1- Compute a null space matrix of a subset of the active jacobian
    SGTELIB::Matrix activeJacobian_mj("activeJacobian_mj",
                                       activeJacobian.get_nb_rows()-1,
                                       activeJacobian.get_nb_cols());
    int k = 0;
    for (int i = 0; i < activeJacobian.get_nb_rows(); i++)
    {
        if (i != cstToDropId)
        {
            for (int j = 0; j < activeJacobian.get_nb_cols(); j++)
            {
                activeJacobian_mj.set(k, j, activeJacobian.get(i,j));
            }
            k++;
        }
    }
    SGTELIB::Matrix Zmj = activeJacobian_mj.null_space();

    // 2- Compute new direction d := sigma_j Zmj Zmj^T activeJacobian[j]
    // where sigma_j = -sign(activeMultipliers[j]) with j = cstToDropId
    SGTELIB::Matrix gradConsToDrop = activeJacobian.get_row(cstToDropId).transpose();
    d = SGTELIB::Matrix::product(Zmj, Zmj.transpose(), gradConsToDrop);
    const double signConsToDrop = activeMultipliers[cstToDropId] > 0 ? 1.0 : -1.0;
    d.multiply(-signConsToDrop);

    // 3- Check if it can be used as a descent direction
    gradConsToDrop.multiply(std::min(signConsToDrop, 0.0));
    gradConsToDrop.add(pseudoGradient);
    const double slope = SGTELIB::Matrix::dot(d, gradConsToDrop);
    if (slope < -0.05)
        return true;
    else
        return false;
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
    const double mu,
    const bool * active,
    const bool * infeasible) const
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
    const Point& X,
    const SGTELIB::Matrix& d,
    const bool* active,
    const bool* feasible,
    const bool* infeasible,
    const bool* active_lb,
    const bool* active_ub,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    const SGTELIB::Matrix& lambda,
    const double mu,
    const double small_gamma, // = 1E-20
    const double gamma_update, // = 1.5
    const double delta /* = 1E-4 // Pk < (P0 - delta) */ ) const
{

    
    
    // Allocations
    // const int nbActive = sum(active, ncon);
    // double slope;
    NOMAD::Point X_k(X);

    SGTELIB::Matrix multiplier_k = get_pseudo_multiplier(active, feasible, infeasible, lambda, mu);
    SGTELIB::Matrix pseudoGradient = getModelLagGradient(X, multiplier_k);

    // The first _nbCons coordinates correspond to the cons constraints,
    // the next _n coordinates correspond to lower bound constraints and
    // the final _n correspond to upper bound constraints
    SGTELIB::Matrix gamma("gamma", _nbCons + 2 * _n, 1);
    gamma.fill(0.0);
    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X);
    SGTELIB::Matrix Jacobian = getModelJacobian(X);
    SGTELIB::Matrix jprod = SGTELIB::Matrix::product(Jacobian, d);

    double ak = SGTELIB::Matrix::dot(d, pseudoGradient); // < 0
    if (ak >= 0)
    {
        return 0.0;
    }

    // Step 1: compute the Ik and gamma sets.
    bool* Ik = new bool [_nbCons + 2 * _n];
    for (int i = 0; i < _nbCons; ++i)
    {
        if (!active[i])
        {
            gamma[i] = -cons[i] / jprod[i];
        }
        Ik[i] = (gamma[i] > 0) && (!active[i]);
    }
    for (int i = _nbCons; i < _nbCons + _n; ++i)
    {
        if (!active_lb[i - _nbCons])
        {
            if (lvar[i - _nbCons] == NOMAD::M_INF)
            {
                gamma[i] = d[i - _nbCons] >= 0 ? NOMAD::M_INF: NOMAD::INF;
            }
            else
            {
                gamma[i] = d[i - _nbCons] == 0 ? 0 // We cannot move towards this direction when the direction is empty
                : - (lvar[i - _nbCons] - X[i - _nbCons].todouble()) / (- d[i - _nbCons]);
            }
        }
        Ik[i] = (gamma[i] > 0) && (!active_lb[i-_nbCons]);
        if (!active_ub[i - _nbCons])
        {
            if (uvar[i - _nbCons] == NOMAD::INF)
            {
                gamma[i + _n] = d[i - _nbCons] >= 0 ? NOMAD::INF: NOMAD::M_INF;
            }
            else
            {
                gamma[i + _n] = d[i - _nbCons] == 0 ? 0 // We cannot move towards this direction when the direction is empty
                : - (-uvar[i - _nbCons] + X[i - _nbCons].todouble()) / (d[i - _nbCons]);
            }
        }
        Ik[i + _n] = (gamma[i + _n] > 0) && (!active_ub[i - _nbCons]);
    }

    // Step 2: if Ik is empty, return a stepsize of length 0 (we are in a nonlinear setting).
    bool unbounded = (sum(Ik, _nbCons + 2 * _n) == 0);
    if (unbounded)
    {
        return 0.0;
    }

    // Step 3: determine stepsize gamma_lk
    double gamma_lk = 1;
    int k = 0;
    const int Iksize = sum(Ik, _nbCons + 2 * _n);
    while (k < Iksize)
    {
        // Find lk such that gamma_lk <= gamma[i], for all i in Ik
        gamma_lk = INF;
        int lk = -1;
        for (int i = 0; i < _nbCons + 2 * _n; ++i)
        {
            if (Ik[i] && gamma[i] <= gamma_lk) {
                lk = i;
                gamma_lk = gamma[i];
            }
        }

        // Update ak
        if (lk < _nbCons) {
            ak += fabs(jprod[lk]);
        } else if (lk < _nbCons + _n) {
            ak += std::abs(d[lk - _nbCons]);
        } else {
            ak += std::abs(d[lk - _nbCons - _n]);
        }

        if (ak >= 0)
            break;
        Ik[lk] = false;

        k++;
    }

    // Step 5:
    for (int i = 0 ; i < d.get_nb_rows() ; i++)
    {
        X_k[i] = X[i] +  gamma_lk * d.get(i, 0);
    }
    
    const double P0 = getPenalizedL1AugLagModelObj(X, cons, lambda, mu);
    getModelCons(&cons, X_k);
    double Pk = getPenalizedL1AugLagModelObj(X_k, cons, lambda, mu);
    bool OK = Pk < (P0 - delta);
    while (!OK)
    {
        // Normally, we should do a cubic interpolation, but we choose
        // to apply an Armijo linesearch instead.
        gamma_lk /= gamma_update;
        for (int i = 0 ; i < d.get_nb_rows() ; i++)
        {
            X_k[i] = X[i] +  gamma_lk * d.get(i, 0);
        }
        getModelCons(&cons, X_k);
        Pk = getPenalizedL1AugLagModelObj(X_k, cons, lambda, mu);
        OK = (Pk < (P0 - delta)) || (gamma_lk <= small_gamma);
    }
    delete [] Ik;

    return gamma_lk;
}

bool NOMAD::QPSolverOptimize::compute_vertical_step(
    const Point& X,
    SGTELIB::Matrix& v_k,
    const SGTELIB::Matrix& activeJacobian_k,
    const SGTELIB::Matrix& cons,
    const bool* active,
    const bool* active_lb,
    const bool* active_ub,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar) const
{
    const int ncon = static_cast<int>(_nbCons);
    const int nbActive = activeJacobian_k.get_nb_rows();
    SGTELIB::Matrix activeCons("activeCons", nbActive, 1);
    int k = 0;
    for (int i = 0; i < ncon; i++)
    {
        if (active[i])
        {
            activeCons.set(k, 0, -cons.get(i,0));
            k++;
        }
    }
    for (int i = ncon; i < ncon + _n; ++i)
    {
        if (active_lb[i - ncon])
        {
            activeCons.set(k, 0, -X[i - ncon].todouble() + lvar.get(i - ncon,0));
            k++;
        }
        if (active_ub[i - ncon])
        {
            activeCons.set(k, 0, X[i - ncon].todouble() - uvar.get(i - ncon,0));
            k++;
        }
    }
    if (nbActive != k)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Active jacobian number of rows do not match active indices.");
    }

    // v_k is computed as: v_k := - A_k (A_kt A_k)^-1 activeCons
    // where A_k = activeJacobian_k.
    //
    // NB: One could also follow the procedure given in Section 4.1 of
    // "Nonlinear programming via an exact penalty function: Global analysis"
    // by T.F. Coleman and A.R. Conn,  Mathematical Programming 24, 137–161 (1982).
    // https://doi.org/10.1007/BF01585101
    v_k = SGTELIB::Matrix::solve_least_squares_SVD(activeJacobian_k, activeCons);
    return true;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::get_pseudo_multiplier(
    const bool* active,
    const bool* feasible,
    const bool* infeasible,
    const SGTELIB::Matrix& lambda,
    const double mu) const
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
    const Point& X,
    SGTELIB::Matrix& h_k,
    const SGTELIB::Matrix& Jacobian_k,
    const bool* active,
    const bool* feasible,
    const bool* infeasible,
    const bool* active_lb,
    const bool* active_ub,
    const SGTELIB::Matrix& lambda_l,
    const double mu_l) const
{
    
    // Compute Z such that activeJacobian_k Z = 0; and Zt Z = I
    SGTELIB::Matrix activeJacobian_k = getModelActiveJacobian(Jacobian_k, active);
    for (int i = 0; i < _n; ++i)
    {
        SGTELIB::Matrix gradActiveBoundCst("grad", 1, _n);
        gradActiveBoundCst.fill(0.0);
        if (active_lb[i])
        {
            gradActiveBoundCst.set(0, i, -1.0);
            activeJacobian_k.add_rows(gradActiveBoundCst);
        }
        if (active_ub[i])
        {
            gradActiveBoundCst.set(0, i, 1.0);
            activeJacobian_k.add_rows(gradActiveBoundCst);
        }
    }
    const SGTELIB::Matrix Z = activeJacobian_k.null_space();
    const SGTELIB::Matrix multiplier_k = get_pseudo_multiplier(active, feasible, infeasible, lambda_l, mu_l);

    // Compute Zt L Z, where L is the Hessian of the penalized augmented Lagrangian.
    const SGTELIB::Matrix Hlag_k = getModelLagHessian(X, multiplier_k);
    const SGTELIB::Matrix ZLZ = SGTELIB::Matrix::product(Z.transpose(), Hlag_k, Z);

    // Compute -Zt g, where g is the gradient of the penalized augmented Lagrangian.
    const SGTELIB::Matrix Glag_k = getModelLagGradient(X, multiplier_k);
    SGTELIB::Matrix ZL = SGTELIB::Matrix::product(Z.transpose(), Glag_k);
    ZL.multiply(-1);

    // Solve (Zt L Z) w = - Zt g
    const SGTELIB::Matrix invZLZ = ZLZ.SVD_inverse();
    SGTELIB::Matrix delta_k = SGTELIB::Matrix::product(invZLZ, ZL);

    // h := Z w*, where w is the solution of the linear system above.
    h_k = SGTELIB::Matrix::product(Z, delta_k);

    // if h is not a descent direction, set h:= -Z Zt g, hence h:= -g.
    const double slope = SGTELIB::Matrix::dot(h_k, Glag_k);
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
    NOMAD::Point& X,
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
    const size_t maxSuccessiveFail // = 3;
   )
{
    // We use slack variables here.
    const int nbVar = _n + _nbCons;

    // Fix tolerance.
    SGTELIB::Matrix Gk("Gk", nbVar, 1);
    getModelGrad(&Gk, X);
    double ng = Gk.norm();
    const double ng0 = ng;
    const double tol = atol + ng0 * rtol;

    // Keep initial f and load constraints values before running the algorithm.
    const double f0 = getModelObj(X);
    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X);

    // Compute initial starting points and bounds.
    SGTELIB::Matrix lvar("lvar", nbVar, 1);
    lvar.fill(0.0);
    SGTELIB::Matrix uvar("uvar", nbVar, 1);
    uvar.fill(INF);
    SGTELIB::Matrix XS("XS", nbVar, 1);
    for (int i=0; i < _n; ++i)
    {
        const double lb = ((_modelLowerBound[i].isDefined())? _modelLowerBound[i].todouble(): NOMAD::M_INF);
        const double ub = ((_modelUpperBound[i].isDefined())? _modelUpperBound[i].todouble(): NOMAD::INF);
        lvar.set(i, 0, lb);
        uvar.set(i, 0, ub);
        XS.set(i, 0, X[i].todouble());
    }

    for (int j=0; j < _nbCons; ++j)
    {
        XS.set(j + _n, 0, -cons.get(j, 0)); // S = -cons
    }


    // Check bound compatibilities and feasibility
    snapToBounds(XS, lvar, uvar);
    for (int i=0; i < nbVar; ++i)
    {
        const bool areBoundsCompatible = (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!areBoundsCompatible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solveAugLag assertion error: Error compatibility lower and upper bound");
        }

        const bool feasible = (XS.get(i, 0) >= lvar.get(i, 0)) && (XS.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            std::cout << XS.get(i, 0) - lvar.get(i, 0) << " " << uvar.get(i, 0) - XS.get(i, 0) << std::endl;
            std::cout << "Error compatibility lower and upper bound" << std::endl;
            throw NOMAD::Exception(__FILE__, __LINE__, "solveAugLag assertion error: Error XS is not feasible");
        }
    }

    // Note: it is not mandatory for the algorithm to start from a feasible point. But it is more efficient
    // in practice. It is a well-known heuristics for augmented Lagrangian algorithms.
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

    double fk = getModelObj(X);
    getModelGrad(&Gk, X);

    // Initialize next iterates
    SGTELIB::Matrix XSp("XSp", nbVar, 1);
    Point Xp(X);

    // Outer loop parameters
    SGTELIB::Matrix lambda_l("lambda_l", _nbCons, 1);
    lambda_l.fill(0.0);

    double Pk = getAugLagModelObj(XS, cons, fk, lambda_l, mu_l);
    SGTELIB::Matrix GradPk("GradAugLag", nbVar, 1);
    getAugLagModelGrad(&GradPk, XS, lambda_l, mu_l);
    // Preallocate matrix hessian of augmented Lagrangian function.
    SGTELIB::Matrix HessPk("HessAugLag", nbVar, nbVar);

    int iterOuterLoop = 0;
    double distXOuterLoop = NOMAD::INF;

    // Compute stopping criterion for outer iterations.
    SGTELIB::Matrix dualFeas("DualFeas", nbVar, 1);
    double ngproj = check_optimality_bounds(XS, GradPk, lvar, uvar, dualFeas);
    bool outerSuccess = ngproj < tol; // nonlinear equality constraints feasibility.
    bool outerFailure = (distXOuterLoop <= tolDistDX) || (iterOuterLoop >= max_iter);

    _verbose && std::cout << "Outer ("<< iterOuterLoop <<"): |G| = " << GradPk.norm();
    _verbose && std::cout << " |Proj(x - grad) - x| = " << ngproj << " P =" << Pk;
    _verbose && std::cout << " f = " << fk << " |c+s| = " << cx;
    _verbose && std::cout << " mu = " << mu_l << " omega = " << omega_l << " eta = " << eta_l << std::endl;

    size_t successiveFailure = 0;
    size_t successiveAcceptable = 0;
    const size_t successiveBeforeUpdate = 2;

    // outer iteration
    while (!outerFailure && !outerSuccess)
    {
        // Solve nonlinear bound-constrained sub-problem
        int innerResult = solveBoundAugLag(XSp, XS, lvar, uvar,
                                           lambda_l, omega_l, mu_l, GradPk, HessPk,
                                           maxIterInner, tolDistDXInner, _verboseFull);
        bool innerSuccess = innerResult > 0;

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
        const double Pkp = getAugLagModelObj(XSp, cons, fk, lambda_l, mu_l);
        
        // Update parameters
        if ((cxp <= eta_l) && innerSuccess)
        {
            // Update lambda multipliers
            cslack.multiply(- 1 / mu_l);
            lambda_l.add(cslack);

            // Update eta
            eta_l = std::max(eta_l * std::pow(mu_l, 0.9), atol);

            // Update omega
            if (innerResult == 1)
            {
                successiveAcceptable = 0;
                omega_l *= mu_l;
            }
            else if (successiveAcceptable >= successiveBeforeUpdate)
            {
                successiveAcceptable += 1;
                omega_l *= std::sqrt(mu_l);
            }
            else
            {
                successiveAcceptable += 1; // Try rerun from new point
            }
            omega_l = std::max(omega_l, std::max(atol * atol, 1E-15));
        }
        else if (innerSuccess) // not feasible
        {
            // Update mu
            if (innerResult == 1)
            {
                successiveAcceptable = 0;
                mu_l /= mu_decrease;
            }
            else if (successiveAcceptable >= successiveBeforeUpdate)
            {
                // The sub-problem has not been solved at optimality, but we still want to keep on
                // iterating. If mu decreases too fast, the algorithm can fail quick.
                successiveAcceptable += 1;
                mu_l /= sqrt(mu_decrease);
            }
            else
            {
                successiveAcceptable += 1; // Try rerun from new point
            }

            // Update eta and omega
            eta_l = std::max(std::pow(mu_l, 0.1), atol);
            omega_l = mu_l;
        }
        else
        { // no success and not feasible
            successiveAcceptable = 0;
            mu_l /= mu_decrease;
            eta_l = std::max(std::pow(mu_l, 0.1), atol);
            omega_l = mu_l;
        }

        // Check optimality and accept new iterate
        getAugLagModelGrad(&GradPk, XSp, lambda_l, mu_l);
        ngproj = check_optimality_bounds(XSp, GradPk, lvar, uvar, dualFeas);
        outerSuccess = (ngproj < tol) && (cxp < tol);

        if (innerSuccess || (Pkp < successRatio * Pk)) // we accept this step
        {
            distXOuterLoop = NOMAD::Point::dist(X, Xp).todouble(); // diff XSp, XS
            XS = XSp;
            X = Xp;
            Pk = Pkp;
            cx = cxp;
            successiveFailure = 0;
        }
        else
        { // for verbose only
            size_t feasibilityPhase = solveLM(Xp, XSp, lvar, uvar, cons, mu_l, omega_l,
                                              30, 1E-15, false, _verbose);
            if (feasibilityPhase > 0)
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
                successiveFailure = 0;
            }
            else
            {
                successiveFailure += 1;
            }
            getAugLagModelGrad(&GradPk, XSp, lambda_l, mu_l);
            ngproj = check_optimality_bounds(XSp, GradPk, lvar, uvar, dualFeas);
        }

        // Increment counter
        iterOuterLoop += 1;

        // Compute stopping criteria
        outerFailure = (iterOuterLoop >= max_iter) || (distXOuterLoop <= tolDistDX);
        outerFailure = outerFailure || (mu_l <= atol / mu_decrease);
        outerFailure = outerFailure || (successiveFailure >= maxSuccessiveFail);

        _verbose && std::cout << "Outer ("<< iterOuterLoop <<") " << innerResult;
        _verbose && std::cout << ": |Proj(x - grad) - x| = " << ngproj << " P = " << Pk << " f = " << fk;
        _verbose && std::cout << " |c+s| = " << cx << " mu = " << mu_l << " omega = " << omega_l << " eta = " << eta_l << std::endl;
    }

    // Ending output
    _verbose && std::cout << "End of solveAugLag" << std::endl;
    _verbose && std::cout << "f(x0) = " << f0 << " f(x*) = " << fk <<std::endl;
    _verbose && std::cout << "|c(x*) + s*| = " << cxp << " tol = " << tol << std::endl;
    if (outerSuccess)
    {
        _verbose && std::cout << "success" << " : " << ngproj << " < " << tol << std::endl;
    }
    else if (outerFailure)
    {
        _verbose && std::cout << " failure" << " iter = " << iterOuterLoop << std::endl;
        _verbose && std::cout << " dist = " << distXOuterLoop << " <= " << tolDistDX << std::endl;
        _verbose && std::cout << " too small parameters? " << mu_l * mu_l << " <= " << atol << std::endl;
        _verbose && std::cout << " successiveFailure = " << successiveFailure << std::endl;
    }
    else
    {
        _verbose && std::cout << "unknown stopping";
    }

    return true;
}

int NOMAD::QPSolverOptimize::solveBoundAugLag(
    SGTELIB::Matrix& XSp,
    const SGTELIB::Matrix& XS,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    const SGTELIB::Matrix& lambda,
    const double omega,
    const double mu,
    SGTELIB::Matrix& GradPk,
    SGTELIB::Matrix& HessPk,
    const size_t maxIterInnerLoop, // = 50
    const double tolerance_distDX, // = 1E-15
    const bool verbose)
{
    const int nvar = XS.get_nb_rows();
    const int ncon = static_cast<int>(_nbCons);

    // Check dimension compatibility ...
    lencheck(nvar, XS);
    lencheck(nvar, XSp);
    lencheck(nvar, lvar);
    lencheck(nvar, uvar);
    lencheck(ncon, lambda);

    // ... and parameters
    if (mu <= 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "solveBoundAugLag assertion error: mu must be positive");
    }
    if (omega < 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "solveBoundAugLag Assertion error: omega must be non-negative");
    }

    // Initialization
    XSp = XS;
    snapToBounds(XSp, lvar, uvar);
    for (int i=0; i < nvar; ++i)
    {
        const bool areBoundsCompatible = (lvar.get(i, 0) <= uvar.get(i, 0));
        if (!areBoundsCompatible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solveBoundAugLag assertion error: Error compatibility lower and upper bound");
        }

        const bool feasible = (XSp.get(i, 0) >= lvar.get(i, 0)) && (XSp.get(i, 0) <= uvar.get(i, 0));
        if (!feasible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solveBoundAugLag assertion error: Error XS is not feasible");
        }
    }

    // Inner loop parameters
    // Trust region update parameters
    const double epsilon_1 = 0.05;
    const double epsilon_2 = 0.9;
    const double gamma_1 = 0.5; // 0.01; // 0.25;
    const double gamma_2 = 2; // 100; // 2.5

    // Trust region radius parameters
    // double delta = 1.0; // trust region initial radius
    const double smallestDelta = 1E-15;
    const double largestDelta = 1E15;

    // BCQP sub-problem parameters
    const size_t maxIterBCQP = 50;
    getAugLagModelHess(&HessPk, XSp, lambda, mu);
    const double atol_BCQP = std::min(1E-8, omega); // By default: omega
    const double rtol_BCQP = 1E-15 * HessPk.norm() / nvar; // 1E-7 by default

    // Others
    const size_t limitUnsuccessful = 40;

    // Compute stopping criteria
    SGTELIB::Matrix dualFeas("dualFeas", nvar, 1);
    check_optimality_bounds(XS, GradPk, lvar, uvar, dualFeas);
    double ngproj = dualFeas.norm_inf();
    const double ngproj0 = ngproj;
    const double tolSuccess = (omega >= 1) ? omega : omega * (1 + ngproj);
    bool innerSuccess = ngproj < tolSuccess; // nonlinear equality constraints feasible

    // Compute the initial trust-region radius
    double delta = 0.1 * ngproj0;

    size_t iterInnerLoop = 0;
    double distXInnerLoop = INF;
    bool innerFailure = (distXInnerLoop <= tolerance_distDX) || (iterInnerLoop >= maxIterInnerLoop);

    // Display initial information
    double Pk = getAugLagModelObj(XSp, lambda, mu);
    verbose && std::cout << "Inner 0 ("<< iterInnerLoop <<") : |Proj(x - grad) - x| = " << ngproj;
    verbose && std::cout << " P = " << Pk << " min tol = " << rtol_BCQP << std::endl;

    // Pre-allocations of trust-region problem variables
    SGTELIB::Matrix Xcan("Xcan", nvar, 1); // x candidate
    SGTELIB::Matrix d("d", nvar, 1); // direction
    SGTELIB::Matrix dlvar("dlvar", nvar, 1); // lower bound of the trust region problem
    SGTELIB::Matrix duvar("duvar", nvar, 1); // upper bound of the trust region problem

    // Initialization of non monotonicity parameters
    const size_t maxNonMonotoneSteps = 10;
    size_t iterNonMonotoneSteps = 0;
    double sigRef = 0.0, sigCan = 0.0;
    double fMin = Pk, fRef = Pk, fCan = Pk;
    double rhoHis = 0;

    bool success = innerSuccess; // overall success
    size_t successiveUnsuccessful = 0;
    bool subPbSuccess = false;
    while (!innerFailure && !innerSuccess)
    {
        // Get quadratic models
        if (successiveUnsuccessful == 0)
        {
            getAugLagModelGrad(&GradPk, XSp, lambda, mu);
            getAugLagModelHess(&HessPk, XSp, lambda, mu);
        }

        // Check bound compatibilities and direction feasibility.
        bool feasible, dfeasible;
        feasible = dfeasible = true;
        for (int j=0; j < nvar; ++j)
        {
            dlvar.set(j, 0, std::max(lvar.get(j, 0) - XSp.get(j, 0), -delta));
            duvar.set(j, 0, std::min(uvar.get(j, 0) - XSp.get(j, 0), delta));
            const bool areBoundsCompatible = dlvar.get(j, 0) <= duvar.get(j, 0);
            if (!areBoundsCompatible)
                throw NOMAD::Exception(__FILE__, __LINE__, "solveBoundAugLag assertion error (d): Error compatibility lower and upper bound");

            feasible = feasible && (dlvar.get(j, 0) <= 0) && (duvar.get(j, 0) >= 0);
            dfeasible = dfeasible && (dlvar.get(j, 0) <= d.get(j, 0)) && (duvar.get(j, 0) >= d.get(j, 0));

        }

        // When x has not been updated (meaning the gradient and the hessian of the trust-region model is still
        // the same) and the last direction is still feasible, we do not need to solve again the BCQP sub-problem,
        // as the solution will be the same.
        if ((successiveUnsuccessful == 0) || !dfeasible)
        {
            d.fill(0);
            if (!feasible)
            {
                return 0;
            }

            subPbSuccess = solveBCQP(d, HessPk, GradPk, 0.0, dlvar, duvar,
                                     maxIterBCQP, atol_BCQP, rtol_BCQP); // you can specify verbose
        }

        // Trust region update
        // Step 1: Compute trust-region ratio = ared / pred
        const double pred = getModelObj(d, HessPk, GradPk);
        if (pred > 0)
        {
            // std::cerr << "Assertion error: prediction " << pred << " > 0" <<std::endl;
            return false;
        }
        Xcan = XSp; Xcan.add(d);
        const double fTrial = getAugLagModelObj(Xcan, lambda, mu);
        const double ared = compute_AugLag_TR_ared(XSp, Xcan, lambda, mu); // (f(x + d) - f(x)

        // Non monotone strategy
        rhoHis = (fRef - fTrial) / (sigRef - pred);
        const double rho = std::max(rhoHis, ared / pred);

        // Step 2 : Update trust-region radius and candidate.
        const double nd = d.norm_inf(); //d.norm();
        double alpha = 1.0;
        bool pointAccepted = false;
        if (subPbSuccess && (rho >= epsilon_1)) // r >= epsilon_1
        {
            // Accept the point
            XSp = Xcan;
            pointAccepted = true;

            // if (ared <= pred * epsilon_2) // r >= epsilon_2
            if (rho >= epsilon_2) // r >= epsilon_2
            {
                delta = std::min(gamma_2 * std::max(nd, delta), largestDelta);
            }

            // Update non monotonicity parameters (taken from npl.py)
            sigRef -= pred;
            sigCan -= pred;
            if (fTrial < fMin)
            {
                fCan = fTrial;
                fMin = fTrial;
                sigCan = 0;
                iterNonMonotoneSteps = 0;
            }
            else
            {
                iterNonMonotoneSteps += 1;
            }

            if (fTrial > fCan)
            {
                fCan = fTrial;
                sigCan = 0;
            }

            if (iterNonMonotoneSteps == maxNonMonotoneSteps)
            {
                fRef = fCan;
                sigRef = sigCan;
            }

            success = true;
            successiveUnsuccessful = 0;
        }
        else
        {
            // First try a backtracking linesearch along the direction d: follow Nocedal and Yuan
            const double slope = SGTELIB::Matrix::dot(d, GradPk);
            const size_t maxBk = 5; // We impose a limit of five backtracking iterations
            size_t bk = 0;
            const double fStart = getAugLagModelObj(XSp, lambda, mu);
            double fTrialBackTrack = fStart;
            bool ArmijoCond = fTrialBackTrack <= fStart + 1e-4 * alpha * slope;
            while ((bk < maxBk) && !ArmijoCond)
            {
                alpha /= 1.2;
                Xcan = d; Xcan.multiply(alpha); Xcan.add(XSp);
                fTrialBackTrack = getAugLagModelObj(Xcan, lambda, mu);
                ArmijoCond = fTrialBackTrack <= fStart + 1e-4 * alpha * slope;
                bk += 1;
            }

            if (ArmijoCond)
            {
                // Accept new candidate
                XSp = Xcan;
                successiveUnsuccessful = 0;

                // Update delta
                delta = std::max(std::min(alpha * nd, delta), smallestDelta);
            }
            else
            {
                // Do not accept the point and reduce trust-region radius delta.
                delta = std::max(gamma_1 * std::min(delta, nd), smallestDelta);
                successiveUnsuccessful += 1;
            }
        }

        // Check optimality
        getAugLagModelGrad(&GradPk, XSp, lambda, mu);
        // ngproj =
        check_optimality_bounds(XSp, GradPk, lvar, uvar, dualFeas);
        ngproj = dualFeas.norm_inf();
        innerSuccess = ngproj < tolSuccess || (ngproj < 1e-7 * ngproj0);

        // Update stopping conditions
        iterInnerLoop += 1;
        distXInnerLoop = sqrt(alpha) * d.norm();
        innerFailure = !subPbSuccess || (iterInnerLoop >= maxIterInnerLoop) || (distXInnerLoop <= tolerance_distDX);
        innerFailure = innerFailure || (successiveUnsuccessful > limitUnsuccessful);
    
        Pk = getAugLagModelObj(XSp, lambda, mu);
        verbose && std::cout << "Inner ("<< iterInnerLoop <<") " << subPbSuccess;
        verbose && std::cout << " " << pointAccepted;
        verbose && std::cout << " |Proj(x - grad) - x| = " << ngproj;
        verbose && std::cout << " P = " << Pk << " |d| = " << distXInnerLoop;
        verbose && std::cout  << " TR radius = " << delta << " TR ratio = ";
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
    const SGTELIB::Matrix& XS,
    const SGTELIB::Matrix& XSp,
    const SGTELIB::Matrix& lambda,
    const double mu) const
{
    const int nbVar = _n + _nbCons;

    // Check dimension compatibility
    lencheck(nbVar, XS);
    lencheck(nbVar, XSp);
    lencheck(_nbCons, lambda);

    const double ared = getAugLagModelObj(XSp, lambda, mu) - getAugLagModelObj(XS, lambda, mu);

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
    const double mu) const
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

    for (int i=0; i < _nbCons; i++)
    {
        const double cpsi = XS.get(i + _n, 0) + cons.get(i, 0); // cons + S
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
    const double mu) const
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
    const double mu) const
{
    const int nbVar = _n + _nbCons;
    SGTELIB::Matrix lagHess("lagHess", nbVar, nbVar);
    getAugLagModelHess(&lagHess, XS, lambda, mu);
    return lagHess;
}

void NOMAD::QPSolverOptimize::getAugLagModelHess(
    SGTELIB::Matrix * lagHess,
    const SGTELIB::Matrix & XS,
    const SGTELIB::Matrix & lambda,
    const double mu) const
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
    NOMAD::Point& X,
    SGTELIB::Matrix& XS,
    SGTELIB::Matrix& lvar,
    SGTELIB::Matrix& uvar,
    const SGTELIB::Matrix& cX) const
{
    for (int i = 0; i < _n; ++i)
    {
        const double lb = ((_modelLowerBound[i].isDefined())? _modelLowerBound[i].todouble(): -NOMAD::INF);
        const double ub = ((_modelUpperBound[i].isDefined())? _modelUpperBound[i].todouble(): NOMAD::INF);
        lvar.set(i, 0, lb);
        uvar.set(i, 0, ub);
        double xi = X[i].todouble();
        if ((xi <= lb) || (xi >= ub))
        {
            if ((_modelLowerBound[i].isDefined()) && !(_modelUpperBound[i].isDefined()))
            {
                xi = lb + 0.5;
            }
            else if (!(_modelLowerBound[i].isDefined()) && (_modelUpperBound[i].isDefined()))
            {
                xi = ub - 0.5;
            }
            else if ((_modelLowerBound[i].isDefined()) && (_modelUpperBound[i].isDefined()))
            {
                const double mid = uvar.get(i, 0) - lvar.get(i, 0);
                xi = lb + mid / 2;
            }
            else
            {
                xi = 0.0;
            }
        }
        XS.set(i, 0, xi);
    }

    for (int j = 0; j < _nbCons; ++j)
    {
        lvar.set(j + _n, 0, 0.0);
        uvar.set(j + _n, 0, INF);
        XS.set(j + _n, 0, std::max(-cX.get(j, 0), 0.5)); // S = -cons
    }

    return check_strict_feasible(XS, lvar, uvar);
}

bool NOMAD::QPSolverOptimize::solveTRIPM(
    NOMAD::Point& X,
    const int max_iter, // 30
    const double tolDistDX, // -1
    const double atol, // 1E-7
    const double rtol, // 1E-7
    const double mu0, // 0.5
    const double muDecrease, // 2
    const size_t maxIterInner, // 40
    const bool verbose, // true
    const bool verbose_SolverBarrier) // true
{
    // We use slack variables here
    const int nbVar = _n + _nbCons;

    if (muDecrease <= 1)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, " muDecrease must be > 1");
    }
    const double mu_decrease = muDecrease;
    const double smallest_tol_mu = atol / 100;
    double mu = mu0;
    double tol_mu = mu0;

    // Compute starting point
    SGTELIB::Matrix cons("cons", _nbCons, 1);
    getModelCons(&cons, X);
    SGTELIB::Matrix lvar("lvar", nbVar, 1);
    SGTELIB::Matrix uvar("uvar", nbVar, 1);
    SGTELIB::Matrix XS("XS", nbVar, 1);

    getStrictFeasiblePoint(X, XS, lvar, uvar, cons);

    // When the difference between lvar and uvar is too small, stop the algorithm
    for (int i = 0; i < nbVar; ++i)
    {
        if (std::abs(uvar.get(i, 0) - lvar.get(i, 0)) <= 1e-8)
        {
            return false;
        }
    }

    solveLM(X, XS, lvar, uvar, cons, mu, tol_mu, 30, 1E-15, true, verbose);

    // Check bounds compatibility
    // bool bound_compat = true;
    for (int i = 0; i < nbVar; ++i)
    {
        const bool areBoundsCompatible = lvar.get(i, 0) <= uvar.get(i, 0);
        if (!areBoundsCompatible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "solveTRIPM assertion error: Error compatibility lower and upper bound");
        }
    }

    // Check feasibility of the starting point
    for (int i = 0; i < nbVar; ++i)
    {
        const bool feasible = (XS.get(i, 0) >= lvar.get(i, 0)) && (XS.get(i, 0) <= uvar.get(i, 0));
        // Levenberg-Marquardt algorithm has failed.
        if (!feasible)
        {
            return false;
        }
    }

    // Compute cslack variables
    SGTELIB::Matrix cslack("c+s", _nbCons, 1);
    for (int j = 0; j < _nbCons; ++j)
    {
        cslack.set(j, 0, XS.get(j + _n, 0) + cons.get(j, 0));
    }
    double cx = cslack.norm();
    double cxp = cx;

    // Compute Lagrange multipliers estimates
    SGTELIB::Matrix Gk("Gk", nbVar, 1);
    getModelGrad(&Gk, X);
    SGTELIB::Matrix Jx = getModelJacobian(X);
    SGTELIB::Matrix lambda("lambda", _nbCons, 1);
    compute_slack_multiplier(lambda, XS, Jx, Gk, 0);

    // Allocate memory for the following vectors
    SGTELIB::Matrix XSp("XSp", nbVar, 1); // Next iterate
    SGTELIB::Matrix p("p", nbVar, 1); // The primal-dual step
    Point Xp(X);

    // TRIPM parameters initialization
    const size_t maxSuccessiveFail = 3;
    size_t successiveFailure = 0;
    size_t successiveAcceptable = 0;
    size_t successiveBeforeUpdate = 2;
    double distXOuterLoop = INF;
    double tolDistDXInner = 1e-15;

    // Compute outer stopping criteria
    double ng = Gk.norm_inf();
    // double ng0 = ng;
    double res = errorTRIPM(XS, lvar, uvar, lambda, cslack, 0);
    double tol = std::max(ng, 1.0) * 1e-6;
    //double tol = atol + std::max(ng0, res) * rtol;
    bool outerSuccess = res <= tol;
    int iterOuterLoop = 0;
    bool outerFailure = (distXOuterLoop <= tolDistDX) || (iterOuterLoop >= max_iter);

    // Initial logs
    double fk = getModelObj(X);
    double f0 = fk;

    verbose && std::cout << "Outer ("<< iterOuterLoop <<"): ";
    verbose && std::cout << " |E(x,s,y;0)| = " << res;
    verbose && std::cout << " f = " << fk << " |c+s| = " << cx;
    verbose && std::cout << " mu = " << mu << " e_mu = " << tol_mu << std::endl;

    // outer iteration
    while (!outerFailure && !outerSuccess)
    {
        const int innerResult = solver_barrier(Xp, XSp, p, cslack, XS, lvar, uvar,
                                               lambda, Gk, cons, Jx,
                                               mu, std::max(mu, 1e-6 - mu), //tol_mu,
                                               maxIterInner,
                                               tolDistDXInner, verbose_SolverBarrier);
        const bool innerSuccess = innerResult > 0;

        for (int i = 0; i < _n; ++i)
        {
            Xp[i] = XSp.get(i, 0);
        }
        getModelCons(&cons, Xp);
        for (int j = 0; j < _nbCons; ++j)
        {
            cslack.set(j, 0, XSp.get(j + _n, 0) + cons.get(j, 0));
        }
        cxp = cslack.norm();
        fk = getModelObj(Xp);

        if (innerSuccess)
        {
            distXOuterLoop = NOMAD::Point::dist(X, Xp).todouble(); // diff XSp, XS
            XS = XSp;
            X = Xp;
            cx = cxp;
            successiveFailure = 0;
            if (innerResult == 1)
            {
                successiveAcceptable = 0;
                mu /= mu_decrease;
                tol_mu /= mu_decrease;
            }
            else if ((successiveFailure > 0) || (successiveAcceptable >= successiveBeforeUpdate))
            { // slower decrease
                successiveAcceptable = 0;
                mu /= sqrt(mu_decrease);
                tol_mu /= sqrt(mu_decrease);
            }
            else
            {
                // Try to re-run from new point.
                successiveAcceptable += 1;
            }
            // Compute stopping criterion
            getModelGrad(&Gk, X);
            ng = Gk.norm_inf();
            tol = std::max(ng, 1.0) * 1e-6;
            res = errorTRIPM(XS, lvar, uvar, lambda, cslack, 0);
            outerSuccess = (res <= tol);
        }
        else
        {
            successiveAcceptable = 0;
            successiveFailure += 1;

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

                    // Compute stopping criterion
                    getModelGrad(&Gk, X);
                    ng = Gk.norm_inf();
                    tol = std::max(ng, 1.0) * 1e-6;
                    res = errorTRIPM(XS, lvar, uvar, lambda, cslack, 0);

                    outerSuccess = (res <= tol);
                    successiveFailure = 0;
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
        outerFailure = (distXOuterLoop <= tolDistDX) || (iterOuterLoop >= max_iter);
        outerFailure = outerFailure || (mu <= atol / mu_decrease);
        outerFailure = outerFailure || (successiveFailure >= maxSuccessiveFail);

        verbose && std::cout << "Outer ("<< iterOuterLoop <<") " << innerResult;
        verbose && std::cout << ": |E(x,s,y;0)| = " << res << " f = " << fk ;
        verbose && std::cout << " |c+s| = " << cx << " mu = " << mu << " e_mu = " << tol_mu << std::endl;
    }

    // Ending output
    verbose && std::cout << "End of solveTR-IPM" << std::endl;
    verbose && std::cout << "f(x0) = " << f0 << " f(x*) = " << fk <<std::endl;
    verbose && std::cout << "|c(x*) + s*| = " << cxp << " tol = " << tol << std::endl;
    if (outerSuccess)
    {
        verbose && std::cout << "success" << " : " << res << " < " << tol << std::endl;
    }
    else if (outerFailure)
    {
        verbose && std::cout << " failure" << " iter = " << iterOuterLoop << std::endl;
        verbose && std::cout << " dist = " << distXOuterLoop << " <= " << tolDistDX << std::endl;
        verbose && std::cout << " too small parameters? " << mu << "<= " << atol << std::endl;
        verbose && std::cout << " successiveFailure = " << successiveFailure << std::endl;
    }
    else
    {
        verbose && std::cout << "unknown stopping";
    }

    const bool success = (fk < f0);
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
    // W = [Jx S] where S = diag(s), s slack variable.
    // bls = [grad f(x) - mu e]^T with e in R^{nbCons X nbCons}
    for (int i = 0; i < _n; i++)
    {
        for (int j = 0; j < _nbCons; j++)
        {
            W.set(i, j, Jx.get(j, i));
        }
        bls.set(i, 0, Gx.get(i, 0));
    }
    for (int i = 0; i < _nbCons; i++)
    {
        for (int j = 0; j < _nbCons; j++)
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
    // Solve least square solution of ||Wy - bls||
    y = SGTELIB::Matrix::solve_least_squares_SVD(W, bls);

    // Enforce sign of y
    for (int i = 0; i < _nbCons; i++)
    {
        if (y.get(i, 0) >= 0)
        {
            const double si = XS.get(_n + i, 0);
            y.set(i, 0, -std::min(std::abs(1E-3), std::abs(mu / si)));
        }
    }
}

double NOMAD::QPSolverOptimize::errorTRIPM(
    const SGTELIB::Matrix& XS,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    const SGTELIB::Matrix& lambda,
    const SGTELIB::Matrix& cslack,
    const double mu)
{
    const int nbVar = _n + _nbCons;

    lencheck(nbVar, XS);
    lencheck(_nbCons, lambda);
    lencheck(_nbCons, cslack);

    SGTELIB::Matrix X("X", _n, 1);

    for (int i = 0; i < _n; ++i)
    {
        X.set(i, 0, XS.get(i, 0));
    }

    auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);

    // Derivative w.r.t. X
    SGTELIB::Matrix lagGradX("tempX", _n, 1);
    SGTELIB::Matrix Mpredict_grad("grad_predict", _nbCons + 1, _n);
    SGTELIB::Matrix Jx("Jx", _nbCons, _n);
    surrogate_prs->getModelLagGrad(&lagGradX, &Mpredict_grad, &Jx, X.transpose(), lambda);

    SGTELIB::Matrix dual_feas = SGTELIB::Matrix("dual_feas", _n, 1);

    // Compute X - P[X - grad L(X)] where P[X] is the projection of X on [lvar, uvar]
    for (int i = 0 ; i < _n ; i++)
    {
        dual_feas.set(i, 0, X.get(i, 0) - lagGradX.get(i, 0));
        if (dual_feas.get(i, 0) < lvar.get(i, 0))
        {
            dual_feas.set(i, 0, lvar.get(i, 0));
        }
        if (uvar.get(i, 0) < dual_feas.get(i, 0))
        {
            dual_feas.set(i, 0, uvar.get(i, 0));
        }
        dual_feas.set(i, 0, dual_feas.get(i, 0) - X.get(i, 0));
    }

    // Compute ||-Sy - mu||_inf
    double lagGradS = 0;
    for (int i=0; i < _nbCons; i++)
    {
        lagGradS = std::max(std::abs(-XS.get(i + _n, 0) * lambda.get(i, 0) - mu), lagGradS);
    }
    return std::max(lagGradS, std::max(cslack.norm_inf(), dual_feas.norm_inf()));
}

size_t NOMAD::QPSolverOptimize::solveLM(
    NOMAD::Point& X,
    SGTELIB::Matrix& XS,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    SGTELIB::Matrix& cX,
    const double Fx,
    const double tol,
    const size_t maxIterInner,
    const double tolDistDXInner,
    const bool checkStrict, // if true, it checks whether the point is strictly feasible.
    const bool verbose)
{
    // We use slack variables here
    const int nbVar = _n + _nbCons;

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
        NOMAD::DoglegTRSolver::solve(vxs, W, cslack, Delta);
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
            if (vxi != 0)
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
            double xcan = XSp.get(i, 0) + vxs.get(i, 0);
            xcan = std::min(xcan, uvar.get(i, 0) - 1e-13);
            xcan = std::max(xcan, lvar.get(i, 0) + 1e-13);
            Xcan.set(i, 0, xcan);
            double xp = X[i].todouble() + vxs.get(i, 0);
            xp = std::min(xp, uvar.get(i, 0) - 1e-13);
            xp = std::max(xp, lvar.get(i, 0) + 1e-13);
            Xp[i] = xp;
        }

        // Magic step: update slack variable s such that s = max(0, c(x)).
        // Allow the function to decrease further.
        getModelCons(&cxp, Xp);
        for (int i = 0; i < _nbCons; ++i)
        {
            const double ci = cxp.get(i, 0);
            if (ci < 0)
            {
                Xcan.set(i + _n, 0, -ci);
            }
            else
            {
                Xcan.set(i + _n, 0, XSp.get(i + _n, 0) + vxs.get(i + _n, 0));
            }
            double si = Xcan.get(i + _n);
            si = std::min(si, uvar.get(i + _n, 0) - 1e-8);
            si = std::max(si, lvar.get(i + _n, 0) + 1e-8);
            Xcan.set(i + _n, 0, si);
        }

        // Make sure the bounds are satisfied strictly
        checkStrict && check_strict_feasible(Xcan, lvar, uvar);

        for (int i = 0; i < nbVar; i++)
        {
            vxs.set(i, 0, Xcan.get(i, 0) - XSp.get(i, 0));
        }

        if (f_normal_model < 0)
        {
            _verbose && std::cout << " solver normal step: |v|=" << vxs.norm() << " f(v)=" << f_normal_model << " |c(xv) + sv|=" << checkslack.norm() << " |c(x) + s|=" << cslack.norm() << std::endl;
        }

        // Compute the residual: r = Jx * vx + vs + (cx + s)
        SGTELIB::Matrix::inplace_product(r, W, vxs);
        r.add(cslack);
        SGTELIB::Matrix::inplace_product(WtWr, W.transpose(), r);

        for (int j = 0; j < _nbCons; ++j)
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
        success = success || (vxs.norm() <= 1e-10);

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
        NOMAD::Point& X,
        SGTELIB::Matrix& XSp,
        SGTELIB::Matrix& p,
        SGTELIB::Matrix& cslack,
        const SGTELIB::Matrix& XS,
        const SGTELIB::Matrix& lvar,
        const SGTELIB::Matrix& uvar,
        SGTELIB::Matrix& lambda,
        SGTELIB::Matrix& Gx,
        SGTELIB::Matrix& cons,
        SGTELIB::Matrix& Jx,
        const double mu,
        const double tol_mu,
        const size_t maxIterInner,
        const double tolDistDXInner,
        const bool verbose,
        const bool verbosePCG)
{
    // We use slack variables here
    const int nbVar = _n + _nbCons;

    // Initialize X
    XSp = XS;
    check_strict_feasible(XSp, lvar, uvar);

    // Allocation of matrices and vectors for solver_barrier
    SGTELIB::Matrix Xcan("Xcan", nbVar, 1);
    Point Xp(X);
    SGTELIB::Matrix r("r", _nbCons, 1); // The residual

    // Allocation of matrices for the normal step
    SGTELIB::Matrix W("W", _nbCons, _nbCons + _n);
    SGTELIB::Matrix Wscal("Wscal", _nbCons, _nbCons + _n);
    SGTELIB::Matrix WtW("WtW", _nbCons + _n, _nbCons + _n);
    SGTELIB::Matrix wq("wq", _nbCons + _n, 1);
    SGTELIB::Matrix vxs("vxs", _nbCons + _n, 1); // The normal step
    SGTELIB::Matrix zer("zer", _nbCons + _n, 1);
    zer.fill(0.0);
    SGTELIB::Matrix temp("-temp-", _nbCons + _n, 1);
    SGTELIB::Matrix vs("vs", _nbCons, 1);
    SGTELIB::Matrix vx("vx", _n, 1);

    // Compute p
    double np = p.norm();
    SGTELIB::Matrix Q("Q", _nbCons + _n, _nbCons + _n);
    Q.fill(0);
    SGTELIB::Matrix qc("qc", _nbCons + _n, 1);
    SGTELIB::Matrix ps("ps", _nbCons, 1);
    SGTELIB::Matrix px("px", _n, 1);
    SGTELIB::Matrix x0PCG("x0PCG", _nbCons + _n, 1); // Starting point for PCG
    SGTELIB::Matrix r0("r0", _nbCons, 1); // Initial residual for PCG (NB: not used, remove ?)
    r0.fill(0.0);

    // Barrier solver parameters:
    const double epsilon_1 = 1e-8; // trust-region successful ratio
    const double epsilon_2 = 0.9; // trust-region very successful ratio
    const double gamma_1 = 0.5; // trust-region decrease factor
    const double gamma_2 = 2; // trust-region increase factor

    // Merit function parameters
    double nu = 1;

    double Delta = 1; // trust-region initial radius
    const double smallestDelta = 1e-15;
    const double largestDelta = 1e15;

    const double tau = 0.995;
    // const double rho = 0.5;
    const double rho = 0.1;

    const double normal_step_regularization = 1e-7; // regularizer of the normal equation in normal step
    const double Delta_normal_step_factor = 0.8; // Factor of Delta used in normal step
    const double small_p = 1e-15; // Below this value `p` is considered 0 (declare success).

    // Stopping criteria for the barrier solver
    // 1- Inner success conditions
    compute_slack_multiplier(lambda, XSp, Jx, Gx, mu);
    double res = errorTRIPM(XSp, lvar, uvar, lambda, cslack, mu);
    bool innerSuccess = (res <= tol_mu);

    // 2- Inner failure condition
    size_t iterInnerLoop = 0;
    double distXInnerLoop = INF;
    bool innerFailure = (distXInnerLoop <= tolDistDXInner) || (iterInnerLoop >= maxIterInner);

    verbose && std::cout << "Inner 0 (-): ";
    verbose && std::cout << " |E(x,s,y; mu)| = " << res << " mu = " << mu;
    verbose && std::cout << " |c+s| = " << cslack.norm();
    verbose && std::cout << " mu = " << mu << " e_mu = " << tol_mu << std::endl;

    size_t successiveUnsuccessful = 0;
    bool success = innerSuccess;
    while (!innerFailure && !innerSuccess)
    {
        // Compute the normal step vxs := (vx, vs), i.e. solve:
        //  min (|| W vxs + clack ||_2)^2
        //  vxs
        //  s.t. || vxs || <= delta_normal_step
        //       vs        >= -tau/2
        //       tau/2 * l <= vx + tau/2 * x <= tau/2 * u
        //       _     _
        // W =  | Jx  S | with S = diag(s), where s are the slack variables of the inequality constraints.
        //      |_     _|
        // NB: the quadratic function is convex, i.e. W^t W is at least semi-definite positive
        if (successiveUnsuccessful == 0)
        {
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
        }

        vxs.fill(0.0);
        NOMAD::DoglegTRSolver::solve(vxs, W, cslack, Delta_normal_step_factor * Delta);
        double f_normal_model = getModelObj(vxs, WtW, wq, 0.5 * std::pow(cslack.norm(), 2));

        ///////////////////////////////////////////////////////////////
        for (int i = 0; i < _n; i++)
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
        for (int j = 0; j < _nbCons; ++j)
        {
            checkslack.set(j, 0, XSp.get(j + _n, 0) + checkcons.get(j, 0));
        }
        ///////////////////////////////////////////////////////////////

        // Backtrack to satisfy vs >= -tau/2
        double backtrack_length = 1;
        for (int i = 0; i < _nbCons; i++)
        {
            const double vsi = vxs.get(_n + i, 0);
            if (vsi < -tau / 2)
            {
                backtrack_length = std::min(backtrack_length, - tau / (2 * vsi));
            }
        }
        // Backtrack to satisfy vx + tau/2 * x >= tau/2 * l and tau/2 * u >= vx + tau/2 * x
        for (int i = 0; i < _n; i++)
        {
            const double vxi = vxs.get(i, 0);
            const double xi = XSp.get(i, 0);
            const double li = lvar.get(i, 0);
            const double ui = uvar.get(i, 0);
            if (vxi != 0)
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

        // Compute p:= (px, ptildes) with projected CG, i.e. solve:
        // min (1/2) p' Q p + qc' p
        //  p
        // s.t. W p = W v
        //      || p || <= Delta
        //      ptildes + tau >= 0,
        //      tau * l <= px + tau * x <= tau * u
        //
        // where:
        //      _                                              _
        // Q = | HLag + mu * (diag(1/(x-l)^2) + diag(1/(u-x)^2) |
        //     |           - S Sigma S                          |
        //     |_                                              _|
        //       _                                           _
        // qc = | grad f - mu (diag(1/(x-l)) - diag(1/(u-x))) |
        //      |         -mu e                               |
        //      |_                                           _|
        // where Sigma = S^-1 Y with Y the estimated Lagrange multipliers (always negative)
        if (successiveUnsuccessful == 0)
        {
            SGTELIB::Matrix HLag = getModelLagHessian(X, lambda);
            for (int i = 0; i < _n; i++)
            {
                const double xi = XSp.get(i, 0);
                const double li = lvar.get(i, 0);
                const double ui = uvar.get(i, 0);
                for (int j = 0; j < _n; j++)
                {
                    if (i == j)
                    {
                        Q.set(i, j, HLag.get(i, j) + mu * (1 / std::pow(xi - li, 2) + 1 / std::pow(ui - xi, 2)));
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
                const double si = XSp.get(_n + i, 0);
                const double lambda_i = lambda.get(i, 0);
                Q.set(i + _n, i + _n, - lambda_i * si); // NB: S Sigma S = S (S^-1 lambda) S = lambda S
                qc.set(i + _n, 0, -mu);
            }
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

        // Compute the residual r:= W v
        SGTELIB::Matrix::inplace_product(r, W, vxs);
        x0PCG.fill(0.0);
        p.fill(0);
        const auto solverStatus = NOMAD::ProjectedConjugateGradientSolver::solve(p, Q, qc, W, r, Delta, verbosePCG);
        if (solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::QUAD_ROOTS_ERROR ||
            solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::TR_PARAM_ERROR ||
            solverStatus == NOMAD::ProjectedConjugateGradientSolverStatus::MATRIX_DIMENSIONS_FAILURE)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "TRIPM: Error with conjugate gradient");
        }

        // Backtrack to satisfy ptildes >= -tau e
        backtrack_length = p.norm() > Delta ? Delta / p.norm() : 1.0;
        for (int i = 0; i < _nbCons; i++)
        {
            const double psi = p.get(_n + i, 0);
            if (psi < -tau)
            {
                backtrack_length = std::min(backtrack_length, - tau / psi);
            }
        }
        // Backtrack to satisfy px + tau * x >= tau * l and tau * u >= px + tau * x
        for (int i = 0; i < _n; i++)
        {
            const double pxi = p.get(i, 0);
            const double xi = XSp.get(i, 0);
            const double li = lvar.get(i, 0);
            const double ui = uvar.get(i, 0);
            if (pxi != 0)
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
        }

        p.multiply(backtrack_length);

        // Compute the residual: r = Jx * px + ps + (cx + s)
        // (as ps = S * pstildes, W p = Jx * px + ps).
        SGTELIB::Matrix::inplace_product(r, W, p);
        r.add(cslack);
        double nr = r.norm();

        // Update nu_l
        double den = cslack.norm() - nr + tol_mu;
        if (den < 0)
        {
            // We try a feasibility step only, in this case
            verbose && std::cout << " p does not increase feasibility " << den <<", we try using v: ";
            nr = checkslack.norm();
            den = cslack.norm() - nr;
            verbose && std::cout << den << std::endl;
            if (den <= 0)
            {
                verbose && std::cout << " normal step: |v| = " << vxs.norm() << " f(v) = " << f_normal_model << " |c(xv) + sv| = " << nr << " |c(x) + s| = " << cslack.norm();
                verbose && std::cout << " Wscalp - Wv = " << SGTELIB::Matrix::sub(SGTELIB::Matrix::product(Wscal, p), SGTELIB::Matrix::product(W, vxs)).norm();
                verbose && std::cout << " Wp = " << SGTELIB::Matrix::product(Wscal, p).norm() << std::endl;
                nu = INF;
            }
            else
            { // den > 0
                p = vxs;
                nu = 1e15;
            }
        }
        else if (den > 0)
        {
            nu = std::max(getModelObj(p, Q, qc) / ((1 - rho) * den) + 1, nu);
        }

        // Compute pred (> 0):
        const double pred = nu * den - getModelObj(p, Q, qc);

        // Compute ared (> 0)
        // a- Perform the following change of variables: ps := S ptildes:
        for (int i = 0; i < _nbCons; i++)
        {
            const double si = XSp.get(_n + i, 0);
            p.set(i + _n, 0, p.get(i + _n, 0) * si);
        }
        np = p.norm();

        // b- Compute ared (> 0):
        for (int i = 0; i < _n; i++)
        {
            Xcan.set(i, 0, XSp.get(i, 0) + p.get(i, 0));
            Xp[i] = X[i] + p.get(i, 0);
        }

        // Make sure the problem is strictly feasible, relatively to the bound constraints
        for (int i = 0; i < _n; ++i)
        {
            const double li = lvar.get(i, 0);
            const double ui = uvar.get(i, 0);
            double xican = std::max(Xcan.get(i, 0), li + 1e-8);
            xican = std::min(xican, ui - 1e-8);
            Xcan.set(i, 0, xican);

            double xi = std::max(Xp[i].todouble(), li + 1e-8);
            xi = std::min(xi, ui - 1e-8);
            Xp[i] = xi;
        }

        // Before the trial point is tested for acceptance by the merit function,
        // set the slack variables for which - ci(Xp) > sip to: sip := - ci(Xp),
        // where sp = s + ps.
        // It enables to decrease further the barrier term.
        // This slack reset strategy is taken from:
        //
        // "An interior algorithm for nonlinear optimization that
        //  combines line search and trust region steps"
        // by R.A. Waltz, J.L. Morales, J. Nocedal and D. Orban
        // Mathematical Programming 107, 391–408 (2006)
        // https://doi.org/10.1007/s10107-004-0560-5
        //
        getModelCons(&cons, Xp);
        for (int i = 0; i < _nbCons; i++)
        {
            const double sip = XSp.get(i + _n, 0) + p.get(i + _n, 0);
            Xcan.set(i + _n, 0, std::max(-cons.get(i, 0), sip));
        }

        const double ared = merit_function_barrier(X, XSp, lvar, uvar, mu, nu)
                            - merit_function_barrier(Xp, Xcan, lvar, uvar, mu, nu);

        // Trust-region update
        if ((ared >= pred * epsilon_1) && (pred > 0))
        { // r >= epsilon_1
            // Accept x and s steps
            XSp = Xcan;
            distXInnerLoop = NOMAD::Point::dist(X, Xp).todouble();
            X = Xp;

            // Increase Delta
            if (ared >= pred * epsilon_2) // r >= epsilon_2
            {
                Delta = std::min(gamma_2 * Delta, std::max(1 / mu, largestDelta));
            }

            success = true;
            successiveUnsuccessful = 0;

            // Check optimality:
            getModelCons(&cons, X);
            for (int j = 0; j < _nbCons; ++j)
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
                verbose && std::cout << "quad = " << getModelObj(p, Q, qc) << " <= " << getModelObj(zer, Q, qc);
                verbose && std::cout << " |Wscalp| = " << SGTELIB::Matrix::product(Wscal, p).norm() << " |Wp| = " << SGTELIB::Matrix::product(W, p).norm() << " |Wvxs| = " << SGTELIB::Matrix::product(W, vxs).norm();
                verbose && std::cout << " |r| = " << nr << " |rv| = " << SGTELIB::Matrix::add(SGTELIB::Matrix::product(W, vxs), cslack).norm();
                verbose && std::cout << " |c+s| = " << cslack.norm() << " nu = " << nu << std::endl;
                verbose && std::cout << "pred = " << nu * cslack.norm();
                verbose && std::cout << " + " << - getModelObj(p, Q, qc);
                verbose && std::cout << " + " << - nu * nr << std::endl;
            }
            // Decrease Delta
            Delta = std::max(gamma_1 * std::min(Delta, np), smallestDelta);
            successiveUnsuccessful += 1;
        }
        iterInnerLoop += 1;

        verbose && std::cout << "Inner " << iterInnerLoop << " (" << successiveUnsuccessful << "):";
        verbose && std::cout << " |E(x,s,y;mu)| = " << res << " |r| = " << nr << " |c+s| = " << cslack.norm();
        verbose && std::cout << " nu = " << nu;
        if (pred != 0)
        {
            verbose && std::cout << " ared/pred = " << ared / pred;
        } else
        {
            verbose && std::cout << " ared/pred = " << ared << "/" << pred;
        }
        verbose && std::cout << " Delta = " << Delta << " |p| = " << np << " |v| = " << vxs.norm() << std::endl;

        innerSuccess = (res <= tol_mu) || (np <= 1e-8) || (Delta < 1e-8); //(np <= small_p);
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
        // At least one step has been made
        resultSolverBarrier = 3;
    }
    else
    {
        resultSolverBarrier = 0;
    }
    return resultSolverBarrier;
}

double NOMAD::QPSolverOptimize::merit_function_barrier(
    const NOMAD::Point& X,
    const SGTELIB::Matrix& XS,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    const double mu,
    const double nu)
{
    double fx = getModelObj(X);

    check_strict_feasible(XS, lvar, uvar);

    double res = 0;
    for (int i = 0; i < _nbCons; i++)
    {
        const double si = XS.get(i + _n, 0);
        res -= mu * std::log(si);
    }

    for (int i = 0; i < _n; i++)
    {
        const double xi = XS.get(i, 0);
        const double ui = uvar.get(i, 0);
        const double li = lvar.get(i, 0);
        res -= mu * std::log(xi - li);
        res -= mu * std::log(ui - xi);
    }

    SGTELIB::Matrix cons("cons", _nbCons, 1);
    SGTELIB::Matrix cslack("cslack", _nbCons, 1);

    getModelCons(&cons, X);
    for (int j = 0; j < _nbCons; ++j)
    {
        cslack.set(j, 0, XS.get(j + _n, 0) + cons.get(j, 0));
    }

    return fx + res + nu * cslack.norm();
}

bool NOMAD::QPSolverOptimize::check_strict_feasible(const SGTELIB::Matrix& X,
                                                    const SGTELIB::Matrix& lvar,
                                                    const SGTELIB::Matrix& uvar) const
{
    bool strict_feasible = true;
    for (int i = 0; i < _n; i++)
    {
        const double xi = X.get(i, 0);
        const double ui = uvar.get(i, 0);
        const double li = lvar.get(i, 0);
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
    const SGTELIB::Matrix& X,
    const SGTELIB::Matrix& gradient,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar)
{
    const int n = X.get_nb_rows();
    SGTELIB::Matrix dual_feas = SGTELIB::Matrix("dual_feas", n, 1);
    return check_optimality_bounds(X, gradient, lvar, uvar, dual_feas);
}

double NOMAD::QPSolverOptimize::check_optimality_bounds(
    const SGTELIB::Matrix& X,
    const SGTELIB::Matrix& gradient,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    SGTELIB::Matrix& dual_feas)
{
    const int n = X.get_nb_rows();
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
    SGTELIB::Matrix* d,
    const SGTELIB::Matrix& g, // gradient
    const SGTELIB::Matrix& gW, // gradient active
    const SGTELIB::Matrix& H, // hessian
    const SGTELIB::Matrix& HW, // hessian active
    int* pp,
    double** D,
    double** L,
    const bool* active, // length n
    const double Delta,
    const bool verbose)
{
    const int n = g.get_nb_rows();
    const int nfree = n - sum(active, n);

    // Check dimension compatibility
    lencheck(n, *d);
    sizecheck(n, n, H);
    sizecheck(nfree, nfree, HW);
    lencheck(n, g);
    lencheck(nfree, gW);

    auto sol = new double[nfree];

    const bool solve_success = computeNewtonDirection(gW, pp, D, L, sol, nfree);
    if (!solve_success)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Error with LDLt solve");
    }

    // Update direction
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

    const double slope = SGTELIB::Matrix::dot(g, *d);
    if (slope > 0)
    {
        verbose && std::cout << "Numerical issue Newton direction is not positive definite, slope= " << slope << std::endl;
    }

    // Check if we solve a "fake" trust-region problem.
    const double nd = d->norm();
    if ((Delta < 1E15) && (nd > Delta))
    {
        verbose && std::cout << " Newton direction is not inside the trust-region: " << nd << " >= " << Delta << std::endl;
        d->multiply(Delta / nd);
    }
    verbose && std::cout << "|d|= " << nd << " slope = " << slope << std::endl;

    delete [] sol;

    return true;
}

bool NOMAD::QPSolverOptimize::computeNewtonDirection(
    const SGTELIB::Matrix& g,
    int* pp,
    double** D,
    double** L, // const SGTELIB::Matrix H,
    double* sol,
    const int n) const
{
    lencheck(n, g);

    auto rhs = new double[n];
    for (int i = 0 ; i < n ; ++i )
    {
        rhs[i] = - g.get(i, 0);
        sol[i] = 0;
    }

    // Solve LDLt x = -g (via a direct method)
    string error_msg;
    bool success = true;
    success = NOMAD::ldl_solve(error_msg, D, L, rhs, sol, pp, n);

    delete [] rhs;

    return success;
}

bool NOMAD::QPSolverOptimize::InverseIteration(
    SGTELIB::Matrix* sol,
    const SGTELIB::Matrix& HW,
    const double eigmin,
    const int nfree,
    const double tol, // > 0
    const bool verbose)
{
    // Check dimension compatibility
    lencheck(nfree, *sol);
    sizecheck(nfree, nfree, HW);

    // We compute the eigenvector corresponding to it
    SGTELIB::Matrix bk("bk", nfree, 1); // depend on X
    bk.fill(1.0 / nfree);
    SGTELIB::Matrix bkp("bkp", nfree, 1);
    bkp.fill(0);

    // Initialization of the inverse iteration algorithm
    // Compute (HW - mu I)^-1 where mu = eigmin + 1e-7
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

    // Main iterations of the inverse iteration method
    SGTELIB::Matrix invHWpbk = SGTELIB::Matrix::product(invHWp, bk);
    double Ck = 1.0;

    bool OK = false;
    size_t count = 0;
    const size_t max_count = 1000;
    while (!OK && count < max_count)
    {
        // b_{k+1} = H^{-1} * b_k * C_k
        SGTELIB::Matrix::inplace_product(bkp, invHWpbk, Ck);
        if (bkp.has_nan())
        {
            return false;
        }

        const double fix_point = SGTELIB::Matrix::sub(bkp, bk).norm();

        bk = bkp;
        SGTELIB::Matrix::inplace_product(invHWpbk, invHWp, bk);

        verbose && std::cout << fix_point << " Ck=" << Ck;
        verbose && std::cout << " |bk|=" << bk.norm() << " |bkp|=" << bkp.norm() << std::endl;

        if (invHWpbk.norm() <= 0.0) {
            return false;
        }

        const double Ckp = 1 / invHWpbk.norm();

        // OK = (fix_point <= tol) || (fabs(Ck - Ckp) <= tol);
        OK = (fix_point <= 1E-7) || (std::abs(Ck - Ckp) <= tol);
        Ck = Ckp;

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
    SGTELIB::Matrix& x,
    const SGTELIB::Matrix& A,
    const SGTELIB::Matrix& b,
    const SGTELIB::Matrix& G,
    const SGTELIB::Matrix& c,
    const double delta,
    const SGTELIB::Matrix& x0,
    const double tol,
    const bool verbose)
{
    // NB: m should be less than n
    const int n = x0.get_nb_rows();
    const int m = b.get_nb_rows();

    const int npm = n + m;

    sizecheck(m, n, A);
    lencheck(m, b);
    sizecheck(n, n, G);
    lencheck(n, c);

    // NB : Other choices are possible for H
    SGTELIB::Matrix H = SGTELIB::Matrix::identity(n);

    ////////// LDLt
    // init matrices for LDLt
    auto M = new double *[npm];
    auto L = new double *[npm];
    auto D = new double *[npm];
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
    auto pp = new int [npm];
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

    auto rhs = new double [npm];
    auto sol = new double [npm];

    // Initial guess
    x = x0;

    // Algorithm PCG may assume that an initial feasible point x0 satisfying A x0 = b
    // is provided. If it is not the case, we try to find one.
    SGTELIB::Matrix Ax = SGTELIB::Matrix::product(A, x);
    SGTELIB::Matrix f = b;
    f.sub(Ax);
    double feas = f.norm();

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

    // Temporary solution
    SGTELIB::Matrix xtmp = x;

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

    // Stopping criteria
    const size_t max_iter = npm * 2;
    double rg = SGTELIB::Matrix::dot(r, g);
    const double tol_CG = 0.01 * std::sqrt(rg);
    double dtGd = SGTELIB::Matrix::dot(d, Gd);

    // Given an iterate x, a direction d and a radius delta, compute theta such
    // that || x + theta d || = delta
    // Requires: || d || ≠ 0, and || x || <= delta.
    auto toBoundary = [](const SGTELIB::Matrix& x,
                         const SGTELIB::Matrix& d,
                         const double delta) -> std::pair<double, double>
    {
        const double xtd = SGTELIB::Matrix::dot(x, d);
        const double nd2 = d.normsquare();
        const double nx2 = x.normsquare();
        const double delta2 = delta * delta;

        // Find the quadratic roots of the following problem
        // q2 theta^2 + q1 theta + q0 = 0,
        // where q2 = d'd, q1 = 2 x'd and q0 = x'x - delta^2.
        // We adopt a numerically stable algorithm (taken from Krylov.jl) but other algorithms could be used.
        const double q2 = nd2;
        const double q1 = 2 * xtd;
        const double q0 = nx2 - delta2;
        double root1, root2;
        const bool hasRealRoots = NOMAD::roots_quadratic(q2, q1, q0, root1, root2);
        if (!hasRealRoots)
            throw NOMAD::Exception(__FILE__, __LINE__, "The quadratic does not have real roots");

        return {root1, root2};
    };

    bool max_iter_reached = true;
    for (size_t iter = 0; iter < max_iter; ++iter)
    {
        verbose && std::cout << "PCG-It " << iter << " :";
        verbose && std::cout << " rg = " << rg << " tol = " << tol_CG;
        verbose && std::cout << " dtGd = " << dtGd;
        verbose && std::cout << " nx = " << x.norm() << " <= delta = " << delta << std::endl;

        // Detection of a negative curvature: in this case, return a solution x
        // such that ||x||_2 = delta.
        if (dtGd <= 0)
        {
            // x := x + theta d
            const auto roots = toBoundary(x, d, delta);
            const double theta = std::max(roots.first, roots.second);
            x.add(SGTELIB::Matrix::product(d, theta));
            verbose && std::cout << "PCG: Detection of a negative curvature : stop" << std::endl;
            max_iter_reached = false;
            break;
        }

        const double alpha = rg / dtGd;

        // x = x + alpha * d
        x.add(SGTELIB::Matrix::product(d, alpha));
        const double xnorm = x.norm();

        // Leave trust-region boundary: return a solution x such that ||x||_2 = delta
        // starting from the temporary solution xtmp
        if (xnorm > delta)
        {
            const auto roots = toBoundary(xtmp, d, delta);
            const double theta = std::max(roots.first, roots.second);
            x = xtmp;
            x.add(SGTELIB::Matrix::product(d, theta));
            verbose && std::cout << "PCG: has reached the trust region boundary : stop" << std::endl;
            max_iter_reached = false;
            break;
        }

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

        // Has reached the maximum tolerance
        const double rgp = SGTELIB::Matrix::dot(rp, gp);
        if (rgp < tol_CG)
        {
            verbose && std::cout << "PCG has reached the minimum tolerance : stop" << std::endl;
            max_iter_reached = false;
            break;
        }

        // beta = r+' g+ / r'g
        const double beta = rgp / rg;

        // d = -(g+) + beta * d
        d.multiply(beta);
        d.sub(gp);

        g = gp;
        r = rp;
        rg = rgp;
        xtmp = x;
        SGTELIB::Matrix::inplace_product(Gd, G, d);
        dtGd = SGTELIB::Matrix::dot(d, Gd);
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

    // Normally, one should retrieve the former direction if the feasibility has not be reached.
    if (max_iter_reached)
    {
        verbose && std::cout << "PCG has reached the maximum number of iterations allowed : stop" << std::endl;
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

    const auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    SGTELIB::Matrix Mpredict_grad( "grad_predict", static_cast<int>(_m), static_cast<int>(_n));
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

    const auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    const SGTELIB::Matrix Hx = surrogate_prs->getModelHessian(XX, j);
    sizecheck(_n, _n, Hx);

    return Hx;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelHessian(const NOMAD::Point& X) const
{   
    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, i, X[i].todouble());
    }

    const auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    const SGTELIB::Matrix Hx = surrogate_prs->getModelHessian(XX);
    sizecheck(_n, _n, Hx);

    return Hx;
}

SGTELIB::Matrix NOMAD::QPSolverOptimize::getModelCons(const NOMAD::Point& X_k) const
{
    SGTELIB::Matrix cons("cons", static_cast<int>(_nbCons), 1);
    getModelCons(&cons, X_k);
    return cons;
}

void NOMAD::QPSolverOptimize::getModelCons(SGTELIB::Matrix* cons, const NOMAD::Point& X_k) const
{
    SGTELIB::Matrix XX("X_k", 1, static_cast<int>(_n));
    for (int i = 0; i < _n; i++)
    {
        XX.set(0, i, X_k[i].todouble());
    }

    const auto surrogate_prs = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
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

double NOMAD::QPSolverOptimize::getModelLag(const NOMAD::Point& x,
                                            const SGTELIB::Matrix& multiplier,
                                            const double sigma) const
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
    const NOMAD::Point& x,
    const SGTELIB::Matrix& multiplier,
    const double sigma) const
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
    const NOMAD::Point& X,
    const SGTELIB::Matrix& multiplier,
    const double sigma) const
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
// Model getters for quadratic functions
//
//*****************************************************************************
double NOMAD::QPSolverOptimize::getModelObj(
    const SGTELIB::Matrix& x,
    const SGTELIB::Matrix& H,
    const SGTELIB::Matrix& g,
    const double g0) const
{
    const int n = x.get_nb_rows();

    lencheck(n, x);
    lencheck(n, g);
    sizecheck(n, n, H);

    // Do not use Sgtelib matrix operations for faster computation
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
    return sum;
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

SGTELIB::Matrix NOMAD::QPSolverOptimize::vector_subset(const SGTELIB::Matrix& X, const bool* active) const
{
    const int n = X.get_nb_rows();
    const int nfree = n - sum(active, n);
    SGTELIB::Matrix Xsub("Xsub", nfree, 1);

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
        throw NOMAD::Exception(__FILE__, __LINE__, "Error dimension");
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

SGTELIB::Matrix NOMAD::QPSolverOptimize::matrix_subset(const SGTELIB::Matrix& X, const bool* active) const
{
    const int n = X.get_nb_rows();
    const int nfree = n - sum(active, n);
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
        throw NOMAD::Exception(__FILE__, __LINE__, "Error dimension");
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
    const SGTELIB::Matrix& X,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    const SGTELIB::Matrix& d)
{
    const int n = X.get_nb_rows();

    // Check dimension compatibility
    lencheck(n, X);
    lencheck(n, lvar);
    lencheck(n, uvar);
    lencheck(n, d);

    // Check bound compatibility and feasibility
    for (int i=0; i < n; ++i)
    {
        const bool areBoundsCompatible = lvar.get(i, 0) <= uvar.get(i, 0);
        if (!areBoundsCompatible)
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Assertion error: Error compatibility lower and upper bound");
        }

        const bool feasible = (X.get(i, 0) >= lvar.get(i, 0)) && (X.get(i, 0) <= uvar.get(i, 0));
    }

    double t_max = INF;
    for (int i = 0; i < n ; ++i)
    {
        const double di = d.get(i, 0);
        double gamma = INF;
        if (di > 0)
        {
            gamma = (uvar.get(i, 0) - X.get(i, 0)) / std::abs(di);
        }
        else if (di < 0)
        {
            gamma = (X.get(i, 0) - lvar.get(i, 0)) / std::abs(di);
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

void NOMAD::QPSolverOptimize::project_bounds(SGTELIB::Matrix& d_k, const bool* active)
{
    for (int i = 0 ; i < _n ; ++i)
    {
        if (active[i])
            d_k.set(i, 0, 0);
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
        if ((X.get(i, 0) == lvar.get(i, 0) && d_k.get(i, 0) < 0 ) ||
            (X.get(i, 0) == uvar.get(i, 0) && d_k.get(i, 0) > 0 ))
        {
            d_k.set(i, 0, 0);
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
    const SGTELIB::Matrix& X,
    const SGTELIB::Matrix& lvar,
    const SGTELIB::Matrix& uvar,
    bool* active_l, // length n
    bool* active_u, // length n
    const double tol)
{
    const int n = X.get_nb_rows();

    // Check dimension compatibility
    lencheck(n, X);
    lencheck(n, lvar);
    lencheck(n, uvar);

    if (tol < 0)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter error tol should be nonnegative");
    }

    for (int i = 0; i < n; ++i)
    {
        active_l[i] = (std::abs(X.get(i, 0) - lvar.get(i, 0)) < tol);
        active_u[i] = (std::abs(X.get(i, 0) - uvar.get(i, 0)) < tol);
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


void NOMAD::QPSolverOptimize::snapToBounds(SGTELIB::Matrix& X,
                                           const SGTELIB::Matrix& lowerBound,
                                           const SGTELIB::Matrix& upperBound) const
{
    const int n = X.get_nb_rows();
    
    if (lowerBound.get_nb_rows() != n || upperBound.get_nb_rows() != n)
    {
        std::string err = "snapToBounds: ";
        err += "Inconsistent dimension for bounds or for X. Expecting a vector";
        err += std::to_string(n);
        err += " but sizes are " + std::to_string(lowerBound.get_nb_rows());
        err += " and " + std::to_string(upperBound.get_nb_rows()) + ".";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
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
    }
}


SGTELIB::Matrix NOMAD::QPSolverOptimize::computeQPModelMatrix()
{
    auto surrogate = std::dynamic_pointer_cast<SGTELIB::Surrogate_PRS>(_model);
    if (surrogate == nullptr)
        return SGTELIB::Matrix("error", 0, 0);

    // The surrogate class returns its coefficient alpha such that each column of alpha corresponds to one
    // model. The alpha matrix is transposed to follow the QPModelMatrix format, i.e. each line of the
    // matrix corresponds to one blackbox output.
    auto alpha = surrogate->get_alpha().transpose();
    alpha.set_name("alpha");

    // First, count the ``real'' number of variables of the model: it may be different from _n.
    const int nvar = _trainingSet->get_nvar();

    // Step 1: the alpha coefficients of the diagonal of the Hessian matrix correspond to the basis monomes
    // xi^2 and not 0.5 xi^2. Fix them.
    for (int ind = 0; ind < _m; ++ind)
    {
        for (int i = 0; i < nvar; ++i)
        {
            const double alphaCoeff = 2.0 * alpha.get(ind, i + nvar + 1);
            alpha.set(ind, i + nvar + 1, alphaCoeff);
        }
    }

    // Step 2: A point passed to the surrogate is firstly scaled: x := D x + b, where D is a diagonal
    // matrix and b a vector translation. We change the alpha matrix to directly take this scaling into account.
    // For a given function, alpha is the vector decomposing into: alpha = [alpha0 alphaL alphaQ]', such that:
    //
    // Q(x) = alpha0 + alphaL' x + 0.5 x' Hq x,
    //
    // where Hq is the Hessian matrix built from alphaQ.
    // We apply the following shift:
    //
    // Q(Dx + b) = alpha0 + alphaL' b + 0.5 b' H b + (D alphaL + 0.5 D H b + 0.5 b' H D) x + 0.5 x' D H D x
    //
    // NB: Some coefficients of the H matrix and the alphaL vector may not be defined, since we may be working
    // on a subset of the _n variables
    SGTELIB::Matrix H("H", _n, _n);
    SGTELIB::Matrix alphaL("alphaL", 1, _n);
    for (int ind = 0; ind < _m; ++ind)
    {
        // Compute alpha0 := alpha0 + alphaL' b + 0.5 b' H b
        double alpha0 = alpha.get(ind, 0);
        int k = 1;
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            const double b = _trainingSet->get_X_scaling_b(i);
            alpha0 += alpha.get(ind, k) * b + 0.5 * alpha.get(ind, k + nvar) * b * b;
            k++;
        }
        k += nvar;
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            const double bi = _trainingSet->get_X_scaling_b(i);
            for (int j = 0; j < i; ++j)
            {
                if (_trainingSet->get_X_nbdiff(j) <= 1)
                    continue;

                const double bj = _trainingSet->get_X_scaling_b(j);
                alpha0 += alpha.get(ind, k) * bi * bj;
                k++;
            }
        }
        alpha.set(ind, 0, alpha0);

        // Compute H
        k = 2 * nvar + 1;
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            for (int j = 0; j < i; ++j)
            {
                if (_trainingSet->get_X_nbdiff(j) <= 1)
                    continue;

                const double alphaQCoeff = alpha.get(ind, k);
                H.set(i, j, alphaQCoeff);
                H.set(j, i, alphaQCoeff);
                k++;
            }
        }
        k = nvar + 1;
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            const double alphaQCoeff = alpha.get(ind, k);
            H.set(i, i, alphaQCoeff);
            ++k;
        }

        // Compute alphaL := D alphaL + 0.5 D H b + 0.5 b' H D
        k = 1;
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            const double alphaLi = alpha.get(ind, k);
            const double a = _trainingSet->get_X_scaling_a(i);
            alphaL.set(0, i, a * alphaLi);
            k++;
        }
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            double alphaLi = alphaL.get(0, i);
            const double a = _trainingSet->get_X_scaling_a(i);
            for (int j = 0; j < _n; ++j)
            {
                if (_trainingSet->get_X_nbdiff(j) <= 1)
                    continue;

                const double b = _trainingSet->get_X_scaling_b(j);
                alphaLi += 0.5 * H.get(i, j) * a * b;
            }
            alphaL.set(0, i, alphaLi);
        }
        for (int j = 0; j < _n; ++j)
        {
            if (_trainingSet->get_X_nbdiff(j) <= 1)
                continue;

            double alphaLj = alphaL.get(0, j);
            const double a = _trainingSet->get_X_scaling_a(j);
            for (int i = 0; i < _n; ++i)
            {
                if (_trainingSet->get_X_nbdiff(i) <= 1)
                    continue;

                const double b = _trainingSet->get_X_scaling_b(i);
                alphaLj += 0.5 * H.get(i, j) * a * b;
            }
            alphaL.set(0, j, alphaLj);
        }
        k = 1;
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            alpha.set(ind, k, alphaL.get(0, i));
            k++;
        }

        // Compute alphaQ corresponding to : Hq := D Hq D
        k = nvar + 1;
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            const double a = _trainingSet->get_X_scaling_a(i);
            const double alphaQcoeff = alpha.get(ind, k);
            alpha.set(ind, k, alphaQcoeff * a * a);
            k++;
        }
        for (int i = 0; i < _n; ++i)
        {
            if (_trainingSet->get_X_nbdiff(i) <= 1)
                continue;

            const double ai = _trainingSet->get_X_scaling_a(i);
            for (int j = 0; j < i ; ++j)
            {
                if (_trainingSet->get_X_nbdiff(j) <= 1)
                    continue;

                const double aj = _trainingSet->get_X_scaling_a(j);
                const double alphaQcoeff = alpha.get(ind, k);
                alpha.set(ind, k, alphaQcoeff * ai * aj);
                k++;
            }
        }
    }

    // Step 3: finally, a scaling is done on the output, i.e.
    // Q(x) := (Q(x) - zb) / za = (alpha0 - zb) / za + alphaL / za x + (0.5 / za) * x' H x.
    // Update the coefficients
    for (int ind = 0; ind < _m; ++ind)
    {
        const double za = _trainingSet->get_Z_scaling_a(ind);
        const double zb = _trainingSet->get_Z_scaling_b(ind);
        int k = 0;

        // Update alpha0
        double alphak = alpha.get(ind, k);
        alphak = (alphak - zb) / za;
        alpha.set(ind, k, alphak);
        k++;

        // Update alphaL
        for (int i = 0; i < nvar; ++i)
        {
            alphak = alpha.get(ind, k);
            alphak /= za;
            alpha.set(ind, k, alphak);
            k++;
        }

        // Update alphaQ
        for (int i = 0; i < nvar; ++i)
        {
            alphak = alpha.get(ind, k);
            alphak /= za;
            alpha.set(ind, k, alphak);
            k++;
        }
        for (int i = 0; i < nvar; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                alphak = alpha.get(ind, k);
                alphak /= za;
                alpha.set(ind, k, alphak);
                k++;
            }
        }
    }

    return alpha;
}

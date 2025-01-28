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
/**
 \file   QPSolverOptimize.hpp
 \brief  Class to create trial points by performing quadratic model optimization using QP solver
 \author Tangi Migot
 \see    QPSolverOptimize.cpp
 */
#ifndef __NOMAD_4_5_QP_SOLVER_OPTIMIZE__
#define __NOMAD_4_5_QP_SOLVER_OPTIMIZE__

#include "../../Algos/Step.hpp"
#include "../../Algos/QuadModel/QuadModelIterationUtils.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to create trial points by performing quadratic model optimization using QP solver
/**
 - Start, run and end tasks are performed.
 - Start: the quadratic model optimization problem is setup and solved by calling startImp. Call ::generateTrialPoints.
 - Run: trial (oracle) points are evaluated with EvalType::BB. Set the stop reason.
 - End: Remove from cache EvalType::MODEL only cache points.
 */
class QPSolverOptimize : public Step, public QuadModelIterationUtils
{
private:

    OutputLevel                         _displayLevel;

    ArrayOfDouble _modelLowerBound; ///> Lower bound: min of trainingSet points
    ArrayOfDouble _modelUpperBound; ///> Upper bound: max of trainingSet points
    Point         _modelFixedVar;   ///> Fixed variables: fixed variables detected from trainingSet

    Double        _modelBoxSizeLimit;  ///> Bounds box size limit to generate trial points
    
    Point _modelCenter;
    
    static EvalPointPtr _prevFeasRefCenter, _prevInfeasRefCenter; ///> Previous refCenter to identify change of pb.
    static Point _prevFeasXopt, _prevInfeasXopt; ///>  Previous QP solver solutions (1 feas, 1 infeas) can be used for initial points

    const std::shared_ptr<PbParameters> _refPbParams; ///< Reference to the original problem parameters.

    std::shared_ptr<PbParameters>       _optPbParams; ///< pb parameters for model optimization

    bool _optWithScaledBounds;
    
    int _n;
    int _m;  // number of bb outputs
    int _nbCons ; // number of constraints (!= _m-1 if multi-obj)
    size_t _quadModelMaxEval;
    
    BBOutputTypeList _bbot;

    bool _verbose, _verboseFull;
    
public:
    /// Constructor
    /* Parent must explicitely be a (pointer to a) QuadModelAlgo.
     * Run parameters will be recomputed for model optimization.
     */
    explicit QPSolverOptimize(const Step* parentStep,
                               const std::shared_ptr<PbParameters>               refPbParams,
                               bool optWithScaledBounds)
      : Step(parentStep),
      QuadModelIterationUtils (parentStep),
        _displayLevel(OutputLevel::LEVEL_DEBUG),
        _modelLowerBound(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelUpperBound(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelFixedVar(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _modelCenter(refPbParams->getAttributeValue<size_t>("DIMENSION"), Double()),
        _refPbParams(refPbParams),
        _optPbParams(nullptr),
        _optWithScaledBounds(optWithScaledBounds),
        _verbose(false),
        _verboseFull(false)
    {
        init();
    }

    /// Generate new points to evaluate
    /**
     - Setup the evaluator control parameters.
     - Manage display of sub-optimization.
     - Setup evaluator (EvalType::MODEL) and success type identification function.
     - Setup the bounds and fixed variables from the trainingSet of the quadratic model.
     - Setup run and pb parameters for Mads
     - Perform start, run and end tasks on Mads.
     - best feasible and best infeasible (if available) are inserted as trial points.
     */
    void generateTrialPointsImp() override;
    
    virtual void startImp() override; ///< The quadratic model optimization problem is setup and solved by calling startImp. Calls ::generateTrialPoints.
    virtual bool runImp() override; ///< Trial (oracle) points are evaluated with EvalType::BB. Set the stop reason.
    virtual void endImp() override {} ///< Remove from cache EvalType::MODEL only cache points.

        
private:
    void init();

    // Helpers
    void setModelBoundsAndFixedVar();

    // X = X + Y
    bool update(NOMAD::Point & X, const SGTELIB::Matrix & Y) { return update(X, Y, 1);}
    // X = X + a * Y
    bool update(NOMAD::Point & X, const SGTELIB::Matrix & Y, const double a) ;
    // return X + a * Y
    bool update(NOMAD::Point & Xup, const NOMAD::Point & X, const SGTELIB::Matrix & Y, const double a) ;
    
    
    // Projected conjugate gradient
    // Solve min (1/2) x' G x + c x
    //        x
    //        s.t  A x = b
    //             || x || <= delta
    bool projected_conjugate_gradient(SGTELIB::Matrix & x ,
                                      const SGTELIB::Matrix & A ,
                                      const SGTELIB::Matrix & b ,
                                      const SGTELIB::Matrix & G ,
                                      const SGTELIB::Matrix & c ,
                                      const double delta ,
                                      const SGTELIB::Matrix & x0 ,
                                      const double tol,
                                      const bool verbose = false);
    
    // Algo 1. Unconstrained solver: Modified Newton method:
    // Solve QP with bound constraints using a basic active-set strategy.
    bool solveBCQP(NOMAD::Point & X, const int max_iter = 10, const double atol = 1E-7, const double rtol = 1E-7, const bool verbose = false);
    bool solveBCQP(SGTELIB::Matrix & X, const SGTELIB::Matrix & H, const SGTELIB::Matrix & g, const double g0, const SGTELIB::Matrix & lvar, const SGTELIB::Matrix & uvar, const int max_iter = 100, const double atol = 1E-12, const double rtol = 1E-12, bool verbose = false);
    // void projectedGradient(SGTELIB::Matrix& X, const SGTELIB::Matrix& H, const SGTELIB::Matrix& g, const double g0,
    //                        const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar,
    //                        bool* active_l, bool* active_u, bool* working, SGTELIB::Matrix& d_k, SGTELIB::Matrix& Temp,
    //                        const double kappa = 0,
    //                        const double tol = 1E-12, const int max_iter = 100, const bool verbose = false);

    void projectedGradient(SGTELIB::Matrix& X, const SGTELIB::Matrix& H, const SGTELIB::Matrix& g, const double g0,
                           const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar,
                           bool* active_l, bool* active_u, SGTELIB::Matrix& d_k, SGTELIB::Matrix& Temp,
                           const double kappa = 0,
                           const double tol = 1E-12, const int max_iter = 100, const bool verbose = false);

    bool conjugateGradient(SGTELIB::Matrix& X, const SGTELIB::Matrix& H, const SGTELIB::Matrix& g,
                           const int maxIter = 100, const double xi = 1e-3,
                           const double atol = 1e-7, const double rtol = 1e-7,
                           const bool verbose = false) const;

    // Levenberg-Marquardt algorithm to find a feasible point
    size_t solveLM(
        NOMAD::Point& X,
        SGTELIB::Matrix& XS,
        const SGTELIB::Matrix& lvar,
        const SGTELIB::Matrix& uvar,
        SGTELIB::Matrix& cX,
        const double Fx,
        const double tol,
        const size_t maxIterInner,
        const double tolDistDXInner,
        const bool checkStrict,
        const bool verbose);
    bool getStrictFeasiblePoint(NOMAD::Point& X,
                                SGTELIB::Matrix& XS,
                                SGTELIB::Matrix& lvar,
                                SGTELIB::Matrix& uvar,
                                const SGTELIB::Matrix& cX) const;

    // Projected Armijo/Wolfe linesearch
    double projected_armijo(const SGTELIB::Matrix & X, const SGTELIB::Matrix & H, const SGTELIB::Matrix & g, const double g0, const SGTELIB::Matrix & lvar, const SGTELIB::Matrix & uvar, const SGTELIB::Matrix & d_k, const double fk, const double slope, SGTELIB::Matrix & Xp, SGTELIB::Matrix & gradientF_kp, const double tmax = 1E20) ;
    // Check whether working set is a subset of binding
    // If false, set the first encountered component to false in working.
    bool check_subset_binding_update(bool* working, const bool* binding, const size_t n);

    // Compute a solution of {d \in argmin_d q(x + d) s.t. d_active = 0, \|d_free\|_2 <= Delta}
    void solve_TR_constrained_QP(SGTELIB::Matrix* d, const SGTELIB::Matrix& X, const SGTELIB::Matrix& H, const SGTELIB::Matrix& g,
                                 SGTELIB::Matrix& grad, const bool* active, const double Delta, const bool verbose = false);
    // Compute a solution of {d \in argmin_d q(x + d) s.t. \|d\|_2 <= Delta}
    void solve_TR_constrained_QP(SGTELIB::Matrix * d, const SGTELIB::Matrix & X, const SGTELIB::Matrix & H, const SGTELIB::Matrix & g, SGTELIB::Matrix & grad, const double Delta);
    // Compute a solution of {d \in argmin_d q(x + d) s.t. d_active = 0}
    void solve_unconstrained_QP(SGTELIB::Matrix* d, const SGTELIB::Matrix& X, const SGTELIB::Matrix& H, const SGTELIB::Matrix& g,
                                SGTELIB::Matrix& grad, const bool* active)
    {
        return solve_TR_constrained_QP(d, X, H, g, grad, active, 1E15);
    }
    // Compute a solution of {d \in argmin_d q(x + d)}
    void solve_unconstrained_QP(SGTELIB::Matrix * d, const SGTELIB::Matrix & X, const SGTELIB::Matrix & H, const SGTELIB::Matrix & g, SGTELIB::Matrix & grad) { return solve_TR_constrained_QP(d, X, H, g, grad, 1E15);}

    // Algo 2.
    // Solve QP with constraint. Method L1 Augmented Lagrangian.
    bool solveL1AugLag(NOMAD::Point & X, const int max_iter = 9,
                       const double atol = 1e-7, const double rtol = 1e-7,
                       const bool verbose = false);
    double getPenalizedL1AugLagModelObj(const NOMAD::Point & X, const SGTELIB::Matrix & cons, const SGTELIB::Matrix & lambda, double mu) const;
    double getPenalizedL1AugLagModelObj(const NOMAD::Point & X, const SGTELIB::Matrix & lambda, double mu) const;
    double piecewise_line_search(const NOMAD::Point& X, const SGTELIB::Matrix& d,
                                 const bool* active, const bool* feasible, const bool* infeasible,
                                 const bool* active_lb, const bool* active_ub,
                                 const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar,
                                 const SGTELIB::Matrix& lambda, const double mu,
                                 const double small_gamma = 1E-20, const double gamma_update = 1.5, const double delta = 1E-4) const;
    SGTELIB::Matrix get_pseudo_multiplier(const bool * active, const bool * feasible, const bool * infeasible, const SGTELIB::Matrix & lambda, const double mu) const;
    // bool compute_horizontal_step(const NOMAD::Point& X, SGTELIB::Matrix& h_k, const SGTELIB::Matrix& Jacobian_k,
    //                              const bool* active, const bool* feasible, const bool* infeasible,
    //                              const SGTELIB::Matrix& lambda_l, const double mu_l) const ;
    bool compute_horizontal_step(const NOMAD::Point& X, SGTELIB::Matrix& h_k, const SGTELIB::Matrix& Jacobian_k,
                                 const bool* active, const bool* feasible, const bool* infeasible, const bool* active_lu, const bool* active_ub,
                                 const SGTELIB::Matrix& lambda_l, const double mu_l) const ;
    bool compute_vertical_step(const NOMAD::Point& X, SGTELIB::Matrix& v_k, const SGTELIB::Matrix& activeJacobian_k, const SGTELIB::Matrix& cons,
                               const bool* active, const bool* active_lb, const bool* active_ub,
                               const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar) const ;
    bool compute_drop_constraint_step(SGTELIB::Matrix& d,
                                      const SGTELIB::Matrix& activeJacobian,
                                      const SGTELIB::Matrix& activeMultipliers,
                                      const SGTELIB::Matrix& pseudoGradient_k,
                                      const double mu) const;
    // Return true if there are less than n active bound constraints, where n is the number of variables
    bool check_active_bound_constraints(const NOMAD::Point& X,
                                        bool* active_lb,
                                        bool* active_ub,
                                        SGTELIB::Matrix& tolBounds,
                                        const SGTELIB::Matrix& lvar,
                                        const SGTELIB::Matrix& uvar) const;
    // Return true if there are less than n - nbActiveBounds constraints, where n is the number of variables
    bool check_active_constraints(const SGTELIB::Matrix& cons,
                                  const int nbActiveBounds,
                                  bool* active,
                                  bool* feasible,
                                  bool* infeasible,
                                  double& innerPrecision) const;
    double check_inner_success(NOMAD::Point & X, const SGTELIB::Matrix & Jacobian_k, SGTELIB::Matrix & multiplier_k, const SGTELIB::Matrix & lambda, const double mu, const bool * active, const bool * infeasible) const ;

    // Algo 3.
    // Augmented Lagrangian
    bool solveAugLag(NOMAD::Point& X,
                     const int max_iter = 30,
                     const double tolDistDX = -1.0,
                     const double atol = 1E-7,
                     const double rtol = 1E-7,
                     const double mu0 = 0.5,
                     const double muDecrease = 2.0,
                     const double eta0 = 1.0,
                     const double omega0 = 1.0,
                     const double successRatio = 0.99,
                     const size_t maxIterInner = 50,
                     const double tolDistDXInner = 1E-15,
                     const size_t maxSuccessiveFail = 3);
    
    int solveBoundAugLag(SGTELIB::Matrix& XSp, const SGTELIB::Matrix& XS,
                         const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar,
                         const SGTELIB::Matrix& lambda, const double omega, const double mu,
                         SGTELIB::Matrix& GradPk, SGTELIB::Matrix& HessPk,
                         const size_t maxIterInnerLoop = 50, const double tolerance_distDX = 1E-15, const bool verbose = false) ;
    double compute_AugLag_TR_ared(const SGTELIB::Matrix & XS, const SGTELIB::Matrix & XSp, const SGTELIB::Matrix & lambda, double mu) const;

// Augmented Lagrangian model getters
    double getAugLagModelObj(const SGTELIB::Matrix & XS, const SGTELIB::Matrix & lambda, const double mu) const;
    double getAugLagModelObj(const SGTELIB::Matrix & XS, const SGTELIB::Matrix & cons, double fx, const SGTELIB::Matrix & lambda, const double mu) const;
    SGTELIB::Matrix getAugLagModelGrad(const SGTELIB::Matrix & XS, const SGTELIB::Matrix & lambda, double mu) const;
    void getAugLagModelGrad(SGTELIB::Matrix * lagGrad, const SGTELIB::Matrix & XS, const SGTELIB::Matrix & lambda, const double mu) const;
    SGTELIB::Matrix getAugLagModelHess(const SGTELIB::Matrix & XS, const SGTELIB::Matrix & lambda, const double mu) const;
    void getAugLagModelHess(SGTELIB::Matrix * lagHess, const SGTELIB::Matrix & XS, const SGTELIB::Matrix & lambda, const double mu) const;

// Algo 4. Trust region interior point
    bool solveTRIPM(NOMAD::Point& X, const int max_iter = 30,
                    const double tolDistDX = -1.0,
                    const double atol = 1E-7, const double rtol = 1E-7, const double mu0 = 0.5,
                    const double muDecrease = 2.0, const size_t maxIterInner = 40,
                    const bool verbose = true,
                    const bool verbose_SolverBarrier = false);
    int solver_barrier(NOMAD::Point& X,
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
                       const bool verbosePCG = false);
    double errorTRIPM(
        const SGTELIB::Matrix & XS,
        const SGTELIB::Matrix & lvar,
        const SGTELIB::Matrix & uvar,
        const SGTELIB::Matrix & lambda,
        const SGTELIB::Matrix & cx, 
        const double mu );
    void compute_slack_multiplier(
        SGTELIB::Matrix & y,
        const SGTELIB::Matrix & XS,
        const SGTELIB::Matrix & Jx,
        const SGTELIB::Matrix & Gx,
        const double mu) ;
    double merit_function_barrier(
        const NOMAD::Point & X,
        const SGTELIB::Matrix & XS,
        const SGTELIB::Matrix & lvar,
        const SGTELIB::Matrix & uvar,
        const double mu,
        const double nu
    ) ;
// Optimality functions

    double check_optimality_bounds(NOMAD::Point & X_k, const SGTELIB::Matrix & gradientF_k) ;

    /// Compute the norm ||X - P_Omega(X - gradient)||, where Omega = {y,lvar <= y <= yvar}.
    double check_optimality_bounds(const SGTELIB::Matrix& X, const SGTELIB::Matrix& gradient,
                                   const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar) ;

    /// Compute the norm ||X - P_Omega(X - gradient)||, where Omega = {y,lvar <= y <= yvar}.
    /// Save X - P_Omega(X - gradient) in dual_feas.
    double check_optimality_bounds(const SGTELIB::Matrix& X, const SGTELIB::Matrix& gradient,
                                   const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar, SGTELIB::Matrix& dual_feas);
    double compute_dual_residual(const SGTELIB::Matrix & Grad_k, const SGTELIB::Matrix & Jacobian_k, const SGTELIB::Matrix & multiplier_k) const ;

// Basic optimization functions

    bool computeNewtonDirection(const SGTELIB::Matrix& g, int* pp, double** D, double** L, double* sol, const int nfree) const;
    // We compute the eigenvector corresponding to it
    bool InverseIteration(SGTELIB::Matrix * sol, const SGTELIB::Matrix & HW, const double eigmin, const int nfree, const double tol = 1E-15, const bool verbose = false) ;
    // Compute in-place a solution of {d \in argmin_d q(x + d) s.t. d_active = 0, \|d_free\|_2 <= Delta} for convex quadratic
    bool Convex_TR_QP(SGTELIB::Matrix* d, const SGTELIB::Matrix& g, const SGTELIB::Matrix& gW, const SGTELIB::Matrix& H,
                      const SGTELIB::Matrix& HW, int* pp, double** D, double** L, const bool* active, const double Delta,
                      const bool verbose) ;
    
// Model getters with NOMAD points
    SGTELIB::Matrix getModelOut(const Point & x) const ;
    // Get model output j at point x
    double getModelOut(const Point & x, int j ) const { return getModelOut(x).get(0,j); }
    // Get gradient at point x of model output j - size nvar x 1
    // SGTELIB::Matrix getModelGrad(const Point & x, int j) const { return getModelJacobian(x).get_row(j).transpose(); }
    // Get model objective at point x
    double getModelObj(const Point & x) const ;
    // Get model objective gradient at point x (size: n x 1)
    SGTELIB::Matrix getModelGrad(const Point & x) const ;
    void getModelGrad(SGTELIB::Matrix * Gx, const Point & x) const ;
    // Get model output j hessian at point x
    SGTELIB::Matrix getModelHessian(const Point& x, int j) const ;
    // Get model objective Hessian at point x (size: n x n)
    SGTELIB::Matrix getModelHessian(const Point& x) const ;
    // Get constraint function at point x (size ncon x 1)
    SGTELIB::Matrix getModelCons(const Point& x) const ;
    void getModelCons(SGTELIB::Matrix* cons, const Point& x) const ;
    // Get jacobian at point x (size ncon x nvar)
    SGTELIB::Matrix getModelJacobian(const Point & x) const ;
    // Get model lagrangian at point x
    double getModelLag(const Point& x, const SGTELIB::Matrix& lambda, const double sigma = 1.0) const;
    // Get gradient of model of lagrangian (given multipliers) at point x
    SGTELIB::Matrix getModelLagGradient(const Point& x, const SGTELIB::Matrix& lambda, const double sigma = 1.0) const;
    // Get hessian of model of lagrangian (given multipliers) at point x
    SGTELIB::Matrix getModelLagHessian(const Point& x, const SGTELIB::Matrix& lambda, const double sigma = 1.0) const;

// Model getters for quadratic functionals
    double getModelObj(const SGTELIB::Matrix & x, const SGTELIB::Matrix & H, const SGTELIB::Matrix & g, const double g0 = 0.0) const ;
    SGTELIB::Matrix getModelGrad(const SGTELIB::Matrix & x, const SGTELIB::Matrix & H, const SGTELIB::Matrix & g) const ;
    void getModelGrad(SGTELIB::Matrix * Gx, const SGTELIB::Matrix & x, const SGTELIB::Matrix & H, const SGTELIB::Matrix & g) const ;

// Utils
    // Number of `true` in `indices` of length `len`
    int sum(const bool * indices, int len) const;
    void lencheck(const int n, const SGTELIB::Matrix & X) const ;
    void sizecheck(const int m, const int n, const SGTELIB::Matrix & X) const ;
    bool check_strict_feasible(const SGTELIB::Matrix & X, const SGTELIB::Matrix & lvar, const SGTELIB::Matrix & uvar) const;
    SGTELIB::Matrix vector_subset(const SGTELIB::Matrix & X, const bool * active) const;
    void vector_broadcast(SGTELIB::Matrix * X, const SGTELIB::Matrix & Xsub, const bool * active) ;
    SGTELIB::Matrix matrix_subset(const SGTELIB::Matrix & X, const bool * active) const;

    // Return the maximum possible step w.r.t to the bounds along the direction d from a point X.
    // X needs to satisfy the bounds.
    double max_step_bounds(const NOMAD::Point & X, const SGTELIB::Matrix & d);
    double max_step_bounds(const SGTELIB::Matrix & X, const SGTELIB::Matrix & lvar, const SGTELIB::Matrix & uvar, const SGTELIB::Matrix & d);
    // Project the direction over bounds
    void project_bounds(NOMAD::Point & X_k, SGTELIB::Matrix & d_k) ;
    void project_bounds(SGTELIB::Matrix& d_k, const bool* active) ;
    void project_bounds(const SGTELIB::Matrix & X, const SGTELIB::Matrix & lvar, const SGTELIB::Matrix & uvar, SGTELIB::Matrix & d_k) ;
    // The vector active_l (resp. active_u) has a `true` i-th component if the lower (resp. upper) bound is active at X.
    void active_bounds(const SGTELIB::Matrix & X, bool * active_l, bool * active_u);
    void active_bounds(const SGTELIB::Matrix & X, const SGTELIB::Matrix & lvar, const SGTELIB::Matrix & uvar, bool * active_l, bool * active_u, double tol = 1E-15);
    // Return the indices of active lower or upper bounds that are binding.
    void binding_bounds(SGTELIB::Matrix & Grad, const bool * active_l, const bool * active_u, bool * binding);
    
    // Return true for indices of active constraints
    void getModelActiveCons(const SGTELIB::Matrix &cons, const double tol, bool * indices) const;
    // Return true for indices of strictly feasible constraints
    void getModelFeasibleCons(const SGTELIB::Matrix &cons, const double tol, bool * indices) const;
    // Return true for indices of strictly feasible constraints
    void getModelInfeasibleCons(const SGTELIB::Matrix &cons, const double tol, bool * indices) const;
    // Return true for indices of strictly feasible constraints
    SGTELIB::Matrix feasibility(const SGTELIB::Matrix &cons) const;
    void feasibility(SGTELIB::Matrix * feas, const SGTELIB::Matrix &cons) const;
    bool isFeasible(const SGTELIB::Matrix &cons, const double tol) const;
    // Get jacobian of active constraints at point x (size ncon x nvar)
    SGTELIB::Matrix getModelActiveJacobian(const SGTELIB::Matrix & Jacobian, const bool * indices) const ;
    void getModelActiveJacobian(SGTELIB::Matrix * activeJacobian, const SGTELIB::Matrix & Jacobian, const bool * indices) const ;
    
    void snapToBounds(SGTELIB::Matrix& x, const SGTELIB::Matrix& lvar, const SGTELIB::Matrix& uvar) const;

    // Return the QPModel matrix corresponding to the QCQP
    SGTELIB::Matrix computeQPModelMatrix();
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_QP_SOLVER_OPTIMIZE__

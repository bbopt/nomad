from ..doe import DOEgenerator, LatinHypercubeDOE, DOEMixin


import numpy as np
from typing import Any, Dict, Optional

from pymoo.core.problem import ElementwiseProblem
from pymoo.optimize import minimize
from pymoo.core.variable import Real, Integer, Binary, Choice
from pymoo.core.mixed import MixedVariableGA
from pymoo.termination import get_termination



class BaseOptimizer(DOEMixin):
    """
    Generic optimizer run parameters.
    """

    def __init__(
        self,
        seed: int = 0,
        initialDOEgenerator: Optional[DOEgenerator] = None,
        quantPollType: str = "ADS",
        isBOSearchUsed: bool = False,   # Attribute passed to a DataManager() in NomadBBO_src.pyx
        modelUpdateHardCap: int = 500,  # Attribute passed to a DataManager() in NomadBBO_src.pyx
        efiStopTol: float = 1e-3,
        efiStopPatience: int = 3
        ):

        # Basic parameter
        self.seed = int(seed)

        self.quantPollType = str(quantPollType).upper().strip()
        allowed_quant_poll = {"ADS", "MADS"}
        if self.quantPollType not in allowed_quant_poll:
            raise ValueError(
                f"Invalid quantPollType='{self.quantPollType}'. "
                f"Expected one of {sorted(allowed_quant_poll)}."
            )


        # BO parameters 
        self.isBOSearchUsed = bool(isBOSearchUsed)
        self.efiStopTol = float(efiStopTol)
        self.efiStopPatience = int(efiStopPatience)
        self.efiStopCount = 0            # Counter initially set to 0: when it reaches efiStopPatience, BO and model update are stopped. Updated from MixedOptimizer()
        self.logEFI_max_found = -np.inf  # Best EFI found , updated each BO step
        
        # Model parameter
        self.modelUpdateHardCap = int(modelUpdateHardCap)


        # DoE
        if initialDOEgenerator is not None and not isinstance(initialDOEgenerator, DOEgenerator):
            raise TypeError(
                f"initialDOEgenerator must be an instance of DOEgenerator, "
                f"got {type(initialDOEgenerator).__name__}."
            )

        if initialDOEgenerator is None:
            initialDOEgenerator = LatinHypercubeDOE(seed=self.seed)

        self.initialDOEgenerator = initialDOEgenerator

        # Initialized later in optimize() from NomadBBO.pyx
        self.problemDefinition = None
        self.dataManager = None


    @property
    def idx_by_type(self):
        return self.problemDefinition.getVariableIndicesByType()

    @property
    def cat_idx(self):
        return self.idx_by_type.get("choice", []) or []

    @property
    def bin_idx(self):
        return self.idx_by_type.get("binary", []) or []

    @property
    def int_idx(self):
        return self.idx_by_type.get("integer", []) or []

    @property
    def real_idx(self):
        return self.idx_by_type.get("real", []) or []

    @property
    def quant_idx(self):
        return self.int_idx + self.real_idx

    @property
    def nbVars(self):
        return self.problemDefinition.getDimension()

    @property
    def nbCat(self):
        return len(self.cat_idx)

    @property
    def nbBin(self):
        return len(self.bin_idx)

    @property
    def nbInt(self):
        return len(self.int_idx)

    @property
    def nbReal(self):
        return len(self.real_idx)

    @property
    def nbQuant(self):
        return len(self.quant_idx)

    @property
    def hasCategoricalVariables(self):
        return self.nbCat > 0

    @property
    def hasBinaryVariables(self):
        return self.nbBin > 0

    @property
    def hasIntegerVariables(self):
        return self.nbInt > 0

    @property
    def hasRealVariables(self):
        return self.nbReal > 0

    @property
    def hasQuantitativeVariables(self):
        return self.nbQuant > 0

    @property
    def nbCategoriesPerVariable(self):
        return self.problemDefinition.getNbCategoriesPerCategoricalVariable() or []

    @property
    def quant_bounds(self):
        return self.problemDefinition.getQuantitativeBounds() or []
    

    # --------------------------------------------------------------- #
    # Search
    # --------------------------------------------------------------- #
    def _best_feasible_objective(self, eps_feas: float = 1e-12) -> float:
        """
        Best (min) objective among feasible cached points.

        If no feasible exists yet (H > eps_feas for all points), return the objective
        value of the point with the smallest H (constraint violation aggregate).
        Ties on H are broken by smaller objective.

        Unconstrained: returns the best overall objective.

        Parameters
        ----------
        eps_feas : float, default=1e-12
            Feasibility threshold for H == 0.
        """

        if self.dataManager.cache.empty or len(self.dataManager.OBJ_cols) != 1:
            return np.inf

        f = self.dataManager.getObjValues(flatten=True)
        if f.size == 0:
            return np.inf

        # Unconstrained => best overall
        if len(self.dataManager.PB_cols) == 0:
            return float(np.min(f))

        H = self.dataManager.getHValues()
        if H.size == 0:
            # Shouldn't happen if PB exists, but be safe
            return float(np.min(f))

        feas = (H <= eps_feas)
        if np.any(feas):
            return float(np.min(f[feas]))

        # No feasible yet => take point with smallest H (tie-break by f)
        Hmin = np.min(H)
        cand = np.where(H <= Hmin + 1e-15)[0]  # tiny tol for numerical ties
        if cand.size == 1:
            return float(f[cand[0]])

        # tie-break by objective
        j = cand[np.argmin(f[cand])]
        return float(f[j])


    def _make_pymoo_var_dict(self) -> Dict[str, Any]:
        """
        Build PyMoo variable dictionary in NOMAD order x0..x{nbVars-1}.

        IMPORTANT:
        - Use bounds=(lb,ub) for Real/Integer/Binary for pymoo compatibility.
        """
        int_idx = self.int_idx
        real_idx = self.real_idx

        quant_bound_map = {j: self.quant_bounds[k] for k, j in enumerate(self.quant_idx)}

        var_dict: Dict[str, Any] = {}
        cat_local = 0

        for j in range(self.nbVars):
            name = f"x{j}"

            if j in self.cat_idx:
                Lj = int(self.nbCategoriesPerVariable[cat_local])
                var_dict[name] = Choice(options=list(range(Lj)))
                cat_local += 1

            elif j in self.bin_idx:
                # explicit bounds for pymoo Variable API
                var_dict[name] = Binary(bounds=(0, 1))

            elif j in int_idx:
                if j not in quant_bound_map:
                    raise ValueError(f"Missing bounds for integer variable x{j}.")
                lb, ub = quant_bound_map[j]
                var_dict[name] = Integer(bounds=(int(lb), int(ub)))

            elif j in real_idx:
                if j not in quant_bound_map:
                    raise ValueError(f"Missing bounds for real variable x{j}.")
                lb, ub = quant_bound_map[j]
                var_dict[name] = Real(bounds=(float(lb), float(ub)))

            else:
                raise ValueError(f"Variable x{j} not classified in choice/binary/integer/real indices.")

        return var_dict


    def _dict_to_numpy_point(self, x_dict: Dict[str, Any]) -> np.ndarray:
        """
        Convert PyMoo mixed decision dict -> numpy point in NOMAD order (float array).
        """
        x = np.zeros((self.nbVars,), dtype=float)
        for j in range(self.nbVars):
            x[j] = float(x_dict[f"x{j}"])
        return x


    def _convert_point_to_smt(self, x: np.ndarray) -> np.ndarray:
        """
        SMT wants discrete columns as integers. Return shape (1, nbVars).
        """
        X = np.atleast_2d(np.asarray(x, dtype=float).copy())

        if self.cat_idx:
            X[:, self.cat_idx] = np.rint(X[:, self.cat_idx]).astype(int)
        if self.bin_idx:
            X[:, self.bin_idx] = np.rint(X[:, self.bin_idx]).astype(int)

        if self.int_idx:
            X[:, self.int_idx] = np.rint(X[:, self.int_idx]).astype(int)

        return X


    def BOSearch(self, incumbentFeas: np.ndarray, incumbentInfeas: np.ndarray, budgetPerVariable: int = 250,   
                         ) -> np.ndarray:
        """
        Bayesian Optimization step using EFI = EI * Π PoF (log-stable),
        solved as an UNCONSTRAINED mixed-variable optimization problem with PyMoo.

        - Mixed variables (Choice/Binary/Integer/Real) handled natively by PyMoo
        - No one-hot, no relaxation of decision variables
        - Surrogate constraints are aggregated via probability of feasibility
        - Objective for GA: maximize log(EI) + sum_i log(PoF_i)
        implemented as minimize -(logEFI)

        Returns
        -------
        np.ndarray
            trialPoints of shape (1, nbVars) best candidate, or empty if not ready.
        """

        trialPoints = np.empty((0, self.nbVars), dtype=np.float64)

        if not self.isBOSearchUsed:
            # Return empty point
            return trialPoints
        
        else:

            # Model updates may have been disabled by constructModel()
            # after reaching modelUpdateHardCap.
            #
            # BO requires repeated model updates while it is active, so
            # reaching the hard cap also stops the BO search.
            if not self.dataManager.isModelUpdated:
                self.isBOSearchUsed = False

                # Recompute whether a surrogate model is still required after disabling BO.
                self.dataManager.isModelRequired = self.isModelRequired

                return trialPoints

            from scipy.stats import norm
            from scipy.special import erfcx, expm1, log1p

            # Need surrogate models ready
            if (not self.dataManager.modelIsReady) or (self.dataManager.gp_models is None):
                return trialPoints

            sm_f = self.dataManager.gp_models.get("objective", None)
            if sm_f is None:
                return trialPoints

            sm_ctrs = self.dataManager.gp_models.get("PB", [])  # PB constraints used for PoF
            nb_ctrs = len(self.dataManager.PB_cols)

            if len(sm_ctrs) != nb_ctrs:
                raise RuntimeError(f"PB GP models mismatch: {len(sm_ctrs)} models vs PB_cols={nb_ctrs}.")

            # ------------------------------------------------------------
            # Reference objective value f_min (minimization): use incumbents first,
            # else fallback to best feasible in cache; if none feasible, fallback to best infeasible.
            # ------------------------------------------------------------
            def _objective_of_point_from_cache(x: np.ndarray) -> float:
                row = self.dataManager._find_cache_index_of_point(x, atol=0.0)
                return float(self.dataManager.cache.loc[row, self.dataManager.OBJ_cols[0]])

            f_min = None

            if incumbentFeas is not None and np.size(incumbentFeas) == self.nbVars:
                try:
                    f_min = _objective_of_point_from_cache(np.asarray(incumbentFeas, dtype=float).ravel())
                except Exception:
                    f_min = None

            if f_min is None and incumbentInfeas is not None and np.size(incumbentInfeas) == self.nbVars:
                try:
                    f_min = _objective_of_point_from_cache(np.asarray(incumbentInfeas, dtype=float).ravel())
                except Exception:
                    f_min = None

            if f_min is None:
                f = self.dataManager.getObjValues(flatten=True)
                if f.size == 0:
                    return trialPoints

                if len(self.dataManager.PB_cols) == 0:
                    f_min = float(np.min(f))
                else:
                    H = self.dataManager.getHValues()
                    feas = (H <= 1e-8)
                    if np.any(feas):
                        f_min = float(np.min(f[feas]))
                    else:
                        infeas = (H > 1e-8)
                        if np.any(infeas):
                            f_min = float(np.min(f[infeas]))
                        else:
                            f_min = float(np.min(f))

            # ------------------------------------------------------------
            # Stable log-EI (based on your log implementation)
            # ------------------------------------------------------------
            def _log_ei(mu: float, var: float, fref: float) -> float:
                var = float(var)
                if (not np.isfinite(var)) or var < 1e-8:
                    return -1e12  # very small

                sigma = np.sqrt(max(var, 1e-16))
                z = (float(fref) - float(mu)) / sigma

                # Stable log(EI) using the "Unexpected Improvements..." style
                # Constants
                cte1 = np.log(2.0 * np.pi) / 2.0
                cte2 = np.log(np.pi / 2.0) / 2.0
                eps = 1e-16

                def log1mexp(x):
                    # Numerically stable log(1 - exp(x)) for x < 0
                    log2 = np.log(2.0)
                    if x >= 0:
                        # outside expected range; degrade gracefully
                        return np.log(max(1e-300, 1.0 - np.exp(-1e-12)))
                    if x > -log2:
                        return np.log(-expm1(x))
                    return log1p(-np.exp(x))

                # log_h part
                if z > -1.0:
                    h = norm.pdf(z) + z * norm.cdf(z)
                    if (not np.isfinite(h)) or h <= 0.0:
                        return -1e10
                    log_h = np.log(h)
                elif -1.0 / np.sqrt(eps) < z <= -1.0:
                    # uses erfcx for stability
                    t = np.log(erfcx(-z / np.sqrt(2.0)) * abs(z)) + cte2
                    log_h = -z * z / 2.0 - cte1 + log1mexp(t)
                else:
                    log_h = -z * z / 2.0 - cte1 - 2.0 * np.log(abs(z))

                log_sigma = 0.5 * np.log(max(var, 1e-16))
                return float(log_h + log_sigma)

            # ------------------------------------------------------------
            # log(PoF): sum_i log Phi( (0 - mu_i)/sigma_i )
            # ------------------------------------------------------------
            def _log_pof(Xsmt: np.ndarray) -> float:
                if nb_ctrs == 0:
                    return 0.0

                logp = 0.0
                for sm_c in sm_ctrs:
                    mu_c = float(sm_c.predict_values(Xsmt).reshape(-1)[0])
                    var_c = float(sm_c.predict_variances(Xsmt).reshape(-1)[0])
                    var_c = max(var_c, 1e-16)
                    sigma_c = float(np.sqrt(var_c))

                    z = (0.0 - mu_c) / sigma_c
                    p = float(norm.cdf(z))

                    # clamp to avoid -inf
                    logp += float(np.log(max(p, 1e-16)))

                return float(logp)

            # ------------------------------------------------------------
            # PyMoo problem: UNCONSTRAINED maximize logEFI = logEI + logPoF
            # ------------------------------------------------------------
            var_dict = self._make_pymoo_var_dict()

            class _EFIProblem(ElementwiseProblem):
                def __init__(self, outer):
                    super().__init__(vars=var_dict, n_obj=1, n_ieq_constr=0)
                    self.outer = outer

                def _evaluate(self, x, out, *args, **kwargs):
                    x_np = self.outer._dict_to_numpy_point(x)
                    Xsmt = self.outer._convert_point_to_smt(x_np)

                    mu = float(sm_f.predict_values(Xsmt).reshape(-1)[0])
                    var = float(sm_f.predict_variances(Xsmt).reshape(-1)[0])

                    logei = _log_ei(mu, var, f_min)
                    logpof = _log_pof(Xsmt)

                    logefi = logei + logpof

                    # maximize logefi  <=>  minimize -logefi
                    out["F"] = np.array([-logefi], dtype=float)

            problem = _EFIProblem(self)

            # ----------------------------
            # Algorithm + termination
            # ----------------------------
            algorithm = MixedVariableGA(pop_size=50)
            n_eval = int(max(500, budgetPerVariable * self.nbVars))

            res = minimize(
                problem,
                algorithm,
                termination=get_termination("n_eval", n_eval),
                seed=self.seed,
                verbose=False
            )

            # No point found
            if res is None or res.X is None:
                return trialPoints

            # res.X is a dict for mixed-variable problems
            x_best = self._dict_to_numpy_point(res.X)

            # Final cleanup: enforce discreteness
            Xsmt = self._convert_point_to_smt(x_best)
            x_best_clean = Xsmt.reshape(-1).astype(float)

            # Update best logEFI_max found
            logefi_current = float(-res.F[0])   # minus sign since we're minizing with Pymoo, but it's a maximization problem

            # ------------------------------------------------------------------ #
            # EFI stopping criterion for BO search stoppage
            previous_logEFI_max = self.logEFI_max_found

            # Compare current logEFI with the best logEFI found before this BO call
            if np.isfinite(previous_logEFI_max):
                relative_decay = logefi_current - previous_logEFI_max
                log_threshold = np.log(max(self.efiStopTol, 1e-300))

                print("LOG EFI VALUE:", logefi_current, " vs previous best ", previous_logEFI_max)

                if relative_decay < log_threshold:
                    self.efiStopCount += 1
                else:
                    self.efiStopCount = 0

            # Update best logEFI after the comparison
            self.logEFI_max_found = max(self.logEFI_max_found, logefi_current)

        # Stop BO search when EFI has been insufficient for
        # efiStopPatience consecutive BO iterations.
        if self.efiStopCount >= self.efiStopPatience:

            self.isBOSearchUsed = False

            # Stop further surrogate-model updates.
            self.dataManager.isModelUpdated = False

            # Recompute whether a surrogate model is still required after disabling BO.
            self.dataManager.isModelRequired = self.isModelRequired

        return np.ascontiguousarray(x_best_clean.reshape(1, -1), dtype=np.float64)
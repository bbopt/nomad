from ..doe import DOEgenerator
from .base_optimizer import BaseOptimizer

import itertools
import warnings
import numpy as np
from scipy.stats import wasserstein_distance
from typing import Any, Dict, Optional



class MixedOptimizer(BaseOptimizer):
    """
    Optimizer configuration for mixed-variable problems.
    """

    def __init__(
        self,
        seed: int = 0,
        initialDOEgenerator: Optional[DOEgenerator] = None,
        quantPollType: str = "ADS",
        catPollType: str = "QUANTILE",
        nbOfNeighbors: Optional[int] = None,     # Default rule based on problem dimension provided in NomadBBO_src.pyx
        catKernelType: Optional[str] = None,
        extendedPollTrigger: float = 0.05,
        updateModel = False,                     # Initial model-update policy passed to a DataManager() in NomadBBO_src.pyx
        isBOSearchUsed: bool = False,            # Passed to a DataManager() in NomadBBO_src.pyx
        modelUpdateHardCap: int = 500,           # Passed to BaseOptimizer(), then to a DataManager() in NomadBBO_src.pyx
        efiStopTol: float = 1e-3,                # Passed to BaseOptimizer()
        efiStopPatience: int = 3):               # Passed to BaseOptimizer()
        super().__init__(
            seed=seed,
            initialDOEgenerator=initialDOEgenerator,
            quantPollType=quantPollType,
            isBOSearchUsed=isBOSearchUsed,
            modelUpdateHardCap=modelUpdateHardCap,
            efiStopTol=efiStopTol,
            efiStopPatience=efiStopPatience)

        # Cat poll type
        self.catPollType = str(catPollType).upper().strip()
        allowed_poll = {"SURROGATE", "QUANTILE", "EXHAUSTIVE"}
        if self.catPollType not in allowed_poll:
            raise ValueError(
                f"Invalid catPollType='{self.catPollType}'. "
                f"Expected one of {sorted(allowed_poll)}."
            )

        # Extended poll is used by default with 5%
        self.extendedPollTrigger = float(extendedPollTrigger) 

        # Model update policy
        if self.isBOSearchUsed:

            # BO requires repeated model updates while it is active.
            if updateModel is False:
                warnings.warn(
                    "updateModel=False is ignored because Bayesian optimization "
                    "requires repeated surrogate-model updates.",
                    UserWarning)

            self.updateModel = True

        elif self.catPollType == "SURROGATE":

            # Without BO, the user decides whether the surrogate categorical poll
            # uses a static model or keeps updating it.
            self.updateModel = False if updateModel is None else bool(updateModel)

        else:

            # Without BO, QUANTILE and EXHAUSTIVE do not use surrogate models.
            if updateModel is True:
                warnings.warn(
                    f"updateModel=True is ignored because Bayesian optimization is not used "
                    f"and catPollType='{self.catPollType}' does not use surrogate models.",
                    UserWarning)

            self.updateModel = False

        # Categorical kernel is only needed when a surrogate model is required.
        if self.isModelRequired:
            if catKernelType is None:
                catKernelType = "GOWER"

            self.catKernelType = str(catKernelType).upper().strip()

            allowed_kernel = {"CONT_RELAX", "HOMO_HSPHERE", "GOWER"}
            if self.catKernelType not in allowed_kernel:
                raise ValueError(
                    f"Invalid catKernelType='{self.catKernelType}'. "
                    f"Expected one of {sorted(allowed_kernel)}.")
        else:
            self.catKernelType = None

        # If the nbOfNeighbors is None, it will be set in the optimize() call in NomadBBO.pyx where it is binded with the problem
        self.nbOfNeighbors = None if nbOfNeighbors is None else int(nbOfNeighbors)


        # Quantile/Wasserstein categorical poll information
        self._quantileObjectiveDistanceMatrices = None
        self._quantileConstraintDistanceMatrices = None
        self._quantileDistanceDataManager = None


    @property
    def isModelRequired(self):
        return self.catPollType == "SURROGATE" or self.isBOSearchUsed

    # --------------------------------------------------------------- #
    #                   Categorical poll utility functions
    #
    # Private Pareto-based ranking for categorical poll
    # --------------------------------------------------------------- #
    def _dominates(self, p1: float, s1: float, p2: float, s2: float) -> bool:
        """
        Check Pareto dominance (minimization) in (p_u, s_u).
        """
        return (p1 <= p2 and s1 <= s2) and (p1 < p2 or s1 < s2)


    def _pareto_front_indices(self, p: np.ndarray, s: np.ndarray, idxs: list) -> list:
        """
        Return indices (subset of idxs) that are non-dominated.
        """
        front = []
        for i in idxs:
            dominated = False
            for j in idxs:
                if j == i:
                    continue
                if self._dominates(p[j], s[j], p[i], s[i]):
                    dominated = True
                    break
            if not dominated:
                front.append(i)
        return front

    def _rank_Pareto(self, p_u: np.ndarray, s_u: np.ndarray, inc_idx: int, eps_zero: float = 1e-12) -> np.ndarray:
        """
        Rank candidates following the desired strategy:

        1) First Pareto front (non-dominated in (p_u, s_u))
        2) Among remaining, points with p_u = 0 ranked by s_u
        3) Remaining points ranked by p_u (then s_u)

        Returns
        -------
        np.ndarray
            Indices sorted from best to worst.
        """
        p = np.asarray(p_u, dtype=float).copy()
        s = np.asarray(s_u, dtype=float).copy()

        if p.shape != s.shape:
            raise ValueError("p_u and s_u must have same shape.")
        if not (0 <= inc_idx < p.size):
            raise ValueError(f"inc_idx out of range: inc_idx={inc_idx}, n={p.size}")

        # Ignore incumbent
        p[inc_idx] = np.inf
        s[inc_idx] = np.inf

        n = p.size
        all_idxs = [
            i for i in range(n)
            if i != inc_idx and np.isfinite(p[i]) and np.isfinite(s[i])]

        # -------------------------
        # Phase 1: Pareto front
        # -------------------------
        remaining = all_idxs.copy()
        pareto_rank = []

        if remaining:
            first_front = self._pareto_front_indices(p, s, remaining)
            first_front = sorted(first_front, key=lambda i: (p[i], s[i]))
            pareto_rank.extend(first_front)
            first_set = set(first_front)
            remaining = [i for i in remaining if i not in first_set]

        # -------------------------
        # Phase 2: p_u == 0 ranked by s_u
        # -------------------------
        zero_constraint = [i for i in remaining if p[i] <= eps_zero]
        zero_constraint = sorted(zero_constraint, key=lambda i: s[i])

        zero_set = set(zero_constraint)
        remaining = [i for i in remaining if i not in zero_set]

        # -------------------------
        # Phase 3: rest by (p_u, s_u)
        # -------------------------
        rest = sorted(remaining, key=lambda i: (p[i], s[i]))

        final_ranking = pareto_rank + zero_constraint + rest
        return np.asarray(final_ranking, dtype=int)


    # --------------------------------------------------------------- #
    # Surrogate-based categorical poll ranking
    # --------------------------------------------------------------- #
    def _surrogateCategoricalPollRanking(self, inc: np.ndarray, wrappedCombinations: list, is_same_cat: np.ndarray, inc_idx: int,
                                            xCat_formatted: list, xBin_for_pt: np.ndarray, xQuant_for_pt: np.ndarray) -> np.ndarray:
        """
        Rank categorical assignments using the active GP surrogate models.

        Parameters
        ----------
        inc : np.ndarray
            Current feasible or infeasible incumbent.

        wrappedCombinations : list
            All possible categorical assignments.

        is_same_cat : np.ndarray
            Boolean mask identifying the incumbent categorical assignment.

        inc_idx : int
            Index of the incumbent categorical assignment.

        xCat_formatted : list
            Categorical part of the incumbent, formatted as integers.

        xBin_for_pt : np.ndarray
            Binary part of the incumbent, kept fixed during the categorical poll.

        xQuant_for_pt : np.ndarray
            Quantitative part of the incumbent, kept fixed during the categorical poll.

        Returns
        -------
        np.ndarray
            Indices of the selected categorical assignments, ordered from best
            to worst.
        """

        if self.dataManager.gp_models is None:
            raise RuntimeError("GP models missing: constructModel(...) was not successful.")

        sm_f = self.dataManager.gp_models.get("objective", None)
        if sm_f is None:
            raise RuntimeError("GP objective model missing in self.dataManager.gp_models['objective'].")

        if self.dataManager._gp_kernel_distance is None:
            raise RuntimeError("GP kernel distance not initialized. Call constructModel() for GP first.")

        nb_combinations = len(wrappedCombinations)

        objective_distances = np.empty(nb_combinations, dtype=float)

        # Build the incumbent mixed point. Binary and quantitative variables
        # remain fixed throughout the categorical poll.
        x_mix = np.zeros(self.nbVars, dtype=float)

        if self.nbBin > 0:
            x_mix[self.bin_idx] = xBin_for_pt

        if self.nbQuant > 0:
            x_mix[self.quant_idx] = xQuant_for_pt

        if self.cat_idx:
            x_mix[self.cat_idx] = np.asarray(xCat_formatted, dtype=float)

        # Compute the GP-kernel distance between the incumbent categorical
        # assignment and every possible categorical assignment.
        for i, ycat in enumerate(wrappedCombinations):
            if is_same_cat[i]:
                objective_distances[i] = np.inf
                continue

            y_mix = x_mix.copy()
            y_mix[self.cat_idx] = np.asarray(ycat, dtype=float)

            objective_distances[i] = self.dataManager._gp_kernel_distance(x_mix, y_mix)

        # Unconstrained problem: rank directly using the objective-model
        # kernel distance.
        if len(self.dataManager.PB_cols) == 0:
            closestIndices = np.argsort(objective_distances)[:self.nbOfNeighbors]
            return np.asarray(closestIndices, dtype=int)

        # Constrained problem.
        parameter = 1.0
        sm_pbs = self.dataManager.gp_models["PB"]

        constraints_values = np.zeros((nb_combinations, len(sm_pbs)), dtype=float)

        x_mix_1xD = np.atleast_2d(x_mix)

        for i, ycat in enumerate(wrappedCombinations):
            y_mix = x_mix_1xD.copy()
            y_mix[0, self.cat_idx] = np.asarray(ycat, dtype=float)

            for j, surrogate in enumerate(sm_pbs):
                mu = float(surrogate.predict_values(y_mix).reshape(-1)[0])
                var = float(surrogate.predict_variances(y_mix).reshape(-1)[0])

                sigma = float(np.sqrt(max(var, 0.0)))

                if mu - parameter * sigma <= 0.0:
                    constraints_values[i, j] = 0.0
                else:
                    constraints_values[i, j] = mu

        # Euclidean distance between the predicted aggregated constraint
        # representations of each category and that of the incumbent.
        diff = constraints_values - constraints_values[inc_idx, :][None, :]
        constraints_distances = np.linalg.norm(diff, axis=1)

        constraints_distances[inc_idx] = np.inf

        inc_is_feasible = self.dataManager.isPointFeasible(inc)

        if inc_is_feasible:
            ranked_indices = self._rank_Pareto(p_u=constraints_distances, s_u=objective_distances, inc_idx=inc_idx)
        else:
            ranked_indices = self._rank_Pareto(p_u=objective_distances, s_u=constraints_distances, inc_idx=inc_idx)

        closestIndices = ranked_indices[:self.nbOfNeighbors]

        return np.asarray(closestIndices, dtype=int)


    # --------------------------------------------------------------- #
    # Pairwise Wasserstein distance matrix
    # --------------------------------------------------------------- #
    @staticmethod
    def _pairwiseWassersteinDistanceMatrix(samples_by_category: list) -> np.ndarray:
        """
        Compute the pairwise 1-Wasserstein distance matrix between
        empirical distributions associated with categorical values.
        """
        nb_categories = len(samples_by_category)

        distances = np.zeros((nb_categories, nb_categories), dtype=float)

        for i in range(nb_categories):
            for j in range(i + 1, nb_categories):
                distance = wasserstein_distance(samples_by_category[i], samples_by_category[j])
                distances[i, j] = distance
                distances[j, i] = distance

        return distances


    # --------------------------------------------------------------- #
    # Fixed quantile categorical distances from the initial DoE
    # --------------------------------------------------------------- #
    def _constructQuantileCategoricalDistances(self) -> None:
        """
        Construct the Wasserstein distance matrices used by the QUANTILE
        categorical poll.

        Distances are learned only from the initial DoE.
        """

        if not self.dataManager.isDoEFull:
            raise RuntimeError("The initial DoE must be fully evaluated before constructing "
                               "QUANTILE categorical distances.")

        initial_doe = self.dataManager.initial_doe

        if initial_doe.empty:
            raise RuntimeError("Cannot construct QUANTILE categorical distances from an empty initial DoE.")

        # The current categorical poll is defined for a single objective.
        if len(self.dataManager.OBJ_cols) != 1:
            raise NotImplementedError("QUANTILE categorical poll currently requires exactly one objective.")

        objective_col = self.dataManager.OBJ_cols[0]
        pb_cols = self.dataManager.PB_cols

        # One matrix for each categorical variable:
        #
        # self._quantileObjectiveDistanceMatrices[j][a, b]
        #
        # gives the Wasserstein objective distance between categories
        # a and b of categorical variable j.
        self._quantileObjectiveDistanceMatrices = []

        # One matrix for each (constraint, categorical variable):
        #
        # self._quantileConstraintDistanceMatrices[k][j][a, b]
        #
        # gives the Wasserstein violation distance for PB constraint k
        # between categories a and b of categorical variable j.
        self._quantileConstraintDistanceMatrices = [[] for _ in pb_cols]

        # Treat each categorical variable marginally.
        for j, variable_idx in enumerate(self.cat_idx):

            nb_categories = self.nbCategoriesPerVariable[j]
            variable_col = f"x{variable_idx}"

            objective_samples_by_category = []

            constraint_samples_by_category = [[] for _ in pb_cols]

            # ----------------------------------------------------------- #
            # Construct empirical distributions for every category
            # ----------------------------------------------------------- #
            for category in range(nb_categories):

                category_mask = (initial_doe[variable_col].to_numpy(dtype=float) == float(category))

                # ---------------- Objective ---------------- #
                objective_samples = initial_doe.loc[category_mask, objective_col].to_numpy(dtype=float)

                # Ignore invalid blackbox outputs
                objective_samples = objective_samples[np.isfinite(objective_samples)]

                if objective_samples.size == 0:
                    raise RuntimeError(
                        f"No valid objective observation in the initial DoE "
                        f"for categorical variable {j}, category {category}."
                    )

                objective_samples_by_category.append(objective_samples)

                # ---------------- Constraints ---------------- #
                for k, pb_col in enumerate(pb_cols):

                    constraint_samples = initial_doe.loc[category_mask, pb_col].to_numpy(dtype=float)

                    # Ignore invalid blackbox outputs
                    constraint_samples = constraint_samples[np.isfinite(constraint_samples)]

                    if constraint_samples.size == 0:
                        raise RuntimeError(
                            f"No valid observation for constraint '{pb_col}' "
                            f"in the initial DoE for categorical variable {j}, "
                            f"category {category}."
                        )

                    # Constraint representation:
                    # feasible values all become zero.
                    constraint_violations = np.maximum(0.0, constraint_samples)

                    constraint_samples_by_category[k].append(constraint_violations)

            # ----------------------------------------------------------- #
            # Pairwise Wasserstein matrices
            # ----------------------------------------------------------- #

            # Objective
            objective_distance_matrix = (self._pairwiseWassersteinDistanceMatrix(objective_samples_by_category))

            self._quantileObjectiveDistanceMatrices.append(objective_distance_matrix)

            # Constraints
            for k in range(len(pb_cols)):

                constraint_distance_matrix = (self._pairwiseWassersteinDistanceMatrix(constraint_samples_by_category[k]))

                self._quantileConstraintDistanceMatrices[k].append(constraint_distance_matrix)

        # Remember which DataManager produced these distances.
        self._quantileDistanceDataManager = self.dataManager


    # --------------------------------------------------------------- #
    # Quantile-based categorical poll ranking
    # --------------------------------------------------------------- #
    def _quantileCategoricalPollRanking(self, inc: np.ndarray, wrappedCombinations: list, is_same_cat: np.ndarray,
                                              inc_idx: int, xCat_formatted: list) -> np.ndarray:
        """
        Rank categorical assignments using Wasserstein distances learned
        from the initial DoE.
        """

        nb_combinations = len(wrappedCombinations)
        nb_cat_variables = len(self.cat_idx)

        objective_distances = np.empty(nb_combinations, dtype=float)

        # ----------------------------------------------------------- #
        # Objective distances
        # ----------------------------------------------------------- #
        for i, ycat in enumerate(wrappedCombinations):

            if is_same_cat[i]:
                objective_distances[i] = np.inf
                continue

            distances_per_variable = np.empty(nb_cat_variables, dtype=float)

            for j in range(nb_cat_variables):

                incumbent_category = xCat_formatted[j]
                candidate_category = ycat[j]

                distances_per_variable[j] = (
                    self._quantileObjectiveDistanceMatrices[j][
                        incumbent_category,
                        candidate_category,
                    ]
                )

            # Euclidean aggregation across categorical variables
            objective_distances[i] = np.linalg.norm(distances_per_variable)

        # ----------------------------------------------------------- #
        # Unconstrained problem
        # ----------------------------------------------------------- #
        if len(self.dataManager.PB_cols) == 0:

            closestIndices = np.argsort(objective_distances)[:self.nbOfNeighbors]

            return np.asarray(closestIndices, dtype=int)

        # ----------------------------------------------------------- #
        # Constrained problem
        # ----------------------------------------------------------- #

        nb_constraints = len(self.dataManager.PB_cols)

        constraints_distances = np.empty(nb_combinations, dtype=float)

        for i, ycat in enumerate(wrappedCombinations):

            if is_same_cat[i]:
                constraints_distances[i] = np.inf
                continue

            distances_per_constraint = np.empty(nb_constraints, dtype=float)

            for k in range(nb_constraints):

                distances_per_variable = np.empty(nb_cat_variables, dtype=float)

                for j in range(nb_cat_variables):

                    incumbent_category = xCat_formatted[j]
                    candidate_category = ycat[j]

                    distances_per_variable[j] = (
                        self._quantileConstraintDistanceMatrices[k][j][
                            incumbent_category,
                            candidate_category,
                        ]
                    )

                # Euclidean aggregation across categorical variables
                # for constraint k.
                distances_per_constraint[k] = np.linalg.norm(distances_per_variable)

            # Euclidean aggregation across PB constraints.
            constraints_distances[i] = np.linalg.norm(distances_per_constraint)

        # ----------------------------------------------------------- #
        # Same Pareto ranking strategy as the surrogate poll
        # ----------------------------------------------------------- #

        inc_is_feasible = self.dataManager.isPointFeasible(inc)

        if inc_is_feasible:

            ranked_indices = self._rank_Pareto(p_u=constraints_distances, s_u=objective_distances, inc_idx=inc_idx)

        else:

            ranked_indices = self._rank_Pareto(p_u=objective_distances, s_u=constraints_distances, inc_idx=inc_idx)

        closestIndices = ranked_indices[:self.nbOfNeighbors]

        return np.asarray(closestIndices, dtype=int)


    # --------------------------------------------------------------- #
    # Categorical poll
    # --------------------------------------------------------------- #
    def categoricalPoll(self, incumbent: np.ndarray) -> np.ndarray:
        """
        Perform categorical polling using the selected categorical poll strategy.

        Returns
        -------
        np.ndarray
            Directions of shape (self.nbOfNeighbors, nbVars).
        """

        # Incumbent provided by NOMAD: it is either the feasible or the infeasible incumbent
        inc = np.asarray(incumbent, dtype=float).ravel()
        if inc.size != self.nbVars:
            raise ValueError(f"incumbent must have length {self.nbVars}, got {inc.size}.")
        
        # Edge-case where categorical poll has no neighbors
        if self.nbOfNeighbors == 0:
            return np.empty((0, self.nbVars), dtype=float)

        # Split incumbent once
        xCat, xBin, xQuant = self.dataManager.splitPoint(inc)
        xCat_formatted = [int(v) for v in xCat.tolist()]

        # Enumerate categorical assignments
        allCombinations = list(itertools.product(*[range(lj) for lj in self.nbCategoriesPerVariable]))
        wrappedCombinations = [list(c) for c in allCombinations]
        if len(wrappedCombinations) <= 1:
            return np.empty((0, self.nbVars), dtype=float)

        # Mask incumbent category
        is_same_cat = np.all(
            np.asarray(wrappedCombinations, dtype=int) == np.asarray(xCat_formatted, dtype=int),
            axis=1)

        # Find the incumbent categorical assignment
        where = np.where(is_same_cat)[0]
        if where.size != 1:
            raise RuntimeError(f"Expected exactly one incumbent categorical match, got {where.size}.")
        inc_idx = int(where[0])

        # Stay fixed for the categorical poll
        xBin_for_pt = np.asarray(xBin, dtype=float) if self.nbBin > 0 else np.zeros((0,), dtype=float)
        xQuant_for_pt = np.asarray(xQuant, dtype=float) if self.nbQuant > 0 else np.zeros((0,), dtype=float)


        # ------------------------- Select categorical vectors --------------------- #
        if self.catPollType == "SURROGATE":
            if not self.dataManager.modelIsReady:
                raise RuntimeError("No active surrogate model. Call constructModel(...) first.")

            closestIndices = self._surrogateCategoricalPollRanking(
                inc=inc,
                wrappedCombinations=wrappedCombinations,
                is_same_cat=is_same_cat,
                inc_idx=inc_idx,
                xCat_formatted=xCat_formatted,
                xBin_for_pt=xBin_for_pt,
                xQuant_for_pt=xQuant_for_pt)

        elif self.catPollType == "QUANTILE":

            # Construct global Wasserstein categorical distances once from the initial DoE
            if self._quantileDistanceDataManager is not self.dataManager:
                self._constructQuantileCategoricalDistances()

            closestIndices = self._quantileCategoricalPollRanking(
                inc=inc,
                wrappedCombinations=wrappedCombinations,
                is_same_cat=is_same_cat,
                inc_idx=inc_idx,
                xCat_formatted=xCat_formatted,
            )

        elif self.catPollType == "EXHAUSTIVE":
                closestIndices = np.asarray([i for i in range(len(wrappedCombinations)) if i != inc_idx], dtype=int)

        else:
            raise RuntimeError(f"Unknown catPollType='{self.catPollType}'.")
        # ---------------------------------------------------------------------------------- #


        # ------------------------- Construct directions-------------- --------------------- #
        closestCatVectors = [wrappedCombinations[i] for i in closestIndices]

        catDirections = []
        for vec in closestCatVectors:
            pt = np.zeros(self.nbVars, dtype=float)

            if self.cat_idx:
                pt[self.cat_idx] = np.asarray(vec, dtype=float)
            if self.nbBin > 0:
                pt[self.bin_idx] = xBin_for_pt
            if self.nbQuant > 0:
                pt[self.quant_idx] = xQuant_for_pt

            direction = pt - inc
            print("Categorical direction:", direction)
            catDirections.append(direction)

        return np.asarray(catDirections, dtype=float)
    # --------------------------------------------- NEW ----------------------------------------------------- #





    # --------------------------------------------------------------- #
    # Old categorical poll
    # --------------------------------------------------------------- #
    def categoricalPollOld(self, incumbent: np.ndarray) -> np.ndarray:
        """
        Perform categorical polling using the active surrogate models.

        Returns
        -------
        np.ndarray
            Directions of shape (self.nbOfNeighbors, nbVars).
        """

        # Incumbent provided by NOMAD: it is either the feasible or the infeasible incumbent
        inc = np.asarray(incumbent, dtype=float).ravel()
        if inc.size != self.nbVars:
            raise ValueError(f"incumbent must have length {self.nbVars}, got {inc.size}.")
        
        # Edge-case where categorical poll has no neighbors
        if self.nbOfNeighbors == 0:
            return np.empty((0, self.nbVars), dtype=float)

        if not self.dataManager.modelIsReady:
            raise RuntimeError("No active surrogate model. Call constructModel(...) first.")

        # Split incumbent once
        xCat, xBin, xQuant = self.dataManager.splitPoint(inc)
        xCat_formatted = [int(v) for v in xCat.tolist()]

        # Enumerate categorical assignments
        allCombinations = list(itertools.product(*[range(lj) for lj in self.nbCategoriesPerVariable]))
        wrappedCombinations = [list(c) for c in allCombinations]
        if len(wrappedCombinations) <= 1:
            return np.empty((0, self.nbVars), dtype=float)

        # Mask incumbent category
        is_same_cat = np.all(
            np.asarray(wrappedCombinations, dtype=int) == np.asarray(xCat_formatted, dtype=int),
            axis=1
        )

        # Mask for GPs
        where = np.where(is_same_cat)[0]
        if where.size != 1:
            raise RuntimeError(f"Expected exactly one incumbent categorical match, got {where.size}.")
        inc_idx = int(where[0])

        # Stay fixed for the categorical poll
        xBin_for_pt = np.asarray(xBin, dtype=float) if self.nbBin > 0 else np.zeros((0,), dtype=float)
        xQuant_for_pt = np.asarray(xQuant, dtype=float) if self.nbQuant > 0 else np.zeros((0,), dtype=float)


        # ======================================================
        # Compute a ranking key
        # ======================================================
        if self.dataManager.gp_models is None:
            raise RuntimeError("GP models missing: constructModel(...) was not successful.")
        sm_f = self.dataManager.gp_models.get("objective", None)
        if sm_f is None:
            raise RuntimeError("GP objective model missing in self.gp_models['objective'].")
        if self.dataManager._gp_kernel_distance is None:
            raise RuntimeError("GP kernel distance not initialized. Call constructModel() for GP first.")

        # Compute kernel-induced distances to all categorical candidates (fixing bin+quant)
        key = np.empty((len(wrappedCombinations),), dtype=float)
        objective_distances = np.empty((len(wrappedCombinations),), dtype=float)
        constraints_distances = np.empty((len(wrappedCombinations),), dtype=float)


        # Build x_mix as a 1D numpy array (length nbVars), with fixed bin+quant
        x_mix = np.zeros((self.nbVars,), dtype=float)

        if self.nbBin > 0:
            x_mix[self.bin_idx] = xBin_for_pt
        if self.nbQuant > 0:
            x_mix[self.quant_idx] = xQuant_for_pt
        if self.cat_idx:
            x_mix[self.cat_idx] = np.asarray(xCat_formatted, dtype=float)


        # Compute objective distance no matter if there are constraints or not
        for i, ycat in enumerate(wrappedCombinations):
            # Distance to itself should be ignored
            if is_same_cat[i]:
                objective_distances[i] = np.inf
            else:
                # y_mix must be a numpy array before .copy()
                y_mix = x_mix.copy()
                y_mix[self.cat_idx] = np.asarray(ycat, dtype=float)
                objective_distances[i] = self.dataManager._gp_kernel_distance(x_mix, y_mix)

        # If there are no constraints keys are directly the objective distances
        if len(self.dataManager.PB_cols) == 0:
            key = objective_distances
        # Otherwise, there are constraints
        else:
            parameter = 1  # value for the pseudo-distance with 2 cases: 1) mu - parameter*sigma<0; 2) otherwise 
            sm_pbs = self.dataManager.gp_models["PB"]
            constraints_values = np.zeros((len(wrappedCombinations), len(sm_pbs)), dtype=float)

            # x_mix should be the incumbent mixed point in original encoding (length nbVars)
            # Make sure it's 2D for SMT when predicting
            x_mix_1xD = np.atleast_2d(np.asarray(x_mix, dtype=float))

            # Compute the constraints distance between the incumbent and all others
            for i, ycat in enumerate(wrappedCombinations):
                # Candidate point = incumbent point, but with categorical part replaced
                # We consider also the incumbent for the constraints_values
                y_mix = x_mix_1xD.copy()  # shape (1, nbVars)
                y_mix[0, self.cat_idx] = np.asarray(ycat, dtype=float)

                for j, surrogate in enumerate(sm_pbs):
                    mu = float(surrogate.predict_values(y_mix).reshape(-1)[0])
                    var = float(surrogate.predict_variances(y_mix).reshape(-1)[0])

                    # Numerical safety: variance should be >= 0, but clamp small negatives if any
                    sigma = float(np.sqrt(max(var, 0.0)))

                    # Rule constraint value is 0 if mu - parameter*sigma <= 0 (predicted to be feasible)
                    if (mu - parameter * sigma) <= 0.0:
                        constraints_values[i, j] = 0.0
                    else:
                        constraints_values[i, j] = mu

            # Euclidean distances to incumbent row
            diff = constraints_values - constraints_values[inc_idx, :][None, :]
            constraints_distances = np.linalg.norm(diff, axis=1)

            # Ignore incumbent itself
            constraints_distances[inc_idx] = np.inf

            # p_u is the primary ranking function, and s_u is the secondary ranking function
            # p_u is prioritized over s_u in case of equality or in case of partial order
            inc_is_feasible = self.dataManager.isPointFeasible(inc)
            #print("point:" , inc, " is feasible ?", inc_is_feasible)
            if inc_is_feasible:
                key = self._rank_Pareto(p_u=constraints_distances, s_u=objective_distances, inc_idx=inc_idx)
            else:
                key = self._rank_Pareto(p_u=objective_distances, s_u=constraints_distances, inc_idx=inc_idx)


        # ======================================================
        # Selection + direction building
        # ======================================================
        if len(self.dataManager.PB_cols) > 0:
            ranked_indices = key  # key is already an ordered list of candidate indices
            closestIndices = ranked_indices[:self.nbOfNeighbors]
        else:
            closestIndices = np.argsort(key)[:self.nbOfNeighbors]

        closestCatVectors = [wrappedCombinations[i] for i in closestIndices]

        catDirections = []
        for vec in closestCatVectors:
            pt = np.zeros(self.nbVars, dtype=float)

            if self.cat_idx:
                pt[self.cat_idx] = np.asarray(vec, dtype=float)
            if self.nbBin > 0:
                pt[self.bin_idx] = xBin_for_pt
            if self.nbQuant > 0:
                pt[self.quant_idx] = xQuant_for_pt

            direction = pt - inc
            catDirections.append(direction)

        return np.asarray(catDirections, dtype=float)
    


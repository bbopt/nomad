"""
Management for CatMads optimization.

This module handles the model management including caching, parameter setting,
and surrogate model operations for CatMads algorithms.
"""

# Needed for `int | None` style annotations to work on Python < 3.10
from __future__ import annotations

# Imported from other local files
from typing import TYPE_CHECKING

# --------------------------------------- #
# Installed packages
import numpy as np
import pandas as pd
# --------------------------------------- #


if TYPE_CHECKING:
    from ..problem_definition import ProblemDefinition

class DataManager:
    """
    Manages data, models and callback related to models.
    
    """
    
    def __init__(self,
                problem_definition: 'ProblemDefinition',
                seed: int = 0,
                isModelRequired: bool = False,  # Whether a surrogate model is currently required
                isModelUpdated: bool = False,   # Whether repeated surrogate-model updates are currently active
                modelUpdateHardCap: int = 500,  # Static, initialized from MixedOptimizer() or QuantitativeOptimizer() in NomadBBO.pyx
                ):
        """
        Initialize the CatMadsDataManager with a problem definition.
        
        Parameters:
        - problem_definition: CatMadsProblemDefinition
            The problem definition containing variables and constraints.
        """
        # ------------------------------------------------ #
        # Basic parameters
        self.seed = int(seed)
        self.modelUpdateHardCap = int(modelUpdateHardCap)
        # ------------------------------------------------ #


        # ------------------------------------------------ #
        # Parameters passed as arguments in the NomadBBO.pyx in optimize()
        self.problem_definition = problem_definition
        if self.problem_definition is None:
            raise ValueError("problem_definition is not set.")
        # ------------------------------------------------ #


        # ------------------------------------------ #
        # Initialize model information
        # States 
        self.modelIsReady = False  # To be set once model is trained successfully
        self.wasModelUpdatedOnce = False
        self.isModelRequired = bool(isModelRequired)
        self.isModelUpdated = bool(isModelUpdated)


        # Common model parameters: fMin and fMax are used for EGO/BO with GPs  
        self.fMin = None 
        self.fMax = None
        # ------------------------------------------ #


        # ------------------------------------------ #
        # Initialize cache information

        # Evaluated points stored in a pandas data structure
        self.cache = pd.DataFrame()
        self.initial_doe = pd.DataFrame()
        self.isDoEFull = False              #

        # We also store the list of the columns to have O(1) operations with the cache in pandas data structure
        self.point_cols = []    # always in format ["x0","x1",..., "xn"]
        self.output_cols = []   # can be any ordering of the outputs, ["OBJ1","PB1","PB2"] or ["PB1","OBJ1","PB2"]
        self.OBJ_cols = []
        self.PB_cols = []
        self.EB_cols = []

        self._init_cache_and_doe()
        # ------------------------------------------ #


    # Private class initializing cache information into the model manager
    def _init_cache_and_doe(self):
        """
        Initialize the cache and initial DOE schemas.

        Both use:
        [x0..x{nbVars-1} | BBOutputType]

        Column order strictly follows BBOutputType.
        """

        # 1) Point columns (fixed NOMAD order)
        self.point_cols = [f"x{j}" for j in range(self.nbVars)]

        # 2) Output columns (nice human-readable naming)
        self.output_cols = []

        obj_count = 0
        pb_count  = 0
        eb_count  = 0

        # For column naming, i.e. distinguish the different OBJ, PB and EB e.g. [OBJ1, OBJ2, PB1, PB2, EB1]
        for t in self.BBOutputType:
            if t == "OBJ":
                obj_count += 1
                col = f"OBJ{obj_count}"
            elif t == "PB":
                pb_count += 1
                col = f"PB{pb_count}"
            elif t == "EB":
                eb_count += 1
                col = f"EB{eb_count}"
            else:
                col = str(t)

            self.output_cols.append(col)

        # 3) Precompute semantic groups (O(1) access later)
        self.OBJ_cols = [c for c in self.output_cols if c.startswith("OBJ")]
        self.PB_cols  = [c for c in self.output_cols if c.startswith("PB")]
        self.EB_cols  = [c for c in self.output_cols if c.startswith("EB")]

        # 4) Final schema
        all_cols = self.point_cols + self.output_cols

        self.cache = pd.DataFrame(columns=all_cols, dtype=float)
        self.initial_doe = pd.DataFrame(columns=all_cols, dtype=float)


    # ------------------------------------------ #
    # Getters from problem definition
    # ------------------------------------------ #
    @property
    def idx_by_type(self):
        return self.problem_definition.getVariableIndicesByType()

    @property
    def cat_idx(self):
        return self.idx_by_type.get("choice", []) or []

    @property
    def bin_idx(self):
        return self.idx_by_type.get("binary", []) or []

    @property
    def quant_idx(self):
        return (self.idx_by_type.get("integer", []) or []) + (self.idx_by_type.get("real", []) or [])

    @property
    def nbBin(self):
        return len(self.bin_idx)

    @property
    def nbQuant(self):
        return len(self.quant_idx)

    @property
    def nbVars(self):
        return self.problem_definition.getDimension()

    @property
    def nbCategoriesPerVariable(self):
        return self.problem_definition.getNbCategoriesPerCategoricalVariable() or []

    @property
    def quant_bounds(self):
        return self.problem_definition.getQuantitativeBounds() or []

    @property
    def lowerBoundsQuant(self):
        if self.nbQuant > 0:
            return np.array([b[0] for b in self.quant_bounds], dtype=float)
        return np.array([], dtype=float)

    @property
    def upperBoundsQuant(self):
        if self.nbQuant > 0:
            return np.array([b[1] for b in self.quant_bounds], dtype=float)
        return np.array([], dtype=float)

    @property
    def BBOutputType(self):
        return self.problem_definition.getBBOutputType()

    # ------------------------------------------ #
    # Getters 
    # ------------------------------------------ #
    def getPoints(self) -> np.ndarray:
        if self.cache.empty:
            return np.empty((0, self.nbVars), dtype=float)
        return self.cache[self.point_cols].to_numpy(dtype=float)

    def getOutputValues(self) -> np.ndarray:
        if self.cache.empty:
            return np.empty((0, len(self.output_cols)), dtype=float)
        return self.cache[self.output_cols].to_numpy(dtype=float)
    
    def getObjValues(self, flatten: bool = True) -> np.ndarray:
        """
        Returns objective values.

        Parameters
        ----------
        flatten : bool, default=True
            If True and there is exactly one objective, returns a 1D array (n,).
            Otherwise returns a 2D array (n, n_obj).

        Returns
        -------
        np.ndarray
            Shape:
            - (n,) if single objective and flatten=True
            - (n, n_obj) otherwise
            - (n, 0) if no objective exists (edge case)
        """
        if self.cache.empty:
            return np.empty((0,), dtype=float) if flatten else np.empty((0, 0), dtype=float)

        if not self.OBJ_cols:
            # Edge case: no objective in BBOutputType
            return np.empty((len(self.cache), 0), dtype=float)

        Y = self.cache[self.OBJ_cols].to_numpy(dtype=float)

        if flatten and Y.shape[1] == 1:
            return Y[:, 0]

        return Y

    def getPBValues(self) -> np.ndarray:
        if self.cache.empty:
            return np.empty((0, len(self.PB_cols)), dtype=float) if self.PB_cols else np.empty((0, 0), dtype=float)
        if not self.PB_cols:
            return np.empty((len(self.cache), 0), dtype=float)
        return self.cache[self.PB_cols].to_numpy(dtype=float)
    
    def getHValues(self) -> np.ndarray:
        if self.cache.empty:
            return np.empty((0,), dtype=float)

        PB = self.getPBValues()
        if PB.size == 0:
            return np.zeros((len(self.cache),), dtype=float)

        aggregation = np.maximum(0.0, PB)
        return np.sum(aggregation * aggregation, axis=1)

    def getModelIsReady(self) -> bool:
        """
        Check if the surrogate model is ready.
        
        Returns:
        - bool: True if the model is ready, False otherwise.
        """
        return self.modelIsReady

    def getProblemDefinition(self) -> 'ProblemDefinition':
        """
        Get the problem definition associated with this manager.
        
        Returns:
        - ProblemDefinition: The problem definition used to initialize this manager.
        """
        return self.problem_definition

    def getSeed(self) -> int:
        """
        Get the random seed used by the manager.
        
        Returns:
        - int: Seed controlling stochastic components (fold shuffling, optimizers, etc.).
        """
        return self.seed


    # ------------------------------------------ #
    # Setter
    # ------------------------------------------ #
    def setInitialDoEPoints(self, x0s) -> None:
        """
        Set the initial DoE points.

        The blackbox outputs are left as NaN until the points are evaluated
        and added through addToCache().
        """
        X = np.asarray(x0s, dtype=float)

        if X.size == 0:
            self.isDoEFull = True
            return

        if X.ndim == 1:
            X = X.reshape(1, -1)

        if X.ndim != 2 or X.shape[1] != self.nbVars:
            raise ValueError(f"Initial DoE must have shape (n, {self.nbVars}), got {X.shape}.")

        self.initial_doe = pd.DataFrame(X, columns=self.point_cols, dtype=float,)

        for col in self.output_cols:
            self.initial_doe[col] = np.nan

        self.isDoEFull = False

    # ------------------------------------------ #
    # Utility methods
    # ------------------------------------------ #
    def splitPoint(self, point: np.ndarray) -> tuple:
        pointArr = np.asarray(point, dtype=float).ravel()

        if pointArr.size != self.nbVars:
            raise ValueError(
                f"splitPoint expected length {self.nbVars}, got {pointArr.size}."
            )

        xCat = pointArr[self.cat_idx]
        xBin = pointArr[self.bin_idx]
        xQuant = pointArr[self.quant_idx]

        return xCat, xBin, xQuant

    def _find_cache_index_of_point(self, point: np.ndarray, atol: float = 0.0) -> int:
        """
        Return the row index in self.cache matching `point` on the decision variables.

        Parameters
        ----------
        point : np.ndarray
            Shape (nbVars,).
        atol : float, default=0.0
            Absolute tolerance for matching (useful if any columns may be float-ish).

        Returns
        -------
        int
            Row index in the cache.

        Raises
        ------
        ValueError
            If point not found or found multiple times.
        """
        if self.cache.empty:
            raise ValueError("Cache is empty; point cannot be found.")

        x = np.asarray(point, dtype=float).ravel()
        if x.size != self.nbVars:
            raise ValueError(f"point must have length {self.nbVars}, got {x.size}.")

        X = self.cache[self.point_cols].to_numpy(dtype=float)

        if atol <= 0.0:
            matches = np.all(X == x[None, :], axis=1)
        else:
            matches = np.all(np.isclose(X, x[None, :], atol=atol, rtol=0.0), axis=1)

        idxs = np.where(matches)[0]
        if idxs.size == 0:
            raise ValueError(f"Point not found in cache: {x.tolist()}")
        if idxs.size > 1:
            raise ValueError(f"Point found multiple times in cache (rows={idxs.tolist()}): {x.tolist()}")

        return int(idxs[0])


    def isPointFeasible(self, point: np.ndarray, atol_lookup: float = 0.0, eps_h: float = 1e-10) -> bool:
        """
        Check feasibility of a cached point using H(point) == 0.

        Feasible means getHValues()[row_idx] is (numerically) zero.

        Parameters
        ----------
        point : np.ndarray
            Shape (nbVars,).
        atol_lookup : float, default=0.0
            Tolerance when searching the point in cache.
        eps_h : float, default=1e-12
            Numerical tolerance for H == 0.

        Returns
        -------
        bool
            True if feasible, False otherwise.
        """
        # No PB constraints => always feasible
        if len(self.PB_cols) == 0:
            return True

        row_idx = self._find_cache_index_of_point(point, atol=atol_lookup)
        H = self.getHValues()

        if H.size == 0:
            # Edge case; but if PB_cols exist, H should not be empty
            return True

        return float(H[row_idx]) <= float(eps_h)


    # Utility function for constructing models with only valid points
    def _get_valid_training_data(self, large_value: float = 1e23):
        """
        Return training data filtered for GP construction.

        Removes rows containing:
        - NaN
        - +/-inf
        - non-real values
        - |value| >= large_value in outputs
        """

        X_train = self.getPoints()
        Y_obj = self.getObjValues(flatten=True)
        Y_pb = self.getPBValues()

        if X_train.shape[0] == 0:
            return X_train, Y_obj, Y_pb

        outputs = [Y_obj.reshape(-1, 1)]
        if Y_pb.size > 0:
            outputs.append(Y_pb)

        Y_all = np.hstack(outputs)

        # Check X and outputs
        X_valid = (
            np.all(np.isreal(X_train), axis=1)
            & np.all(np.isfinite(X_train), axis=1)
        )

        Y_valid = (
            np.all(np.isreal(Y_all), axis=1)
            & np.all(np.isfinite(Y_all), axis=1)
            & np.all(np.abs(Y_all) < large_value, axis=1)
        )

        valid_mask = X_valid & Y_valid

        return (
            X_train[valid_mask],
            Y_obj[valid_mask],
            Y_pb[valid_mask, :] if Y_pb.size > 0 else Y_pb,
        )


    def _find_initial_doe_index_of_point(self, point: np.ndarray) -> int | None:
        """
        Return the row index in self.initial_doe matching `point`.

        Returns None if the point is not part of the initial DoE.
        """
        x = np.asarray(point, dtype=float).ravel()

        if x.size != self.nbVars:
            raise ValueError(
                f"point must have length {self.nbVars}, got {x.size}."
            )

        X_doe = self.initial_doe[self.point_cols].to_numpy(dtype=float)

        matches = np.all(X_doe == x[None, :], axis=1)
        idxs = np.where(matches)[0]

        if idxs.size == 0:
            return None

        if idxs.size > 1:
            raise ValueError(
                f"Point found multiple times in initial DoE: {x.tolist()}"
            )

        return int(idxs[0])


    def _is_initial_doe_full(self) -> bool:
        if self.initial_doe.empty:
            return True

        X_doe = self.initial_doe[self.point_cols].to_numpy(dtype=float)
        X_cache = self.cache[self.point_cols].to_numpy(dtype=float)

        doe_points = {tuple(x) for x in X_doe}
        cache_points = {tuple(x) for x in X_cache}

        return doe_points.issubset(cache_points)


    # ---------------------------------------------------------- #
    # Callbacks passing data info between Nomad and Python side
    # ---------------------------------------------------------- #
    def addToCache(self, evaluatedPoint: np.ndarray) -> bool:
        arr = np.asarray(evaluatedPoint, dtype=float).ravel()

        expected = self.nbVars + len(self.output_cols)
        if arr.size != expected:
            raise ValueError(f"Cache row must have length {expected}, got {arr.size}.")

        doe_idx = None
        already_evaluated = False

        # Check initial DoE information before adding the point to the cache
        if not self.isDoEFull:
            x = arr[:self.nbVars]
            doe_idx = self._find_initial_doe_index_of_point(x)

            if doe_idx is not None:
                X_cache = self.cache[self.point_cols].to_numpy(dtype=float)

                already_evaluated = np.any(np.all(X_cache == x[None, :], axis=1))

        # Add evaluated point to complete cache
        self.cache.loc[len(self.cache), self.point_cols + self.output_cols] = arr

        # Complete initial DoE information
        if doe_idx is not None:
            # Preserve the first evaluation of each initial DoE point
            if not already_evaluated:
                self.initial_doe.loc[doe_idx, self.output_cols] = arr[self.nbVars:]

            # Check whether every initial DoE point has been evaluated
            self.isDoEFull = self._is_initial_doe_full()

        return True


    def completeTrialPointsInformation(self, trialPoints: np.ndarray) -> np.ndarray:
        """
        Complete trial points information using the active surrogate model.

        Uses:
        - objective and PB constraints surrogates (self.gp_models):  EB is not modeled/considered

        Returns
        -------
        np.ndarray
            Shape (n_points, len(self.output_cols)), ordered as self.output_cols.
        """

        tpMat = np.asarray(trialPoints, dtype=float)
        if tpMat.ndim != 2:
            raise ValueError(f"trialPoints must be a 2D array, got ndim={tpMat.ndim}.")
        if tpMat.shape[1] != self.nbVars:
            raise ValueError(f"trialPoints must have {self.nbVars} columns, got {tpMat.shape[1]}.")

        if not self.modelIsReady:
            raise RuntimeError(
                "No active surrogate model. Call constructModel(...) first."
            )

        # Enforce single objective (current design)
        if len(self.OBJ_cols) != 1:
            raise RuntimeError(
                f"CatMadsDataManager enforces exactly 1 objective, got {len(self.OBJ_cols)}: {self.OBJ_cols}"
            )

        obj_col = self.OBJ_cols[0]
        obj_pos = self.output_cols.index(obj_col)

        n_pts = tpMat.shape[0]
        results = np.zeros((n_pts, len(self.output_cols)), dtype=np.float64)

        if self.gp_models is None:
            raise RuntimeError("GP models missing: constructModel(...) was not successful.")

        sm_f = self.gp_models.get("objective", None)
        if sm_f is None:
            raise RuntimeError("GP objective model missing in self.gp_models['objective'].")

        sm_pbs = self.gp_models.get("PB", [])
        if len(sm_pbs) != len(self.PB_cols):
            raise RuntimeError(
                f"GP PB models mismatch: have {len(sm_pbs)} models but PB_cols={len(self.PB_cols)}."
            )

        # SMT expects discrete variables as integers
        Xsmt = np.array(tpMat, copy=True)

        if self.cat_idx:
            Xsmt[:, self.cat_idx] = np.rint(Xsmt[:, self.cat_idx]).astype(int)
        if self.bin_idx:
            Xsmt[:, self.bin_idx] = np.rint(Xsmt[:, self.bin_idx]).astype(int)

        int_idx = (self.idx_by_type.get("integer", []) or [])
        if int_idx:
            Xsmt[:, int_idx] = np.rint(Xsmt[:, int_idx]).astype(int)

        # Objective prediction
        y_obj = sm_f.predict_values(Xsmt).reshape(-1)
        results[:, obj_pos] = y_obj

        # PB prediction(s)
        for k, pb_col in enumerate(self.PB_cols):
            pb_pos = self.output_cols.index(pb_col)
            y_pb = sm_pbs[k].predict_values(Xsmt).reshape(-1)
            results[:, pb_pos] = y_pb

        return np.ascontiguousarray(results, dtype=np.float64)


   
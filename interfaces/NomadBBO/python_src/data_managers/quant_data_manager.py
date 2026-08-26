import numpy as np
import io
import re
import copy
import contextlib

from .data_manager import DataManager
from typing import TYPE_CHECKING

from smt.applications.mixed_integer import MixedIntegerKrigingModel
from smt.surrogate_models import KRG
from smt.design_space import DesignSpace, FloatVariable, IntegerVariable

if TYPE_CHECKING:
    from ..problem_definition import ProblemDefinition


class QuantitativeDataManager(DataManager):
    """
    Child class of DataManager for quantitative problems.

    Handles real, integer, and binary variables, but no categorical variables.
    """

    def __init__(
        self,
        problem_definition: 'ProblemDefinition',
        seed: int = 0,
        isModelRequired: bool = False,
        isModelUpdated: bool = False,
        modelUpdateHardCap: int = 500,
    ):
        super().__init__(
            problem_definition=problem_definition,
            seed=seed,
            isModelRequired=isModelRequired,
            isModelUpdated=isModelUpdated,
            modelUpdateHardCap=modelUpdateHardCap
        )

        # ------------------------------------------ #
        # Model parameters for GP/SMT (to reuse after training GPs)
        self.gp_models = None
        self.gp_design_space = None
        # ------------------------------------------ #


    def constructModel(self, **kwargs) -> bool:
        """
        Construct the active surrogate model from the cache.
        """

        # No component of the optimizer currently requires a surrogate model.
        if not self.isModelRequired:
            return False


        # ------------ Model update stopping rules -------- #
        cache_size = len(self.cache)
        isHardCapReached = cache_size >= self.modelUpdateHardCap

        # Keep using the last successfully trained model once repeated
        # model updates have been disabled.
        if self.wasModelUpdatedOnce and not self.isModelUpdated:
            return self.modelIsReady

        # Stop repeated model updates once the hard cap is reached.
        # Since quantitative surrogate models are only used by BO,
        # BOSearch() will consequently disable BO.
        if self.wasModelUpdatedOnce and isHardCapReached:
            self.isModelUpdated = False
            return self.modelIsReady
        # ------------------------------------------------- #

        self.modelIsReady = False

        if len(self.OBJ_cols) != 1:
            raise ValueError(
                f"DataManager enforces exactly 1 objective, got {len(self.OBJ_cols)}: {self.OBJ_cols}"
            )

        X_train, Y_obj, Y_pb = self._get_valid_training_data(large_value=1e23)

        if X_train.size == 0 or X_train.shape[0] < 2:
            return False  # Not enough data to train a model

        n, d = X_train.shape

        if d != self.nbVars:
            raise RuntimeError(
                f"Training X has {d} variables but expected nbVars={self.nbVars}."
            )

        if Y_obj.shape[0] != n:
            raise RuntimeError(
                f"Y_obj has {Y_obj.shape[0]} rows but X_train has {n}."
            )

        if Y_pb.shape[0] != n:
            raise RuntimeError(
                f"Y_pb has {Y_pb.shape[0]} rows but X_train has {n}."
            )

        self.fMin = float(np.min(Y_obj))
        self.fMax = float(np.max(Y_obj))

        print("Constructing quantitative GP surrogate models (SMT) from cache...")

        n_start_mult = int(kwargs.get("n_start_mult", 15))

        self.gp_models = None
        self.gp_design_space = None

        def _train_with_warnings(model, X_train, y_train) -> bool:
            error_pattern = r"exception :\s+\d+-th leading minor of the array is not positive definite"
            captured_output = io.StringIO()

            model.set_training_values(X_train, y_train)

            with contextlib.redirect_stderr(captured_output), contextlib.redirect_stdout(captured_output):
                try:
                    model.train()
                except Exception as e:
                    print(f"An exception occurred during training: {e}")

            return re.search(error_pattern, captured_output.getvalue()) is not None

        # SMT expects binary and integer columns to be integers.
        Xsmt = np.array(X_train, copy=True)

        if self.cat_idx:
            raise ValueError(
                "QuantitativeDataManager does not support categorical variables. "
                "Use MixedDataManager instead."
            )

        if self.bin_idx:
            Xsmt[:, self.bin_idx] = np.rint(Xsmt[:, self.bin_idx]).astype(int)

        int_idx = self.idx_by_type.get("integer", []) or []
        real_idx = self.idx_by_type.get("real", []) or []

        if int_idx:
            Xsmt[:, int_idx] = np.rint(Xsmt[:, int_idx]).astype(int)

        quant_bound_map = {j: self.quant_bounds[k] for k, j in enumerate(self.quant_idx)}

        var_ds = []
        for j in range(self.nbVars):
            if j in self.bin_idx:
                var_ds.append(IntegerVariable(0, 1))

            elif j in int_idx:
                if j not in quant_bound_map:
                    raise ValueError(f"Missing bounds for integer variable x{j}.")
                lb, ub = quant_bound_map[j]
                var_ds.append(IntegerVariable(int(lb), int(ub)))

            elif j in real_idx:
                if j not in quant_bound_map:
                    raise ValueError(f"Missing bounds for real variable x{j}.")
                lb, ub = quant_bound_map[j]
                var_ds.append(FloatVariable(float(lb), float(ub)))

            else:
                raise ValueError(
                    f"Variable x{j} not classified in binary/integer/real indices."
                )

        ds = DesignSpace(var_ds)
        self.gp_design_space = ds

        def make_gp():
            return MixedIntegerKrigingModel(
                surrogate=KRG(
                    design_space=copy.deepcopy(ds),
                    hyper_opt="Cobyla",
                    corr="squar_exp",
                    n_start=max(1, int(n_start_mult * self.nbVars)),
                    print_prediction=False,
                ),
            )

        sm_f = make_gp()
        warning_triggered = _train_with_warnings(sm_f, copy.deepcopy(Xsmt), Y_obj)

        if warning_triggered:
            print("\n>Warning triggered during training of objective. GP construction aborted.\n")
            self.modelIsReady = False
            return False

        sm_pbs = []
        for i in range(Y_pb.shape[1]):
            sm_c = make_gp()
            warning_triggered = _train_with_warnings(sm_c, copy.deepcopy(Xsmt), Y_pb[:, i])

            if warning_triggered:
                print(f"\n>Warning triggered during training of PB constraint {i}. GP construction aborted.\n")
                self.modelIsReady = False
                return False

            sm_pbs.append(sm_c)

        self.gp_models = {"objective": sm_f, "PB": sm_pbs}
        self.modelIsReady = True
        self.wasModelUpdatedOnce = True

        print("GP objective and PB-constraint surrogates trained successfully.")
        return self.modelIsReady
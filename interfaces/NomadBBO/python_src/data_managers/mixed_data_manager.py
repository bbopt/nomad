
import numpy as np
import io
import re
import copy
import contextlib

from .data_manager import DataManager
from typing import TYPE_CHECKING

from smt.applications.mixed_integer import MixedIntegerKrigingModel
from smt.surrogate_models import KRG, MixIntKernelType
from smt.design_space import (CategoricalVariable, DesignSpace, FloatVariable, IntegerVariable)
from smt.utils.kriging import (MixHrcKernelType, cross_distances, gower_componentwise_distances, cross_levels)

if TYPE_CHECKING:
    from ..problem_definition import ProblemDefinition    

class MixedDataManager(DataManager):
    """
    Child class of DataManager specifically for mixed-variable problems with categorical variables. 
    
    """

    def __init__(
        self,
        # DataManager parent class attributes
        problem_definition: 'ProblemDefinition',
        seed=0,
        isModelUpdated=False,
        modelUpdateHardCap=500,
        isModelRequired=False,
        # Specific to MixedDataManager
        catKernelType=None):
        #
        super().__init__(
            problem_definition=problem_definition,
            seed=seed,
            isModelRequired=isModelRequired,
            isModelUpdated=isModelUpdated,
            modelUpdateHardCap=modelUpdateHardCap,
        )


        # ------------------------------------------------ #
        # Categorical kernel initialization
        if self.isModelRequired:
            if catKernelType is None:
                catKernelType = "GOWER"

            catKernelType = str(catKernelType).upper().strip()

            allowed = {"CONT_RELAX", "HOMO_HSPHERE", "GOWER"}
            if catKernelType not in allowed:
                raise ValueError(
                    f"Invalid catKernelType='{catKernelType}'. "
                    f"Expected one of {sorted(allowed)}."
                )

            self.catKernelType = catKernelType

        else:
            self.catKernelType = None
        # ------------------------------------------------ #


        # ------------------------------------------ #
        # Model parameters for GP/SMT (to reuse after training GPs)
        self.gp_models = None
        self.gp_design_space = None
        self._gp_kernel_distance = None  # callable set after GP training
        # ------------------------------------------ #


    def getCatKernelType(self) -> str:
        """
        Get the categorical kernel type used by the GP surrogate.

        Returns:
        - str: One of {"GOWER", "CONT_RELAX", "HOMO_HSPHERE"}.
        """
        return self.catKernelType
    

    # Utility function creating the kernel used for comparing points
    def _init_gp_kernel_distance(self, trained_surrogate: MixedIntegerKrigingModel) -> None:
        """
        Prepare a callable self._gp_kernel_distance(x_mix, y_mix) that returns
        the kernel-induced distance a+b-2c using SMT internals.

        Must be called AFTER:
        - self.gp_design_space is set
        - trained_surrogate is trained (optimal_theta exists)
        """

        if self.gp_design_space is None:
            raise RuntimeError("self.gp_design_space is None; cannot init GP kernel distance.")

        design_space = self.gp_design_space
        cat_kernel_type = self.catKernelType

        allowed = {"CONT_RELAX", "EXP_HOMO_HSPHERE", "HOMO_HSPHERE", "GOWER"}
        if cat_kernel_type not in allowed:
            raise ValueError(f"Invalid cat_kernel_type='{cat_kernel_type}'. Allowed: {sorted(allowed)}")

        def mixed_kernel(x1, x2):
            x1 = np.atleast_2d(x1)
            x2 = np.atleast_2d(x2)

            x1x2, is_acting, n_eval, ij, Lij, dx = trained_surrogate._surrogate._predict_init(
                np.concatenate((x1, x2)), None
            )

            if not trained_surrogate._surrogate.is_continuous:
                x1c, is_acting = design_space.correct_get_acting(x1)
                _, ij = cross_distances(x1c, x2)

                dx = gower_componentwise_distances(
                    x1c,
                    x_is_acting=is_acting,
                    design_space=design_space,
                    hierarchical_kernel=MixHrcKernelType.ARC_KERNEL,
                    y=x2,
                    y_is_acting=design_space.correct_get_acting(x2)[1],
                )[0]

                if cat_kernel_type == "CONT_RELAX":
                    X_train2 = trained_surrogate._surrogate.design_space.unfold_x(x1x2)[0]
                    dx = cross_distances(X_train2)[0]

                listcatdecreed = design_space.is_conditionally_acting[trained_surrogate._surrogate.cat_features]
                if np.any(listcatdecreed):
                    dx = trained_surrogate._surrogate._correct_distances_cat_decreed(
                        dx,
                        is_acting,
                        listcatdecreed,
                        ij,
                        is_acting_y=design_space.correct_get_acting(x2)[1],
                        mixint_type=MixIntKernelType.GOWER,
                    )

                Lij, _ = cross_levels(X=x1x2, ij=ij, design_space=design_space, y=x2)

            cat_kernel_enum = getattr(MixIntKernelType, cat_kernel_type)

            corr = trained_surrogate._surrogate._matrix_data_corr(
                corr="abs_exp",
                design_space=design_space,
                power=1.9,
                theta=trained_surrogate._surrogate.optimal_theta,
                theta_bounds=trained_surrogate._surrogate.options["theta_bounds"],
                dx=dx,
                Lij=Lij,
                n_levels=trained_surrogate._surrogate.n_levels,
                cat_features=trained_surrogate._surrogate.cat_features,
                cat_kernel=cat_kernel_enum,
                x=x1x2,
                kplsk_second_loop=False,
            )
            return corr

        def kernel_distance(x_mix, y_mix) -> float:
            a = mixed_kernel(x_mix, x_mix)
            b = mixed_kernel(y_mix, y_mix)
            c = mixed_kernel(x_mix, y_mix)
            return float((a + b - 2.0 * c).squeeze())

        self._gp_kernel_distance = kernel_distance



    # CALLBACK for constructing the surrogate model from the cache
    def constructModel(self, **kwargs) -> bool:
        """
        Construct the active surrogate model from the cache.

        Parameters
        ----------
        seed : int, default=0
            Seed for stochastic parts.

        Returns
        -------
        bool
            True if training succeeded (GP may return False on warning-trigger case).
            # bool wasUpdated
            # bool isModelReady
        """

        # No component of the optimizer currently requires a surrogate model.
        if not self.isModelRequired:
            return False

        # ------------ Model update stopping rules -------- #
        cache_size = len(self.cache)
        isHardCapReached = cache_size >= self.modelUpdateHardCap

        # The model must always be constructed at least once when required.
        #
        # If repeated model updates are disabled, keep using the last
        # successfully trained model.
        if self.wasModelUpdatedOnce and not self.isModelUpdated:
            return self.modelIsReady

        # Repeated model updates stop once the hard cap is reached.
        #
        # If BO is active, BOSearch() will detect that model updates have
        # been disabled and will consequently stop BO.
        #
        # If BO is not active, the last trained model remains available
        # for the SURROGATE categorical poll.
        if self.wasModelUpdatedOnce and isHardCapReached:
            self.isModelUpdated = False
            return self.modelIsReady
        # ------------------------------------------------- #


        # ------------------------- Model construction below ------------- #
        self.modelIsReady = False

        # Enforce single-objective (current design choice)
        if len(self.OBJ_cols) != 1:
            raise ValueError(
                f"DataManager enforces exactly 1 objective, got {len(self.OBJ_cols)}: {self.OBJ_cols}"
            )

        # ------------------------------------------------------------ #
        # filtered training data 
        X_train, Y_obj, Y_pb = self._get_valid_training_data(large_value=1e23)

        if X_train.size == 0 or X_train.shape[0] < 2:
            raise ValueError(
                f"Not enough valid cached points to train GP models. "
                f"Got {X_train.shape[0]} valid points. "
                f"Try changing the random seed or manually providing valid points in the problem object."
            )

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
        # ------------------------------------------------------------ #


        # ----------------------------------------------------- #
        # SMT
        print("Constructing GP surrogate models (SMT) from cache...")
        
        n_start_mult = int(kwargs.get("n_start_mult", 15))

        # Kernel configuration
        categorical_kernel = self.catKernelType

        if not hasattr(MixIntKernelType, categorical_kernel):
            raise ValueError(
                f"Unknown categorical_kernel='{categorical_kernel}'. "
                f"Expected MixIntKernelType attribute name."
            )

        cat_kernel_enum = getattr(MixIntKernelType, categorical_kernel)

        # Reset GP-specific state
        self.gp_models = None
        self.gp_design_space = None

        # Local helper (detect your specific SMT warning pattern while capturing output)
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


        # SMT expects discrete columns to be integers
        Xsmt = np.array(X_train, copy=True)

        if self.cat_idx:
            Xsmt[:, self.cat_idx] = np.rint(Xsmt[:, self.cat_idx]).astype(int)
        if self.bin_idx:
            Xsmt[:, self.bin_idx] = np.rint(Xsmt[:, self.bin_idx]).astype(int)

        int_idx = (self.idx_by_type.get("integer", []) or [])
        real_idx = (self.idx_by_type.get("real", []) or [])

        if int_idx:
            Xsmt[:, int_idx] = np.rint(Xsmt[:, int_idx]).astype(int)

        # ------------------ DesignSpace in NOMAD variable order ------------------ #
        quant_bound_map = {j: self.quant_bounds[k] for k, j in enumerate(self.quant_idx)}

        var_ds = []
        cat_local = 0
        for j in range(self.nbVars):
            if j in self.cat_idx:
                Lj = int(self.nbCategoriesPerVariable[cat_local])
                var_ds.append(CategoricalVariable(list(range(Lj))))
                cat_local += 1

            elif j in self.bin_idx:
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
                raise ValueError(f"Variable x{j} not classified in choice/binary/integer/real indices.")

        ds = DesignSpace(var_ds)
        self.gp_design_space = ds

        # Common SMT model factory (forced settings)
        def make_gp():
            return MixedIntegerKrigingModel(
                surrogate=KRG(
                    design_space=copy.deepcopy(ds),
                    categorical_kernel=cat_kernel_enum,
                    hyper_opt="Cobyla",
                    corr="squar_exp",
                    n_start=max(1, int(n_start_mult * self.nbVars)),
                    print_prediction=False,
                ),
            )

        # ------------------ objective GP ------------------ #
        sm_f = make_gp()
        warning_triggered = _train_with_warnings(sm_f, copy.deepcopy(Xsmt), Y_obj)
        if warning_triggered:
            print("\n>Warning triggered during training of objective. GP construction aborted.\n")
            self.modelIsReady = False
            return False

        # ------------------ PB constraint GPs ------------------ #
        sm_pbs = []
        for i in range(Y_pb.shape[1]):
            sm_c = make_gp()
            warning_triggered = _train_with_warnings(sm_c, copy.deepcopy(Xsmt), Y_pb[:, i])
            if warning_triggered:
                print(f"\n>Warning triggered during training of PB constraint {i}. GP construction aborted.\n")
                self.modelIsReady = False
                return False
            sm_pbs.append(sm_c)
        # ----------------------------------------------------- #

        # Store
        self.gp_models = {"objective": sm_f, "PB": sm_pbs}
        self.modelIsReady = True
        self.wasModelUpdatedOnce = True

        self._init_gp_kernel_distance(sm_f)

        print("GP objective and PB-constraint surrogates trained successfully.")
        return self.modelIsReady
            
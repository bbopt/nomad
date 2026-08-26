from typing import Optional

from ..doe import DOEgenerator
from .base_optimizer import BaseOptimizer


class QuantitativeOptimizer(BaseOptimizer):
    """
    Optimizer configuration for quantitative problems.
    """
    
    def __init__(
        self,
        seed: int = 0,
        initialDOEgenerator: Optional[DOEgenerator] = None,
        quantPollType: str = "ADS",
        isBOSearchUsed: bool = False,   # Attribute passed to a DataManager() in NomadBBO_src.pyx 
        modelUpdateHardCap: int = 500,  # Attribute passed to BaseOptimizer(), then to a DataManager() in NomadBBO_src.pyx
        efiStopTol: float = 1e-3,       # Attribute passed to BaseOptimizer()
        efiStopPatience: int = 3        # Attribute passed to BaseOptimizer()
    ):
        super().__init__(
            seed=seed,
            initialDOEgenerator=initialDOEgenerator,
            quantPollType=quantPollType,
            isBOSearchUsed=isBOSearchUsed,
            modelUpdateHardCap=modelUpdateHardCap,
            efiStopTol=efiStopTol,
            efiStopPatience=efiStopPatience
        )

    
        # Model update policy
        # Quantitative surrogate models are only used by BO.
        self.updateModel = self.isBOSearchUsed

    @property
    def isModelRequired(self):
        return self.isBOSearchUsed
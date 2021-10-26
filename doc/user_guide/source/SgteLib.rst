.. _sgtelib:

Surrogate Library
========================

The *SGTELIB* library is a dynamic surrogate modeling library. It is used in the *Search* step of Mads to dynamically construct models from the previous evaluations.
During a *Search* step that uses *SGTELIB*, the selected models of the objective and the constraints are contrusted and a surrogate subproblem involving these modes is optimized.
The resulting solutions are the next candidate for evaluation by the true problem.

| Models from the *SGTELIB* library can be used by setting ``SGTELIB_MODEL_SEARCH`` to ``yes``.


Models
-------------------

Models in sgtelib are defined by using a succession of field names and field values.
To choose a model, the ``SGTELIB_MODEL_DEFINITION`` parameter must be used followed by the field name ``TYPE``, and then by the model type.
The subsequent fields enable to define the settings of the model.
Each field name is made of one single word and each field value is made of one single word or numerical value.

Example : ``SGTELIB_MODEL_DEFINITION TYPE <model type> FIELD1 <field 1 value> FIELD2 <field 2 value>``

The section below describes the models and settings availables.


Types of models
""""""""""""""""""""""

.. _prs:

``PRS``
""""""""
| PRS (Polynomial Response Surface) is a type of model.
| Authorized fields:

* :ref:`degree` (Can be optimized)

* :ref:`ridge` (Can be optimized)

* :ref:`budget`: Defines the budget allocated for parameter optimization.

* :ref:`output`: Defines the output text file.

| Examples:
| ``TYPE PRS DEGREE 2``
| ``TYPE PRS DEGREE OPTIM RIDGE OPTIM``


.. _prs_edge:

``PRS_EDGE``
""""""""""""""

| PRS_EDGE (Polynomial Response Surface EDGE) is a type of model that allows to model discontinuities at :math:`0` by using additional basis functions.
| Authorized fields:

* :ref:`degree` (Can be optimized)

* :ref:`ridge` (Can be optimized)

* :ref:`budget`: Defines the budget allocated for parameter optimization.

* :ref:`output`: Defines the output text file.

| Examples:
| ``TYPE PRS_EDGE DEGREE 2``
| ``TYPE PRS_EDGE DEGREE OPTIM RIDGE OPTIM``


.. _prs_cat:

``PRS_CAT``
""""""""""""""
| PRS_CAT (Categorical Polynomial Response Surface) is a type of model that allows to build one PRS model for each dierent value of the first component of :math:`x`.
| Authorized fields:

* :ref:`degree` (Can be optimized)
* :ref:`ridge` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE PRS_CAT DEGREE 2``
| ``TYPE PRS_CAT DEGREE OPTIM RIDGE OPTIM``


.. _rbf:

``RBF``
""""""""""""""
| RBF (Radial Basis Function) is a type of model.
| Authorized fields:

* :ref:`kernel_type` (Can be optimized)
* :ref:`kernel_shape` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`ridge` (Can be optimized)
* :ref:`preset`: Defines the type of RBF model used.
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE RBF KERNEL_TYPE D1 KERNEL_SHAPE OPTIM DISTANCE TYPE NORM2``


.. _ks:

``KS``
""""""""""""""
| KS (Kernel Smoothing) is a type of model.
| Authorized fields:

* :ref:`kernel_type` (Can be optimized)
* :ref:`kernel_shape` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM`` 


.. _kriging:

``KRIGING``
""""""""""""""
| KRIGING is a type of model.
| Authorized fields:

* :ref:`ridge` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE KRIGING``


.. _lowess:

``LOWESS``
""""""""""""""
| LOWESS (Locally Weighted Regression) is a type of model [TaAuKoLed2016]_.
| Authorized fields:

* :ref:`degree`: Must be 1 (default) or 2 (Can be optimized).
* :ref:`ridge` (Can be optimized)
* :ref:`kernel_type` (Can be optimized)
* :ref:`kernel_shape` (Can be optimized)
* :ref:`distance_type` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE LOWESS DEGREE 1``
| ``TYPE LOWESS DEGREE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE D1``
| ``TYPE LOWESS DEGREE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE TYPE OPTIM``


.. _cn:

``CN``
""""""""""""""
| CN (Closest Neighbours) is a type of model.
| Authorized fields:

* :ref:`distance_type` (Can be optimized)
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE CN``


.. _ensemble:

``ENSEMBLE``
""""""""""""""
| ENSEMBLE is a type of model.
| Authorized fields:

* :ref:`weight`: Defines how the ensemble weights are computed.
* :ref:`metric`: Defines which metric is used to compute the weights.
* :ref:`budget`: Defines the budget allocated for parameter optimization.
* :ref:`distance_type`: This parameter is transfered to the models contained in the Ensemble.
* :ref:`output`: Defines the output text file.

| Example:
| ``TYPE ENSEMBLE WEIGHT SELECT METRIC OECV``
| ``TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV DISTANCE TYPE NORM2 BUDGET 100``


.. _ensemble_stat:

``ENSEMBLE_STAT``
""""""""""""""""""

[AuLedSa2021]_

TODO


The following table summarizes the possible fields for every model.

.. csv-table:: Model authorized fields
   :header: "Model type", :ref:`degree`, :ref:`ridge`, :ref:`kernel_type`, :ref:`kernel_shape`, :ref:`distance_type`, :ref:`preset`, :ref:`weight`, :ref:`metric`, :ref:`uncertainty`,:ref:`budget`, :ref:`output`

   :ref:`prs`,          ✔,  ✔,  ,    ,    ,   ,  ,  ,   , ✔, ✔
   :ref:`prs_edge`,     ✔,  ✔,  ,    ,    ,   ,  ,  ,  , ✔, ✔
   :ref:`prs_cat`,      ✔,  ✔,  ,    ,    ,   ,  ,  ,  , ✔, ✔
   :ref:`rbf`,           ,  ✔,  ✔,  ✔,  ✔, ✔,   ,  ,  , ✔, ✔
   :ref:`ks`,            ,   ,  ✔,  ✔,  ✔,   ,   ,  ,  , ✔, ✔
   :ref:`kriging`,       ,  ✔,  ,    ,   ✔,  ,   ,  ,  , ✔, ✔
   :ref:`lowess`,       ✔, ✔, ✔,  ✔,   ✔,   ,   ,  ,  , ✔, ✔
   :ref:`cn`,            ,   ,  ,    ,   ✔,  ,    ,  ,  , ✔, ✔
   :ref:`ensemble`,      ,   ,  ,    ,   ✔,  ,  ✔, ✔,  , ✔, ✔
   :ref:`ensemble_stat`, ,   ,  ,    ,   ✔,  ,  ✔, ✔, ✔, ✔, ✔


Main model parameters
""""""""""""""""""""""""""

.. _degree:

``DEGREE``
""""""""""""""
| The field name DEGREE defines the degree of a polynomial response surface. The value must be an integer :math:`\geq 1`.
| Allowed for models of type: :ref:`prs`, :ref:`prs_edge`, :ref:`prs_cat`, :ref:`lowess`.
| Default value: 5

* For PRS models, the default degree is 2.
* For LOWESS models, the degree must be 1 (default) or 2.

| Example:
| ``TYPE PRS DEGREE 3 defines a PRS model of degree 3.``
| ``TYPE PRS_EDGE DEGREE 2 defines a PRS_EDGE model of degree 2.``
| ``TYPE LOWESS DEGREE OPTIM defines a LOWESS model where the degree is optimized.``


.. _ridge:

``RIDGE``
""""""""""""""
| The field name RIDGE defines the regularization parameter of the model.
| Allowed for models of type: :ref:`prs`, :ref:`prs_edge`, :ref:`prs_cat`, :ref:`lowess`, :ref:`rbf`.
| Possible values: Real value :math:`\geq 0`. Recommended values are :math:`0` and :math:`0.001`.
| Default value: :math:`0.001`.

| Example:
| ``TYPE PRS DEGREE 3 RIDGE 0`` defines a PRS model of degree 3 with no ridge.
| ``TYPE PRS DEGREE OPTIM RIDGE OPTIM`` defines a PRS model where the degree and ridge coefficient are optimized.


.. _kernel_type:

``KERNEL_TYPE``
""""""""""""""""
| The field name KERNEL_TYPE defines the type of kernel used in the model. The field name ``KERNEL`` is equivalent.
| Allowed for models of type: :ref:`rbf`, :ref:`kriging`, :ref:`lowess` and :ref:`ks`.
| Possible values:

* ``D1``: Gaussian kernel (default)
* ``D2``: Inverse Quadratic Kernel
* ``D3``: Inverse Multiquadratic Kernel
* ``D4``: Bi-quadratic Kernel
* ``D5``: Tri-cubic Kernel
* ``D6``: Exponential Sqrt Kernel
* ``D7``: Epanechnikov Kernel
* ``I0``: Multiquadratic Kernel
* ``I1``: Polyharmonic splines, degree 1
* ``I2``: Polyharmonic splines, degree 2
* ``I3``: Polyharmonic splines, degree 3
* ``I4``: Polyharmonic splines, degree 4
* ``OPTIM``: The type of kernel is optimized

| Example:
| ``TYPE KS KERNEL_TYPE D2`` defines a KS model with Inverse Quadratic Kernel.
| ``TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM`` defines a KS model with optimized kernel shape and type.


.. _kernel_shape:

``KERNEL_SHAPE``
""""""""""""""""""
| The field name KERNEL_SHAPE defines the shape coefficient of the kernel function. Note that this field name has no impact for kernel types ``I1``, ``I2``, ``I3`` and ``I4`` because these kernels do not include a shape parameter.
| Allowed for models of type: :ref:`rbf`, :ref:`ks`, :ref:`kriging`, :ref:`lowess`.
| Possible values: Real value :math:`\geq 0`. Recommended range is :math:`[0:1; 10]`. For KS and LOWESS model, small values lead to smoother models.
| Default value: By default, the kernel coefficient is optimized.

| Example:
| ``TYPE RBF KERNEL_SHAPE 10`` defines a RBF model with a shape coefficient of :math:`10`.
| ``TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM`` defines a KS model with optimized kernel shape and type.


.. _distance_type:

``DISTANCE_TYPE``
""""""""""""""""""
| The field name DISTANCE TYPE defines the distance function used in the model.
| Allowed for models of type: :ref:`rbf`, :ref:`ks`, :ref:`lowess`, :ref:`cn`.
| Possible values:

* ``NORM1``: Euclidian distance
* ``NORM2``: Distance based on norm :math:`1`
* ``NORMINF``: Distance based on norm :math:`1`
* ``NORM2_IS0``: Tailored distance for discontinuity in :math:`0`
* ``NORM2_CAT``: Tailored distance for categorical models

| Default value: NORM2.

| Example:
| ``TYPE KS DISTANCE NORM2 IS0`` defines a KS model tailored for VAN optimization.


.. _preset:

``PRESET``
""""""""""""""
| The field name PRESET defines the type of RBF model used.
| Allowed for models of type: :ref:`rbf`.
| Possible values:

* ``O``: RBF with linear terms and orthogonal constraints
* ``R``: RBF with linear terms and regularization term
* ``I``: RBF with incomplete set of basis functions (see [AuKoLedTa2016]_ for RBFI models)

| Default value: I.

| Example:
| ``TYPE RBF PRESET O``

.. _weight:

``WEIGHT``
""""""""""""""
| The field name WEIGHT defines the method used to compute the weights :math:`\boldsymbol{w}` of the ensemble of models. The keyword ``WEIGHT_TYPE`` is equivalent.
| Allowed for models of type: :ref:`ensemble`, :ref:`ensemble_stat`.
| Possible values:

* ``WTA1``: :math:`w_k \propto \mathcal{E}_{sum} - \mathcal{E}_k`  (default)
* ``WTA3``: :math:`w_k \propto (\mathcal{E}_k + \alpha\mathcal{E}_{mean})^{\beta}`
* ``SELECT``: :math:`w_k \propto 1` if :math:`\mathcal{E}_k = \mathcal{E}_{min}`
* ``OPTIM``: :math:`\boldsymbol{w}` minimizes :math:`\mathcal{E}(\boldsymbol{w})`
* TODO

| Default value: TODO

| Example:
| ``TYPE ENSEMBLE WEIGHT SELECT METRIC RMSECV`` defines an ensemble of models which selects the model that has the best RMSECV.
| ``TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV`` defines an ensemble of models where the weights :math:`\boldsymbol{w}` are computed to minimize the RMSECV of the model.


.. _uncertainty:

``UNCERTAINTY``
""""""""""""""""
| The field name UNCERTAINTY defines the type of uncertainty used in ENSEMBLE_STAT models. 
| Allowed for models of type: :ref:`ensemble_stat`.
| Possible values:

* ``SMOOTH``: Smooth alternative of the uncertainty (default)
* ``NONSMOOTH``: Nonmooth alternative of the uncertainty

| Example:
| ``TYPE ENSEMBLE_STAT UNCERTAINTY NONSMOOTH``


.. _output:

``OUTPUT``
""""""""""""""
Defines a text file in which model information are recorded. Allowed for ALL types of model.


(specific to ENSEMBLE_STAT models)

``SIZE_PARAM``
""""""""""""""""

TODO

``SIGMA_MULT``
""""""""""""""""

``LAMBDA_P``
""""""""""""""

``LAMBDA_PI``
""""""""""""""


Parameter optimization and selection
""""""""""""""""""""""""""""""""""""""""

.. _optim:

``OPTIM``
""""""""""""""
| The field value OPTIM indicate that the model parameter must be optimized. The default optimization criteria is the AOECV error metric (except for ENSEMBLE_STAT models where it is OECV).
| Parameters that can be optimized:

* :ref:`degree`
* :ref:`ridge`
* :ref:`kernel_type`
* :ref:`kernel_shape`
* :ref:`distance_type`

| Example:
| ``TYPE PRS DEGREE OPTIM``
| ``TYPE LOWESS DEGREE OPTIM KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM METRIC ARMSECV``


.. _metric:

``METRIC``
""""""""""""""
| The field name METRIC defines the metric used to select the parameters of the model (including the weights of Ensemble models).
| Allowed for ALL types of model.
| Possible values:

* ``EMAX``: Error Max
* ``EMAXCV``: Error Max with Cross-Validation
* ``RMSE``: Root Mean Square Error
* ``RMSECV``: RMSE with Cross-Validation
* ``OE``: Order Error
* ``OECV``: Order Error with Cross-Validation [AuKoLedTa2016]_
* ``LINV``: Invert of the Likelihood
* ``AOE``: Aggregate Order Error
* ``AOECV``: Aggregate Order Error with Cross-Validation

| Default value: ``AOECV``. TODO

| Example:
| ``TYPE ENSEMBLE WEIGHT SELECT METRIC RMSECV`` defines an ensemble of models which selects the model that has the best RMSECV.


.. _budget:

``BUDGET``
""""""""""""""
| Budget for model parameter optimization. The number of sets of model parameters that are tested is equal to the optimization budget multiplied by the the number of parameters to optimize.
| Allowed for ALL types of model.
| Default value: :math:`20`

| Example:
| ``TYPE LOWESS KERNEL_SHAPE OPTIM METRIC AOECV BUDGET 100``
| ``TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV BUDGET 50``




Surrogate subproblem formulations
-------------------------------------

The *SGTELIB* library offers different formulations of the surrogate subproblem to be optimized at the *Search* step [TaLeDKo2014]_.
The ``SGTELIB_MODEL_FORMULATION`` parameter enables to choose a formulation, and the ``SGTELIB_MODEL_DIVERSIFICATION`` parameter enables to adjust a diversification parameter.


``SGTELIB_MODEL_FORMULATION``
""""""""""""""""""""""""""""""

| The formulations of the surrogate subproblem involve various quantities.
| :math:`\hat f` denotes a model of the objective :math:`\hat f` and :math:`\hat c_j` a model of the constraint :math:`c_j`, :math:`j=1,2,\dots,m`. For :math:`x\in X`, :math:`\sigma_f(x)` denotes the uncertainty associated with the prediction :math:`\hat f(x)`, and :math:`\sigma_j(x)` the uncertainty associated with the prediction :math:`\hat c_j(x)`, :math:`j=1,2,\dots,m`. This uncertainty depends on the model chosen.

| For a :ref:`kriging` model, :math:`\sigma_f(x)` (or :math:`\sigma_j(x)`) is readily available through the standard deviation that the model natively produces.
| For an :ref:`ensemble_stat` model, the uncertainty is constructed by comparing the predictions of the ensemble models as in [AuLedSa2021]_.
| For any other model except ENSEMBLE, :math:`\sigma_f(x)` (or :math:`\sigma_j(x)`) is computed with the distance from :math:`x` to previously evaluated points.
| Finally, for an :ref:`ensemble` model, the uncertainty is computed through a weighted sum of the squared uncertainties of the ensemble models.

| There are eight different formulations that can be chosen with the parameter ``SGTELIB_MODEL_FORMULATION``. Some formulations involve a :math:`\lambda` parameter that is described later.

* ``FS`` (default):

.. math::

      \min_{x\in X}&\ \ \hat f(x)-\lambda\hat\sigma_f(x) \\
      \mathrm{s.t.}&\ \ \hat c_j(x)-\lambda\hat\sigma_j(x)\leq0,\ \ j=1,2,\dots,m

* ``FSP``:

.. math::

      \min_{x\in X}&\ \ \hat f(x)-\lambda\hat\sigma_f(x) \\
      \mathrm{s.t.}&\ \ \mathrm{P}(x)\geq 0.5

where :math:`\mathrm{P}` is the *probability of feasibility* which is the probability that a given point is feasible.

* ``EIS``:

.. math::

      \min_{x\in X}&\ -\mathrm{EI}(x)-\lambda\hat\sigma_f(x) \\
      \mathrm{s.t.}&\ \ \hat c_j(x)-\lambda\hat\sigma_j(x)\leq0,\ \ j=1,2,\dots,m

where :math:`\mathrm{EI}` is the *expected improvement* that takes into account the probability of improvement and
the expected amplitude thereof.

* ``EFI``:

.. math::
 
      \min_{x\in X}\ -\mathrm{EFI}(x)

where :math:`\mathrm{EFI}` is the *expected feasible improvement* : :math:`\mathrm{EFI} = \mathrm{EI}\times\mathrm{P}`

* ``EFIS``:

.. math::
  
      \min_{x\in X}\ -\mathrm{EFI}(x)-\lambda\hat\sigma_f(x)

* ``EFIM``:

.. math::
  
      \min_{x\in X}\ -\mathrm{EFI}(x)-\lambda\hat\sigma_f(x)\mu(x)

where :math:`\mu` is the *uncertainty in the feasibility* : :math:`\mu = 4\mathrm{P}\times(1-\mathrm{P})`

* ``EFIC``:

.. math::

      \min_{x\in X}\ -\mathrm{EFI}(x)-\lambda(\mathrm{EI}(x)\mu(x)
      +\mathrm{P}(x)\hat\sigma_f(x))

* ``PFI``:

.. math::
  
      \min_{x\in X}\ -\mathrm{PFI}(x)

where :math:`\mathrm{PFI}` is the *probability of improvement* : :math:`\mathrm{PFI} = \mathrm{PI}\times\mathrm{P}`,
with :math:`\mathrm{PI}` being the *probability of improvement* which is the probability that the objective decreases from the best known value at a given point.


| Example:
| ``SGTELIB_MODEL_DEFINITION TYPE KRIGING``
| ``SGTELIB_MODEL_FORMULATION EFIC``
| The two lines above define a surrogate subproblem based on the EFIC formulation that will involve kriging models.


``SGTELIB_MODEL_DIVERSIFICATION``
""""""""""""""""""""""""""""""""""

| The exploration parameter :math:`\lambda` enables to control the exploration of the search space against the intensification in the most promising areas. A higher :math:`\lambda` favors exploration whereas a lower :math:`\lambda` favors intensification.

| :math:`\lambda` is a real value in :math:`[0,1]` defined by the parameter ``SGTELIB_MODEL_DIVERSIFICATION``.
| Default value : :math:`0.01`.

| Example:
| ``SGTELIB_MODEL_DEFINITION TYPE ENSEMBLE``
| ``SGTELIB_MODEL_FORMULATION FSP``
| ``SGTELIB_MODEL_DIVERSIFICATION 0.1``
| The three lines above define a surrogate subproblem based on the FSP formulation with an exploration parameter equals to :math:`0.1` that will involve ensemble models.



.. topic:: References


  .. [TaAuKoLed2016] B.Talgorn, C.Audet, M.Kokkolaras and S.Le Digabel.
    Locally weighted regression models for surrogate-assisted design optimization.
    *Optimization and Engineering*, 19(1):213–238, 2018.
  
  .. [TaLeDKo2014] B.Talgorn, S.Le Digabel and M.Kokkolaras.
    Statistical Surrogate Formulations for Simulation-Based Design Optimization.
    *Journal of Mechanical Design*, 137(2):021405–1–021405–18, 2015
  
  .. [AuKoLedTa2016] C.Audet, M.Kokkolaras, S.Le Digabel and B.Talgorn.
    Order-based error for managing ensembles of surrogates in mesh adaptive direct search
    *Journal of Global Optimization*, 70(3):645–675, 2018.

  .. [AuLedSa2021] C.Audet, S.Le Digabel and R.Saltet.
    Quantifying uncertainty with ensembles of surrogates for blackbox optimization.
    Rapport technique G-2020-58, Les cahiers du GERAD, 2020.
    http://www.optimization-online.org/DB_HTML/2021/07/8489.html
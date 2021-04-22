.. _tricks_of_the_trade:

Tricks of the trade
===================

NOMAD has default values for all algorithmic parameters. These values represent a compromise between robustness and performance obtained by developers on sets of problems used for benchmarking. But you might want to improve NOMAD performance for your problem by tuning the parameters or use advanced functionalities. The following sections provide tricks that may work for you.

Here are a few suggestions for tuning NOMAD when facing different symptoms. The suggestions can be tested one by one or all together.

.. csv-table:: Suggestions for tuning NOMAD
   :header: "Symptom", "Suggestion", "Ref."
   :widths: 20,7,20

   I want to see more display ,	Increase display degree , :ref:`display_degree`
   Quantifiable constraints	,	Try PB  EB or combinations , :ref:`bb_output_type`
   Difficult constraint , Try PB instead of EB , :ref:`bb_output_type`
   No initial point		,	 	Add a LH search , :ref:`LH Search and X0 <x0>`
   Variables of different magnitudes , 	Change blackbox input scaling	, :ref:`create_blackbox_program`
   " ", Change :math:`\Delta_0` per variable , :ref:`initial_mesh_size`
   " ",	Tighten bounds , :ref:`bounds`
   Many variables ,	Fix some variables  , :ref:`fixed_variable`
   " ",	Use *PSD-MADS*  , :ref:`psd_mads`
   Unsatisfactory solution ,	Change direction type to ``N+1 UNI`` , :ref:`direction_type`
   " ", 	Change initial point , :ref:`LH Search and X0 <x0>`
   " ", 	Add a LH search , :ref:`LH Search and X0 <x0>`
   " ", 	Tighten bounds , :ref:`bounds`
   " ", 	Change :math:`\Delta_0` , :ref:`initial_mesh_size`
   " ", 	Modify seeds that affect algorithms , :ref:`seed`
   " ", 	Disable quadratic models , set ``QUAD_MODEL_SEARCH no``
   " ", 	Unable *SGTELIB* models , set ``SGTELIB_MODEL_SEARCH yes``
   " ",   Disable opportunistic evaluations, set ``EVAL_OPPORTUNISTIC no``
   " ",  Disable anisotropic mesh , set ``ANISOTROPIC_MESH no``
   " ",  Change anisotropy factor , set ``ANISOTROPY_FACTOR 0.05``
   Improvements get negligible ,	Change stopping criteria	, Type ``nomad -h stop``
   " ",	Disable quadratic models  , set ``QUAD_MODEL_SEARCH no``
   It takes long to improve :math:`f` ,	Decrease :math:`\Delta_0` , :ref:`initial_mesh_size`
   Optimization is time consuming	, 	Perform parallel blackbox evaluations , :ref:`block_evaluations` and :ref:`parallel_evaluations`
   Blackbox is not that expensive		, Setup maximum wall-clock time , remove ``MAX_BB_EVAL`` and set ``MAX_TIME``
   " ", 	Add a LH search , :ref:`LH Search and X0 <x0>`

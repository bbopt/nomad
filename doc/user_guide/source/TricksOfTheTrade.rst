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
   No initial point		,	 	Add a LH search , :ref:`x0`
   Variables of widely , 	Change blackbox input scaling	, :ref:`create_blackbox_program`
   Different magnitudes	, Change :math:`\Delta_0` per variable , :ref:`initial_mesh_size`

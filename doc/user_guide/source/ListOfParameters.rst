Detailed information
^^^^^^^^^^^^^^^^^^^^


BB_INPUT_TYPE

::

  Type: NOMAD::BBInputTypeList

  Default: * R

  Description:

  . Blackbox input types

  . List of types for each variable

  . Available types:
    . B: binary
    . I: integer
    . R: continuous

  . Examples:
    . BB_INPUT_TYPE * I       # all variables are integers
    . BB_INPUT_TYPE ( R I B ) # for all 3 variables
    . BB_INPUT_TYPE 1-3 B     # NOT YET SUPPORTED ( variables 1 to 3 are binary )
    . BB_INPUT_TYPE 0 I       # NOT YET SUPPORTED ( first variable is integer )

.. _dimension_2:

DIMENSION

::

  Type: size_t

  Default: 0

  Description :

  . Number of variables

  . Argument: one positive integer

  . Example: DIMENSION 3

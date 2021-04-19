.. _library_mode:

Optimization in library mode
----------------------------

The library mode allows to tailor the evaluation of the objectives and constraints within a
specialized executable that contains NOMAD shared object library.

For example, it is possible to link your own code with the NOMAD library (provided during installation)
in a single executable that can define and run optimization for your problem. Contrary to the batch
mode, this has the disadvantage that a crash within the executable (for example during the evaluation of a point)
will end the optimization unless a special treatment of exception is provided by the user.
But, as a counterpart, it offers more options and flexibility for blackbox integration and
optimization management (display, pre- and post-processing, multiple optimizations, user search, etc.).

The library mode requires additional coding and compilation before conducting optimization.
First, we will briefly review the compilation of source code to obtain NOMAD binaries
(executable and shared object libraries) and how to use library.
Then, details on how to interface your own code are presented.

Compilation of the source code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOMAD source code files are located in ``$NOMAD_HOME/src``.
Examples are provided in ``$NOMAD_HOME/examples/basic/library`` and ``$NOMAD_HOME/examples/advanced/library``.

The compilation procedure uses the provided ``CMake`` files along with the source code.

In what follows it is supposed that you have a write access to the source codes directory.
If it is not the case, please consider making a copy in a more convenient location.

Using NOMAD libraries
^^^^^^^^^^^^^^^^^^^^^

Using the routines that are in the pre-compiled NOMAD shared object libraries (so) or dll
(not yet available for Windows) with a ``C++`` program requires building an executable
(:ref:`installation` describes how to build the libraries and the examples). This is illustrated on the example located in the directory::

  $NOMAD_HOME/examples/basic/library/example1

It is supposed that the environment variable ``NOMAD_HOME`` is defined and NOMAD shared
object libraries are built. A basic knowledge of object oriented programming with ``C++`` is assumed.
For this example, just one ``C++`` source file is used, but there could be a lot more.

Test the basic example 1
""""""""""""""""""""""""

Let us first test the basic example to check that libraries are working fine and accessible.
Library mode examples are built during the installation procedure (unless the flag ``BUILD_LIBMODE_EXAMPLES`` is set to ``OFF``)::

  > cd $NOMAD_HOME/examples/basic/library/example1
  > ls
  CMakeLists.txt		example1_lib.cpp	example1_lib.exe
  > ./example1_lib.exe
  All variables are granular. MAX_EVAL is set to 1000000 to prevent algorithm from circling around best solution indefinetely
  BBE OBJ
  1 -28247.525326  (Phase One)
  12   -398.076167  (Phase One)
  55   -413.531262
  81  -1084.90725
  136  -1632.176507
  188  -1754.758402
  201  -1787.5835
  260  -1967.858372
  764  -1967.860497
  765  -1967.866871
  766  -1967.885992
  767  -1967.943355
  768  -1968.115446
  769  -1968.63171
  770  -1970.180451
  771  -1974.828012
  772  -1988.763221
  862  -1989.292074
  863  -1990.878576
  864  -1995.637558
  895  -1999.023493
  940  -1999.116474
  942  -1999.395208
  957  -1999.556452
  961  -1999.835309
  A termination criterion is reached: Max number of blackbox evaluations (Eval Global) No more points to evaluate 1001

  Best feasible solution:     #44485 ( 1.81678 5.21878 4.40965 8.1 15.1 8.8 5 10.1 1.6 5.5 )	Evaluation OK	 f = -1999.835309000000052     	 h =   0

  Best infeasible solution:   #44525 ( -26570.1 0 -5895.58 -3.58722e+06 -8.60934e+06 -5.73956e+06 -1.36315e+07 5.73957e+06 1.14791e+07 2.86978e+06 )	Evaluation OK	 f = -1151.1639900000000125    	 h =   0.5625

  Blackbox evaluations:        1001
  Total model evaluations:     41241
  Cache hits:                  104
  Total number of evaluations: 1105

Modify ``CMake`` files
""""""""""""""""""""""

As a first task, you can create a ``CMakeLists.txt`` for your source code(s) based on the one for the basic example 1.


.. TODO add the CMake procedure for an example out of Nomad subdirectories.

.. code-block:: cmake

  add_executable(example1_lib.exe example1_lib.cpp )
  target_include_directories(example1_lib.exe PRIVATE ${CMAKE_SOURCE_DIR}/src)
  set_target_properties(example1_lib.exe PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  if(OpenMP_CXX_FOUND)
    target_link_libraries(example1_lib.exe PUBLIC nomadAlgos nomadUtils nomadEval OpenMP::OpenMP_CXX)
  else()
    target_link_libraries(example1_lib.exe PUBLIC nomadAlgos nomadUtils nomadEval)
  endif()

  # installing executables and libraries
  install(TARGETS example1_lib.exe  RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR} )

  # Add a test for this example
  if(BUILD_TESTS MATCHES ON)
     message(STATUS "    Add example library test 1")

     # Can run this test after install
     add_test(NAME Example1BasicLib COMMAND ${CMAKE_BINARY_DIR}/examples/runExampleTest.sh ./example1_lib.exe WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
  endif()

If you include your problem into the ``$NOMAD_HOME/examples`` directories, you just need to copy
the example ``CMakeLists.txt`` into your own problem directory (for example ``$NOMAD_HOME/examples/basic/library/myPb``),
change the name ``example1_lib`` with your choice and add the subdirectory into ``$NOMAD_HOME/examples/CMakeLists.txt``::

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/basic/library/myPb)


Modify ``C++`` files
""""""""""""""""""""

We now describe the other steps required for the creation of the source file (let us use ``example1.cpp``)
which is divided into two parts: a class for the description of the problem, and the main function.

The use of standard ``C++`` types for reals and vectors is of course allowed within your code, but it
is suggested that you use the NOMAD types as much as  possible. For reals, NOMAD uses the class ``NOMAD::Double``,
and for vectors, the classes ``NOMAD::Point`` or ``NOMAD::ArrayOfDouble``.
A lot of functionalities have been coded for theses classes, which are visible  in files ``$NOMAD_HOME/src/Math/*.hpp``.

The namespace \sComp{NOMAD} is used for all NOMAD types, and you must type ``NOMAD::``
in front of all types  unless you type ``using namespace NOMAD;``  at the beginning of your program.

Providing the blackbox evaluation of objective and constraints directly in the code avoids
the use of temporary files and system calls by the algorithm. This is achieved by defining a derived
class (let us call it ``My_Evaluator``) that inherits from the class ``NOMAD::Evaluator``.
The blackbox evaluation is programmed in a user-defined class that will  be automatically called by the algorithm.}

.. code-block:: c++

  /**
   \file   example1_lib.cpp
   \brief  Library example for nomad
   \author Viviane Rochon Montplaisir
   \date   2017
   */

  #include "Nomad/nomad.hpp"

  /*----------------------------------------*/
  /*               The problem              */
  /*----------------------------------------*/
  class My_Evaluator : public NOMAD::Evaluator
  {
  public:
      My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
      : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
      {}

      ~My_Evaluator() {}

      bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override
      {
          bool eval_ok = false;
          // Based on G2.
          NOMAD::Double f = 1e+20, g1 = 1e+20, g2 = 1e+20;
          NOMAD::Double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, prod1 = 1.0, prod2 = 1.0;
          size_t n = x.size();

          try
          {
              for (size_t i = 0; i < n ; i++)
              {
                  sum1  += pow(cos(x[i].todouble()), 4);
                  sum2  += x[i];
                  sum3  += (i+1)*x[i]*x[i];
                  prod1 *= pow(cos(x[i].todouble()), 2);
                  if (prod2 != 0.0)
                  {
                      if (x[i] == 0.0)
                      {
                          prod2 = 0.0;
                      }
                      else
                      {
                          prod2 *= x[i];
                      }
                  }
              }

              g1 = -prod2 + 0.75;
              g2 = sum2 -7.5 * n;

              f = 10*g1 + 10*g2;
              if (0.0 != sum3)
              {
                  f -= ((sum1 -2*prod1) / sum3.sqrt()).abs();
              }
              // Scale
              if (f.isDefined())
              {
                  f *= 1e-5;
              }

              NOMAD::Double c2000 = -f-2000;
              auto bbOutputType = _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
              std::string bbo = g1.tostring();
              bbo += " " + g2.tostring();
              bbo += " " + f.tostring();
              bbo += " " + c2000.tostring();

              x.setBBO(bbo);

              eval_ok = true;
          }
          catch (std::exception &e)
          {
              std::string err("Exception: ");
              err += e.what();
              throw std::logic_error(err);
          }

          countEval = true;
          return eval_ok;
      }
    };


The argument ``x`` (in/out in ``eval_x()``) corresponds to an evaluation point, i.e. a vector containing the
coordinates of the point to be evaluated, and also the result of the evaluation.
The coordinates are accessed with the operator ``[]`` (``x[0]`` for the first coordinate) and outputs are set with ``x.setBBO(bbo);``.
The outputs are returned as a string that will be interpreted by NOMAD based on the ``BB_OUTPUT_TYPE`` defined by the user.
We recall that constraints must be represented by values :math:`c_j` for a constraint :math:`c_j \leq 0`.

The second argument, the real ``h_max`` (in), corresponds to the current value of the barrier :math:`h_{max}` parameter.
It is not used in this example but it may be used to interrupt an expensive evaluation if the constraint violation value :math:`h` grows larger than :math:`h_{max}`.
See [AuDe09a]_ for the definition of :math:`h` and :math:`h_{max}` and of the *Progressive Barrier* method for handling constraints.

The third argument, ``countEval`` (out), needs to be set to ``true`` if the evaluation counts as a blackbox
evaluation, and ``false`` otherwise (for example, if the user interrupts an evaluation with the :math:`h_{max}`
criterion before it costs some expensive computations, then set ``countEval`` to ``false``).

Finally, note that the call to ``eval_x()`` inside the NOMAD code  is inserted into a ``try`` block.
This means that if an error is detected inside the ``eval_x()`` function,  an exception should be thrown.
The choice for the type of this exception is left to the user, but  ``NOMAD::Exception`` is available.
If an exception is thrown by the user-defined function, then the associated evaluation  is tagged as a failure
and not counted unless the user explicitely set the flag ``countEval`` to ``true``.


Setting parameters
""""""""""""""""""

Once your problem has been defined, the main function can be written. NOMAD routines may throw ``C++`` exceptions,
so it is recommended that you put your code into a ``try`` block.

.. code-block:: c++

  /*------------------------------------------*/
  /*            NOMAD main function           */
  /*------------------------------------------*/
  int main (int argc, char **argv)
  {

      NOMAD::MainStep TheMainStep;

      auto params = std::make_shared<NOMAD::AllParameters>();
      initAllParams(params);
      TheMainStep.setAllParameters(params);

      std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
      TheMainStep.setEvaluator(std::move(ev));

      try
      {
          TheMainStep.start();
          TheMainStep.run();
          TheMainStep.end();
      }

      catch(std::exception &e)
      {
          std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
      }

      return 0;
  }

The execution of NOMAD is controlled by the ``NOMAD::MainStep`` class using the ``start``, ``run`` and ``end`` functions.
The user defined ``NOMAD::Evaluator`` is set into the ``NOMAD::MainStep``.

The base evaluator constructor takes an ``NOMAD::EvalParameters`` as input.
The evaluation parameters are included into a ``NOMAD::AllParameters``.

Hence, in library mode, the main function must declare a ``NOMAD::AllParameters`` object to set all types of parameters.
Parameter names are the same as in batch mode but may be defined programmatically.

A parameter ``PNAME`` is set with the method ``AllParameters::setAttributeValue( "PNAME", PNameValue)``.
The ``PNameValue`` must be of a type registered for the ``PNAME`` parameter.

.. warning:: If the ``PNameValue`` has not the type associated to the ``PName`` parameters, the compilation
   will succeed but execution will be stopped when setting or getting the value.

.. note:: A brief description (including the ``NOMAD::`` type) of all parameters is given :ref:`appendix_parameters`.
   More information on parameters can be obtained by running ``$NOMAD_HOME/bin/nomad -h KEYWORD``.

For the example, the parameters are set in

.. code-block:: c++

  void initAllParams(std::shared_ptr<NOMAD::AllParameters> allParams)
  {
      // Parameters creation
      // Number of variables
      size_t n = 10;
      allParams->setAttributeValue( "DIMENSION", n);
      // The algorithm terminates after
      // this number of black-box evaluations
      allParams->setAttributeValue( "MAX_BB_EVAL", 1000);
      // Starting point
      allParams->setAttributeValue( "X0", NOMAD::Point(n, 7.0) );

      allParams->getPbParams()->setAttributeValue("GRANULARITY", NOMAD::ArrayOfDouble(n, 0.0000001));

      // Constraints and objective
      NOMAD::BBOutputTypeList bbOutputTypes;
      bbOutputTypes.push_back(NOMAD::BBOutputType::PB);     // g1
      bbOutputTypes.push_back(NOMAD::BBOutputType::PB);     // g2
      bbOutputTypes.push_back(NOMAD::BBOutputType::OBJ);    // f
      bbOutputTypes.push_back(NOMAD::BBOutputType::EB);     // c2000
      allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

      allParams->setAttributeValue("DISPLAY_DEGREE", 2);
      allParams->setAttributeValue("DISPLAY_ALL_EVAL", false);
      allParams->setAttributeValue("DISPLAY_UNSUCCESSFUL", false);

      allParams->getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
      allParams->getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);


      // Parameters validation
      allParams->checkAndComply();

  }

The ``checkAndComply`` function must be called to ensure that parameters are compatible.
Otherwise an exception is triggered.

Access to solution and optimization data
""""""""""""""""""""""""""""""""""""""""

**TODO**

.. In the basic example 1, final information is displayed at the end of an algorithm. More specialized access to solution and optimization data is allowed.

.. To access the best feasible and infeasible points, use the methods \sComp{NOMAD::Mads::get\_best}\-\sComp{\_feasible()} and \sComp{NOMAD::Mads::get\_best\_infeasible()}. To access optimization data or statistics, call the method \sComp{NOMAD::Mads::get\_stats()} which returns access to  a \sComp{NOMAD::Stats} object. Then, use the access methods defined in \sComp{Stats.hpp}. For example, to display the number of blackbox evaluations, write:

.. NOMAD::CacheBase::getInstance()->findBestFeas(bf, NOMAD::Point(n), NOMAD::EvalType::BB,NOMAD::ComputeType::STANDARD, nullptr);
.. NOMAD::CacheBase::getInstance()->findBestInf(bi, NOMAD::INF, NOMAD::Point(n), NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD,nullptr);


PyNomad interface compilation
-----------------------------

.. note:: The Python interface requires Python 3.6 and Cython 0.24.

A Python interface for NOMAD is provided for Mac OS X and Linux.
Some examples and source codes are provided in ``$NOMAD_HOME/interfaces/PyNomad``.
To enable the building of the Python interface, option ``-DBUILD_INTERFACES=ON`` must be
set when building NOMAD, as such: ``cmake -DBUILD_TESTS=ON -S . -B build/release``.
The build procedure relies on Python 3.6 and Cython 0.24 or higher.
A simple way to make it work is to first install the `Anaconda <http://www.anaconda.org/>`_ package.
The command ``cmake --install build/release`` must be run before using the PyNomad module.

All functionalities of NOMAD are available in PyNomad.
NOMAD parameters are provided in a list of strings using the same syntax as used in the NOMAD parameter
files.
Several tests and examples are proposed in the ``PyNomad`` directory to check that everything is up and
running.

C interface
-----------

**TODO**

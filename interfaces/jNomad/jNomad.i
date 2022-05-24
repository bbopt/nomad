%module(directors="1") jNomad

%include <std_vector.i>
%include <std_string.i>
%include <std_set.i>
%include <std_list.i>
%include <std_shared_ptr.i>
%include <typemaps.i>
%include <various.i>
%include <enums.swg>

#if SWIG_VERSION >= 0x040000
              %extend std::vector {
                vector(size_type count) { return new std::vector< T >(count); }
              }
#endif

// force the generated Java code to use the C enum values rather than making a JNI call
%javaconst(1);

// convert references bool & to boolean[]
%apply bool & INOUT { bool & count_eval };
// %apply char **STRING_ARRAY { char **argv }



// generate directors for class Evaluator
%feature("director") Evaluator;

%include "../../src/nomad_version.hpp"
%include "../../src/nomad_platform.hpp"

%{
#include "Algos/Step.hpp"
#include "Algos/MainStep.hpp"
#include "Eval/Evaluator.hpp"
#include "Param/AllParameters.hpp"
#include "Math/ArrayOfDouble.hpp"
#include "Math/Double.hpp"
#include "Math/Point.hpp"
#include "Type/BBInputType.hpp"
#include "Type/BBOutputType.hpp"
#include "Type/DirectionType.hpp"
#include "Type/EvalType.hpp"
#include <locale.h>
%}


%shared_ptr(NOMAD::AllParameters)
%shared_ptr(NOMAD::EvalParameters)
%shared_ptr(NOMAD::Evaluator)
namespace NOMAD{

  class DLL_UTIL_API Double {
    public:
      Double ( double v );
      const double & todouble ( void ) const;
  };

  class DLL_UTIL_API ArrayOfDouble {
    public:
      explicit ArrayOfDouble ( const size_t n = 0 , const NOMAD::Double & d = NOMAD::Double() );
      explicit ArrayOfDouble( const std::vector<double> & v) ;
      void set ( size_t j , const NOMAD::Double & v, bool relative = false, const NOMAD::Double & lb = NOMAD::Double(), const NOMAD::Double & ub = NOMAD::Double() );
      NOMAD::Double& operator[](size_t i) const;
  };


  class DLL_UTIL_API Point : public ArrayOfDouble {
    public:
      explicit Point ( const size_t n = 0 , const NOMAD::Double & d = NOMAD::Double() ) : ArrayOfDouble (n, d);
      explicit Point ( const std::vector<double> & v) : ArrayOfDouble (v);
      void set ( size_t j , const NOMAD::Double & v, bool relative = false, const NOMAD::Double & lb = NOMAD::Double(), const NOMAD::Double & ub = NOMAD::Double() );
  };

//  enum class DLL_UTIL_API BBOutputType
//  {
//    OBJ,        ///< Objective value
//    EB,         ///< Extreme barrier constraint
//    PB,         ///< Progressive barrier constraint
//    CNT_EVAL,   ///< Output set to 0 or 1 to count the blackbox evaluation or not
//    //STAT_AVG, ///< Stat (average)
//    //STAT_SUM, ///< Stat (sum)
//    BBO_UNDEFINED ///< Output ignored
//  };
  // %template(NomadBBOutputTypeList) std::list<BBOutputType>;


//  enum class DLL_UTIL_API BBInputType
//  {
//    CONTINUOUS  ,     ///< Continuous variable (default) (R)
//    ALL_CONTINUOUS  , ///< All variables are continuous variable (default, *R). Need a checkAndComply to set the BBInputTypeList.
//    INTEGER     ,     ///< Integer variable (I)
//    ALL_INTEGER ,     ///< All variables are integer (*I). Need a checkAndComply to set the BBInputTypeList.
//    //CATEGORICAL ,   ///< Categorical variable          (C)
//    BINARY         ,  ///< Binary variable               (B)
//    ALL_BINARY       ///< All variables are binary (*B). Need a checkAndComply to set the BBInputTypeList.
//  };
  // %template(NomadBBInputTypeList) std::list<BBInputType>;


  // Direction type
//  enum class DLL_UTIL_API DirectionType
//  {
//    ORTHO_2N,
//    CS,
//    ORTHO_NP1_NEG,
//    ORTHO_NP1_QUAD,
//    NP1_UNI,
//    SINGLE,
//    DOUBLE,
//    LT_2N,
//    LT_1,
//    LT_2,
//    LT_NP1,
//    GPS_2N_STATIC,
//    GPS_2N_RAND,
//    GPS_BINARY,
//    GPS_NP1_STATIC,
//    GPS_NP1_STATIC_UNIFORM,
//    GPS_NP1_RAND,
//    GPS_NP1_RAND_UNIFORM,
//    GPS_1_STATIC,
//    UNDEFINED_DIRECTION
//    ///< DirectionType is mandatory
//  };
  // %template(NomadDirectionTypeList) std::list<DirectionType>;


  // Evaluator type
  enum DLL_EVAL_API class EvalType
  {
      BB,                 ///< The evaluator is a blackbox.
      MODEL,              ///< The evaluator is a model function,
                          /// potentially much faster to run than a blackbox.
      SURROGATE,          ///< The evaluator is a static surrogate,
                          /// potentially much faster to run than a blackbox.
      LAST,               ///< For iterations; Note: UNDEFINED evals are ignored.
      UNDEFINED           ///< Undefined: This value may be used when the
                          ///< EvalType is not mandatory
  };

  /// Enum for the type of Evaluator.
  enum class EvalXDefined
  {
    EVAL_BLOCK_DEFINED_BY_USER, ///< User redefined eval_block() in library mode; Default value
    EVAL_X_DEFINED_BY_USER,     ///< User redefined eval_x() in library mode
    USE_BB_EVAL                 ///< Neither eval_x() nor eval_block() were redefined by library mode. An external executable is provided.
  };

  class DLL_EVAL_API EvalPoint : public NOMAD::Point {
    public:

     explicit EvalPoint() : Point();

      void setBBO(const std::string &bbo,
                  const std::string &sBBOutputTypes = "",
                  NOMAD::EvalType evalType = NOMAD::EvalType::BB,
                  const bool evalOk = true);

  };
  
  class DLL_UTIL_API EvalParameters : public NOMAD::Parameters {
    public:

      explicit EvalParameters () : Parameters() ;
  };

  class DLL_EVAL_API Evaluator {
    public:
      Evaluator ( const std::shared_ptr<NOMAD::EvalParameters> & p,
                  NOMAD::EvalType evalType = EvalType::BB,
                  NOMAD::EvalXDefined evalXDefined = EvalXDefined::EVAL_BLOCK_DEFINED_BY_USER );

      virtual ~Evaluator();
      virtual bool eval_x ( NOMAD::EvalPoint & x , const NOMAD::Double & h_max , bool & count_eval ) const override {};
      // virtual bool eval_x ( std::list<NOMAD::EvalPoint *> &x , const NOMAD::Double & h_max, std::list<bool> & count_eval ) {};
  };

  class DLL_UTIL_API AllParameters  {
    public:

      explicit AllParameters ();

      void readParamLine(const std::string& line);

      template<typename T> void setAttributeValue(std::string name, T value);
      %template(setAttributeValueSizeT) setAttributeValue<size_t>;
      %template(setAttributeValueAOD) setAttributeValue<ArrayOfDouble>;
      %template(setAttributeValuePoint) setAttributeValue<Point>;
      %template(setAttributeValueBool) setAttributeValue<bool>;
      %template(setAttributeValueInt) setAttributeValue<int>;
      // %template(setAttributeValueBBOutputTypeList) setAttributeValue<NomadBBOutputTypeList>;
      // %template(setAttributeValueBBInputTypeList) setAttributeValue<NomadBBInputTypeList>;
      //%template(setAttributeValueDirectionTypeList) setAttributeValue<NomadDirectionTypeList>;

      const std::shared_ptr<NOMAD::EvalParameters>& getEvalParams() const;

      std::string getSetAttributeAsString() const;

      void checkAndComply();

  };

  class DLL_ALGO_API Step
  {
    public:
       explicit Step();

       void start();
       bool run();
       void end();

     protected:
      virtual void startImp() = 0 ;
   };

  class DLL_ALGO_API MainStep : public Step
  {
    public:
      explicit MainStep() ;

      void setAllParameters(const std::shared_ptr<AllParameters> & allParams);

      void setEvaluator(std::shared_ptr<Evaluator> ev);

      void displayHelp(const std::string& helpSubject = "all", bool devHelp = false);

      static void resetComponentsBetweenOptimization();

    protected:
      virtual void startImp() override;

  };


}


%extend NOMAD::Evaluator {
  std::vector<bool> eval_block(std::vector<std::shared_ptr<NOMAD::EvalPoint>> &block,
                         const NOMAD::Double& hMax,
                         std::vector<bool> &countEval) const override
  {
    std::vector<bool> evalOk(block.size(), false);
    countEval.resize(block.size(), false);

    for (size_t index = 0; index < block.size(); index++)
    {
      bool countEval1 = false;
      evalOk[index] = (*($self)).eval_x(*block[index], hMax, countEval1);
      countEval[index] = countEval1;
      }

      return evalOk;
  }
}

%extend NOMAD::EvalPoint {
    NOMAD::Double & get(size_t i) {
        return (*($self))[i];
    }
}

%extend NOMAD::MainStep {
    void init() {
       setlocale(LC_ALL,"C");
       NOMAD::MainStep::resetComponentsBetweenOptimization();
    }
}

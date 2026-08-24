%module(directors="1") jNomad

%include <std_vector.i>
%include <std_string.i>
%include <std_set.i>
%include <std_list.i>
%include <std_shared_ptr.i>
%include <typemaps.i>
%include <various.i>
%include <enums.swg>

// force the generated Java code to use the C enum values rather than making a JNI call
%javaconst(1);

// convert references bool & to boolean[]
%apply bool & INOUT { bool & count_eval };
%apply bool & INOUT { bool & rel };
%apply bool & INOUT { bool & stop };

// generate directors for class Evaluator
%feature("director") Evaluator;
%feature("director") EvalStopCheckCallback;

%{
#include "Algos/Step.hpp"
#include "Algos/MainStep.hpp"
#include "Eval/Evaluator.hpp"
#include "Eval/EvalQueuePoint.hpp"
#include "Param/Parameters.hpp"
#include "Param/AllParameters.hpp"
#include "Param/EvalParameters.hpp"
#include "Math/Double.hpp"
#include "Math/ArrayOfDouble.hpp"
#include "Math/Point.hpp"
#include "Type/BBInputType.hpp"
#include "Type/BBOutputType.hpp"
#include "Type/DirectionType.hpp"
#include "Type/EvalType.hpp"
#include "Type/EvalSortType.hpp"
#include "Cache/CacheBase.hpp"
#include "Util/ArrayOfString.hpp"
#include "nomad_platform.hpp"
#include <locale.h>
%}

%include "../../src/nomad_version.hpp"
%include "../../src/nomad_platform.hpp"

%shared_ptr(NOMAD::Parameters)
%shared_ptr(NOMAD::AllParameters)
%shared_ptr(NOMAD::EvalParameters)
%shared_ptr(NOMAD::Evaluator)
%shared_ptr(NOMAD::EvalStopCheckCallback)

%rename("equals") *::operator==;
%rename("notEquals") *::operator!=;
%rename("set") *::operator=;
%rename("add") *::operator+;
%rename("sub") *::operator-;
%rename("mul") *::operator*;
%rename("div") *::operator/;
%rename("incr") *::operator++;
%rename("decr") *::operator--;

%ignore *::operator<<;
%ignore *::operator>>;
%ignore *::operator-=;
%ignore *::operator<;
%ignore *::operator>;
%ignore *::operator<=;
%ignore *::operator>=;
%ignore *::operator[];
%ignore *::operator*=;
%ignore *::operator/=;

namespace NOMAD{
  %ignore Double::NotDefined;
  %ignore Double::InvalidValue;
  %ignore Double::operator+=;
  %ignore Double::operator-=;
  %ignore Double::operator-=;
  %include "../../src/Math/Double.hpp"

  %include "../../src/Util/ArrayOfString.hpp"

  %ignore ArrayOfDouble::set(size_t, Double const *); // duplicated translation.
  %include "../../src/Math/ArrayOfDouble.hpp"

  %ignore Point::operator<;
  %ignore Point::operator+;
  %ignore Point::vectorize; // needs Direction. Ignored to avoid too big wrapper
  %include "../../src/Math/Point.hpp"
  %include "../../src/Type/BBOutputType.hpp"
  %include "../../src/Type/BBInputType.hpp"
  %include "../../src/Type/DirectionType.hpp"
  %include "../../src/Type/EvalType.hpp"
  %include "../../src/Type/EvalSortType.hpp"

  // Class required for creating a fully functional EvalQueuePointPtr.java class
  class EvalQueuePointPtr{
    public:
      NOMAD::EvalQueuePoint* get();
  }; 
  typedef std::shared_ptr<EvalQueuePoint> EvalQueuePointPtr;

  class DLL_UTIL_API Parameters {
    protected:
      explicit Parameters() : _typeName("Unknown"), _toBeChecked(true),  _attributes();
  };

  class DLL_UTIL_API AllParameters  {
    public:
      explicit AllParameters ();

      void readParamLine(const std::string& line);

      template<typename T> void setAttributeValue(std::string name, T value);
      %template(setAttributeValueSizeT) setAttributeValue<size_t>;
      %template(setAttributeValueAOD) setAttributeValue<ArrayOfDouble>;
      %template(setAttributeValueAOS) setAttributeValue<ArrayOfString>;
      %template(setAttributeValuePoint) setAttributeValue<Point>;
      %template(setAttributeValueBool) setAttributeValue<bool>;
      %template(setAttributeValueInt) setAttributeValue<int>;
      %template(setAttributeValueDT) setAttributeValue<DirectionType>;
      %template(setAttributeValueEST) setAttributeValue<EvalSortType>;
      %template(setAttributeValueBBOutputTypeList) setAttributeValue<NOMAD::BBOutputTypeList>;
      %template(setAttributeValueBBInputTypeList) setAttributeValue<NOMAD::BBInputTypeList>;

      const std::shared_ptr<NOMAD::EvalParameters>& getEvalParams() const;

      std::string getSetAttributeAsString() const;

      void checkAndComply();
  };

  class DLL_UTIL_API EvalParameters : public NOMAD::Parameters {
    public:
      explicit EvalParameters () : Parameters() ;
  };

  class DLL_EVAL_API Eval {
    public:
      explicit Eval();
  };

  class DLL_EVAL_API EvalPoint : public NOMAD::Point {
    public:
      explicit EvalPoint() : Point();

      void setBBO(const std::string &bbo,
                  const std::string &sBBOutputTypes = "",
                  NOMAD::EvalType evalType = NOMAD::EvalType::BB,
                  const bool evalOk = true);

      Eval* getEval(EvalType evalType) const;

      const Point* getX() const;
    };

    class EvalQueuePoint : public EvalPoint{
      public:
        explicit EvalQueuePoint(const EvalPoint& evalPoint, EvalType evalType): EvalPoint(evalPoint);
        const EvalType& getEvalType() const;
        const NOMAD::SuccessType& getSuccess() const;
    };

  // Enum for the type of Evaluator.
  enum class EvalXDefined
  {
    EVAL_BLOCK_DEFINED_BY_USER, ///< User redefined eval_block() in library mode; Default value
    EVAL_X_DEFINED_BY_USER,     ///< User redefined eval_x() in library mode
    USE_BB_EVAL                 ///< Neither eval_x() nor eval_block() were redefined by library mode. An external executable is provided.
  };

  class DLL_EVAL_API Evaluator {
    public:
      Evaluator ( const std::shared_ptr<NOMAD::EvalParameters> & p,
                  NOMAD::EvalType evalType,
                  NOMAD::EvalXDefined evalXDefined = EvalXDefined::EVAL_BLOCK_DEFINED_BY_USER );

      virtual ~Evaluator();
      virtual bool eval_x ( NOMAD::EvalPoint & x , const NOMAD::Double & h_max , bool & count_eval ) const override {};
      // virtual bool eval_x ( std::list<NOMAD::EvalPoint *> &x , const NOMAD::Double & h_max, std::list<bool> & count_eval ) {};
  };

  class DLL_EVAL_API EvalStopCheckCallback
   {
    public:
      EvalStopCheckCallback() = default;
      virtual ~EvalStopCheckCallback() = default;

      EvalStopCheckCallback(const EvalStopCheckCallback&) = delete;
      EvalStopCheckCallback& operator=(const EvalStopCheckCallback&) = delete;

      virtual void call(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool& stop) const = 0;
  };

  // Success type of a step.
  enum class SuccessType {
    UNDEFINED,          ///< Default type set at start
    NO_TRIALS,          ///< No trial points produced
    UNSUCCESSFUL,       ///< Trial point is not a success
    PARTIAL_SUCCESS,    ///< Partial success (improving). Found an infeasible
    ///< solution with a better h. f is worse.
    FULL_SUCCESS        ///< Full success (dominating)
  };

  class DLL_ALGO_API Step {
    public:
       explicit Step();

       void start();
       bool run();
       void end();

     protected:
      virtual void startImp() = 0 ;
   };

  class DLL_ALGO_API MainStep : public Step {
    public:
      explicit MainStep() ;

      void setAllParameters(const std::shared_ptr<AllParameters> & allParams);

      void addEvaluator(std::shared_ptr<Evaluator> ev);

      void addEvalStopCheckCallback(std::shared_ptr<EvalStopCheckCallback> cbEval);

      void displayHelp(const std::string& helpSubject = "all", bool devHelp = false);

      static void resetComponentsBetweenOptimization();

    protected:
      virtual void startImp() override;
  };

  class DLL_ALGO_API CacheBase {
    private:
      CacheBase();
      ~CacheBase();
  };

  DLL_UTIL_API EvalSortType stringToEvalSortType(const std::string &s);

  DLL_UTIL_API std::string evalSortTypeToString(const EvalSortType& evalSortType);

  // End of NOMAD namespace
}

%extend NOMAD::Evaluator {
  std::vector<bool> eval_block(std::vector<std::shared_ptr<NOMAD::EvalPoint>> &block,
                         const NOMAD::Double& hMax,
                         std::vector<bool> &countEval) const override{
    std::vector<bool> evalOk(block.size(), false);
    countEval.resize(block.size(), false);

    for (size_t index = 0; index < block.size(); index++) {
      bool countEval1 = false;
      evalOk[index] = (*($self)).eval_x(*block[index], hMax, countEval1);
      countEval[index] = countEval1;
    }
    return evalOk;
  }

}

%extend NOMAD::Eval {
  Double getObjective() const  {
    return (*($self)).getBBOutput().getObjective((*($self)).getBBOutputTypeList());
  }

  ArrayOfDouble getObjectives() const  {
    return (*($self)).getBBOutput().getObjectives((*($self)).getBBOutputTypeList());
  }

  ArrayOfDouble getConstraints() const  {
    return (*($self)).getBBOutput().getConstraints((*($self)).getBBOutputTypeList());
  }
}

%extend NOMAD::ArrayOfDouble {
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

%extend NOMAD::CacheBase {
  static const CacheBase& getInstance() {
    return *NOMAD::CacheBase::getInstance();
  }

  size_t findBestFeas(std::vector<EvalPoint> &evalPointList) const {
    return $self->findBestFeas(evalPointList);
  }

  size_t findBestInf(std::vector<EvalPoint> &evalPointList) const {
    return $self->findBestInf(evalPointList);
  }
}

namespace std {
  %template(EvalPointVector) vector<NOMAD::EvalPoint>;
  %template(BBOutputTypeList) vector<NOMAD::BBOutputType>;
  %template(BBInputTypeList) vector<NOMAD::BBInputType>;
  %template(DirectionTypeList) vector<NOMAD::DirectionType>;
}

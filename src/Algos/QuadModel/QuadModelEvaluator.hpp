#ifndef __NOMAD400_QUAD_MODEL_EVALUATION__
#define __NOMAD400_QUAD_MODEL_EVALUATION__

#include "../../Eval/Evaluator.hpp"
#include "../../Output/OutputInfo.hpp"
#include "../../../ext/sgtelib/src/Surrogate.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for evaluating trial points as EvalType::SGTE.
class QuadModelEvaluator : public Evaluator
{
private:
    const std::shared_ptr<SGTELIB::Surrogate> _model;
    std::string                _modelDisplay;
    OutputLevel                _displayLevel;
    Point                      _fixedVariable;  ///< Points are sent to evaluator in full space. Evaluator works in its own dimension. This member is used for conversions.


public:
    /// Constructor
    /**
     Usually, Evaluators work in full dimension. In this case, the models may work better in subdimension. This is why a fixed variable is used.
     */
    explicit QuadModelEvaluator(const std::shared_ptr<EvalParameters>& evalParams,
                                const std::shared_ptr<SGTELIB::Surrogate>& model,
                                const std::string& modelDisplay,
                                const Point& fixedVariable)
      : Evaluator(evalParams, EvalType::SGTE),
        _model(model),
        _modelDisplay(modelDisplay),
        _displayLevel(OutputLevel::LEVEL_INFO),
        _fixedVariable(fixedVariable)
    {
        init();
    }

    virtual ~QuadModelEvaluator();

    /**
     Points for evaluations are given in a block. Sgtelib models handle the points as a matrix and return a matrix for outputs.
     */
    std::vector<bool> eval_block(Block &block,
                                 const Double &hMax  __attribute__((unused)),
                                 std::vector<bool> &countEval) const override;

    static void evalH(const ArrayOfDouble& bbo,
                      const BBOutputTypeList& bbot,
                      Double &h);


private:
    void init();


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_EVALUATION__

#ifndef __NOMAD400_EVALUATORCONTROLGLOBALPARAMETERS__
#define __NOMAD400_EVALUATORCONTROLGLOBALPARAMETERS__

#include "../Param/Parameters.hpp"

#include "../nomad_nsbegin.hpp"

/// The class for global EvaluatorControl parameters.
/**
- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class EvaluatorControlGlobalParameters final : public Parameters
{
public:

    explicit EvaluatorControlGlobalParameters()
      : Parameters()
    {
        init();
    }

    /// Check the sanity of parameters.
    void checkAndComply( );

private:

    /// Helper for constructor
    /**
     Register and set default values for all evaluator control attributes. The information to register all the attributes is contained in evaluatorControlAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable, called automatically by makefile when the evaluatorControlAttributeDefinition.txt file is modified.
     */
    void init() override;

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVALUATORCONTROLGLOBALPARAMETERS__


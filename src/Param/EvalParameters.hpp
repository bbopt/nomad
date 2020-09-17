#ifndef __NOMAD400_EVALPARAMETERS__
#define __NOMAD400_EVALPARAMETERS__

#include "../Param/Parameters.hpp"
#include "../Param/RunParameters.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for Evaluator parameters
/**
- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class EvalParameters final : public Parameters
{
public:

    explicit EvalParameters()
    : Parameters()
    {
        init();
    }

    /// Check the sanity of parameters.
    void checkAndComply( const std::shared_ptr<RunParameters> & runParams );

private:
    /// Helper for constructor
    /**
     Register and set default values for all evaluation attributes. The information to register all the attributes is contained in evalAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable, called automatically by makefile when the evalAttributeDefinition.txt file is modified.
     */
    void init() override ;

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_EVALPARAMETERS__


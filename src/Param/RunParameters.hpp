#ifndef __NOMAD400_RUNPARAMETERS__
#define __NOMAD400_RUNPARAMETERS__

#include "../Param/EvaluatorControlGlobalParameters.hpp"
#include "../Param/EvaluatorControlParameters.hpp"
#include "../Param/Parameters.hpp"
#include "../Param/PbParameters.hpp"

#include "../nomad_nsbegin.hpp"

/// The class for the parameters defining the type of optimization/task to perform.
/**
The RunParameter are used by other parameters to update their value during sanity check.

- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class RunParameters final : public Parameters
{
private:
    static bool _warningUnknownParamShown;

public:
    // Constructor: init() will be called.
    // This will register and set default values to all attributes.
    explicit RunParameters()
    : Parameters()
    {
        init();
    }

    /// The copy constructor is not implemented in the parent class
    RunParameters& operator=(const RunParameters& params) { copyParameters(params) ; return *this; }
    RunParameters(const RunParameters& params) : RunParameters() { copyParameters(params); }

    /// Check the sanity of parameters.
    /**
     Register and set default values for all run attributes. The information to register all the attributes is contained in runAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable, called automatically by makefile when the runAttributeDefinition.txt file is modified.
    */
    void checkAndComply(const std::shared_ptr<EvaluatorControlGlobalParameters>& evaluatorControlGlobalParams,
                        const std::shared_ptr<PbParameters>& pbParams);

private:
    /// Initialization
    /**
     * This will register and set default values to all attributes.
     */
    void init() override ;

    /// Helper for checkAndComply()
    void setStaticParameters();

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_RUNPARAMETERS__

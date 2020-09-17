#ifndef __NOMAD400_PBPARAMETERS__
#define __NOMAD400_PBPARAMETERS__

#include "../Param/Parameters.hpp"

#include "../nomad_nsbegin.hpp"

/// The class for the parameters defining the optimization problem.
/**
- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class PbParameters final : public Parameters
{
private:
    bool _showWarningMeshSizeRedefined;

public:

    explicit PbParameters()
      : Parameters(),
        _showWarningMeshSizeRedefined(true)
    {
        init();
    }

    /**
     The copy constructor is not implemented in the parent class to allow some control over what parameters can be copied or not. Use the deep copy function of parameters: Parameters::copyParameters.
     */
    PbParameters& operator=(const PbParameters& params) { copyParameters(params) ; return *this; }

    /**
     The copy constructor is not implemented in the parent class to allow some control over what parameters can be copied or not. Use the deep copy function of parameters: Parameters::copyParameters.
     */
    PbParameters(const PbParameters& params) : PbParameters() { copyParameters(params); }

    /// Check the sanity of parameters.
    void checkAndComply( );

    /// Do not show certain warnings
    void doNotShowWarnings() { _showWarningMeshSizeRedefined = false; }

private:
    /// Helper for constructor
    /**
     Register and set default values for all problem attributes. The information to register all the attributes is contained in pbAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable, called automatically by makefile when the pbAttributeDefinition.txt file is modified.
     */
    void init() override ;

    ///  Helper for checkAndComply()
    void setGranularityAndBBInputType();
    ///  Helper for checkAndComply()
    void setVariableGroups();
    ///  Helper for checkAndComply()
    void setFixedVariables();
    ///  Helper for checkAndComply()
    void setMinMeshParameters(const std::string &paramName);
    ///  Helper for checkAndComply()
    void setInitialMeshParameters();
    ///  Helper for checkAndComply()
    void checkX0ForGranularity() const;
    ///  Helper for checkAndComply()
    void checkX0AgainstBounds() const;
    ///  Helper for checkAndComply()
    void checkForGranularity(const std::string &paramName) const;
    ///  Helper for checkAndComply()
    void checkForGranularity(const std::string &paramName,
                             const ArrayOfDouble &arrayToCheck) const;

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_PBPARAMETERS__


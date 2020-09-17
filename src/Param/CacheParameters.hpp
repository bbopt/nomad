#ifndef __NOMAD400_CACHEPARAMETERS__
#define __NOMAD400_CACHEPARAMETERS__

#include "../Param/Parameters.hpp"
#include "../Param/RunParameters.hpp"

#include "../nomad_nsbegin.hpp"

/// The class for all Cache parameters.
/**
- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class CacheParameters final : public Parameters
{
public:
    /// Constructor
    explicit CacheParameters()
    : Parameters()
    {
        init();
    }

    /// Check the sanity of parameters.
    /**
     RunParameters is needed to obtain the value of PROBLEM_DIR parameter.
     */
    void checkAndComply( std::shared_ptr<RunParameters> runParams );

private:
    /// Helper for constructor
    /**
     Register and set default values for all cache attributes. The information to register all the attributes is contained in cacheAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable from the cacheAttributeDefinition.txt.
     */
    void init() override ;

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_CACHEPARAMETERS__


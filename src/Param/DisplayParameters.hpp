#ifndef __NOMAD400_DISPLAYPARAMETERS__
#define __NOMAD400_DISPLAYPARAMETERS__


#include "../Param/Parameters.hpp"
#include "../Param/RunParameters.hpp"
#include "../Param/PbParameters.hpp"


#include "../nomad_nsbegin.hpp"

/// The class for Display parameters
/**
- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class DisplayParameters final : public Parameters
{
public:

    explicit DisplayParameters()
    : Parameters()
    {
        init();
    }


    /// Check the sanity of parameters.
    void checkAndComply( const std::shared_ptr<RunParameters> & runParams ,
                        const std::shared_ptr<PbParameters> & pbParams );


private:

    /// Helper for constructor
    /**
     Register and set default values for all display attributes. The information to register all the attributes is contained in displayAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable, called automatically by makefile when the displayAttributeDefinition.txt file is modified.
     */
    void init() override ;

    /// Helper for checkAndComply
    ArrayOfDouble setFormatFromGranularity(const ArrayOfDouble & aod );

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_DISPLAYPARAMETERS__


#ifndef __NOMAD400_SGTELIB_MODEL_INITIALIZATION__
#define __NOMAD400_SGTELIB_MODEL_INITIALIZATION__

#include "../../Algos/Initialization.hpp"

#include "../../nomad_nsbegin.hpp"


class SgtelibModelInitialization: public Initialization
{
public:
    /// Constructor
    explicit SgtelibModelInitialization(const Step* parentStep)
      : Initialization(parentStep)
    {
        init();
    }

    virtual ~SgtelibModelInitialization();


private:
    void init();

    // Step methods
    void startImp() override;
    bool runImp() override;
    void endImp() override;

    void validateX0s() const;
    bool eval_x0s();
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIB_MODEL_INITIALIZATION__


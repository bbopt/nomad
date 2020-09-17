#ifndef __NOMAD400_MADSINITIALIZATION__
#define __NOMAD400_MADSINITIALIZATION__

#include "../../Algos/Initialization.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for Mads initialization (step 0)
/**
 The run function of this step validates and evaluates X0(s).
 */
class MadsInitialization final: public Initialization
{
private:
    std::shared_ptr<MeshBase> _initialMesh;

public:
    /// Constructor
    /*
     \param parentStep      The parent of this step -- \b IN.
     */
    explicit MadsInitialization(const Step* parentStep)
      : Initialization(parentStep),
        _initialMesh(nullptr)
    {
        init();
    }

    virtual ~MadsInitialization() {}

    std::shared_ptr<MeshBase> getMesh() const { return _initialMesh; }

private:
    void init();

    void validateX0s() const;
    bool eval_x0s();

    virtual bool runImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_MADSINITIALIZATION__

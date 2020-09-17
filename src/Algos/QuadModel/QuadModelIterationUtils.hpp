
#ifndef __NOMAD400_QUADMODELITERATIONUTILS__
#define __NOMAD400_QUADMODELITERATIONUTILS__

#include "../../Algos/IterationUtils.hpp"
#include "../../../ext/sgtelib/src/Surrogate.hpp"
#include "../../../ext/sgtelib/src/TrainingSet.hpp"

#include "../../nomad_nsbegin.hpp"


/// Class of utils for QuadModel iterations.
class QuadModelIterationUtils : public IterationUtils
{
private:

    void init();

protected:
    std::shared_ptr<SGTELIB::TrainingSet>   _trainingSet;
    std::shared_ptr<SGTELIB::Surrogate>     _model;

public:
    /// Constructor
    /**
     The model and training set are obtained from QuadModelIteration.

     \param parentStep      The calling iteration Step.
     */
    explicit QuadModelIterationUtils(const Step* parentStep)
      : IterationUtils(parentStep),
        _trainingSet(nullptr),
        _model(nullptr)
    {
        init();
    }

    void displayModelInfo() const;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUADMODELITERATIONUTILS__

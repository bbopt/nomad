#ifndef __NOMAD400_SPECULATIVESEARCHMETHOD__
#define __NOMAD400_SPECULATIVESEARCHMETHOD__

#include "../../Algos/Mads/SearchMethodSimple.hpp"

#include "../../nomad_nsbegin.hpp"


/// Speculative search.
/**
 The speculative search consists in looking further away along
 the successful direction after an improvement.
 */
class SpeculativeSearchMethod  final : public SearchMethodSimple
{
public:
    /// Constructor
    /**
     \param parentStep      The parent of this search step -- \b IN.
     */
    explicit SpeculativeSearchMethod(const Step* parentStep )
    : SearchMethodSimple( parentStep )
    {
        init();
    }

private:
    void init();

    /// Generate new points to evaluate
    /**
     \copydoc SearchMethodSimple::generateTrialPointsImp \n
     The speculative search generates points in the direction of success.
     */
    void generateTrialPointsImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SPECULATIVESEARCHMETHOD__

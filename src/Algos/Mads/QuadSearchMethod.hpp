#ifndef __NOMAD400_QUADSEARCHMETHOD__
#define __NOMAD400_QUADSEARCHMETHOD__

#include "../../Algos/Mads/SearchMethodSimple.hpp"
#ifdef USE_SGTELIB
#include "../../Algos/QuadModel/QuadModelSinglePass.hpp"
#endif

#include "../../nomad_nsbegin.hpp"

/// Class to perform a Search method using the quadratic model optimization algorithm.
/**
 If Quad Search is enabled (check is done in QuadSearchMethod::init), the QuadSearchMethod::run function manages the execution (start, run, end) of the algorithm based on Quad Model. \n
 The new trial points can be generated during a single pass of Quad model construction and optimization (generateTrialPoints) or as repeated Quad model optimizations (multiple constructions, optimizations and evaluations -> run).
 */
class QuadSearchMethod final : public SearchMethodSimple
{
private:
    OutputLevel _displayLevel;

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit QuadSearchMethod(const Step* parentStep)
      : SearchMethodSimple(parentStep),
        _displayLevel(OutputLevel::LEVEL_NORMAL)
    {
        init();
    }

private:

    /// Helper for constructor.
    /**
     Test if the Quad Search is enabled or not. Test if the Sgtelib library has been linked. Manage displays.
     */
    void init();

    ///Generate new points (no evaluation)
    /**
     \copydoc SearchMethodSimple::generateTrialPointsImp \n
     A quadratic model (built with Sgtelib) is constructed around the best feasible and around the best infeasible point. For each model a  single sub-optimization is performed to obtain a best feasible and a best infeasible point. At most 4 trial points can be generated.
     */
    void generateTrialPointsImp() override ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUADSEARCHMETHOD__


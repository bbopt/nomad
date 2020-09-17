#ifndef __NOMAD400_SEARCHMETHODALGO__
#define __NOMAD400_SEARCHMETHODALGO__

#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for generic search method of MADS. Run by Search.
class SearchMethodAlgo: public SearchMethodBase
{


public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SearchMethodAlgo( const Step* parentStep )
    : SearchMethodBase( parentStep ) {}

    /**
         An empty (disabled) startImp is required for a search method that launches an iterative algorithm during the run.

     */
    void startImp() override {}

    /**
     - Pure virtual function.
     - This function must be implemented for algo based search methods that can perform several iterations.
     - The derived runImp implementation of this function launches the sequence of start, run and end on an algorithm.
     - When running the algorithm, evaluations must be performed.
     - Calling this function implies that ::generateTrialPointsImp is not called.
     - This function is used only when the option to generate all points before evaluation is disabled, that is the ::generateTrialPointsImp is not called.
     */
    virtual bool runImp() override = 0 ;

    /**
     - Pure virtual function.
     - This function must be implemented for algo based search methods that can perform a single iteration for generating points.
     - This function is used only when the option to generate all points before evaluation is enabled, that is the runImp is not called.
     - Evaluations are automatically performed when running SearchMethodSimple::runImp.
     */
    virtual void generateTrialPointsImp() override = 0 ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SEARCHMETHODALGO__


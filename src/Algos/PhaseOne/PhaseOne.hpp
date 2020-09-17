#ifndef __NOMAD400_PHASE_ONE__
#define __NOMAD400_PHASE_ONE__

#include "../../Eval/EvalPoint.hpp"
#include "../../Algos/Algorithm.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/Mads.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for phase one search of MADS to satisfy Extreme Barrier (EB) constraints.
/**
 This algorithm is called as a SearchMethod for Mads that modifies the problem to minimize the infeasibility of EB constraints. This is done when initial point has infeasible EB constraint. Then a sub-optimization using Mads solves the modified problem. Once completed, the parent Mads continues on the regular problem. Points for which an EB constraint is infeasible have their objective set to infinity.
 */
class PhaseOne: public Algorithm
{
private:

    std::shared_ptr<Mads>    _mads;
    std::shared_ptr<AlgoStopReasons<MadsStopType>>    _madsStopReasons;


    /**
      The list of ::BBOutputType parameters used for this Phase One.
      Used to recompute h values at the end of Phase One.
      Since the recompute methods are static, this member
      needs to be static.
     */
    static BBOutputTypeList _bboutputtypes;

public:
    /// Constructor
    /**
     \param parentStep    The parent of this step -- \b IN.
     \param stopReasons   The Phase One stop reasons -- \b IN/OUT.
     \param runParams     Parameters for algorithm -- \b IN.
     \param refPbParams   Parameters for original optimization problem. Phase One use its own copy -- \b IN.
     */
    explicit PhaseOne(const Step* parentStep,
                      std::shared_ptr<AlgoStopReasons<PhaseOneStopType>> stopReasons,
                      const std::shared_ptr<RunParameters>& runParams,
                      const std::shared_ptr<PbParameters>& refPbParams)
      : Algorithm(parentStep, stopReasons, runParams, std::make_shared<PbParameters>(*refPbParams)),
        _mads(nullptr)
    {
        init();
    }
    virtual ~PhaseOne() {}

    static void setBBOutputTypes(const BBOutputTypeList& bboutputtypes) { _bboutputtypes = bboutputtypes; }

    /**
     - Setup EvalPoint success computation to be based on h rather than f.
     - Recompute points in cache.
     - Setup stop if feasible criterion.
     - Setup Mads

     */
    virtual void    startImp() override;
    virtual bool    runImp()   override;
    virtual void    endImp()   override;

    virtual void readInformationForHotRestart() override;

private:
    /// Helper for constructor
    void init();

    /*------------------------*/
    /* Private helper methods */
    /*------------------------*/
    /**
     Static function called by Cache::processOnAllPoints().
     */
    static void recomputeH(EvalPoint& evalPoint);

    /**
     Static function called by Cache::processOnAllPoints().
     */
    static void recomputeHPB(EvalPoint& evalPoint);

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_PHASE_ONE__

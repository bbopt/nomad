#ifndef __NOMAD400_SSDMADS__
#define __NOMAD400_SSDMADS__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Algorithm.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for SSD-Mads
class SSDMads: public Algorithm
{

public:
    /// Constructor
    /**
     \param parentStep    The parent of this step -- \b IN.
     \param stopReasons   The SSD Mads stop reasons -- \b IN/OUT.
     \param runParams     Parameters for algorithm -- \b IN.
     \param refPbParams   Parameters for original optimization problem. SSD-Mads use its own copy -- \b IN.
     */
    explicit SSDMads(const Step* parentStep,
                     std::shared_ptr<AlgoStopReasons<MadsStopType>> stopReasons,
                     const std::shared_ptr<RunParameters>& runParams,
                     const std::shared_ptr<PbParameters>& refPbParams)
      : Algorithm(parentStep, stopReasons, runParams, std::make_shared<PbParameters>(*refPbParams))
    {
        init();
    }
    virtual ~SSDMads() {}

    virtual bool runImp()   override;

    virtual void readInformationForHotRestart() override;

private:
    /// Helper for constructor
    void init();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SSDMADS__

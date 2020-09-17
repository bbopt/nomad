#ifndef __NOMAD400_QUAD_MODEL_ALGO__
#define __NOMAD400_QUAD_MODEL_ALGO__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Algorithm.hpp"
#include "../../Algos/EvcInterface.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for implementation of quadratic model optimization algorithm using Bastien Talgorn's sgtelib.
/**
 * Use the start, run and end tasks. Iterate on the following sequence:
 *
 * 1- Points provided as X0s and points in cache are put in a training set.
 * 2- These points are used to build a surrogate model.
 * 3- The model is optimized. This gives oracle points.
 * 4- The oracle points are evaluated by the blackbox.
 * 5- As long as new oracle points are found, the process is repeated.
 *
 * When used by Mads SearchMethod (QuadSearchMethod):
 * - Steps 1, 2, 3 and 4 are the same.
 * - The oracle points are send back to QuadSearchMethod, which takes care
 *   of projecting them to mesh and evaluate them.
 *
 * Training set and model are stored here to allow access to other Quad classes.
 *
 */

class QuadModelAlgo: public Algorithm
{
public:
    /// Constructor
    explicit QuadModelAlgo(const Step* parentStep,
                           std::shared_ptr<AlgoStopReasons<ModelStopType>> stopReasons,
                           const std::shared_ptr<RunParameters>& runParams,
                           const std::shared_ptr<PbParameters>& pbParams)
      : Algorithm(parentStep, stopReasons, runParams, pbParams)
    {
        init();
    }

    virtual ~QuadModelAlgo();

    // Utility function to get BB_OUTPUT_TYPE parameter, which is buried in Evaluator.
    static BBOutputTypeList getBBOutputType()
    {
        if (nullptr == EvcInterface::getEvaluatorControl()
            || nullptr == EvcInterface::getEvaluatorControl()->getEvalParams())
        {
            throw Exception(__FILE__, __LINE__, "Error in QuadModel::getBBOutputType()");
        }
        return EvcInterface::getEvaluatorControl()->getEvalParams()->getAttributeValue<BBOutputTypeList>("BB_OUTPUT_TYPE");
    }

    void readInformationForHotRestart() override {}

private:
    void init();

    void startImp() override;
    bool runImp() override;
    void endImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_ALGO__


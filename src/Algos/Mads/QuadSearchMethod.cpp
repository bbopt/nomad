
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/QuadSearchMethod.hpp"
#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Output/OutputQueue.hpp"

//
// Reference: File Sgtelib_Model_Search.cpp in NOMAD 3.9.1
// Author: Bastien Talgorn

void NOMAD::QuadSearchMethod::init()
{
    setName("Quad Model Search Method");
    //setComment("(QuadSearch)");
    verifyParentNotNull();

    const auto parentSearch = getParentStep()->getParentOfType<NOMAD::QuadSearchMethod*>(false);

    setEnabled((nullptr == parentSearch) && _runParams->getAttributeValue<bool>("QUAD_MODEL_SEARCH"));
#ifndef USE_SGTELIB
    if (isEnabled())
    {
        OUTPUT_INFO_START
        AddOutputInfo(_name + " cannot be performed because NOMAD is compiled without sgtelib library");
        OUTPUT_INFO_END
        setEnabled(false);
    }
#endif

#ifdef USE_SGTELIB
    // Check that there is exactly one objective
    if (isEnabled())
    {
        const auto bbot = NOMAD::QuadModelAlgo::getBBOutputType();
        auto nbObj = NOMAD::getNbObj(bbot);
        if (0 == nbObj)
        {
            OUTPUT_INFO_START
            AddOutputInfo(_name + " not performed when there is no objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }
        else if (nbObj > 1)
        {
            OUTPUT_INFO_START
            AddOutputInfo(_name + " not performed on multi-objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }

        auto modelDisplay = _runParams->getAttributeValue<std::string>("MODEL_DISPLAY");
        _displayLevel = modelDisplay.empty()
                            ? NOMAD::OutputLevel::LEVEL_DEBUGDEBUG
                            : NOMAD::OutputLevel::LEVEL_INFO;
    }
#endif
}


void NOMAD::QuadSearchMethod::generateTrialPointsImp()
{
#ifdef USE_SGTELIB
    // The trial points are generated for a feasible frame center and an infeasible one.
    if ( ! _stopReasons->checkTerminate() )
    {
        auto madsIteration = getParentOfType<MadsIteration*>();

        // MegaIteration's barrier member is already in sub dimension.
        auto bestXFeas = madsIteration->getMegaIterationBarrier()->getFirstXFeas();
        auto bestXInf  = madsIteration->getMegaIterationBarrier()->getFirstXInf();

        if (nullptr != bestXFeas)
        {
            NOMAD::QuadModelSinglePass singlePassFeas(this, bestXFeas, madsIteration->getMesh());

            // Generate the trial points
            singlePassFeas.generateTrialPoints();

            // Pass the generated trial pts to this
            auto trialPtsSinglePassFeas = singlePassFeas.getTrialPoints();
            for (auto point : trialPtsSinglePassFeas)
            {
                insertTrialPoint(point);
            }
        }
        if (nullptr != bestXInf)
        {
            NOMAD::QuadModelSinglePass singlePassInf(this, bestXInf, madsIteration->getMesh());

            // Generate the trial points
            singlePassInf.generateTrialPoints();

            // Pass the generated trial pts to this
            auto trialPtsSinglePassInf = singlePassInf.getTrialPoints();
            for (auto point : trialPtsSinglePassInf)
            {
                insertTrialPoint(point);
            }
        }
    }
#endif
}
 // end generateTrialPoints

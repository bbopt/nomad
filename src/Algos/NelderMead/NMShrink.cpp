
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/NelderMead/NMShrink.hpp"
#include "../../Algos/NelderMead/NMUpdate.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::NMShrink::init()
{
    _name = getAlgoName() + "Shrink";

    _currentStepType = NOMAD::NMStepType::SHRINK;

    _gamma = _runParams->getAttributeValue<NOMAD::Double>("NM_GAMMA");

    if ( _gamma <= 0.0 || _gamma > 1 )
        throw NOMAD::Exception(__FILE__,__LINE__,"Gamma value not compatible with shrink");

    verifyParentNotNull();

}



void NOMAD::NMShrink::startImp()
{

    // Update main barrier.
    NOMAD::NMUpdate update( this );
    update.start();
    update.run();
    update.end();

    // Create EvalPoints
    generateTrialPoints();

}


bool NOMAD::NMShrink::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }
    if ( getNbEvalPointsThatNeededEval() == 0 )
        setStopReason( );

    return foundBetter;
}


void NOMAD::NMShrink::endImp()
{
    // From IterationUtils
    postProcessing(NOMAD::EvcInterface::getEvaluatorControl()->getEvalType());

}


void NOMAD::NMShrink::generateTrialPoints ()
{

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    OUTPUT_INFO_START
    AddOutputInfo("Shrink simplex with " + _name +" (gamma=" + _gamma.tostring() +") with " + std::to_string(_nmY->size()) + " points.");
    OUTPUT_INFO_END


    // Shrink simplex Y
    std::set<NOMAD::EvalPoint>::const_iterator it = _nmY->begin() ;
    const NOMAD::EvalPoint & y0 = (*it);
    int i=0;
    for ( ; it !=_nmY->end(); ++it, ++i )
    {
        OUTPUT_INFO_START
        AddOutputInfo("y" + std::to_string(i) + ": " + (*it).display() );
        OUTPUT_INFO_END

        NOMAD::Point yi(n,0);
        for (size_t k = 0 ; k < n ; ++k )
        {
            yi[k] = (*it)[k];
        }
        NOMAD::Point y(n,0);
        for (size_t k = 0 ; k < n ; ++k )
        {
            y[k] = y0[k] + _gamma*(yi[k]-y0[k]);
        }


        // Shrink should not generate points outside the bounds
        // Shrink is not used with mesh
        // ----> No need to use snapPointToBoundsAndProjectOnMesh

        // New EvalPoint to be evaluated.
        // Add it to the list.
        bool inserted = insertTrialPoint(NOMAD::EvalPoint(y));

        OUTPUT_INFO_START
        std::string s = "xr:";
        s += (inserted) ? " " : " not inserted: ";
        s += y.display();
        AddOutputInfo(s);
        OUTPUT_INFO_END

        // Test for too small shrink (not the first point)
        if (i > 0 && y == yi )
        {
            OUTPUT_INFO_START
            AddOutputInfo("Shrink point to close to simplex point.");
            OUTPUT_INFO_END
            setStopReason( );
            clearTrialPoints();
            return;
        }
    }
}


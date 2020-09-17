
#include "../Util/AllStopReasons.hpp"

NOMAD::StopReason<NOMAD::BaseStopType> NOMAD::AllStopReasons::_baseStopReason = NOMAD::StopReason<NOMAD::BaseStopType>();
NOMAD::StopReason<NOMAD::EvalStopType> NOMAD::AllStopReasons::_evalStopReason = NOMAD::StopReason<NOMAD::EvalStopType>();

void NOMAD::AllStopReasons::setStarted()
{
    _baseStopReason.setStarted();
    _evalStopReason.setStarted();
    _iterStopReason.setStarted();
}


std::string NOMAD::AllStopReasons::getStopReasonAsString() const
{
    std::string stopReason="";
    bool flagTerminate = false;

    if (_baseStopReason.checkTerminate())
    {
        stopReason += getBaseStopReasonAsString();
        flagTerminate = true;
    }

    if (_evalStopReason.checkTerminate())
    {
        stopReason += getEvalStopReasonAsString();
        flagTerminate = true;
    }

    if (_iterStopReason.checkTerminate())
    {
        stopReason += _iterStopReason.getStopReasonAsString() + " (Iter) ";
        flagTerminate = true;
    }


    if (!flagTerminate)
    {
        stopReason = "No termination (all). ";
    }

    return stopReason;

}


std::string NOMAD::AllStopReasons::getEvalStopReasonAsString()
{
    return _evalStopReason.getStopReasonAsString() + " (Eval) ";
}


std::string NOMAD::AllStopReasons::getBaseStopReasonAsString()
{
    return _evalStopReason.getStopReasonAsString() + " (Base) ";
}


bool NOMAD::AllStopReasons::checkTerminate() const
{
    return ( _baseStopReason.checkTerminate()
            || _evalStopReason.checkTerminate()
            || _iterStopReason.checkTerminate());
}



#include "../Output/OutputInfo.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Math/Point.hpp"

void NOMAD::OutputInfo::addMsgAndSol(const std::string& msg, const NOMAD::Point & point)
{
    // Note: Currently displaying solution in a formated mode.
    auto solFormat = NOMAD::OutputQueue::getInstance()->getSolFormat();
    addMsg(msg + point.display(solFormat));
}

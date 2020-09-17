
#include "../../Algos/QuadModel/QuadModelIterationUtils.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::QuadModelIterationUtils::init()
{


    auto iter = dynamic_cast<const QuadModelIteration*>(_iterAncestor);
    if ( nullptr != iter )
    {
        _model = iter->getModel();
        _trainingSet = iter->getTrainingSet();
    }
}

void NOMAD::QuadModelIterationUtils::displayModelInfo() const
{
    if ( nullptr == _model || nullptr == _trainingSet )
        NOMAD::Exception(__FILE__, __LINE__, "The iteration utils must have a model and a training set to work with");

    OUTPUT_DEBUG_START
    NOMAD::OutputInfo dbgInfo("Quad Model iteration utils", "", NOMAD::OutputLevel::LEVEL_DEBUG );

//    dbgInfo.addMsg("The quad model from sgtelib: ");
//        dbgInfo.addMsg( _model->display()) ;

    NOMAD::OutputQueue::Add(std::move(dbgInfo));
    NOMAD::OutputQueue::Flush();
    OUTPUT_DEBUG_END
}

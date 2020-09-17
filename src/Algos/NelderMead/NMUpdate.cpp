
#include "../../Algos/NelderMead/NMUpdate.hpp"

void NOMAD::NMUpdate::init()
{
    _name = getAlgoName() + "Update";
    verifyParentNotNull();
}


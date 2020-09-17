/**
 \file   BBInput.cpp
 \brief  Input to a Blackbox evaluation
 \author Viviane Rochon Montplaisir
 \date   March 2018
 \see    BBInput.hpp
 */
#include "../Eval/BBInput.hpp"

/*-------------------------------------*/
/*            Constructor              */
/*-------------------------------------*/
NOMAD::BBInput::BBInput(const NOMAD::BBInputTypeList& bbInputTypeList,
                        const NOMAD::Point& point)
{
    // TODO
}


void NOMAD::BBInput::display(std::ostream &out) const
{
    // _rawBBI not implemented. Currently nothing to display.
    //out << _rawBBI;
}




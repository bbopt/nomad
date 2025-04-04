/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/
/**
 \file   BBInput.hpp
 \brief  Input of a Blackbox evaluation
 \author Viviane Rochon Montplaisir
 \date   March 2018
 \see    BBInput.cpp
 */

// Manage input to blackbox:
// Variable types - Continuous, integer, binary
// Scaling (future work)
//
// Note: As of September 2019, this class is not used by NOMAD.

#ifndef __NOMAD_4_5_BB_INPUT__
#define __NOMAD_4_5_BB_INPUT__


#include "../Math/Point.hpp"
#include "../Type/BBInputType.hpp"

#include "../nomad_nsbegin.hpp"


/// Class for the representation of the input to a blackbox evaluation.
class DLL_EVAL_API BBInput
{
    //std::string             _rawBBI;        // Actual input string (currently not implemented).

public:

    /*---------------*/
    /* Class Methods */
    /*---------------*/
    /// Constructor
    /**
      Currently does nothing.
     \param bbInputTypeList     The list of blackbox input types  -- \b IN.
     \param point               The point -- \b IN.
     */
    explicit BBInput(const BBInputTypeList& bbInputTypeList,
                     const Point& point);

    /*---------*/
    /* Get/Set */
    /*---------*/
    // To implement as needed.
    // Maybe something like this:
    //void setBBInput(const std::string bbInputString);
    //void setBBInput(const Point point);
    //void getBBInput(Point& point) const;

    // Currently does nothing.
    void display(std::ostream &out) const;

};


inline std::ostream& operator<<(std::ostream& out, const BBInput &bbinput)
{
    bbinput.display(out);
    return out;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD_4_5_BB_INPUT__

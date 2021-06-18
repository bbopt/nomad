/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
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
 \file   Subproblem.hpp
 \brief  Subproblem of lesser dimension than the original problem
 \author Viviane Rochon Montplaisir
 \date   February 2019
 */
#ifndef __NOMAD_4_0_SUBPROBLEM__
#define __NOMAD_4_0_SUBPROBLEM__

#include "../Math/Point.hpp"
#include "../Param/PbParameters.hpp"

#include "../nomad_nsbegin.hpp"

/// Class to define an optimization subproblem
/**
*  Subproblem of lesser dimension than the original problem
*
* - Sets up the new parameters
* - Keeps the necessary information to bridge the gap between subproblem and
    original problem
*/
class Subproblem
{
private:
    /**
      * The elements of this point that have defined values are fixed value
      * "variables". The elements that are undefined are for true variables.
      * This Point is always in full dimension.
     */
    const Point  _fixedVariable;
    size_t       _dimension;   ///< Dimension of the subproblem.

    /**
     Reference to the original problem's PbParameters.
     */
    const std::shared_ptr<PbParameters>  _refPbParams;

    /**
     PbParameters converted to subdimension
     */
    std::shared_ptr<PbParameters>        _subPbParams;

public:
    /// Constructor
    /**
     Pb parameters will be recomputed as dimension has changed.
     */
    explicit Subproblem(const std::shared_ptr<PbParameters> refPbParams,
                        const Point& fullFixedVariable)
      : _fixedVariable(fullFixedVariable),
        _dimension(refPbParams->getAttributeValue<size_t>("DIMENSION")),
        _refPbParams(refPbParams),
        _subPbParams(nullptr)
    {
        init();
    }

    /// Destructor
    virtual ~Subproblem();

    // Get/Set

    const Point& getFixedVariable() const { return _fixedVariable; }
    const std::shared_ptr<PbParameters>& getPbParams() const { return _subPbParams; }

private:
    /// Helper for constructor calls to Subproblem::setupProblemParameters
    void init();

    /// Helper for constructor
    /**
     Construct the subproblem parameters (X0, LB, UB, mesh sizes, variable groups...) based on Subproblem::_fixedVariable
     \note If a new parameter with dimension (ex. a parameter of type ArrayOfDouble, Point, or Dimension)
     is added to the class PbParameters, this method will break.
     Currently supported parameters: X0 LOWER_BOUND UPPER_BOUND BB_INPUT_TYPE INITIAL_MESH_SIZE INITIAL_FRAME_SIZE MIN_MESH_SIZE MIN_FRAME_SIZE GRANULARITY VARIABLE_GROUP
    */
    void setupProblemParameters();


    ///  Helper for setupProblemParameters()
    void resetVariableGroupsAgainstFixedVariables(ListOfVariableGroup& lvg, const Point& fixedVar) const;
};


#include "../nomad_nsend.hpp"


#endif  // __NOMAD_4_0_SUBPROBLEM__

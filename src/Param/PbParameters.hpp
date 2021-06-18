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
#ifndef __NOMAD_4_0_PBPARAMETERS__
#define __NOMAD_4_0_PBPARAMETERS__

#include "../Param/Parameters.hpp"

#include "../nomad_nsbegin.hpp"

/// The class for the parameters defining the optimization problem.
/**
- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class PbParameters final : public Parameters
{
private:
    bool _showWarningMeshSizeRedefined;

public:

    explicit PbParameters()
      : Parameters(),
        _showWarningMeshSizeRedefined(true)
    {
        init();
    }

    /**
     The copy constructor is not implemented in the parent class to allow some control over what parameters can be copied or not. Use the deep copy function of parameters: Parameters::copyParameters.
     */
    PbParameters& operator=(const PbParameters& params) { copyParameters(params) ; return *this; }

    /**
     The copy constructor is not implemented in the parent class to allow some control over what parameters can be copied or not. Use the deep copy function of parameters: Parameters::copyParameters.
     */
    PbParameters(const PbParameters& params) : PbParameters() { copyParameters(params); }

    /// Check the sanity of parameters.
    void checkAndComply( );

    /// Do not show certain warnings
    void doNotShowWarnings() { _showWarningMeshSizeRedefined = false; }

private:
    /// Helper for constructor
    /**
     Register and set default values for all problem attributes. The information to register all the attributes is contained in pbAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable, called automatically by makefile when the pbAttributeDefinition.txt file is modified.
     */
    void init() override ;

    ///  Helper for checkAndComply()
    void setGranularityAndBBInputType();
    ///  Helper for checkAndComply()
    void setVariableGroups();
    ///  Helper for checkAndComply()
    void setFixedVariables();
    ///  Helper for checkAndComply()
    void setMinMeshParameters(const std::string &paramName);
    ///  Helper for checkAndComply()
    void setInitialMeshParameters();
    ///  Helper for checkAndComply()
    void checkX0ForGranularity() const;
    ///  Helper for checkAndComply()
    void checkX0AgainstBounds() const;
    ///  Helper for checkAndComply()
    void checkForGranularity(const std::string &paramName) const;
    ///  Helper for checkAndComply()
    void checkForGranularity(const std::string &paramName,
                             const ArrayOfDouble &arrayToCheck) const;

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_PBPARAMETERS__


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
#ifndef __NOMAD_4_0_RUNPARAMETERS__
#define __NOMAD_4_0_RUNPARAMETERS__

#include "../Param/EvaluatorControlGlobalParameters.hpp"
#include "../Param/Parameters.hpp"
#include "../Param/PbParameters.hpp"

#include "../nomad_nsbegin.hpp"

/// The class for the parameters defining the type of optimization/task to perform.
/**
The RunParameter are used by other parameters to update their value during sanity check.

- Register all parameters during construction.
- Implement the checkAndComply function for sanity check.
*/
class RunParameters final : public Parameters
{
private:
    static bool _warningUnknownParamShown;

public:
    // Constructor: init() will be called.
    // This will register and set default values to all attributes.
    explicit RunParameters()
    : Parameters()
    {
        init();
    }

    /// The copy constructor is not implemented in the parent class
    RunParameters& operator=(const RunParameters& params) { copyParameters(params) ; return *this; }
    RunParameters(const RunParameters& params) : RunParameters() { copyParameters(params); }

    /// Check the sanity of parameters.
    /**
     Register and set default values for all run attributes. The information to register all the attributes is contained in runAttributesDefinition.hpp as a set of strings to be interpreted. This file is created by the writeAttributeDefinition executable, called automatically by makefile when the runAttributeDefinition.txt file is modified.
    */
    void checkAndComply(const std::shared_ptr<EvaluatorControlGlobalParameters>& evaluatorControlGlobalParams,
                        const std::shared_ptr<PbParameters>& pbParams);

private:
    /// Initialization
    /**
     * This will register and set default values to all attributes.
     */
    void init() override ;

    /// Helper for checkAndComply()
    void setStaticParameters();

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_RUNPARAMETERS__

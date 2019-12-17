/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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
#ifndef __NOMAD400_RUNPARAMETERS__
#define __NOMAD400_RUNPARAMETERS__

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <typeindex>
#include <typeinfo>

#include "../Param/Parameters.hpp"
#include "../Param/EvaluatorControlParameters.hpp"
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
    void checkAndComply(const std::shared_ptr<EvaluatorControlParameters>& evaluatorControlParams,
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

#endif // __NOMAD400_RUNPARAMETERS__

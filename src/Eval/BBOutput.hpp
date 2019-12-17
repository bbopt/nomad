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
/**
 \file   BBOutput.hpp
 \brief  Output of a Blackbox evaluation
 \author Viviane Rochon Montplaisir
 \date   January 2018
 \see    BBOutput.cpp
 */



#ifndef __NOMAD400_BB_OUTPUT__
#define __NOMAD400_BB_OUTPUT__

#include "../Type/BBOutputType.hpp"
#include "../Math/ArrayOfDouble.hpp"
#include "../Math/Double.hpp"
#include "../Util/ArrayOfString.hpp"

#include "../nomad_nsbegin.hpp"


/// Class for the representation of the output of a blackbox evaluation.
/**
 *
 * Manage output from blackbox:
 *  - Raw output (string)
 *  - Is eval ok. This is a boolean indicating that there were no problem during evaluation.
 *  - Scaling (future work)
 */
class BBOutput {
public:
    static const std::string bboStart; ///< Static variable used for field delimitation.
    static const std::string bboEnd; ///< Static variable used for field delimitation.


private:
    std::string             _rawBBO;    ///< Actual output string
    bool                    _evalOk;    ///< Flag for evaluation

public:

    /*---------------*/
    /* Class Methods */
    /*---------------*/

    /// Constructor
    /**
     Usually if we have a rawBBO we can assume that the evaluation went OK.
     \param rawBBO  The outputs of the blackbox as a string -- \b IN.
     \param evalOk  The eval ok flag -- \b IN.
     */
    explicit BBOutput(const std::string &rawBBO, const bool evalOk = true);

    /*---------*/
    /* Get/Set */
    /*---------*/

    /// Get the objective from raw blackbox evaluation
    /**
     \param bbOutputType    The list of blackbox output types -- \b IN.
     \return                The objective value (single objective).
     */
    Double getObjective(const BBOutputTypeList &bbOutputType) const;
    
    /// Get the constraints from raw blackbox evaluation
    /**
     \param bbOutputType    The list of blackbox output types -- \b IN.
     \return                The constraints as a an array of values.
     */
    ArrayOfDouble getConstraints(const BBOutputTypeList &bbOutputType) const;
    
    /// Set each blackbox output separately from a string.
    /**
     \param bbOutputString    The string returned by blackbox evaluation -- \b IN.
     \param evalOk            The evaluation status -- \b IN.
     */
    void setBBO(const std::string &bbOutputString, const bool evalOk = true);

    /// Get if this evaluation proceeded properly
    /**
     \return \c True if the evaluation ended normally; \c False if there was an error.
     */
    bool getEvalOk() const { return _evalOk; }

    /// Does this evaluation count?
    /**
     \return \c True if it does
     */
    bool getCountEval(const BBOutputTypeList &bbOutputType) const;

    /// Get the raw blackbox outputs
    /**
     \return    A single string containing the raw blackbox outputs.
     */
    const std::string& getBBO() const { return _rawBBO; }
    
    
    /// Get the blackbox outputs separately.
    /**
     \return    An array of double of the blackbox outputs.
     */
    ArrayOfDouble getBBOAsArrayOfDouble() const;

    /// Display
    void display (std::ostream & out) const;

private:
    /// Helper functions that can trigger exception.
    /**
     Exception is triggered if the BBOutputList and
     the ArrayOfString have inconsistent size.
     \param bbOutputType    The list of blackbox output type -- \b IN.
     \param array           The array of string -- \b IN.
     
     */
    void checkSizeMatch(const BBOutputTypeList &bbOutputType,
                        const ArrayOfString &array) const;

};


/// Display, using field delimitors BBOutput::bboStart and BBOutput::bboEnd.
std::ostream& operator<<(std::ostream& os, const BBOutput &bbo);

/// Read bbo that has been written by BBOutput display operator.
std::istream& operator>>(std::istream& is, BBOutput &bbo);


#include "../nomad_nsend.hpp"
#endif // __NOMAD400_BB_OUTPUT__

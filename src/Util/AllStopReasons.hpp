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
#ifndef __NOMAD_4_0_ALLSTOPREASONS__
#define __NOMAD_4_0_ALLSTOPREASONS__

#include "../nomad_platform.hpp"
#include "../Util/StopReason.hpp"

#include "../nomad_nsbegin.hpp"


/// Class combining all stop reasons that are not algorithmic stop reasons.
/**

 Several stop reasons are members of this class. The stop reasons are templated on stop type. Several stop types are available in this class:
 - a ::BaseStopType for high level stop reasons.
 - an ::EvalGlobalStopType for global evaluation stop reasons.
 - an ::IterStopType for stop reasons during iteration of an algorithm (for example, maximum iteration number reached).

 The static stop reasons ::BaseStopType and ::EvalGlobalStopType are shared.
 */
class AllStopReasons
{
public:
    /// Constructor
    explicit AllStopReasons ()
    {
    }

    /// Destructor
    virtual ~AllStopReasons()
    {}

private:
    DLL_UTIL_API static StopReason<BaseStopType> _baseStopReason; ///< A single base stop reason is considered for NOMAD.
    DLL_UTIL_API static StopReason<EvalGlobalStopType> _evalGlobalStopReason; ///< An eval stop reason valable for the whole of NOMAD.
    StopReason<IterStopType> _iterStopReason; ///< An iteration stop reason.

public:
    /*---------*/
    /* Get/Set */
    /*---------*/

    static const StopReason<BaseStopType>& getBaseStopReason() { return _baseStopReason; }
    const StopReason<EvalGlobalStopType>& getEvalGlobalStopReason() { return _evalGlobalStopReason; }
    const StopReason<IterStopType>& getIterStopReason() const { return _iterStopReason; }

    static void set(const BaseStopType& s)
    {
        _baseStopReason.set(s);
    }

    static void set(const EvalGlobalStopType& s)
    {
        _evalGlobalStopReason.set(s);
    }

    void set(const IterStopType& s)
    {
        _iterStopReason.set(s);
    }

    /*---------*/
    /* Other   */
    /*---------*/

    /// Test static BaseStopType
    static bool testIf(const BaseStopType& s)
    {
        return (_baseStopReason.get() == s);
    }

    /// Test static EvalGlobalStopType
    static bool testIf (const EvalGlobalStopType& s)
    {
        return (_evalGlobalStopReason.get() == s);
    }

    /// Test IterStopType
    bool testIf (IterStopType s)
    {
        return (_iterStopReason.get() == s);
    }

    /// Reset all stop reasons to their default STARTED state
    virtual void setStarted();

    /// Get the stop reason that requires termination as a string.
    /**
     If no termination is required, an empty string is returned.
     */
    virtual std::string getStopReasonAsString() const;


    /// Get the global eval stop reason as a string.
    /**
    \return An empty string is in STARTED state, the stop reason otherwise.
     */
    static std::string getEvalGlobalStopReasonAsString();

    /// Get the base stop reason as a string.
    /**
     \return An empty string is in STARTED state, the stop reason otherwise.
     */
    static std::string getBaseStopReasonAsString();

    /// Check if among all stop reasons, one requires a termination.
    /**
     \see StopReason::checkTerminate()

     \return \c true if a termination is required, \c false otherwise.
     */
    virtual bool checkTerminate() const;

    static bool checkBaseTerminate()
    {
        return _baseStopReason.checkTerminate();
    }

    static bool checkEvalGlobalTerminate()
    {
        return _evalGlobalStopReason.checkTerminate();
    }
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_ALLSTOPREASONS__

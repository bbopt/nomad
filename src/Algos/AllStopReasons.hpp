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
#ifndef __NOMAD400_ALLSTOPREASONS__
#define __NOMAD400_ALLSTOPREASONS__

#include "../Util/Exception.hpp"
#include "../Util/StopReason.hpp"

#include "../nomad_nsbegin.hpp"


/// Class combining all stop reasons that are not algorithmic stop reasons.
/**

 Several stop reasons are members of this class. The stop reasons are templated on stop type. Several stop types are available in this class:
 - a ::BaseStopType for high level stop reasons.
 - a ::EvalStopType for evaluation stop reasons.
 - an ::IterStopType for stop reasons during iteration of an algorithm (for example, maximum iteration number reached).
 
 The static stop reasons ::BaseStopType and ::EvalStopType are shared.
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
    static StopReason<BaseStopType> _baseStopReason ; ///< A single base stop reason is considered for NOMAD.
    static StopReason<EvalStopType> _evalStopReason ; ///< A single evaluation stop reason is considered for NOMAD.
    StopReason<IterStopType> _iterStopReason; ///< An iteration stop reason.


public:
    
    
    /*---------*/
    /* Get/Set */
    /*---------*/

    StopReason<BaseStopType> & getBaseStopReason() { return _baseStopReason; }
    StopReason<EvalStopType> & getEvalStopReason() { return _evalStopReason; }
    StopReason<IterStopType> & getIterStopReason() { return _iterStopReason; }
    
    static void set( BaseStopType s )
    {
        _baseStopReason.set(s);
    }
    
    static void set( EvalStopType s )
    {
        _evalStopReason.set(s);
    }
    
    void set( IterStopType s )
    {
        _iterStopReason.set(s);
    }
    
    /*---------*/
    /* Other   */
    /*---------*/
    
    /// Test static BaseStopType
    static bool testIf ( BaseStopType s )
    {
        return ( _baseStopReason.get() == s );
    }
    
    /// Test static EvalStopType
    static bool testIf ( EvalStopType s )
    {
        return ( _evalStopReason.get() == s );
    }
    
    /// Test IterStopType
    bool testIf ( IterStopType s )
    {
        return ( _iterStopReason.get() == s );
    }
    
    /// Reset all stop reasons to their default STARTED state
    virtual void setStarted ()
    {
        _baseStopReason.setStarted();
        _evalStopReason.setStarted();
        _iterStopReason.setStarted();
    }

    /// Get the stop reason that requires termination as a string.
    /**
     If no termination is required, an empty string is returned.
     */
    virtual std::string getStopReasonAsString() const
    {
        std::string stopReason="";
        bool flagTerminate = false;
        
        if ( _baseStopReason.checkTerminate() )
        {
            stopReason += _baseStopReason.getStopReasonAsString() + " (Base stop reason) ";
            flagTerminate = true;
        }
        if ( _evalStopReason.checkTerminate() )
        {
            stopReason += _evalStopReason.getStopReasonAsString() + " (Eval stop reason) ";
            flagTerminate = true;
        }
        
        if ( _iterStopReason.checkTerminate() )
        {
            stopReason += _iterStopReason.getStopReasonAsString() + " (Iteration stop reason) ";
            flagTerminate = true;
        }
        
        if ( ! flagTerminate )
            stopReason = "No termination (all). ";
        
        return stopReason;
        
    }
    
    
    /// Get the eval stop reason as a string.
    /**
    \return An empty string is in STARTED state, the stop reason otherwise.
     */
    static std::string getEvalStopReasonAsString()
    {
        std::string stopReason="";
        
        if ( ! _evalStopReason.isStarted() )
            stopReason += _evalStopReason.getStopReasonAsString() + " (Eval) ";
        
        return stopReason;
        
    }
    
    /// Get the base stop reason as a string.
    /**
     \return An empty string is in STARTED state, the stop reason otherwise.
     */
    static std::string getBaseStopReasonAsString()
    {
        std::string stopReason="";
        
        if ( ! _baseStopReason.isStarted() )
            stopReason += _baseStopReason.getStopReasonAsString() + " (Eval) ";
        
        return stopReason;
        
    }

    /// Check if among all stop reasons, one requires a termination.
    /**
     \see StopReason::checkTerminate()
     
     \return \c true if a termination is required, \c false otherwise.
     */
    virtual bool checkTerminate () const
    {
        return ( _baseStopReason.checkTerminate()
                || _evalStopReason.checkTerminate()
                || _iterStopReason.checkTerminate() );
    }
    
    static bool checkBaseTerminate ()
    {
        return _baseStopReason.checkTerminate () ;
    }
    
    static bool checkEvalTerminate ()
    {
        return _evalStopReason.checkTerminate () ;
    }
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_ALLSTOPREASONS__

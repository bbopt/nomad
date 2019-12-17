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
#ifndef __NOMAD400_ALGOSTOPREASONS__
#define __NOMAD400_ALGOSTOPREASONS__

#include "../Algos/AllStopReasons.hpp"

#include "../Util/Exception.hpp"
#include "../Util/StopReason.hpp"

#include "../nomad_nsbegin.hpp"


/// Template class for algorithm stop reasons.
/**
 
 The class is templated with a StopType defined according to which algorithm is considered. For example, we have ::MadsStopType, ::LHStopType and ::NMStopType. \n
 At some point during an algorithm a stop reason is set. It can be specific to the algorithm or generic (that is be an AllStopReasons stop type). \n
 The stop reasons in AllStopReasons are private and not directly accessible. But the AlgoStopReasons::checkTerminate() function, checks both AllStopReasons and AlgoStopReasons.
 */
template <typename StopType>
class AlgoStopReasons : public AllStopReasons
{
public:
    /// Constructor
    /*
     */
    explicit AlgoStopReasons () : AllStopReasons()
    {
    }
    
    /// Destructor
    virtual ~AlgoStopReasons()
    {}
    
    
private:
    StopReason<StopType> _algoStopReason;
    
public:
    
    /// Access to the algo stop reason (no the other generic stop reasons).
    StopReason<StopType> & getAlgoStopReason() { return _algoStopReason; }

    /// Set the algo stop reason to a specific stop type.
    void set( StopType s )
    {
        _algoStopReason.set(s);
    }
    
    std::string getStopReasonAsString() const override
    {
        std::string stopReason= AllStopReasons::getStopReasonAsString();
        
        if ( ! _algoStopReason.isStarted() )
            stopReason += _algoStopReason.getStopReasonAsString() + " (Algo) ";
        
        return stopReason;
        
    }
    
    /// Check among generic stop reasons and algo stop reason if the algorithm must terminate
    bool checkTerminate () const override
    {
        return ( AllStopReasons::checkTerminate()
                || _algoStopReason.checkTerminate() );
    }
    
    /// Access to the AlgoStopReasons
    static std::shared_ptr<AlgoStopReasons<StopType>> get ( std::shared_ptr<AllStopReasons> allStopReasons )
    {
        std::shared_ptr<AlgoStopReasons<StopType>> stopReasons = std::dynamic_pointer_cast<AlgoStopReasons<StopType>>( allStopReasons );
        
        if ( stopReasons == nullptr )
            throw Exception(__FILE__, __LINE__, "Invalid shared pointer cast");
        return stopReasons;
    }
    
    /// Test for a specific algorithm stop type.
    /**
     Used to pass a sub-algorithm stop reason to a parent algorithm stop reason.
     */
    bool testIf ( StopType s )
    {
        return ( _algoStopReason.get() == s );
    }

    /// Reset stop reasons to their default STARTED state.
    void setStarted () override
    {
        _algoStopReason.setStarted();
        AllStopReasons::setStarted();
    }
    
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_ALGOSTOPREASONS__

#ifndef __NOMAD400_ALGOSTOPREASONS__
#define __NOMAD400_ALGOSTOPREASONS__

#include <memory>   // for shared_ptr
#include "../Util/Exception.hpp"
#include "../Util/AllStopReasons.hpp"

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

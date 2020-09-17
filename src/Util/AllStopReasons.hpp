#ifndef __NOMAD400_ALLSTOPREASONS__
#define __NOMAD400_ALLSTOPREASONS__

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
    static StopReason<BaseStopType> _baseStopReason; ///< A single base stop reason is considered for NOMAD.
    static StopReason<EvalStopType> _evalStopReason; ///< A single eval stop reason is considered for NOMAD.
    StopReason<IterStopType> _iterStopReason; ///< An iteration stop reason.

public:
    /*---------*/
    /* Get/Set */
    /*---------*/

    static const StopReason<BaseStopType>& getBaseStopReason() { return _baseStopReason; }
    static const StopReason<EvalStopType>& getEvalStopReason() { return _evalStopReason; }
    const StopReason<IterStopType>& getIterStopReason() const { return _iterStopReason; }

    static void set(const BaseStopType& s)
    {
        _baseStopReason.set(s);
    }

    static void set(const EvalStopType& s)
    {
        _evalStopReason.set(s);
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

    /// Test static EvalStopType
    static bool testIf (const EvalStopType& s)
    {
        return (_evalStopReason.get() == s);
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


    /// Get the eval stop reason as a string.
    /**
    \return An empty string is in STARTED state, the stop reason otherwise.
     */
    static std::string getEvalStopReasonAsString();

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

    static bool checkEvalTerminate()
    {
        return _evalStopReason.checkTerminate();
    }
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_ALLSTOPREASONS__

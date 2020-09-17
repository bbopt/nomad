#ifndef __NOMAD400_USERSEARCHMETHOD__
#define __NOMAD400_USERSEARCHMETHOD__

#include "../../Algos/Mads/SearchMethodSimple.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class UserSearchMethod: Search method defined by user.
class UserSearchMethod: public SearchMethodSimple
{
public:
    /// Constructor
    /**
     \param parentStep      The parent of this search step -- \b IN.
     */
    explicit UserSearchMethod(const Step* parentStep)
      : SearchMethodSimple(parentStep)
    {
        init();
    }

private:

    /// Helper for constructor
    void init();

    /**
     \copydoc SearchMethodSimple::generateTrialPointsImp \n
     A user can implement this function.
     */
    virtual void generateTrialPointsImp() override
    {
        throw Exception(__FILE__, __LINE__, "User search not implemented.");
    };

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_USERSEARCHMETHOD__


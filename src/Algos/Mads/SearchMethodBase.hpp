#ifndef __NOMAD400_SEARCHMETHODBASE__
#define __NOMAD400_SEARCHMETHODBASE__

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for generic search method of MADS. Run by Search.
/**
 Pure virtual class from which SearchMethodSimple and SearchMethodAlgo derive.
 */
class SearchMethodBase: public Step, public IterationUtils
{
private:

    bool _enabled; ///< Should this search method be used? Modified by parameters.

    std::string _comment; ///<  Comment shown when a search method is used

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SearchMethodBase( const Step* parentStep )
      : Step( parentStep ),
        IterationUtils ( parentStep ),
        _enabled(true),
        _comment("")
    {
        init();
    }

    bool isEnabled() const { return _enabled; }
    void setEnabled(const bool enabled) { _enabled = enabled; }

    const std::string& getComment() const { return _comment; }
    bool hasComment() const { return (!_comment.empty()); }
    void setComment(const std::string& comment) { _comment = comment; }

    /**
     - Pure virtual function.
     - The implementation of startImp function in the derived class generates trial points  (in SearchMethodSimple) OR does nothing (in SearchMethodAlgo).
     */
    virtual void startImp() override =0 ;

    /**
     - Pure virtual function.
     - The implementation of runImp function in the derived class evaluates the trial points (in SearchMethodSimple) OR launches an algo (in SearchMethodAlgo).
     */
    virtual bool runImp() override = 0 ;

    /// Implementation of endImp (not virtual)
    /**
        Call to the postProcessing function to update the Barrier
    */
    void endImp() override ;

    /// Intermediate function (not yet implementation that can generate the trial points)
    /**
     - Display before and after generation comments.
     - Launches the implementation of the search method to generate the trial points (::generateTrialPointsImp).
     - Snap the points to bounds and mesh.
     */
    void generateTrialPoints() override;

    /**
     - Pure virtual function.
     - See derived classes (SearchMethodSimple, SearchMethodAlgo) for implementations.
     */
    virtual void generateTrialPointsImp() = 0 ;


protected:
    void init();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SEARCHMETHODBASE__


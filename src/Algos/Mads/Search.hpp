#ifndef __NOMAD400_SEARCH__
#define __NOMAD400_SEARCH__

#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to manage the SearchMethods used by MADS algorithm during its search step.
class Search final : public Step , public IterationUtils
{
private:
    std::vector<std::shared_ptr<SearchMethodBase>> _searchMethods;
#ifdef TIME_STATS
    static std::vector<double>  _searchTime;        ///< Total time spent running each search
    static std::vector<double>  _searchEvalTime;    ///< Total time spent evaluating search points
#endif // TIME_STATS

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit Search(const Step* parentStep )
      : Step( parentStep ),
        IterationUtils( parentStep ),
        _searchMethods()
    {
        init();
    }

    virtual ~Search() {}

    /**
     - Generate new points to evaluate. Use all enabled search methods.
     - To be used only when parameter GENERATE_ALL_POINTS_BEFORE_EVAL is true.
     */
    void generateTrialPoints() override;

#ifdef TIME_STATS
    /// Time stats
    static std::vector<double> getSearchTime()       { return _searchTime; }
    static std::vector<double> getSearchEvalTime()   { return _searchEvalTime; }
#endif // TIME_STATS

private:

    void init();

    /// Implementation of the start task.
    /**
     Just perform a sanity check on GENERATE_ALL_POINTS_BEFORE_EVAL that must be false.
     */
    virtual void startImp() override;

    /// The implementation of run tasks.
    /**
      Perform start+run+end for all search methods in the vector _searchMethods.
     */
    virtual bool runImp() override;

    /// Implementation of the end tasks
    /**
      If a sub optimization is used during search we probably set a stop reason to terminate. The parent optimization must go on. The stop reason is set to started if sub optimization reached its evaluation budget.
     */
    virtual void endImp() override;

    /**
     Identify if there is at least one search enabled. If there are none, do not print Search step at all.
     */
    bool isEnabled() const;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SEARCH__



#ifndef __NOMAD400_SUBPROBLEMMANAGER__
#define __NOMAD400_SUBPROBLEMMANAGER__

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include "../Algos/Algorithm.hpp"
#include "../Algos/Subproblem.hpp"

#include "../nomad_nsbegin.hpp"


/// Class to associate Algorithms with Subproblems
/**
 Ease the passage between sub-dimension and full dimension. Algorithm works in
 a sub dimentsion and does not know the full dimension.
 */
class SubproblemManager
{
private:
    static std::map<const Algorithm*, const Subproblem> _map;
#ifdef _OPENMP
    static omp_lock_t _mapLock;
#endif // _OPENMP

public:
    /// Constructor
    explicit SubproblemManager()
    {
        init();
    }

    /// Destructor
    virtual ~SubproblemManager()
    {
        destroy();
    }

    static const Subproblem& getSubproblem(const Step* step);

    static const Point& getSubFixedVariable(const Step* step);

    static void addSubproblem(const Algorithm* algo, const Subproblem& subproblem);

    static void removeSubproblem(const Algorithm* algo);

    static void reset();

private:
    /// Helper for constructor
    void init();

    /// Helper for destructor
    void destroy();
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_SUBPROBLEMMANAGER__

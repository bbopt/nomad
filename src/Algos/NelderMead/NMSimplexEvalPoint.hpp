#ifndef __NOMAD400_NMSIMPLEX_EVALPOINT__
#define __NOMAD400_NMSIMPLEX_EVALPOINT__

#include "../../Eval/EvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Compare EvalPoints for NelderMead algorithm
class NMSimplexEvalPointCompare
{
public:
    bool operator() (const EvalPoint& lhs, const EvalPoint& rhs) const;
};

// Set of EvalPoints for NelderMead algorithm
typedef std::set<EvalPoint, NMSimplexEvalPointCompare> NMSimplexEvalPointSet;
typedef NMSimplexEvalPointSet::iterator NMSimplexEvalPointSetIterator;

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_NMSIMPLEX_EVALPOINT__

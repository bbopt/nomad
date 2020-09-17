
#include "../Algos/SubproblemManager.hpp"

// Initialize static variables
std::map<const NOMAD::Algorithm*, const NOMAD::Subproblem> NOMAD::SubproblemManager::_map = std::map<const NOMAD::Algorithm*, const NOMAD::Subproblem>();
#ifdef _OPENMP
omp_lock_t NOMAD::SubproblemManager::_mapLock;
#endif // _OPENMP

void NOMAD::SubproblemManager::init()
{
#ifdef _OPENMP
    omp_init_lock(&_mapLock);
#endif // _OPENMP
}


void NOMAD::SubproblemManager::destroy()
{
#ifdef _OPENMP
    omp_destroy_lock(&_mapLock);
#endif // _OPENMP
}


void NOMAD::SubproblemManager::addSubproblem(const NOMAD::Algorithm* algo, const NOMAD::Subproblem& subproblem)
{
    auto algoSubPair = std::pair<const NOMAD::Algorithm*, const NOMAD::Subproblem&>(algo, subproblem);
#ifdef _OPENMP
    omp_set_lock(&_mapLock);
#endif // _OPENMP
    _map.insert(algoSubPair);
#ifdef _OPENMP
    omp_unset_lock(&_mapLock);
#endif // _OPENMP
}


void NOMAD::SubproblemManager::removeSubproblem(const Algorithm* algo)
{
#ifdef _OPENMP
    omp_set_lock(&_mapLock);
#endif // _OPENMP
    int nbErased = _map.erase(algo);
#ifdef _OPENMP
    omp_unset_lock(&_mapLock);
#endif // _OPENMP
    if (0 == nbErased)
    {
        std::cerr << "Warning: SubproblemManager could not remove subproblem for Algorithm " << algo->getName() << std::endl;
    }

}


void NOMAD::SubproblemManager::reset()
{
    if (_map.size() > 0)
    {
        // Shoud not happen. Warn the user.
        std::cerr << "Warning: SubproblemManager::clear() called on non-empty SubproblemManager" << std::endl;
    }
#ifdef _OPENMP
    omp_set_lock(&_mapLock);
#endif // _OPENMP
    _map.clear();
#ifdef _OPENMP
    omp_unset_lock(&_mapLock);
#endif // _OPENMP
}


const NOMAD::Subproblem& NOMAD::SubproblemManager::getSubproblem(const NOMAD::Step* step)
{
    NOMAD::Algorithm* algo;
    std::string s;

    if (step->isAnAlgorithm())
    {
        algo = dynamic_cast<NOMAD::Algorithm*>(const_cast<NOMAD::Step*>(step));
    }
    else
    {
        algo = step->getParentOfType<NOMAD::Algorithm*>();
    }

    if (nullptr == algo)
    {
        s = "Algorithm not found for step " + step->getName();
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    try
    {
        return _map.at(algo);
    }
    catch (const std::out_of_range& oor)
    {
        std::cerr << "Subproblem not found for Algorithm " << algo->getName() << std::endl;
    }

    s = "SubproblemManager could not get Subproblem for step " + step->getName();
    throw NOMAD::Exception(__FILE__,__LINE__,s);

}


const NOMAD::Point& NOMAD::SubproblemManager::getSubFixedVariable(const NOMAD::Step* step)
{
    return getSubproblem(step).getFixedVariable();
}

/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
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

#include "../Algos/SubproblemManager.hpp"

// Initialize singleton
std::unique_ptr<NOMAD::SubproblemManager> NOMAD::SubproblemManager::_single=nullptr;


#ifdef _OPENMP
omp_lock_t NOMAD::SubproblemManager::_mapLock;
#endif

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
    auto retPair = _map.insert(algoSubPair);
    if (false == retPair.second)
    {
        std::string err = "Error: SubproblemManager: could not add subproblem for Algorithm ";
        err += algo->getName();
        throw NOMAD::StepException(__FILE__,__LINE__, err, algo);
    }
#ifdef _OPENMP
    omp_unset_lock(&_mapLock);
#endif // _OPENMP
}


void NOMAD::SubproblemManager::removeSubproblem(const Algorithm* algo)
{
#ifdef _OPENMP
    omp_set_lock(&_mapLock);
#endif // _OPENMP
    size_t nbErased = _map.erase(algo);
#ifdef _OPENMP
    omp_unset_lock(&_mapLock);
#endif // _OPENMP
    if (0 == nbErased)
    {
        std::string err = "Warning: SubproblemManager could not remove subproblem for Algorithm " + algo->getName();
        throw NOMAD::StepException(__FILE__,__LINE__, err, algo);
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
    std::string err;

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
        err = "Algorithm not found for step " + step->getName();
        throw NOMAD::StepException(__FILE__,__LINE__,err,step);
    }

    try
    {
        // Lock is needed even if map is not modified by this method,
        // because it could be modified by another main thread at the same moment.
#ifdef _OPENMP
        omp_set_lock(&_mapLock);
#endif // _OPENMP
        const NOMAD::Subproblem& sub = _map.at(algo);
#ifdef _OPENMP
        omp_unset_lock(&_mapLock);
#endif // _OPENMP
        return sub;
    }
    catch (const std::out_of_range&)
    {
        std::cerr << "Error: Subproblem not found for Algorithm " << algo->getName() << std::endl;
    }

    err = "SubproblemManager could not get Subproblem for step " + step->getName();
    throw NOMAD::StepException(__FILE__,__LINE__,err, step);

}


const NOMAD::Point& NOMAD::SubproblemManager::getSubFixedVariable(const NOMAD::Step* step)
{
    return getSubproblem(step).getFixedVariable();
}

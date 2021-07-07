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
/**
 \file   ComparePriority.hpp
 \brief  Compare EvalQueuePoints for sorting
 \author Viviane Rochon Montplaisir
 \date   November 2020
 \see    ComparePriority.cpp
 */

#ifndef __NOMAD_4_0_COMPAREPRIORITY__
#define __NOMAD_4_0_COMPAREPRIORITY__

#include "../Eval/EvalQueuePoint.hpp"
#include "../Math/Direction.hpp"
#include "../Math/RandomPickup.hpp"

#include "../nomad_nsbegin.hpp"

//typedef std::function<bool(EvalQueuePointPtr &p1, EvalQueuePointPtr &p2)> ComparePriorityFunction;

/// Definition for compare priority method.
/**
 An instance of this class has a comp() method that compares two EvalQueuePoints for ordering.
*/
class ComparePriorityMethod
{
private:
    std::string _name;  ///< Method name, useful for information or debugging

public:
    virtual bool comp(EvalQueuePointPtr& NOMAD_UNUSED(point1), EvalQueuePointPtr& NOMAD_UNUSED(point2)) const
    {
        return false;
    }

    void setName(const std::string& name) { _name = name; }
    const std::string& getName() const { return _name; }
};


class BasicComp : public ComparePriorityMethod
{
public:
    /// Constructor
    explicit BasicComp()
    {
        setName("BasicComp");
    }

    bool comp(EvalQueuePointPtr& point1, EvalQueuePointPtr& point2) const override;
};


// Class for comparison using a direction.
class OrderByDirection : public ComparePriorityMethod
{
private:
    /** Vector of directions: One per main thread; one list for feasible and
      * infeasible points. Makes it possible
      * to compare points from different algorithms.
     **/
    std::vector<std::shared_ptr<Direction>> _lastSuccessfulFeasDirs;
    std::vector<std::shared_ptr<Direction>> _lastSuccessfulInfDirs;

public:
    /// Constructor
    explicit OrderByDirection(const std::vector<std::shared_ptr<Direction>>& feasDirs,
                              const std::vector<std::shared_ptr<Direction>>& infDirs)
      : _lastSuccessfulFeasDirs(feasDirs),
        _lastSuccessfulInfDirs(infDirs)
    {
        setName("OrderByDirection");
    }

    bool comp(EvalQueuePointPtr& point1, EvalQueuePointPtr& point2) const override;
};


// Class for mixing points randomly
class RandomComp : public ComparePriorityMethod
{
private:
    mutable RandomPickup                _randomPickup;
    mutable std::map<size_t, size_t>    _tagToRank;

public:
    /// Constructor
    explicit RandomComp(const size_t n);

    bool comp(EvalQueuePointPtr& point1, EvalQueuePointPtr& point2) const override;
};


// Class for comparison using static surrogate evaluations.
class OrderBySurrogate : public ComparePriorityMethod
{
public:
    /// Constructor
    explicit OrderBySurrogate()
    {
        setName("OrderBySurrogate");
    }

    bool comp(EvalQueuePointPtr& point1, EvalQueuePointPtr& point2) const override;
};

/// Class to compare priority of two EvalQueuePoint
class ComparePriority
{
private:
    std::shared_ptr<ComparePriorityMethod>  _compMethod; ///< Comparison method to be used to sort eval queue points

public:
    /// Constructor
    explicit ComparePriority(const std::shared_ptr<ComparePriorityMethod>& compMethod)
      : _compMethod(compMethod)
    {}

    ///  Function call operator
    /**
     \param p1  First eval queue point -- \b IN.
     \param p2  Second eval queue point -- \b IN.
     \return    \c true if p1 has a lower priority than p2. \c false otherwise.
     */
    bool operator()(EvalQueuePointPtr& p1, EvalQueuePointPtr& p2);
};


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_COMPAREPRIORITY__



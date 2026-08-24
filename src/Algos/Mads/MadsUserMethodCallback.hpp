/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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

#ifndef __NOMAD_4_6_MADSUSERMETHODCALLBACK__
#define __NOMAD_4_6_MADSUSERMETHODCALLBACK__

#include "../../Type/CallbackType.hpp"

#include <memory>
#include <functional>

#include "../../nomad_nsbegin.hpp"

class Step;

// Base class for algo callbacks (can be stored into a vector)
class MadsCallbackBase
{
public:
    virtual ~MadsCallbackBase() = default;
    
    virtual NOMAD::MadsCallbackType getType() const = 0;
    
    virtual std::unique_ptr<MadsCallbackBase> clone() const = 0;
    
};

/// Templated general class for algo callbacks
template<NOMAD::MadsCallbackType CT>
class  MadsCallback : public MadsCallbackBase
{
public:
    MadsCallback() = default;
    virtual ~MadsCallback() = default;
    
    NOMAD::MadsCallbackType getType() const override { return CT;}
    
    std::unique_ptr<MadsCallbackBase> clone() const override
    {
        return std::make_unique<MadsCallback<CT>>(*this);
    }
};

//
// Define types for all the types of mads callback functions available
//
typedef std::function<bool(const Step& step, std::list<Direction> & dir, const size_t n)> UserPollMethodCbFunc;  ///< Type definitions for callback functions for user Poll method.
typedef std::function<bool(const Step& step, EvalPointSet & trialPoint)> UserSearchMethodCbFunc;  ///< Type definitions for callback functions for user Search method.
typedef std::function<bool(const Step& step)> UserMethodEndCbFunc;  ///< Type definitions for callback functions used after evaluations of trial points proposed by user Search and Poll methods.


///
/// Template specializations (for Mads only) according to MadsCallbackType
///
template<>
class  MadsCallback<NOMAD::MadsCallbackType::USER_METHOD_POLL> : public MadsCallbackBase
{
private:
    UserPollMethodCbFunc _fn;
public:
    explicit MadsCallback(UserPollMethodCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~MadsCallback() = default;
    
    // Intermediate public callback function
    bool call(const NOMAD::Step& step,
               std::list<Direction> & dirs,
              const size_t n) const {return _fn(step, dirs, n);};
   
    NOMAD::MadsCallbackType getType() const override { return MadsCallbackType::USER_METHOD_POLL;}
    
    std::unique_ptr<MadsCallbackBase> clone() const override
    {
        return std::make_unique<MadsCallback<MadsCallbackType::USER_METHOD_POLL>>(*this);
    }
};
template<>
class  MadsCallback<NOMAD::MadsCallbackType::USER_METHOD_FREE_POLL> : public MadsCallbackBase
{
private:
    UserPollMethodCbFunc _fn;
public:
    explicit MadsCallback(UserPollMethodCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~MadsCallback() = default;
    
    // Intermediate public callback function
    bool call(const NOMAD::Step& step,
               std::list<Direction> & dirs,
              const size_t n) const {return _fn(step, dirs, n);};
    
    NOMAD::MadsCallbackType getType() const override { return MadsCallbackType::USER_METHOD_FREE_POLL;}
    
    std::unique_ptr<MadsCallbackBase> clone() const override
    {
        return std::make_unique<MadsCallback<MadsCallbackType::USER_METHOD_FREE_POLL>>(*this);
    }
    
};
template<>
class  MadsCallback<NOMAD::MadsCallbackType::USER_METHOD_FREE_POLL_END> : public MadsCallbackBase
{
private:
    UserMethodEndCbFunc _fn;
public:
    explicit MadsCallback(UserMethodEndCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~MadsCallback() = default;

    // Intermediate public callback function
    bool call(const NOMAD::Step& step) const {return _fn(step);};
    
    NOMAD::MadsCallbackType getType() const override { return MadsCallbackType::USER_METHOD_FREE_POLL_END;}
    
    std::unique_ptr<MadsCallbackBase> clone() const override
    {
        return std::make_unique<MadsCallback<MadsCallbackType::USER_METHOD_FREE_POLL_END>>(*this);
    }
    
};

template<>
class  MadsCallback<NOMAD::MadsCallbackType::USER_METHOD_SEARCH> : public MadsCallbackBase
{
private:
    UserSearchMethodCbFunc _fn;
public:
    explicit MadsCallback(UserSearchMethodCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~MadsCallback() = default;

    // Intermediate public callback function
    bool call(const NOMAD::Step& step,
              NOMAD::EvalPointSet & trialPoints) const {return _fn(step, trialPoints);};
    
    NOMAD::MadsCallbackType getType() const override { return MadsCallbackType::USER_METHOD_SEARCH;}
    
    std::unique_ptr<MadsCallbackBase> clone() const override
    {
        return std::make_unique<MadsCallback<MadsCallbackType::USER_METHOD_SEARCH>>(*this);
    }
    
};
template<>
class  MadsCallback<NOMAD::MadsCallbackType::USER_METHOD_SEARCH_END> : public MadsCallbackBase
{
private:
    UserMethodEndCbFunc _fn;
public:
    explicit MadsCallback(UserMethodEndCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~MadsCallback() = default;
    
    // Intermediate public callback function
    bool call(const NOMAD::Step& step) const {return _fn(step);};
    
    NOMAD::MadsCallbackType getType() const override { return MadsCallbackType::USER_METHOD_SEARCH_END;}
    
    std::unique_ptr<MadsCallbackBase> clone() const override
    {
        return std::make_unique<MadsCallback<MadsCallbackType::USER_METHOD_SEARCH_END>>(*this);
    }
    
};

// Factory function to create base class for MadsCallback
// Template based on callback type and function
template<NOMAD::MadsCallbackType CT, typename Fn>
std::unique_ptr<MadsCallbackBase> makeMadsCallbackBase(Fn&& fn)
{
    return std::make_unique<MadsCallback<CT>>(std::forward<Fn>(fn));
}



#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_6_MADSUSERMETHODCALLBACK__



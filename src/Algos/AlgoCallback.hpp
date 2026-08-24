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


#ifndef __NOMAD_4_6_ALGOCALLBACK__
#define __NOMAD_4_6_ALGOCALLBACK__

#include "../Type/CallbackType.hpp"

#include <memory>
#include <functional>

#include "../nomad_nsbegin.hpp"


class Step;

// Base class for algo callbacks (can be stored into a vector)
class AlgoCallbackBase
{
public:
    virtual ~AlgoCallbackBase() = default;
    
    virtual NOMAD::AlgoCallbackType getType() const = 0 ;
    
    virtual std::unique_ptr<AlgoCallbackBase> clone() const = 0;
    
};

/// Templated general class for algo callbacks
template<NOMAD::AlgoCallbackType ACT>
class  AlgoCallback : public AlgoCallbackBase
{
public:
    AlgoCallback() = default;
    virtual ~AlgoCallback() = default;
    
    NOMAD::AlgoCallbackType getType() const override { return ACT;}
    
    std::unique_ptr<AlgoCallbackBase> clone() const override
    {
        return std::make_unique<AlgoCallback<ACT>>(*this);
    }
};


//
// Define type for the iteration control
//
typedef std::function<void(const Step& step, bool & stop)> IterCbFunc;




///
/// Specialized classes for algo callbacks. Mads has its own specialized callbacks for poll/search methods.
///
template<>
class  AlgoCallback<NOMAD::AlgoCallbackType::ITERATION_END> : public AlgoCallbackBase
{
private:
    IterCbFunc _fn;
public:
    explicit AlgoCallback(IterCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~AlgoCallback() = default;
    
    // Intermediate public callback function
    void call(const NOMAD::Step& step,
              bool & stop ) const {return _fn(step, stop);};
    
    NOMAD::AlgoCallbackType getType() const override { return NOMAD::AlgoCallbackType::ITERATION_END;}
    
    std::unique_ptr<AlgoCallbackBase> clone() const override
    {
        return std::make_unique<AlgoCallback<NOMAD::AlgoCallbackType::ITERATION_END>>(*this);
    }
    
};
template<>
class  AlgoCallback<NOMAD::AlgoCallbackType::MEGA_ITERATION_START> : public AlgoCallbackBase
{
private:
    IterCbFunc _fn;
public:
    explicit AlgoCallback(IterCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~AlgoCallback() = default;
    
    // Intermediate public callback function
    void call(const NOMAD::Step& step,
              bool & stop ) const {return _fn(step, stop);};
    
    NOMAD::AlgoCallbackType getType() const override { return NOMAD::AlgoCallbackType::MEGA_ITERATION_START;}
    
    std::unique_ptr<AlgoCallbackBase> clone() const override
    {
        return std::make_unique<AlgoCallback<NOMAD::AlgoCallbackType::MEGA_ITERATION_START>>(*this);
    }
    
};
template<>
class  AlgoCallback<NOMAD::AlgoCallbackType::MEGA_ITERATION_END> : public AlgoCallbackBase
{
private:
    IterCbFunc _fn;
public:
    explicit AlgoCallback(IterCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~AlgoCallback() = default;
    
    // Intermediate public callback function
    void call(const NOMAD::Step& step,
              bool & stop ) const {return _fn(step, stop);};
    
    NOMAD::AlgoCallbackType getType() const override { return NOMAD::AlgoCallbackType::MEGA_ITERATION_END;}
    
    std::unique_ptr<AlgoCallbackBase> clone() const override
    {
        return std::make_unique<AlgoCallback<NOMAD::AlgoCallbackType::MEGA_ITERATION_END>>(*this);
    }
    
};
template<>
class  AlgoCallback<NOMAD::AlgoCallbackType::POSTPROCESSING_CHECK> : public AlgoCallbackBase
{
private:
    IterCbFunc _fn;
public:
    explicit AlgoCallback(IterCbFunc fn) : _fn(std::move(fn)) {}
    virtual ~AlgoCallback() = default;
    
    // Intermediate public callback function
    void call(const NOMAD::Step& step,
              bool & stop ) const {return _fn(step, stop);};
    
    NOMAD::AlgoCallbackType getType() const override { return NOMAD::AlgoCallbackType::POSTPROCESSING_CHECK;}
    
    std::unique_ptr<AlgoCallbackBase> clone() const override
    {
        return std::make_unique<AlgoCallback<NOMAD::AlgoCallbackType::POSTPROCESSING_CHECK>>(*this);
    }
    
};


// Factory function to create base class for AlgoCallback
// Template based on callback type and function
template<NOMAD::AlgoCallbackType CT, typename Fn>
std::unique_ptr<AlgoCallbackBase> makeAlgoCallbackBase(Fn&& fn)
{
    return std::make_unique<AlgoCallback<CT>>(std::forward<Fn>(fn));
}


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_6_ALGOCALLBACK__



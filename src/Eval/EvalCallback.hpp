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


#ifndef __NOMAD_4_6_EVALCALLBACK__
#define __NOMAD_4_6_EVALCALLBACK__

#include <memory>
#include <functional>

#include "../nomad_nsbegin.hpp"

#include "../Eval/EvalQueuePoint.hpp"

/// Abstract base class for eval callbacks
/// Derive from EvalCallback to make sure to have default call functions.
class EvalCallbackBase
{
protected:
    NOMAD::EvalCallbackType _ect = NOMAD::EvalCallbackType::COUNT;
    bool _custom = true;
public:
    EvalCallbackBase() = default;
    virtual ~EvalCallbackBase() = default;
    
    bool isCustom() const { return _custom;}
    
    NOMAD::EvalCallbackType getType() const { return _ect;} ;
    
    virtual std::unique_ptr<EvalCallbackBase> clone() const = 0;

};

template<NOMAD::EvalCallbackType CT>
class  EvalCallback : public EvalCallbackBase
{

public:
    EvalCallback() = default;
    virtual ~EvalCallback() = default;

    std::unique_ptr<EvalCallbackBase> clone() const override
    {
        return std::make_unique<EvalCallback<CT>>(*this);
    }
};

//
// Define type for the stopping control / eval check / opportunism
//
typedef std::function<void(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool& globalStop)> EvalStopCheckCbFunc;
typedef std::function<void(NOMAD::EvalQueuePointPtr & evalQueuePoint, const NOMAD::Double & hMax, bool & countEval)> PreEvalUpdateCbFunc;
typedef std::function<void(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool & opportunisticEvalStop, bool &opportunisticIterStop)> EvalOpportunisticCheckCbFunc;
typedef std::function<void(NOMAD::BlockForEval & block)> PreBlocEvalCbFunc;
typedef std::function<void(NOMAD::EvalQueuePointPtr & evalQueuePoint)> PostEvalCbFunc;


///
/// Specialized classes for evak callbacks.
///
template<>
class  EvalCallback<NOMAD::EvalCallbackType::EVAL_STOP_CHECK> : public EvalCallbackBase
{
private:
    EvalStopCheckCbFunc _fn;
public:
    explicit EvalCallback(EvalStopCheckCbFunc fn, bool custom): _fn(std::move(fn)) { _custom = custom; _ect = NOMAD::EvalCallbackType::EVAL_STOP_CHECK; }
    virtual ~EvalCallback() = default;
    
    // Intermediate public callback function
    void call(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool& globalStop ) const {return _fn(evalQueuePoint, globalStop);};

    std::unique_ptr<EvalCallbackBase> clone() const override
    {
        return std::make_unique<EvalCallback<NOMAD::EvalCallbackType::EVAL_STOP_CHECK>>(*this);
    }
};
template<>
class  EvalCallback<NOMAD::EvalCallbackType::PRE_EVAL_UPDATE> : public EvalCallbackBase
{
private:
    PreEvalUpdateCbFunc _fn;
public:
    explicit EvalCallback(PreEvalUpdateCbFunc fn, bool custom) : _fn(std::move(fn)) {_custom = custom; _ect = NOMAD::EvalCallbackType::PRE_EVAL_UPDATE; }
    virtual ~EvalCallback() = default;
    
    // Intermediate public callback function
    void call(NOMAD::EvalQueuePointPtr & evalQueuePoint, const NOMAD::Double & hMax, bool & countEval) const {return _fn(evalQueuePoint, hMax, countEval);};

    std::unique_ptr<EvalCallbackBase> clone() const override
    {
        return std::make_unique<EvalCallback<NOMAD::EvalCallbackType::PRE_EVAL_UPDATE>>(*this);
    }
};
template<>
class  EvalCallback<NOMAD::EvalCallbackType::EVAL_OPPORTUNISTIC_CHECK> : public EvalCallbackBase
{
private:
    EvalOpportunisticCheckCbFunc _fn;
public:
    explicit EvalCallback(EvalOpportunisticCheckCbFunc fn, bool custom):_fn(std::move(fn)){_custom = custom; _ect = NOMAD::EvalCallbackType::EVAL_OPPORTUNISTIC_CHECK; }
    virtual ~EvalCallback() = default;
    
    // Intermediate public callback function
    void call(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool & opportunisticEvalStop, bool &opportunisticIterStop ) const {return _fn(evalQueuePoint, opportunisticEvalStop, opportunisticIterStop);};

    std::unique_ptr<EvalCallbackBase> clone() const override
    {
        return std::make_unique<EvalCallback<NOMAD::EvalCallbackType::EVAL_OPPORTUNISTIC_CHECK>>(*this);
    }
};
template<>
class  EvalCallback<NOMAD::EvalCallbackType::PRE_EVAL_BLOCK_UPDATE> : public EvalCallbackBase
{
private:
    PreBlocEvalCbFunc _fn;
public:
    explicit EvalCallback(PreBlocEvalCbFunc fn, bool custom) : _fn(std::move(fn)) {_custom = custom ; _ect = NOMAD::EvalCallbackType::PRE_EVAL_BLOCK_UPDATE; }
    virtual ~EvalCallback() = default;
    
    // Intermediate public callback function
    void call(NOMAD::BlockForEval & block) const {return _fn(block);};
    
    std::unique_ptr<EvalCallbackBase> clone() const override
    {
        return std::make_unique<EvalCallback<NOMAD::EvalCallbackType::PRE_EVAL_BLOCK_UPDATE>>(*this);
    }
};
template<>
class  EvalCallback<NOMAD::EvalCallbackType::POST_EVAL_UPDATE> : public EvalCallbackBase
{
private:
    PostEvalCbFunc _fn;
public:
    explicit EvalCallback(PostEvalCbFunc fn, bool custom) : _fn(std::move(fn)) {_custom = custom ; _ect = NOMAD::EvalCallbackType::POST_EVAL_UPDATE ; }
    virtual ~EvalCallback() = default;
    
    // Intermediate public callback function
    void call(NOMAD::EvalQueuePointPtr & evalQueuePoint) const {return _fn(evalQueuePoint);};
    std::unique_ptr<EvalCallbackBase> clone() const override
    {
        return std::make_unique<EvalCallback<NOMAD::EvalCallbackType::POST_EVAL_UPDATE>>(*this);
    }
    
};


template<>
class  EvalCallback<NOMAD::EvalCallbackType::EVAL_FAIL_CHECK> : public EvalCallbackBase
{
private:
    PostEvalCbFunc _fn;
public:
    explicit EvalCallback(PostEvalCbFunc fn, bool custom) : _fn(std::move(fn)) {_custom = custom ; _ect = NOMAD::EvalCallbackType::EVAL_FAIL_CHECK; }
    virtual ~EvalCallback() = default;
    
    // Intermediate public callback function
    void call(NOMAD::EvalQueuePointPtr & evalQueuePoint) const {return _fn(evalQueuePoint);};
    
    std::unique_ptr<EvalCallbackBase> clone() const override
    {
        return std::make_unique<EvalCallback<NOMAD::EvalCallbackType::EVAL_FAIL_CHECK>>(*this);
    }
        
};

// Factory function to create base class for EvalCallback
// Template based on callback type and function
template<NOMAD::EvalCallbackType CT, typename Fn>
std::unique_ptr<EvalCallbackBase> makeEvalCallbackBase(Fn&& fn, bool custom)
{
    return std::make_unique<EvalCallback<CT>>(std::forward<Fn>(fn), custom);
}


// Some classes for Java interfaces. Those duplicate the template classes
class  EvalStopCheckCallback
{
public:
    explicit EvalStopCheckCallback() = default;
    virtual ~EvalStopCheckCallback() = default;
    
    // Intermediate public callback function
    virtual void call(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool& globalStop ) const = 0;

};


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_6_EVALCALLBACK__



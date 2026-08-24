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


#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/MainStep.hpp"
#include "../../Algos/MadsMixedPenaltyLogBarrier/MadsPIP.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Eval/ProgressiveBarrier.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/fileutils.hpp"
#ifdef TIME_STATS
#include "../../Util/Clock.hpp"
#endif

#include <cmath>


void NOMAD::MadsPIP::init()
{
    // Define the single-objective merit function with inequality constraints as penalty.
    std::function<NOMAD::Double(const NOMAD::BBOutputTypeList& bbOutputTypeList,
                                const NOMAD::BBOutput& bbOutput)> meritZCompute = [&](const NOMAD::BBOutputTypeList& bbOutputTypeList,
                                                                                      const NOMAD::BBOutput& bbOutput) -> NOMAD::Double
    {
        if (!bbOutput.getEvalOk() || bbOutputTypeList.empty())
        {
            return NOMAD::INF;
        }
        
        if (!bbOutput.checkSizeMatch(bbOutputTypeList))
        {
            return NOMAD::INF;
        }
        
        // IMPORTANT: the _Glog, _GextI and _GextE are indices based on ALL outputs.
        NOMAD::Double meritZ = 0;
        const auto & allO = bbOutput.getBBOAsArrayOfDouble();
        
        if (allO.size() == 1)
        {
            // This is in case the algo is used for an unconstrained pb!
            meritZ += bbOutput.getObjective(bbOutputTypeList);
            return meritZ;
        }
        
        // Initialization of Glog and Gext
        if (_Glog.empty() && _GextI.empty() && _GextE.empty())
        {
            
            // Possible improvments: use EB and handle it properly.
            // EB: constraint that should never goes into ext. If initially
            // infeasible -> exception
            
            // Initialization of indices for Glog and Gext (done once)
            for (size_t i=0 ; i < allO.size(); i++)
            {
                if (bbOutputTypeList[i] == NOMAD::BBOutputType::Type::PB)
                {
                    if (allO[i].todouble() < 0)
                    {
                        _Glog.push_back(i);
                    }
                    else
                    {
                        _GextI.push_back(i);
                    }
                }
                else if (bbOutputTypeList[i] == NOMAD::BBOutputType::Type::EQPB)
                {
                    _GextE.push_back(i);
                }
                else if (bbOutputTypeList[i] == NOMAD::BBOutputType::Type::OBJ)
                {
                    continue;
                }
                else
                {
                    throw NOMAD::Exception(__FILE__,__LINE__,"BB Output Type not handled: "+bbOutputTypeList[i].display());
                }
                
            }
        }
        
        // Initialization of rhoExt
        if (!_b_rho_cext.isDefined())
        {
            NOMAD::ArrayOfDouble Fs = bbOutput.getObjectives(bbOutputTypeList);
            if (Fs.size() > 1)
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"Cannot handle multi-objective in MadsPIP");
            }
            if (Fs[0]==0)
            {
                _b_rho_cext = 1.0;
            }
            else
            {
                // Use order of Fs instead of Fs
                int expFx = static_cast<int>(std::floor(std::log10(std::abs(Fs[0].todouble()))));
                _b_rho_cext =std::max(pow(10.0,expFx),1.0);
            }
        }
        
        // Compute the log barrier contribution to the merit function
        double logSum=0;
        for (const auto i: _Glog)
        {
            // Strict feasibility
            if (allO[i].todouble() < 0)
            {
                logSum += log(std::min(1.0,-allO[i].todouble()));
            }
            else
            {
                meritZ = NOMAD::INF;
                return meritZ;
            }
        }
        meritZ = -_b_rho_cint * _rho * logSum;
        
        for (const auto i: _GextI)
        {
            // Inequality constraint
            if (allO[i].todouble() > 0.0)
            {
                meritZ += _b_rho_cext * pow(allO[i].todouble(), 2.0)/ _rho;
            }
        }
        for (const auto i: _GextE)
        {
            // Equality constraint, NO THRESHOLD
            meritZ += _b_rho_cext * pow(allO[i].todouble(), 2.0)/_rho;
        }
        
        auto obj = bbOutput.getObjective(bbOutputTypeList);
        if (! obj.isDefined())
        {
            meritZ = NOMAD::INF;
        }
        else
        {
            meritZ += obj;
        }
        
        return meritZ;
    };
    
    // Feasibility threshold for log barrier
    _logSwitchFeasThres = _runParams->getAttributeValue<NOMAD::Double>("MADSPIP_OPTIMIZATION_SWITCH_FEASIBILITY_THRESHOLD").todouble();
    if (_logSwitchFeasThres < 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"MADSPIP_OPTIMIZATION_SWITCH_FEASIBILITY_THRESHOLD must be positive");
    }

    // Constant on rho^beta for the rho penalty parameter update
    _b_rhobeta_rho = _runParams->getAttributeValue<NOMAD::Double>("MADSPIP_OPTIMIZATION_BRHOBETA_RHO").todouble();
    
    // Make the penalized merit function compute available to the evaluator control.
    // The infeasibility H function must always return 0 (feasible)
    // Must be placed after main step start
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
    evc->setComputeType(NOMAD::ComputeType::USER, meritZCompute, [](const NOMAD::BBOutputTypeList& bbOutputTypeList,
                                                                    const NOMAD::BBOutput& bbOutput){return NOMAD::Double(0);});
    
}

///*-----------------------------------------------------*/
///* Mega Iteration start callback                       */
///* Update objective penalization weights               */
///*-----------------------------------------------------*/
void NOMAD::MadsPIP::megaIterationCallback(const NOMAD::Step& step,
                                            bool &stop)
{
    // Important: by default USER_CALLS are disabled when doing quad model optimization.

    // Several NOMAD::Algorithm are used by NOMAD.
    // We are interested only on the main Mads (Mega) Iteration.
    // Use a dynamic cast to make sure with have the Mads (Mega) Iteration.
    auto iter = dynamic_cast<const NOMAD::MadsMegaIteration*>(&step);

    // Fetch the iteration success, not the mega iteration success (because it is reset at start)
    // to get previous iteration success type.
    auto success_type = iter->getMadsIterationSuccess();

    if (nullptr != iter)
    {
        const auto & barrier = iter->getBarrier();
        const auto & bf = barrier->getCurrentIncumbentFeas();
        const auto & evalType = barrier->getEvalType();
        const auto & computeType = barrier->getFHComputeType().Short() ;

        bool rhoHasChanged = false, typeHasChanged = false;

        // Get all outputs but only inequality constraints are needed.
        // Need all output because _G indices are defined with respect to all outputs.
        const auto & allO = bf->getEval(evalType)->getBBOutput().getBBOAsArrayOfDouble();

        OUTPUT_DEBUG_START
        std::string s = "Best feasible before update: BBO (" + allO.display() + "); meritZ=" + bf->getEval(evalType)->getFs(computeType).display() + "\n" ;
        AddOutputDebug(s);
        OUTPUT_DEBUG_END

        if (success_type == SuccessType::UNSUCCESSFUL)
        {
            double comp;
            double gmax = NOMAD::M_INF; // for bf, g_l(xk)<=0 for all l in Glog 
            for (const auto i: _Glog)
            {
                double g = allO[i].todouble();
                if (g > gmax)
                {
                    gmax = g;
                }
            }

            if (gmax > NOMAD::M_INF)
            {
                comp = std::min(_b_rhobeta_rho*pow(_rho.todouble(),_beta),_b_cint_rho*pow(gmax,2.0) );
            }
            else
            {
                comp = pow(_rho.todouble(),_beta);
            }

            const auto & maxFrameSize = iter->getMesh()->getDeltaFrameSize().max();

            if ( maxFrameSize.todouble() < comp)
            {
                _rho *= _theta_rho;
                rhoHasChanged = true;
            }
        }
        {
            // Implement switch of type GextI -> Glog
            if (_extToLogSwitch)
            {
                auto it=_GextI.begin();
                while (it != _GextI.end())
                {

                    // We can switch from _Gext to _Glog if constraints becomes strictly feasible
                    if (allO[*it].todouble() < -_logSwitchFeasThres )
                    {
                        _Glog.push_back(*it);
                        it = _GextI.erase(it);
                        typeHasChanged = true;
                    }
                    else
                    {
                        it++;
                    }
                }
            }
        }

        // Reset bf
        if ( _updateBestPoint && (rhoHasChanged || typeHasChanged))
        {
            bf->resetFValues();
            auto evc = NOMAD::EvcInterface::getEvaluatorControl();
            evc->getBestIncumbent(-1)->resetFValues();

            // Update barrier with the point from the cache having the best merit
            std::function<bool(const EvalPoint&)> func = [&, evalType, computeType](const EvalPoint& evalPoint) -> bool
            {
                if(evalPoint.isEvalOk(evalType) && evalPoint.getEval(evalType)->isFeasible(computeType) && evalPoint.getEval(evalType)->getFs(computeType)[0] < bf->getEval(evalType)->getFs(computeType)[0])
                {
                    return true;
                }
                return false;
            };
            std::vector<NOMAD::EvalPoint> evalPointList;
            NOMAD::CacheInterface cacheInterface(this);
            cacheInterface.find(func, evalPointList);
            barrier->updateWithPoints(evalPointList,false, true /*update incumbent*/);
        }
    }

};

void NOMAD::MadsPIP::startImp()
{
    // Default start to perform Mads initialization (evalX0,....)
    NOMAD::Mads::startImp();

    // Set callback for mega iteration start is after barrier update and mesh update
    IterCbFunc cb = [&](const NOMAD::Step& step,
                        bool &stop)
    {
        // megaIterationCallback is a private member function that has access to all
        // member attributes.
        megaIterationCallback(step, stop);
    };
    NOMAD::Algorithm::addCallback<NOMAD::AlgoCallbackType::MEGA_ITERATION_START>(cb);
}

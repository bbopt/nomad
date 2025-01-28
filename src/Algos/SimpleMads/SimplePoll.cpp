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

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/DoublePollMethod.hpp"
#include "../../Algos/Mads/Ortho2NPollMethod.hpp"
#include "../../Algos/SimpleMads/SimplePoll.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Type/DirectionType.hpp"


/// <#Description#>
void NOMAD::SimplePoll::init()
{
    setStepType(NOMAD::StepType::SIMPLE_POLL);

    // Create the model evaluator for a block of points
    if (nullptr != _model)
    {
        _eval_x = [&](std::vector<NOMAD::SimpleEvalPoint>& trialPoints)->bool
        {
            size_t m = _trialPoints.size();
            
            // Init the matrices for prediction
            // Creation of matrix for input / output of SGTELIB model
            SGTELIB::Matrix M_predict (  "M_predict", static_cast<int>(m), static_cast<int>(_nbOutputs));
            SGTELIB::Matrix X_predict("X_predict", static_cast<int>(m), static_cast<int>(_nSimple));
            
            int j = 0;
            
            
            for (const auto & ep : _trialPoints)
            {
                
                // Set the input matrix
                size_t k = 0;
                for (size_t i = 0; i < _nSimple; i++)
                {
                    if ( _fixedVariable[i].isDefined() )
                    {
                        X_predict.set(j, static_cast<int>(i), _fixedVariable[i].todouble());
                    }
                    else
                    {
                        X_predict.set(j, static_cast<int>(i), ep[k].todouble());
                        k++;
                    }
                }
                j++;
            }

            
            _model->check_ready(__FILE__,__FUNCTION__,__LINE__);
            
            if (!_model->is_ready())
            {
                return false;
            }
            
            _model->predict(X_predict, &M_predict);
            
            j = 0;
            
            for (auto & ep : _trialPoints)
            {
                
                // ----------------------------------------- //
                //   Compute F and H from model outputs      //
                // ----------------------------------------- //
                NOMAD::ArrayOfDouble out(_bbot.size(),NOMAD::INF);
                
                for (size_t i = 0; i < _bbot.size(); i++)
                {
                    out[i] = M_predict.get(j,static_cast<int>(i));
                }
                
                ep.setF(getF(out));
                ep.setH(getH(out));
                
                j++;
            }
            return true;
        };
    }
    
        
    // The direction types for primary and secondary poll centers
    _primaryDirectionType = NOMAD::DirectionType::ORTHO_2N;
    _secondaryDirectionType = NOMAD::DirectionType::DOUBLE;
    
    // Rho parameter of the progressive barrier. Used to choose if the primary frame center is the feasible or infeasible incumbent.
    _rho = _runParams->getAttributeValue<NOMAD::Double>("RHO");
        
    _n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    _fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(_parentStep);
    _nSimple = _fixedVariable.size();
    _nbOutputs = _bbot.size();
    
    // Evaluated X0
    const auto X0s = _pbParams->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
    if (X0s.size() != 1 && !X0s[0].isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Simple Mads needs a single valid X0.");
    }
    
    NOMAD::SimpleEvalPoint evalPointX0(X0s[0]);
    _trialPoints.push_back(evalPointX0);
    evalTrialPoints(); // Compute f and h according to standard method
    
    // Create a barrier from X0. Two cases
    // 1- PhaseOne (X0 has EB constraints not feasible -> h(X0)=INF) -> first step: minimize infeasibility
    // 2- Regular (h(X0) != INF)
    if ( NOMAD::INF == _trialPoints[0].getH() )
    {

        _phaseOneSearch = true;
        
        // Recompute f and h for trial points using PhaseOne
        _trialPoints.clear();
        _trialPoints.push_back(evalPointX0);
        evalTrialPoints();
        
    }
    _barrier = std::make_unique<NOMAD::SimpleProgressiveBarrier>(NOMAD::INF,
                                                                 NOMAD::Point(_n) /*undefined fixed variables */,
                                                                 _trialPoints);

    // Create the mesh
    _mesh = std::make_shared<NOMAD::GMesh>(_pbParams, _runParams);
    
    // No need to enforce sanity check
    _mesh->setEnforceSanityChecks(false);
}

void NOMAD::SimplePoll::startImp()
{
    // Manage the case when PhaseOne is done
    // Switch to PhaseOne barrier to Regular barrier
    if ( _phaseOneSearch )
    {
        // Test if F=0
        auto xIncFeas = _barrier->getCurrentIncumbentFeas();
        
        if ( !xIncFeas.isDefined())
        {
            throw NOMAD::Exception(__FILE__,__LINE__,"Barrier for PhaseOneSearch has no feasible incumbent.");
        }
        
        // F target is reached. Switch to standard compute type
        if ( xIncFeas.getF() <= 0 )
        {
            // Recompute f and h for trial points using standard method
            _trialPoints.clear();
            _trialPoints.push_back(xIncFeas);
            _phaseOneSearch = false;
            evalTrialPoints();
            
            // Put incumbent in a fresh barrier.
            _barrier = std::make_unique<NOMAD::SimpleProgressiveBarrier>(NOMAD::INF,
                                                                         NOMAD::Point(_n) /*undefined fixed variables */,
                                                                         _trialPoints);
        }
        
    }
}


bool NOMAD::SimplePoll::runImp()
{

    // 1- Create poll methods and generate points (just one pass)
    generateTrialPoints();
    
    if (!_trialPoints.empty())
    {
        // 2- Evaluate points
        evalTrialPoints();
        
        // 3- Update barrier
        _barrier->updateWithPoints(_trialPoints);
        
        return true;
    }
    else
    {
        return false;
    }

    
}


void NOMAD::SimplePoll::endImp()
{

    // Update Mesh. Compute hMax and update Barrier.
    
    // Barrier is already updated during run.
    // Get ref best feasible and infeasible, and then update
    // reference values at the end.
    const auto & refBestFeas = _barrier->getRefBestFeas();
    const auto & refBestInf  = _barrier->getRefBestInf();

    const auto & newBestFeas = _barrier->getCurrentIncumbentFeas();
    const auto & newBestInf  = _barrier->getCurrentIncumbentInf();
    std::string s;


    if (refBestFeas.isDefined() || refBestInf.isDefined())
    {
        // Compute success
        // Get which of newBestFeas and newBestInf is improving
        // the solution. Check newBestFeas first.
        // F and H of points are available.
        NOMAD::SuccessType success = _barrier->computeSuccessType(newBestFeas, refBestFeas);
        
        NOMAD::SimpleEvalPoint newBest;
        // NOMAD::SuccessType success = computeSuccess(newBestFeas, refBestFeas, _barrier->getHMax());
        if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
        {
            // newBestFeas is the improving point.
            newBest = newBestFeas;

        }
        else
        {
            // Check newBestInf
            NOMAD::SuccessType success2 = _barrier->computeSuccessType(newBestInf, refBestInf);
            if (success2 > success)
            {
                success = success2;
            }
            if (success >= NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                // newBestInf is the improving point.
                newBest = newBestInf;
            }
        }


        // Update Mesh
        const NOMAD::SuccessType trigger = NOMAD::SuccessType::PARTIAL_SUCCESS;
        
        if (success >= trigger)
        {
            OUTPUT_INFO_START
            if (success == NOMAD::SuccessType::PARTIAL_SUCCESS)
            {
                AddOutputInfo("Last Iteration Improving. Delta remains the same.");
            }
            OUTPUT_INFO_END

            if (success >= NOMAD::SuccessType::FULL_SUCCESS)
            {
                // Use empty direction. It is not needed because the mesh is not anisotropic.
                if (_mesh->enlargeDeltaFrameSize(NOMAD::Direction()))
                {
                    OUTPUT_INFO_START
                    AddOutputInfo("Last Iteration Successful. Delta is enlarged.");
                    OUTPUT_INFO_END
                }
                else
                {
                    OUTPUT_INFO_START
                    AddOutputInfo("Last Iteration Successful. Delta remains the same.");
                    OUTPUT_INFO_END
                }
            }
        }
        else
        {
            OUTPUT_INFO_START
            AddOutputInfo("Last Iteration Unsuccessful. Delta is refined.");
            OUTPUT_INFO_END
            _mesh->refineDeltaFrameSize();
        }
    }

    OUTPUT_INFO_START
    AddOutputInfo("delta mesh  size = " + _mesh->getdeltaMeshSize().display());
    AddOutputInfo("Delta frame size = " + _mesh->getDeltaFrameSize().display());
    OUTPUT_INFO_END

    
    // Make new best the ref bests
    _barrier->updateRefBests();
    
    // Reset success
    _success = NOMAD::SuccessType::UNSUCCESSFUL;

}


// Generate new points to evaluate
void NOMAD::SimplePoll::generateTrialPoints()
{
    _trialPoints.clear();
    
    // Poll methods depend on poll centers
    createPollMethodsForPollCenters();

    // Use poll methods to create trial points
    for (const auto& pollMethod : _pollMethods)
    {
        // Generate trial points. Snap to bounds and project on mesh is also done.
        pollMethod->generateTrialPoints();

        // Pass the points from Poll method to Poll for later evaluation
        const auto & pollMethodPoints = pollMethod->getTrialPoints();
        for (const auto & point : pollMethodPoints)
        {
            _trialPoints.push_back(static_cast<SimpleEvalPoint>(point));
        }
    }
}




void NOMAD::SimplePoll::computePrimarySecondaryPollCenters(NOMAD::SimpleEvalPoint & primaryCenter,
                                                           NOMAD::SimpleEvalPoint  &secondaryCenter) const
{
    
    auto firstXIncFeas = _barrier->getCurrentIncumbentFeas();
    auto firstXIncInf  = _barrier->getCurrentIncumbentInf();
    bool primaryIsInf = false;
    
    // Negative rho means make no distinction between primary and secondary polls.
    const bool usePrimarySecondary = (_rho >= 0) && (firstXIncFeas.isDefined()) && (firstXIncInf.isDefined());
    if (usePrimarySecondary)
    {
        
        const NOMAD::Double & fFeas = firstXIncFeas.getF();
        const NOMAD::Double & fInf  = firstXIncFeas.getF();
        
        if (fFeas.isDefined() && fInf.isDefined())
        {
            // Single objective case
            if ((fFeas - _rho) > fInf)
            {
                // xFeas' f is too large, use xInf as primary poll instead.
                primaryIsInf = true;
            }
        }
    }
    
    if (usePrimarySecondary)
    {
        if (primaryIsInf)
        {
            primaryCenter=firstXIncInf;
            secondaryCenter=firstXIncFeas;
        }
        else
        {
            primaryCenter=firstXIncFeas;
            secondaryCenter=firstXIncInf;
        }
    }
    else
    {
        if (firstXIncFeas.isDefined())
        {
            primaryCenter=firstXIncFeas;
        }
        else
        {
            if (firstXIncInf.isDefined())
            {
                primaryCenter=firstXIncInf;
            }
        }
    }
}

void NOMAD::SimplePoll::createPollMethod(const bool isPrimary, const SimpleEvalPoint & frameCenter)
{

    if ( !frameCenter.isDefined())
    {
        return;
    }
    
    // Select the poll methods to be executed
    std::shared_ptr<NOMAD::PollMethodBase> pollMethod;
    auto fc = std::make_shared<NOMAD::EvalPoint>(frameCenter);
    if (isPrimary)
    {
        pollMethod = std::make_shared<NOMAD::Ortho2NPollMethod>(this, fc);
    }
    else
    {
        pollMethod = std::make_shared<NOMAD::DoublePollMethod>(this, fc);
    }
    _frameCenters.push_back(frameCenter);
    _pollMethods.push_back(pollMethod);

}


void NOMAD::SimplePoll::createPollMethodsForPollCenters()
{
    // Compute primary and secondary poll centers
    NOMAD::SimpleEvalPoint primaryCenter, secondaryCenter;
    computePrimarySecondaryPollCenters(primaryCenter, secondaryCenter);

    // Add poll methods for primary polls
    _pollMethods.clear();
    _frameCenters.clear();

    createPollMethod(true, primaryCenter);

    // Add poll methods for secondary polls
    createPollMethod(false, secondaryCenter);
}


void NOMAD::SimplePoll::evalTrialPoints()
{
    
    bool evalOk = _eval_x(_trialPoints);
    
    OUTPUT_INFO_START
    if (!evalOk)
    {
        AddOutputInfo("Eval not ok for one of the trial point ");
    }
    OUTPUT_INFO_END
    
    
    // Count eval
    _nbEval += _trialPoints.size();
}



NOMAD::Double NOMAD::SimplePoll::getF(const NOMAD::ArrayOfDouble & out) const
{
    NOMAD::Double f=0;
    
    if (out.size() != _bbot.size())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Inconsistent bbot and model output.");
    }
    
    size_t i=0;
    if (_phaseOneSearch)
    {
        bool fPos=false;
        for (const auto & bbo : _bbot)
        {
            if (bbo != NOMAD::BBOutputType::Type::EB)
            {
                i++;
                continue;
            }
            else if (out[i] > 0.0)
            {
                fPos = true;
                f += out[i]*out[i];
            }
            i++;
        }
        // Failsafe: If at least one EB constraint is positive, f is set
        // to at least epsilon.
        if (fPos && f.isDefined() && (0 == f))
        {
            f = NOMAD::Double::getEpsilon();
        }
    }
    else
    {
        NOMAD::BBOutput bboutput(out);
        f = _singleObjCompute(_bbot, bboutput);
    }
    return f;
}


/*-------------------------------------*/
/*      Get h. .      */
/*-------------------------------------*/
NOMAD::Double NOMAD::SimplePoll::getH(const NOMAD::ArrayOfDouble & out) const
{
    NOMAD::Double h=0.0;
    
    bool hPos=false;
    
    if (out.size() != _bbot.size())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Inconsistent bbot and model output.");
    }
    
    if (!_phaseOneSearch)
    {
        size_t i=0;
        for (const auto & bbo : _bbot)
        {
            if (!bbo.isConstraint())
            {
                i++;
                continue;
            }
            else if (!out[i].isDefined())
            {
                h = NOMAD::Double();    // h is undefined
                break;
            }
            else if (out[i] > 0.0)
            {
                hPos = true;
                NOMAD::Double hTemp = 0.0;
                if (bbo == NOMAD::BBOutputType::Type::EB)
                {
                    hTemp = NOMAD::INF;
                }
                else if (bbo == NOMAD::BBOutputType::Type::PB || bbo == NOMAD::BBOutputType::Type::RPB )
                {
                    hTemp = out[i] * out[i];
                }
                
                // Violated Extreme Barrier constraint:
                // Set h to infinity and break.
                if (NOMAD::INF == hTemp)
                {
                    h = NOMAD::INF;
                    break;
                }
                h += hTemp;
            }
            i++;
        }
    }
    
    // Failsafe: If at least one PB constraint is positive, h must be set
    // to at least epsilon so that the Eval is recognized as infeasible.
    // Catch cases such as constraint violated by 1e-8, which gives h = 1e-16
    // which is considered as 0.
    if (hPos && (0 == h))
    {
        h = NOMAD::Double::getEpsilon();
    }
    
    return h;
    
}



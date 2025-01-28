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
/**
 \file   ComputeSuccessType.cpp
 \brief  Comparison methods for EvalPoints (implementation)
 \author Viviane Rochon Montplaisir
 \date   March 2017 / September 2020
 \see    ComputeSuccessType.hpp
 */
#include "../Eval/ComputeSuccessType.hpp"

/*--------------------------*/
/* Class ComputeSuccessType */
/*--------------------------*/


NOMAD::SuccessType NOMAD::ComputeSuccessType::operator()(const NOMAD::EvalPointCstPtr &p1,
                                                         const NOMAD::EvalPointCstPtr &p2,
                                                         const NOMAD::Double& hMax)
{
    auto compactComputeType = _computeType.Short();
    auto computeType = compactComputeType.computeType;
    auto evalType = _computeType.evalType;
    
    if (NOMAD::EvalType::UNDEFINED == evalType)
    {
        std::string err = "Cannot compute success for undefined eval type";
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }
    
    NOMAD::SuccessType success = NOMAD::SuccessType::UNDEFINED;
    
    if (nullptr != p1)
    {
        if (nullptr == p2)
        {
            if (NOMAD::ComputeType::PHASE_ONE == computeType)
            {
                success = NOMAD::SuccessType::FULL_SUCCESS;
            }
            else if (NOMAD::ComputeType::STANDARD == computeType || NOMAD::ComputeType::DMULTI_COMBINE_F == computeType || NOMAD::ComputeType::USER == computeType)
            {
                
                NOMAD::Double h = p1->getH(_computeType);
                if (!h.isDefined() || h > hMax || h == NOMAD::INF)
                {
                    // Even if evalPoint2 is NULL, this case is still
                    // not a success.
                    success = NOMAD::SuccessType::UNSUCCESSFUL;
                }
                else if (p1->isFeasible(_computeType))
                {
                    // New feasible point: full success
                    success = NOMAD::SuccessType::FULL_SUCCESS;
                }
                else
                {
                    // New infeasible makes for partial success, not full success
                    success = NOMAD::SuccessType::PARTIAL_SUCCESS;
                }
            }
            else
            {
                std::string err = "Compute success type function for " + NOMAD::computeTypeToString(computeType) + " not available";
                throw NOMAD::Exception(__FILE__,__LINE__,err);
            }
        }
        else
        {
            success = NOMAD::Eval::computeSuccessType(p1->getEval(evalType),
                                                      p2->getEval(evalType),
                                                      compactComputeType,
                                                      hMax);
        }
    }
    
    return success;
    
    
}



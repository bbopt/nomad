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
 \file   SgtelibModelFormulationType.hpp
 \brief  types for parameter SGTELIB_MODEL_FORMULATION
 \author Viviane Rochon Montplaisir
 \date   July 2019
 \see    SgtelibModel.hpp
 */
#ifndef __NOMAD_4_0_SGTELIB_MODEL_FORMULATION_TYPE__
#define __NOMAD_4_0_SGTELIB_MODEL_FORMULATION_TYPE__

#include <string>
#include <sstream>

#include "../nomad_nsbegin.hpp"


/// Formulations for sgtelib model Search
enum class SgtelibModelFormulationType
{
    FS    ,  /// min f-lambda*sigma, st c-lambda*sigma < 0
    FSP   ,  /// min f-lambda*sigma, st P(x) > 1/2
    EIS   ,  /// min -EI-lambda*sigma, st c-lambda*sigma < 0
    EFI   ,  /// min -EFI
    EFIS  ,  /// min -EFI-lambda*sigma
    EFIM  ,  /// min -EFI-lambda*sigma*mu
    EFIC  ,  /// min -EFI-lambda*(EI*sigma+P*mu)
    PFI   ,  /// min -PFI
    D     ,  /// min -distance_to_closest
    EXTERN,  /// min f, st c, with extern executable model
    UNDEFINED /// Undefined
};


// Convert a string (ex "FS", "EIS", "FSP"...)
// to a SgtelibModelFormulationType.
SgtelibModelFormulationType stringToSgtelibModelFormulationType(const std::string &s);

std::string SgtelibModelFormulationTypeToString(const SgtelibModelFormulationType &smft);


inline std::ostream& operator<<(std::ostream& os, const SgtelibModelFormulationType &smft)
{
    switch (smft)
    {
        case SgtelibModelFormulationType::FS:
            os << "FS";
            break;
        case SgtelibModelFormulationType::FSP:
            os << "FSP";
            break;
        case SgtelibModelFormulationType::EIS:
            os << "EIS";
            break;
        case SgtelibModelFormulationType::EFI:
            os << "EFI";
            break;
        case SgtelibModelFormulationType::EFIS:
            os << "EFIS";
            break;
        case SgtelibModelFormulationType::EFIM:
            os << "EFIM";
            break;
        case SgtelibModelFormulationType::EFIC:
            os << "EFIC";
            break;
        case SgtelibModelFormulationType::PFI:
            os << "PFI";
            break;
        case SgtelibModelFormulationType::D:
            os << "D";
            break;
        case SgtelibModelFormulationType::EXTERN:
            os << "EXTERN";
            break;
        default:
            return os << "UNDEFINED";
            break;
    }

    return os;
}



#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_SGTELIB_MODEL_FORMULATION_TYPE__

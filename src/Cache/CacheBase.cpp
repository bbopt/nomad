/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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
 \file   CacheBase.cpp
 \brief  Code for base class CacheBase
 \author Viviane Rochon Montplaisir
 \date   February 2019
 \see    CacheBase.hpp
 */

#include "../Cache/CacheBase.hpp"
#include "../Util/fileutils.hpp"


// Initialize CacheBase class.
// To be called by the Constructor.
void NOMAD::CacheBase::init()
{
    // Default cache parameters are considered if cacheParams is not set
    if (nullptr == _cacheParams)
    {
        _cacheParams = std::shared_ptr<CacheParameters>(new NOMAD::CacheParameters());
    }

    _maxSize  = _cacheParams->getAttributeValue<size_t>("CACHE_SIZE_MAX") ;
    _filename = _cacheParams->getAttributeValue<std::string>("CACHE_FILE");
    // Verify filename has full path, otherwise, confusion will arise
    if (!_filename.empty() && !NOMAD::isAbsolute(_filename))
    {
        std::string err = "Error: Cache file name should have been converted to full path: ";
        err += _filename;
        std::cerr << err;
        // Be kind. Do not throw an exception from the constructor. Yet.
        //throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
}


bool isTrue(const NOMAD::EvalPoint& NOMAD_UNUSED(evalPoint))
{
    return true;
}


size_t NOMAD::CacheBase::getAllPoints(std::vector<NOMAD::EvalPoint> &evalPointList) const
{
    evalPointList.clear();
    return find(isTrue, evalPointList);
}


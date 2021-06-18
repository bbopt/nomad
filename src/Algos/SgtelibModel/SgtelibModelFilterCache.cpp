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

#include "../../Algos/CacheInterface.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelFilterCache.hpp"
#include "../../Type/SgtelibModelFormulationType.hpp"
#include "../../Output/OutputQueue.hpp"


NOMAD::SgtelibModelFilterCache::~SgtelibModelFilterCache()
{
    freeSpace();
}


void NOMAD::SgtelibModelFilterCache::init()
{
    //_name = getAlgoName() + "Search Filter";
    verifyParentNotNull();

    // Find cache points with model evaluation
    NOMAD::CacheInterface cacheInterface(this);
    size_t nbModelEval = cacheInterface.find(NOMAD::EvalPoint::hasModelEval, _cacheModelEval);

    // Initialize structures
    // Objective function (prediction)
    _f.resize(nbModelEval);
    // Aggregate constraint (prediction)
    _h.resize(nbModelEval);
    // Feasibility value (max of cj)
    _hmax.resize(nbModelEval);
    // Distance to main cache.
    _DX.resize(nbModelEval);
    // Distance between each pair of points
    _DSS.resize(nbModelEval);
    for (size_t i = 0; i < nbModelEval; i++)
    {
        _DSS[i].resize(nbModelEval);
    }
    // Initial isolation distances
    _distIsolation.resize(nbModelEval);

    _keep.resize(nbModelEval);
    _DT.resize(nbModelEval);
    _DTX.resize(nbModelEval);
    _nIsolation.resize(nbModelEval);
    _nDensity.resize(nbModelEval);
    for (size_t i = 0; i < nbModelEval; i++)
    {
        _keep[i] = false;
        _DT[i] = NOMAD::INF;
        _nIsolation[i] = -1;
        _nDensity[i] = -1;
    }

}


void NOMAD::SgtelibModelFilterCache::startImp()
{
    computeInitialValues();
}


bool NOMAD::SgtelibModelFilterCache::runImp()
{
    std::string s;  // Used for output
    size_t nbModelEval = _cacheModelEval.size();
    auto modelFormulation = _runParams->getAttributeValue<NOMAD::SgtelibModelFormulationType>("SGTELIB_MODEL_FORMULATION");

    // Get used methods from parameter
    const size_t methodMax = 10;
    bool useMethod[methodMax] = { false };
    size_t nbMethods = 0;
    std::string filterParam = _runParams->getAttributeValue<std::string>("SGTELIB_MODEL_SEARCH_FILTER");
    for (size_t i = 0; i < methodMax; i++)
    {
        if (std::string::npos != filterParam.find(NOMAD::itos(i)))
        {
            useMethod[i] = true;
            nbMethods++;
        }
    }

    // Display info
    OUTPUT_INFO_START
    s = "Used method: ";
    if (NOMAD::SgtelibModelFormulationType::D == modelFormulation)
    {
        s += "Method override. Use method 1";
    }
    else
    {
        for (size_t i = 0; i < methodMax; i++)
        {
            if (useMethod[i])
            {
                s += NOMAD::itos(i) + " ";
            }
        }
        s += "(total nb methods = " + NOMAD::itos(nbMethods) + ")";
    }
    NOMAD::OutputQueue::Add(s, _displayLevel);
    OUTPUT_INFO_END

    if (0 == nbMethods)
    {
        s = "No filter method selected, parameter SGTELIB_MODEL_SEARCH_FILTER";
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    OUTPUT_INFO_START
    s = "Filter: Start greedy selection";
    NOMAD::OutputQueue::Add(s, _displayLevel);
    OUTPUT_INFO_END

    // Find _nbCandidates points
    size_t nbKeep = 0;
    size_t nbFail = 0;  // Number of consecutive failures
    NOMAD::FilterSelectionMethod method = NOMAD::FilterSelectionMethod(0);
    while ( (nbKeep < _nbCandidates) && (nbFail < 2 * nbMethods))
    {
        // Select method.
        if (NOMAD::SgtelibModelFormulationType::D == modelFormulation)
        {
            // SGTELIB_MODEL_FORMULATION_D: Distance to closest
            method = NOMAD::FilterSelectionMethod::METHOD_MOST_DISTANT;
        }
        else
        {
            // Otherwise, cycle through all the methods
            while (true)
            {
                // Next method
                method = NOMAD::FilterSelectionMethod(size_t(method)+1);
                if (NOMAD::FilterSelectionMethod::NB_METHODS == method)
                {
                    method = NOMAD::FilterSelectionMethod(0);
                }
                if (useMethod[size_t(method)])
                {
                    break;
                }
            }
        }
        OUTPUT_INFO_START
        s = "Method " + NOMAD::FilterSelectionMethodDict.at(method);
        NOMAD::OutputQueue::Add(s, _displayLevel);
        OUTPUT_INFO_END

        int iSelect = applyMethod(method);

        // iSelect has a value >= 0 if it is the index to an interesting point.
        // _keep[i] (meaning the point is already selected) is updated.

        if ( (iSelect >= 0) && (!_keep[iSelect]))
        {
            OUTPUT_INFO_START
            s = "--> Selection of search point " + _cacheModelEval[iSelect].displayAll();
            NOMAD::OutputQueue::Add(s, _displayLevel);
            OUTPUT_INFO_END

            _keep[iSelect] = true;
            nbKeep++;
            nbFail = 0;

            // _DT and _distIsolation are updated.
            // _nIsolation is reset if needed.
            for (size_t i = 0; i < nbModelEval; i++)
            {
                if (_DT[i] > _DSS[i][iSelect])
                {
                    // Update _DT
                    _DT[i] = _DSS[i][iSelect];
                    _DTX[i] = std::min(_DTX[i], _DT[i]);
                    _nDensity[i] = -1;
                    // Update _distIsolation
                    if (_distIsolation[i] > _DT[i])
                    {
                        // If _distIsolation is updated, then _nIsolation is reset
                        _distIsolation[i] = _DT[i];
                        _nIsolation[i] = 0;
                    }
                }
            }
        }
        else
        {
            OUTPUT_INFO_START
            s = "Method " + NOMAD::FilterSelectionMethodDict.at(method) + " did not return a point";
            NOMAD::OutputQueue::Add(s, _displayLevel);
            OUTPUT_INFO_END

            nbFail++;
        }
    }

    // Update oracle points
    for (size_t i = 0; i < nbModelEval; i++)
    {
        if (_keep[i])
        {
            _oraclePoints.insert(_cacheModelEval[i]);
        }
        if (_oraclePoints.size() >= _nbCandidates)
        {
            break;
        }
    }

    OUTPUT_INFO_START
    s = "Cache filter found " + std::to_string(_oraclePoints.size()) + " points";
    AddOutputInfo(s);
    OUTPUT_INFO_END

    return true;
}


void NOMAD::SgtelibModelFilterCache::endImp()
{
}


void NOMAD::SgtelibModelFilterCache::computeInitialValues()
{
    auto modelDisplay = _runParams->getAttributeValue<std::string>("QUAD_MODEL_DISPLAY");
    _displayLevel = (std::string::npos != modelDisplay.find("F"))
                                            ? NOMAD::OutputLevel::LEVEL_INFO
                                            : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

    size_t nbModelEval = _cacheModelEval.size();
    std::string s;



    // Compute values for _f, _h, _hmax, _DX, _DTX.
    for (size_t i = 0; i < nbModelEval; i++)
    {
        NOMAD::EvalPoint x(_cacheModelEval[i]);
        _f[i] = x.getF(NOMAD::EvalType::MODEL, NOMAD::ComputeType::STANDARD).todouble();
        _h[i] = x.getH(NOMAD::EvalType::MODEL, NOMAD::ComputeType::STANDARD).todouble();

        // Compute hmax = max_j c_j for x_i
        // NOTE: this computation looks cumbersome, there may be some
        // simplifications that could be done.
        _hmax[i] = -NOMAD::INF;
        NOMAD::ArrayOfDouble bbo = x.getEval(NOMAD::EvalType::MODEL)->getBBOutput().getBBOAsArrayOfDouble();
        auto evalParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvalParams();
        const auto bbot = evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
        for (size_t j = 0; j < bbo.size(); j++)
        {
            if (NOMAD::isConstraint(bbot[j]))
            {
                _hmax[i] = std::max(_hmax[i], bbo[j].todouble());
            }
        }

        // Compute distance to main cache.
        // In cache, find a point that has a bb (non-model) eval, and which is at
        // a minimal euclidian distance from x.
        // Keep this distance DX.
        double d = NOMAD::INF;
        std::vector<NOMAD::EvalPoint> cacheBB;
        NOMAD::CacheInterface cacheInterface(this);
        cacheInterface.find(NOMAD::EvalPoint::hasBbEval, cacheBB);
        for (auto epWithBBEval : cacheBB)
        {
            d = std::min(d, NOMAD::Point::dist(x, epWithBBEval).todouble());
        }
        _DX[i] = d;
        _DTX[i] = _DX[i];
    }

    //  Compute _DSS - distance between each pair of points of S
    OUTPUT_INFO_START
    s = "Compute distances";
    NOMAD::OutputQueue::Add(s, _displayLevel);
    OUTPUT_INFO_END
    for (size_t i = 0; i < nbModelEval; i++)
    {
        _DSS[i][i] = 0;
        for (size_t j = i + 1; j < nbModelEval; j++)
        {
            _DSS[i][j] = NOMAD::Point::dist(_cacheModelEval[i], _cacheModelEval[j]).todouble();
            _DSS[j][i] = _DSS[i][j];
        }
    }

    // Compute initial isolation distances.
    // The isolation of a point i of the cache
    // is the distance to the closest point that is better than i.
    OUTPUT_INFO_START
    s = "Compute isolations";
    NOMAD::OutputQueue::Add(s, _displayLevel);
    OUTPUT_INFO_END

    for (size_t i = 0; i < nbModelEval; i++)
    {
        double d = INF;
        for (size_t j = 0; j < nbModelEval; j++)
        {
            // If the point j is better than i
            if ( (_h[j] < _h[i]) || ((_h[j] == _h[i]) && (_f[j] < _f[i])) )
            {
                d = std::min(d, _DSS[i][j]);
            }
        }
        _distIsolation[i] = d;
    }

    // Compute _hmaxThreshold
    for (size_t i = 0; i < nbModelEval; i++)
    {
        if (_hmax[i] < 0)
        {
            _hmaxThreshold = std::max(_hmaxThreshold, _hmax[i]);
        }
    }

}


void NOMAD::SgtelibModelFilterCache::freeSpace()
{
    _distIsolation.clear();
    _nIsolation.clear();
    _nDensity.clear();
    _DX.clear();
    _DT.clear();
    _DTX.clear();
    for (size_t i = 0; i < _DSS.size(); i++)
    {
        _DSS[i].clear();
    }
    _DSS.clear();

    // Delete f and h
    _f.clear();
    _h.clear();
    _hmax.clear();

    _keep.clear();

    _cacheModelEval.clear();
}


int NOMAD::SgtelibModelFilterCache::applyMethod(NOMAD::FilterSelectionMethod method)
{
    std::string s;
    size_t nbModelEval = _cacheModelEval.size();

    int iSelect = -1;
    double fmin = NOMAD::INF;
    double hmin = NOMAD::INF;
    double dmin = 0;
    double dmax = 0;
    int nmax = 0;
    double deltaMNorm = 0;
    if (_modelAlgo->getDeltaMNorm().isDefined())
    {
        deltaMNorm = _modelAlgo->getDeltaMNorm().todouble();
    }


    switch(method)
    {
        case NOMAD::FilterSelectionMethod::METHOD_BEST:
            // Select the best points
            for (size_t i = 0; i < nbModelEval; i++)
            {
                if ((!_keep[i]) && (_DTX[i] > 0))
                {
                    // Check if i is better than iSelect.
                    if ( (_h[i] < hmin) || ( (_h[i] == hmin) && (_f[i] < fmin)) )
                    {
                        hmin = _h[i];
                        fmin = _f[i];
                        iSelect = static_cast<int>(i);
                    }
                }
            }
            break;

        case NOMAD::FilterSelectionMethod::METHOD_MOST_DISTANT:
            // Special case for formulation D
            // Select the most distant point
            for (size_t i = 0; i < nbModelEval; i++)
            {
                if ( (!_keep[i]) && (_DTX[i] >= dmax) )
                {
                    dmax = _DTX[i];
                    iSelect = static_cast<int>(i);
                }
            }
            break;

        case NOMAD::FilterSelectionMethod::METHOD_BEST_MIN_DIST:
            // Select the best point but with a minimum distance to points already selected
            // dmin is 0 at this point
            OUTPUT_INFO_START
            s = "dmin = " + NOMAD::Double(dmin).display();
            NOMAD::OutputQueue::Add(s, _displayLevel);
            OUTPUT_INFO_END
            for (size_t i = 0; i < nbModelEval; i++)
            {
                if ( (!_keep[i]) && (_DTX[i] >= dmin) )
                {
                    // Check if i is better than iSelect
                    if ( (_h[i] < hmin) || ( (_h[i] == hmin) && (_f[i] < fmin)) )
                    {
                        hmin = _h[i];
                        fmin = _f[i];
                        iSelect = static_cast<int>(i);
                    }
                }
                if (-1 != iSelect)
                {
                    OUTPUT_INFO_START
                    s = "d = " + NOMAD::Double(_DTX[iSelect]).display();
                    NOMAD::OutputQueue::Add(s, _displayLevel);
                    s = "h select = " + NOMAD::Double(hmin).display();
                    NOMAD::OutputQueue::Add(s, _displayLevel);
                    OUTPUT_INFO_END

                    dmin = std::max(dmin, _DTX[iSelect] + deltaMNorm);
                }
            }
            break;

        case NOMAD::FilterSelectionMethod::METHOD_BEST_GOOD_HMAX:
            // Select the best points with a good enough value of hmax
            for (size_t i = 0; i < nbModelEval; i++)
            {
                if ( (!_keep[i]) && (_hmax[i] <= _hmaxThreshold)
                    && (_f[i] < fmin) && (_DTX[i] > deltaMNorm) )
                {
                    fmin = _f[i];
                    iSelect = static_cast<int>(i);
                }
                if (-1 == iSelect)
                {
                    _hmaxThreshold *= 0.5;
                }
                else
                {
                    _hmaxThreshold = 2.0 * _hmax[iSelect];
                }
            }

            break;

        case NOMAD::FilterSelectionMethod::METHOD_HIGHEST_ISOLATION:
            // Select point with highest isolation number

            for (size_t i = 0; i < nbModelEval; i++)
            {
                if ( (!_keep[i]) && (_distIsolation[i] > 0) )
                {
                    int ni = _nIsolation[i];
                    // If criteria undef, then compute.
                    if (-1 == ni)
                    {
                        ni = 0;
                        for (size_t j = 0; j < nbModelEval; j++)
                        {
                            if (_DSS[i][j] <= _distIsolation[i])
                            {
                                ni++;
                            }
                        }
                        _nIsolation[i] = ni;
                    }
                    // Keep biggest
                    if (ni > nmax)
                    {
                        nmax = ni;
                        iSelect = static_cast<int>(i);
                    }
                }

            }// End for
            break;

        case NOMAD::FilterSelectionMethod::METHOD_HIGHEST_DENSITY:
            // Select point with highest density number

            nmax = 0;
            for (size_t i = 0; i < nbModelEval; i++)
            {
                if ( (!_keep[i]) && (_DTX[i] > 0) )
                {
                    int ni = _nDensity[i];
                    // If criteria undef, then compute.
                    if (-1 == ni)
                    {
                        ni = 0;
                        for (size_t j = 0; j < nbModelEval; j++)
                        {
                            if ( _DSS[i][j] <= _DTX[i] )
                            {
                                ni++;
                            }
                        }
                        _nDensity[i] = ni;
                    }
                    // Keep biggest
                    if (ni > nmax)
                    {
                        nmax = ni;
                        iSelect = static_cast<int>(i);
                    }
                }
            } // End for
            break;

        default:
            s = "Method is not valid: " + NOMAD::itos(int(method));
            throw NOMAD::Exception(__FILE__, __LINE__, s);
            break;
    }

    return iSelect;
}

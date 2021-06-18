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
#ifndef __NOMAD_4_0_SGTELIB_MODEL_FILTER_CACHE__
#define __NOMAD_4_0_SGTELIB_MODEL_FILTER_CACHE__

#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Output/OutputInfo.hpp"  // for OutputLevel

#include "../../nomad_nsbegin.hpp"

// Methods available for greedy selection
enum class FilterSelectionMethod
{
    METHOD_BEST = 0,            // Method 0: Select the best point
    METHOD_MOST_DISTANT,        // Method 1: Select the most distant point
    METHOD_BEST_MIN_DIST,       // Method 2: Select the best point but with a minimum distance to points already selected
    METHOD_BEST_GOOD_HMAX,      // Method 3: Select the best points with a good enough value of hmax
    METHOD_HIGHEST_ISOLATION,   // Method 4: Select point with highest isolation number
    METHOD_HIGHEST_DENSITY,     // Method 5: Select point with highest density number
    NB_METHODS
};

// Dictionary for outputs
const std::map<FilterSelectionMethod, std::string> FilterSelectionMethodDict =
{
    {FilterSelectionMethod::METHOD_BEST,                "Select the best point"},
    {FilterSelectionMethod::METHOD_MOST_DISTANT,        "Select the most distant point"},
    {FilterSelectionMethod::METHOD_BEST_MIN_DIST,       "Select the best point but with a minimum distance to points already selected"},
    {FilterSelectionMethod::METHOD_BEST_GOOD_HMAX,      "Select the best points with a good enough value of hmax"},
    {FilterSelectionMethod::METHOD_HIGHEST_ISOLATION,   "Select point with highest isolation number"},
    {FilterSelectionMethod::METHOD_HIGHEST_DENSITY,     "Select point with highest density number"}
};



class SgtelibModelFilterCache : public Step
{
private:
    const SgtelibModel*     _modelAlgo;
    const size_t            _nbCandidates;
    EvalPointSet            _oraclePoints;  // Out - points generated by selection methods
    OutputLevel             _displayLevel;

    // Vector of EvalPoints which have a model eval
    std::vector<EvalPoint>  _cacheModelEval;

    // Structures used for filtering computations

    // Objective function (prediction)
    std::vector<double> _f;
    // Aggregate constraint (prediction)
    std::vector<double> _h;
    // Feasibility value (max of cj)
    std::vector<double> _hmax;
    // Distance to main cache.
    std::vector<double> _DX;
    // Distance between each pair of points
    std::vector<std::vector<double>> _DSS;
    // Initial isolation distances
    std::vector<double> _distIsolation;

    // Values for greedy selection
    std::vector<bool>   _keep;
    std::vector<double> _DT;
    std::vector<double> _DTX;
    std::vector<int>    _nIsolation;
    std::vector<int>    _nDensity;
    double              _hmaxThreshold;


public:
    explicit SgtelibModelFilterCache(const SgtelibModel* modelAlgo,
                                     const size_t nbCandidates)
      : Step(modelAlgo),
        _modelAlgo(modelAlgo),
        _nbCandidates(nbCandidates),
        _oraclePoints(),
        _displayLevel(OutputLevel::LEVEL_INFO),
        _cacheModelEval(0),
        _f(),
        _h(),
        _hmax(),
        _DX(),
        _DSS(),
        _distIsolation(),
        _keep(),
        _DT(),
        _DTX(),
        _nIsolation(),
        _nDensity(),
        _hmaxThreshold(-INF)
    {
        init();
    }

    virtual ~SgtelibModelFilterCache();

    // Get / Set
    EvalPointSet getOraclePoints() const { return _oraclePoints; }

private:
    void init();

    virtual void startImp() override;
    virtual bool runImp() override;
    virtual void endImp() override;

    void computeInitialValues();
    int applyMethod(FilterSelectionMethod method);
    void freeSpace();


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_SGTELIB_MODEL_FILTER_CACHE__

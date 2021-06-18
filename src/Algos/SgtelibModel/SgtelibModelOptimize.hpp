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
#ifndef __NOMAD_4_0_SGTELIB_MODEL_OPTIMIZE__
#define __NOMAD_4_0_SGTELIB_MODEL_OPTIMIZE__

#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Step.hpp"
#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Output/OutputInfo.hpp"  // for OutputLevel

#include "../../nomad_nsbegin.hpp"

class SgtelibModelOptimize : public Step
{
private:
    const SgtelibModel*                     _modelAlgo;
    OutputLevel                             _displayLevel;
    //const Point &                        _incumbent;
    //std::vector<std::shared_ptr<EvalPoint>> _x0s;
    std::shared_ptr<Mads>                   _mads;
    EvalPointSet                            _oraclePoints;

    // Reference to the original problem's RunParameters and PbParameters.
    const std::shared_ptr<RunParameters>    _refRunParams;
    const std::shared_ptr<PbParameters>     _refPbParams;

    // RunParameters and PbParameters converted for model optimization
    std::shared_ptr<RunParameters>          _optRunParams;
    std::shared_ptr<PbParameters>           _optPbParams;


public:
    /// Constructor
    // Parent must explicitely be a (pointer to a) SgtelibModel.
    // Run parameters will be recomputed for model optimization.
    explicit SgtelibModelOptimize(const SgtelibModel* modelAlgo,
                                  const std::shared_ptr<RunParameters> refRunParams,
                                  const std::shared_ptr<PbParameters>  refPbParams)
                                  //const Point& incumbent,
                                  //std::vector<std::shared_ptr<EvalPoint>> x0s,
      : Step(modelAlgo),
        _modelAlgo(modelAlgo),
        _displayLevel(OutputLevel::LEVEL_INFO),
        //_incumbent(incumbent),
        //_x0s(x0s),
        _mads(nullptr),
        _oraclePoints(),
        _refRunParams(refRunParams),
        _refPbParams(refPbParams),
        _optRunParams(nullptr),
        _optPbParams(nullptr)
    {
        init();
    }

    /*-----------*/
    /* Get / Set */
    /*-----------*/
    void setupPbParameters(const ArrayOfDouble& lowerBound,
                           const ArrayOfDouble& upperBound);

    const EvalPointSet& getOraclePoints() const { return _oraclePoints; }


private:
    void init();

    virtual void startImp() override;
    virtual bool runImp() override;
    virtual void endImp() override;

    // Helpers
    void setupRunParameters();
    void updateOraclePoints();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_SGTELIB_MODEL_OPTIMIZE__

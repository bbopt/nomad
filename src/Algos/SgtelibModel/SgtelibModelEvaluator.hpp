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

#ifndef __NOMAD400_SGTELIB_MODEL_EVALUATION__
#define __NOMAD400_SGTELIB_MODEL_EVALUATION__

#include "../../nomad_platform.hpp"
#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Eval/Evaluator.hpp"
#include "../../Type/SgtelibModelFeasibilityType.hpp"
#include "../../Output/OutputInfo.hpp"  // for OutputLevel

#include "../../nomad_nsbegin.hpp"

class SgtelibModelEvaluator : public Evaluator
{
private:
    const SgtelibModel*         _modelAlgo;
    std::string                 _modelDisplay;
    Double                      _diversification;
    SgtelibModelFeasibilityType _modelFeasibility;
    double                      _tc;
    OutputLevel                 _displayLevel;
    Point                       _fixedVariable;

public:
    /// Constructor
    explicit SgtelibModelEvaluator(const std::shared_ptr<EvalParameters>& evalParams,
                                   const SgtelibModel* modelAlgo,
                                   const std::string& modelDisplay,
                                   const Double& diversification,
                                   const SgtelibModelFeasibilityType& modelFeasibility,
                                   const double tc,
                                   const Point& fixedVariable);

    virtual ~SgtelibModelEvaluator();

    bool eval_x(EvalPoint &x,
                const Double &NOMAD_UNUSED(hMax),
                bool &countEval) const override;

private:
    void init();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIB_MODEL_EVALUATION__

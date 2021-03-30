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
#ifndef __NOMAD400_QUAD_MODEL_EVALUATION__
#define __NOMAD400_QUAD_MODEL_EVALUATION__

#include "../../Eval/Evaluator.hpp"
#include "../../Output/OutputInfo.hpp"
#include "../../../ext/sgtelib/src/Surrogate.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for evaluating trial points as EvalType::SGTE.
class QuadModelEvaluator : public Evaluator
{
private:
    const std::shared_ptr<SGTELIB::Surrogate> _model;
    std::string                _modelDisplay;
    OutputLevel                _displayLevel;
    Point                      _fixedVariable;  ///< Points are sent to evaluator in full space. Evaluator works in its own dimension. This member is used for conversions.


public:
    /// Constructor
    /**
     Usually, Evaluators work in full dimension. In this case, the models may work better in subdimension. This is why a fixed variable is used.
     */
    explicit QuadModelEvaluator(const std::shared_ptr<EvalParameters>& evalParams,
                                const std::shared_ptr<SGTELIB::Surrogate>& model,
                                const std::string& modelDisplay,
                                const Point& fixedVariable)
      : Evaluator(evalParams, EvalType::SGTE),
        _model(model),
        _modelDisplay(modelDisplay),
        _displayLevel(OutputLevel::LEVEL_INFO),
        _fixedVariable(fixedVariable)
    {
        init();
    }

    virtual ~QuadModelEvaluator();

    /**
     Points for evaluations are given in a block. Sgtelib models handle the points as a matrix and return a matrix for outputs.
     */
    std::vector<bool> eval_block(Block &block,
                                 const Double &hMax  __attribute__((unused)),
                                 std::vector<bool> &countEval) const override;

private:
    void init();


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_QUAD_MODEL_EVALUATION__

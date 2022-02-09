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
#ifndef __NOMAD_4_2_QUAD_MODEL_SLD_EVALUATION__
#define __NOMAD_4_2_QUAD_MODEL_SLD_EVALUATION__

#include "../../Eval/Evaluator.hpp"
#include "../../Output/OutputInfo.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSld.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for evaluating trial points as EvalType::MODEL.
class QuadModelSldEvaluator : public Evaluator
{
private:
    const std::shared_ptr<QuadModelSld> _model;
    std::string                _modelDisplay;
    OutputLevel                _displayLevel;
    
    int       _n;           ///< Number of variables.
    int       _nm1;         ///< Number of variables minus one.
    int       _m;           ///< Number of blackbox outputs.
    double  * _x;           ///< An evaluation point.
    double ** _alpha;       ///< Model parameters.
    bool      _model_ready; ///< \c true if model ready to evaluate.
    

public:
    /// Constructor
    /**
     Quad model evaluators work in the local full space. No need to pass the fixed variables
     */
    explicit QuadModelSldEvaluator(const std::shared_ptr<EvalParameters>& evalParams,
                                const std::shared_ptr<QuadModelSld>& model,
                                const std::string& modelDisplay)
      : Evaluator(evalParams, EvalType::MODEL),
        _model(model),
        _modelDisplay(modelDisplay),
        _displayLevel(OutputLevel::LEVEL_INFO),
        _n                  ( model->get_n()         ) ,
        _nm1                ( _n-1                  ) ,
        _m                  ( 0                     ) ,
        _x                  ( NULL                  ) ,
        _alpha              ( NULL                  ) ,
        _model_ready        ( model->check()         )
    {
        init();
    }

    virtual ~QuadModelSldEvaluator();

    /**
        Evaluation of given eval point
     */
    bool eval_x(EvalPoint &x,
                const Double& hMax,
                bool &countEval) const override;

private:
    void init();


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_2_QUAD_MODEL_SLD_EVALUATION__

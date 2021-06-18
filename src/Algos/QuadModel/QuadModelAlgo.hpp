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
#ifndef __NOMAD_4_0_QUAD_MODEL_ALGO__
#define __NOMAD_4_0_QUAD_MODEL_ALGO__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Algorithm.hpp"
#include "../../Algos/EvcInterface.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for implementation of quadratic model optimization algorithm using Bastien Talgorn's sgtelib.
/**
 * Use the start, run and end tasks. Iterate on the following sequence:
 *
 * 1- Points provided as X0s and points in cache are put in a training set.
 * 2- These points are used to build a dynamic model.
 * 3- The model is optimized. This gives oracle points.
 * 4- The oracle points are evaluated by the blackbox.
 * 5- As long as new oracle points are found, the process is repeated.
 *
 * When used by Mads SearchMethod (QuadSearchMethod):
 * - Steps 1, 2, 3 and 4 are the same.
 * - The oracle points are send back to QuadSearchMethod, which takes care
 *   of projecting them to mesh and evaluate them.
 *
 * Training set and model are stored here to allow access to other Quad classes.
 *
 */

class QuadModelAlgo: public Algorithm
{
public:
    /// Constructor
    explicit QuadModelAlgo(const Step* parentStep,
                           std::shared_ptr<AlgoStopReasons<ModelStopType>> stopReasons,
                           const std::shared_ptr<RunParameters>& runParams,
                           const std::shared_ptr<PbParameters>& pbParams)
      : Algorithm(parentStep, stopReasons, runParams, pbParams)
    {
        init();
    }

    virtual ~QuadModelAlgo();

    // Utility function to get BB_OUTPUT_TYPE parameter, which is buried in Evaluator.
    static BBOutputTypeList getBBOutputType()
    {
        if (nullptr == EvcInterface::getEvaluatorControl()
            || nullptr == EvcInterface::getEvaluatorControl()->getEvalParams())
        {
            throw Exception(__FILE__, __LINE__, "Error in QuadModel::getBBOutputType()");
        }
        return EvcInterface::getEvaluatorControl()->getEvalParams()->getAttributeValue<BBOutputTypeList>("BB_OUTPUT_TYPE");
    }

    void readInformationForHotRestart() override {}

private:
    void init();

    void startImp() override;
    bool runImp() override;
    void endImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_QUAD_MODEL_ALGO__


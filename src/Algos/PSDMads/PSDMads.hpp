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
#ifndef __NOMAD400_PSDMADS__
#define __NOMAD400_PSDMADS__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Algorithm.hpp"
#include "../../Eval/Evaluator.hpp"
#include "../../Math/RandomPickup.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for PSD-Mads
class PSDMads: public Algorithm
{
private:
    RandomPickup                _randomPickup;          ///< To manage selection of variables for subproblems. Reset only by pollster.
    std::shared_ptr<MeshBase>   _psdMainMesh;           ///< Base Mesh to create subproblem Mads. Updated only by pollster.
    std::shared_ptr<Barrier>    _barrier;               ///< Barrier with the latest successful values. Updated by all Mads.

    std::atomic<bool>           _lastMadsSuccessful;    ///< Used as indication to enlarge or refine the mesh. Updated by all Mads.

    static omp_lock_t           _psdMadsLock;           ///< Lock access to the previous elements when they are updated.

public:
    /// Constructor
    /**
     \param parentStep    The parent of this step -- \b IN.
     \param evaluator     Evaluator to initialize all main threads -- \b IN.
     \param evalContParams Parameters to initialize all main threads -- \b IN.
     \param stopReasons   The PSD Mads stop reasons -- \b IN/OUT.
     \param runParams     Parameters for algorithm -- \b IN.
     \param refPbParams   Parameters for original optimization problem. PSD-Mads use its own copy -- \b IN.
     */
    explicit PSDMads(const Step* parentStep,
                     const std::shared_ptr<Evaluator>& evaluator,
                     const std::shared_ptr<EvaluatorControlParameters>& evalContParams,
                     std::shared_ptr<AlgoStopReasons<MadsStopType>> stopReasons,
                     const std::shared_ptr<RunParameters>& runParams,
                     const std::shared_ptr<PbParameters>& refPbParams)
      : Algorithm(parentStep, stopReasons, runParams, std::make_shared<PbParameters>(*refPbParams)),
        _randomPickup(_pbParams->getAttributeValue<size_t>("DIMENSION")),
        _psdMainMesh(nullptr),
        _barrier(nullptr),
        _lastMadsSuccessful(false)
    {
        init(evaluator, evalContParams);
    }

    virtual ~PSDMads()
    {
        destroy();
    }

    virtual void startImp() override;
    virtual bool runImp() override;
    virtual void endImp() override;

    void readInformationForHotRestart() override {}

    void setupSubproblemParams(std::shared_ptr<PbParameters> &subProblemPbParams,
                               std::shared_ptr<RunParameters> &subProblemRunParams,
                               const Point& bestPoint,
                               const bool isPollster);

private:
    /// Helper for constructor
    void init(const std::shared_ptr<Evaluator>& evaluator,
              const std::shared_ptr<EvaluatorControlParameters>& evalContParams);
    /// Helper for destructor
    void destroy();

    /// Wait for _barrier to be initialized by pollster before running a worker.
    void waitForBarrier() const;

    /// Assess if it is time to update _mainMesh
    bool doUpdateMesh() const;

    /// Setup fixedVariable to define subproblem
    void generateSubproblem(Point &fixedVariable);
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_PSDMADS__

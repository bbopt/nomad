/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
#ifndef __NOMAD_4_4_DISCOMADSMEGAITERATION__
#define __NOMAD_4_4_DISCOMADSMEGAITERATION__

#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/DiscoMads/DiscoMadsIteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the iterations of DiscoMADS
/**
Manager for Mads iterations.

*/
class DiscoMadsMegaIteration: public MadsMegaIteration
{

private:
    // vector of indices of revealing outputs
    std::vector<int> _idxRevealingOutput;

    // DiscoMads parameters  
    NOMAD::Double  _exclusionRadius;   // used to penalize revealed regions 
    NOMAD::Double  _detectionRadius;   // detectionRadius and limitRate are used to reveal discontinuities
    NOMAD::Double  _limitRate;   

    bool _detectHiddConst;               // if true, discoMads used to reveal hidden constraints instead of discontinuities      
    NOMAD::Double _hiddConstOutputValue;  // only used to reveal hidden constraints regions

    bool _isRevealing ;  // Flag for indicating if MegaIteration is revealing. Reset at start, after DiscoMadsUpdate because it uses the flag.

public:
    /// Constructor
    /**
     \param parentStep      The parent step of this step -- \b IN.
     \param k               The main iteration counter -- \b IN.
     \param barrier         The barrier for constraints handling -- \b IN.
     \param mesh            Mesh on which other Iteration meshes are based -- \b IN.
     \param success         Success type of the previous MegaIteration. -- \b IN.
     */
    explicit DiscoMadsMegaIteration(const Step* parentStep,
                                  size_t k,
                                  std::shared_ptr<BarrierBase> barrier,
                                  MeshBasePtr mesh,
                                  SuccessType success)
      : MadsMegaIteration(parentStep, k, barrier, mesh, success),
        _isRevealing(false)
    {
        // Replace MadsIteration by DiscoMadsIteration
        _madsIteration = std::make_unique<NOMAD::DiscoMadsIteration>(this, k, mesh);
        init();
    }
    // No Destructor needed - keep defaults.

    // Access to the revealing status of the MegaIteration run.
    bool isRevealing() const { return _isRevealing; }


private:
    /// Implementation of the start tasks for DiscoMads mega iteration.
    void startImp() override ;

    /// Implementation of the run tasks for DiscoMads mega iteration.
    /**
     Manages the generation of points: either all poll and search points are generated all together before starting evaluation using the MegaSearchPoll or they are generated using a MadsIteration with search and poll separately. A run parameter controls the behavior.
     */
    virtual bool runImp() override;

    /// Implementation of the end tasks for DiscoMads mega iteration.
    /**
    Increment number of revealing iterations if required and calls MegaIteration::endImp()
     */
    virtual void endImp() override;

    void init();

    // Test presence of weak discontinuity between x1 and x2
    bool discontinuityTest(const NOMAD::EvalPoint & x1, const NOMAD::EvalPoint & x2);

    // Test if x1 and x2 are at distance less than exclusion radius and that x2 is revealing
    bool proximityTestOnRevealingPoint(const NOMAD::Point & x1, const NOMAD::EvalPoint & x2);



    // Callback attached to evaluator: test if evalQueuePoint is a revealing point (for discontinuities or hidden constraints) and update its revealed constraint
    void callbackCheckIfRevealingAndUpdate(EvalQueuePointPtr & evalQueuePoint);

    // Callback attached to evaluator: put high value of objective function f for failed eval
    /**
    Only used if DiscoMads is used to reveal hidden constraints regions)
    */
    void callbackFailedEval(EvalQueuePointPtr & evalQueuePoint);

    // Callback attached to evaluator: check after each evaluation if there was a revealation and set opportunisticStop to True
    void callbackEvalOpportStop(bool &opportunisticStop, EvalQueuePointPtr & evalQueuePoint );

    // Callback attached to postProcessing : if there has been a revelation, stop current step and stop megaiteration without doing remaining eval
    void callbackPostProcessing(const NOMAD::Step & step, bool &stop);

    // Only use for debug: save a special text file with cache information at the end of each megaiteration
    void exportCache(std::string cacheFile);

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_DISCOMADSMEGAITERATION__

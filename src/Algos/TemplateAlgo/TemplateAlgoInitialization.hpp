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
#ifndef __NOMAD_4_4_TEMPLATEALGOINITIALIZATION__
#define __NOMAD_4_4_TEMPLATEALGOINITIALIZATION__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Initialization.hpp"
#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgoInitialization.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for template algo initialization
/**
 * Evaluate initial trial points (x0).
 * If x0 is not provided simply pass.
 */
class TemplateAlgoInitialization: public Initialization, public IterationUtils
{
private:

    std::shared_ptr<AlgoStopReasons<RandomAlgoStopType>> _templateAlgoStopReason;

public:
    /// Constructor
    /*
     \param parentStep      The parent of this step -- \b IN.
     */
    explicit TemplateAlgoInitialization(const Step* parentStep)
      : Initialization(parentStep),
        IterationUtils(parentStep)
    {
        init();
    }

    /// Destructor
    virtual ~TemplateAlgoInitialization() {}


private:
    /// Helper for constructor
    void init();

    /// Implementation of start task
    /**
     If needed, randomly generate trial points and put them in cache.
     For a standalone optimization (RANDOM_ALGO_OPTIMIZATION true), initial trial point can be provided in x0 or we simply pass.
     */
    virtual void startImp() override ;

    /// Implementation of run task
    /**
     For a standalone template algorithm, evaluate the trial points generated during start.
     Otherwise, there are no trial points available and a failed stop reason is set.
     */
    virtual bool runImp() override ;

    // Update _evalPointList member with evaluated trial points for future use
    void endImp() override;

    /// Generate initial trial point randomly
    void generateTrialPointsImp() override;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_TEMPLATEALGOINITIALIZATION__

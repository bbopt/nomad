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

#ifndef __NOMAD_4_0_ALGORITHM__
#define __NOMAD_4_0_ALGORITHM__

#include "../Algos/Initialization.hpp"
#include "../Algos/MegaIteration.hpp"
#include "../Algos/Step.hpp"
#include "../Algos/Termination.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../nomad_nsbegin.hpp"


/// Generic class for any direct search optimizer algorithm
/**
  \note: AllParameters and EvaluatorControl are held by MainStep.
 \note: Cache is a singleton all by itself.
 \note MegaIteration holds the algorithm-related structures.
 */
class Algorithm: public Step
{
protected:

    std::unique_ptr<Initialization>  _initialization;   ///< To initialize the algorithm (X0)
    std::unique_ptr<Termination>     _termination;      ///< To verify termination conditions
    std::shared_ptr<MegaIteration>   _megaIteration;    ///< MegaIteration used to keep information between steps

    bool                             _endDisplay;

    SuccessType                      _algoBestSuccess ; ///< The best succes type of the algorithm (cannot always get this information from _megaIteration).

    bool                             _algoSuccessful;

#ifdef TIME_STATS
    size_t _totalRealAlgoTime;
    double _startTime;
    double _totalCPUAlgoTime;
#endif // TIME_STATS

public:
    /// Constructor
    /**
     \param parentStep          The parent of this Step -- \b IN.
     \param stopReasons         The stop reasons of this algo -- \b IN.
     \param runParams           The run parameters that control the algorithm -- \b IN.
     \param pbParams            The problem parameters that control the algorithm -- \b IN.
     */
    explicit Algorithm(const Step* parentStep,
                       std::shared_ptr<AllStopReasons> stopReasons,
                       const std::shared_ptr<RunParameters>& runParams,
                       const std::shared_ptr<PbParameters>& pbParams )
      : Step(parentStep, stopReasons, runParams, pbParams),
        _initialization(nullptr),
        _termination(nullptr),
        _megaIteration(nullptr),
        _endDisplay(true)
#ifdef TIME_STATS
        ,_totalRealAlgoTime(0),
        _startTime(0.0),
        _totalCPUAlgoTime(0.0)
#endif // TIME_STATS
    {
        init();
    }

    /// Destructor
    virtual ~Algorithm();

    /*---------*/
    /* Get/Set */
    /*---------*/
    const std::shared_ptr<MegaIteration>& getMegaIteration() const { return _megaIteration; }
    void setMegaIteration(const std::shared_ptr<MegaIteration> megaIteration) { _megaIteration = megaIteration; }

    void setEndDisplay( bool endDisplay ) { _endDisplay = endDisplay; }

protected:


    /// Default implementation of the start tasks of an algorithm
    /**
     If doing a hot restart get the algorithm ready to continue. \n
     If starting a new algorithm, reset the stop reason, the lap evaluation counter, and perform initialization.
     */
    virtual void startImp() override;

    /// Default implementation of the end tasks of an algorithm
    /**
     Display some information, reset the lap counters, set success type for a search method and save information for a potential hot restart.
     */
    virtual void endImp() override;

    /// Each algorithm must implement its run tasks.
    /**
     Run algorithm execution for single-objective.
     \return \c true
     */
    virtual bool runImp() override = 0;

    /// Helper for start() when doing a hot restart.
    virtual void readInformationForHotRestart() = 0;

    /// Helper for end()
    void saveInformationForHotRestart() const;
    /// Helper for end()
    void displayBestSolutions() const;
    /// Helper for end()
    void displayEvalCounts() const;

    /// Helper for hot restart
    void hotRestartOnUserInterrupt() override;

public:
    /**
     Sub-algo: an algorithm can be part of an algorithm.
     */
    bool isSubAlgo() const;
    bool isRootAlgo() const { return !isSubAlgo(); }

    /*---------*/
    /* Others  */
    /*---------*/
    /// Verify if this Algorithm is ready to be terminated
    bool terminate(size_t iteration);

    virtual void read(std::istream& is);
    virtual void display(std::ostream& os) const;

private:

    ///  Helper for Constructor.
    void init();

};

/// Operator to write parameters used for hot restart.
std::ostream& operator<<(std::ostream& os, const Algorithm& algo);

/// Operator to read parameters used for hot restart.
std::istream& operator>>(std::istream& is, Algorithm& algo);


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_ALGORITHM__

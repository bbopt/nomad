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
/**
  \file   MainStep.hpp
  \brief  Main Step to hold MADS, or other Algorithms
  \author Viviane Rochon Montplaisir
  \date   June 2018
*/
#ifndef __NOMAD_4_0_MAINSTEP__
#define __NOMAD_4_0_MAINSTEP__

#include "../Algos/Algorithm.hpp"
#include "../Eval/Evaluator.hpp"
#include "../Param/AllParameters.hpp"
#include "../Type/LHSearchType.hpp"

#include "../nomad_nsbegin.hpp"

/// Main step to manage an algorithm execution.
/**
* MainStep is the root step of an executable running an algorithm.
* Think of it as a wrapper around the algorithm.
* MainStep manages all algorithm needs: Parameters, Evaluator
* and has utility functions (ex. displayUsage).
* MainStep takes care of the OpenMP parallelism.
* An algorithm can call other algorithms during its execution.
*/
class MainStep: public Step
{
protected:
    std::string                         _paramFileName;  ///< Name of the file containing the parameters.
    std::shared_ptr<AllParameters>      _allParams;
    std::shared_ptr<Evaluator>          _evaluator; ///< Used in library running mode (not batch mode)
    std::vector<std::shared_ptr<Algorithm>>  _algos;

#ifdef TIME_STATS
    size_t _totalRealTime;
    double _startTime;
    double _totalCPUTime;
#endif // TIME_STATS

public:
    /// Constructor
    explicit MainStep()
    : Step(),
        _paramFileName(""),
        _evaluator(nullptr),
        _algos()
#ifdef TIME_STATS
        ,_totalRealTime(0),
        _startTime(0.0),
        _totalCPUTime(0.0)
#endif // TIME_STATS
    {
        init();
    }

    /// Destructor
    virtual ~MainStep();

    /*---------*/
    /* Get/Set */
    /*---------*/

    /**
     In batch mode: Set the parameter file name, which will be read in start().
     In library mode: Set the parameters directly using set_PARAM_NAME(...). In library mode, it is also possible to read a parameter file.
     */
    void setParamFileName(const std::string& paramFileName) { _paramFileName = paramFileName;}

    void setAllParameters(const std::shared_ptr<AllParameters> &allParams);

    /**
     The evaluator may be shared between main threads.
     */
    void setEvaluator(std::shared_ptr<Evaluator> ev) { _evaluator = ev;}


    /*---------*/
    /* Others  */
    /*---------*/

    /**
     Once all algos have been added are executed (call start, run and end for each algo) in the MainStep::run.
     */
    void addAlgo(const std::shared_ptr<Algorithm> algo) { _algos.push_back(algo); }

    void clearAlgos() { _algos.clear(); }

    /// Helper function called by the code main function if necessary.
    void displayUsage(const char* exeName);

    /// Helper function called by the code main function if necessary.
    void displayVersion();

    /// Helper function called by the code main function if necessary.
    void displayInfo();

    /// Helper function called by the code main function if necessary.
    void displayHelp(const std::string& helpSubject = "all", bool devHelp = false);

    /**
     The user has requested a hot restart. Update the parameters with the changes requested by the user (read file or set inline).
     */
    void hotRestartOnUserInterrupt() override;

    /// Helper to reset some components (used by the runner when running multiple optimization)
    static void resetComponentsBetweenOptimization();

    /// Helper to reset the cache
    static void resetCache();

    NOMAD::ArrayOfPoint suggest() override;

    /**
      Observe method updates cache, computes new mesh size and new hMax.
      */
    void observe(const std::vector<NOMAD::EvalPoint>& evalPointList) override;
    /**
      Observe version to be called by the Python interface.
      \return new values of key parameters.
      */
    std::vector<std::string> observe(const NOMAD::ArrayOfPoint & xs, const std::vector<NOMAD::ArrayOfDouble> & fxs, const std::string & destinationCacheFileName="");

protected:
    /// Specific implementation to start NOMAD
    /**
     Implementation called by Step::start.
     During the main step start, the parameters are read (if parameter file is available), the cache, the evaluator and the evaluator control are created. If an algorithm is set for run in the parameter file, it is added (MainStep::addAlgo).
     */
    virtual void startImp() override;

    /**
     Implementation called by Step::run.
     Once all algos have been added are executed (call start, run and end for each algo) in the MainStep::run. \n

     If a stop reason (not the default STARTED) is propagated to the MainStep, the sequence of algos is stopped.
     */
    virtual bool runImp() override;

    /**
     Implementation called by Step::end.
     */
    virtual void endImp() override;

    /// Helper function when creating the evaluator
    int getNumThreads() const;

    /// Set the number of threads for the evaluator
    void setNumThreads() const;

    void printNumThreads() const;

    /// Detect if a Phase One search is required
    /**
     A phase one search is required if an EB type constraint
     is not feasible for X0.
     */
    bool detectPhaseOne();

    /// Helper for start
    void createCache() const;

    /// Helper for start
    void updateX0sFromCache() const;

    /// Helper for start
    ArrayOfPoint suggestFromLH(const size_t nbPoints) const;

private:
    /// Helper for constructor
    void init();

    ///  Detailed stats
    void displayDetailedStats() const;

};


#include "../nomad_nsend.hpp"


#endif // __NOMAD_4_0_MAINSTEP__

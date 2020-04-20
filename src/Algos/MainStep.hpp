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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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
#ifndef __NOMAD400_MAINSTEP__
#define __NOMAD400_MAINSTEP__

#include "Step.hpp"
#include "Subproblem.hpp"

#include "../Cache/CacheBase.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Param/AllParameters.hpp"

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
private:
    std::string                         _paramFileName;  ///< Name of the file containing the parameters.
    std::shared_ptr<AllParameters>      _allParams;
    std::unique_ptr<Evaluator>          _evaluator; ///< Used in library running mode (not batch mode)
    std::vector<std::shared_ptr<Step>>  _algos;
    static std::string                  _algoComment;   ///< Comment to appear in the stats, e.g. "Phase One"
    static std::vector<std::string>     _prevAlgoComment; ///< Pile of previous comments, used when going back to the main algo after running a sub-algo.  
    static bool                         _forceAlgoComment; ///< When true, do not change comment until reset is called

    /// Subproblem definition
    /**
     \todo Generalize / Make an array to handle multiple subproblems
     */
    std::vector<Subproblem>             _subproblems;


public:
    /// Constructor
    explicit MainStep()
    : Step( ),
        _paramFileName(""),
        _evaluator(nullptr),
        _algos(),
        _subproblems()
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
     The evaluator is unique. All algorithms must use the same evaluator.
     */
    void setEvaluator(std::unique_ptr<Evaluator> ev) { _evaluator = std::move(ev);}

    /// Set comment to be added at the end of the display stats, e.g., "Phase One"
    static void setAlgoComment(const std::string& algoComment, const bool force = false);
    /// Reset comment to the previous in the stack
    static void resetPreviousAlgoComment(const bool force = false);
    /// Get current comment on the top of the stack
    static const std::string& getAlgoComment() { return _algoComment; }

    /// Get the subproblem that is currently treated.
    /**
     Currently, there is only one.
     \todo Generalize subproblem management. Make it possible to have more than one subproblem defined.
     */
    std::shared_ptr<Subproblem> getCurrentSubproblem() const;


    /*---------*/
    /* Others  */
    /*---------*/

    /**
     Once all algos have been added are executed (call start, run and end for each algo) in the MainStep::run.
     */
    void addAlgo(const std::shared_ptr<Step> algo) { _algos.push_back(algo); }

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


private:

    /// Helper for constructor
    void init();
    
    /// Specific implementation to start NOMAD
    /**
     Implementation called by Step::start.
     During the main step start, the parameters are read (if parameter file is available), the cache, the evaluator and the evaluator control are created. If an algorithm is set for run in the parameter file, it is added ( MainStep::addAlgo ).
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
    virtual void endImp() override {}

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

};


#include "../nomad_nsend.hpp"


#endif // __NOMAD400_MAINSTEP__

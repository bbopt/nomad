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

#ifndef __NOMAD_4_5_MADS__
#define __NOMAD_4_5_MADS__

#include "../../Algos/Algorithm.hpp"
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

typedef std::function<bool(const Step& step, std::list<Direction> & dir, const size_t n)> UserPollMethodCbFunc;  ///< Type definitions for callback functions for user Poll method.
typedef std::function<bool(const Step& step, EvalPointSet & trialPoint)> UserSearchMethodCbFunc;  ///< Type definitions for callback functions for user Search method.
typedef std::function<bool(const Step& step)> UserMethodEndCbFunc;  ///< Type definitions for callback functions used after evaluations of trial points proposed by user Search and Poll methods.


/// The (M)esh (A)daptive (D)irect (S)earch algorithm.
/**
\note AllParameters and EvaluatorControl are held by MainStep.
Cache is a singleton all by itself.
MegaIteration holds the algorithm-related structures: Mesh, Barrier.
 */
class DLL_ALGO_API Mads: public Algorithm
{
private:

    static UserSearchMethodCbFunc       _cbUserSearchMethod;
    static UserSearchMethodCbFunc       _cbUserSearchMethod_2;
    static UserMethodEndCbFunc          _cbUserSearchMethodEnd;
    static UserPollMethodCbFunc         _cbUserPollMethod;
    static UserPollMethodCbFunc         _cbUserFreePollMethod;
    static UserMethodEndCbFunc          _cbUserFreePollMethodEnd;

    // Flags for user method callbacks.
    // Flags are set to true when adding callback. This is done only if USER_CALLS_ENABLED==true.
    bool _hasUserSearchMethod, _hasUserPollMethod, _hasUserFreePollMethod;


private:
    std::vector<std::pair<std::size_t,std::shared_ptr<SearchMethodBase>>> _extraSearchMethods;

public:
    /// Constructor
    /**
     \param parentStep          The parent of this step -- \b IN.
     \param stopReasons         The stop reasons for MADS -- \b IN.
     \param runParams           The run parameters that control MADS -- \b IN.
     \param pbParams            The problem parameters that control MADS -- \b IN.
     \param barrierInitializedFromCache  Flag to initialize barrier from cache or not -- \b IN.
     \param useOnlyLocalFixedVariables   Flag to use only local fixed variables or not (not the ones from the original problem) -- \b IN.
     */
    explicit Mads(const Step* parentStep,
                  std::shared_ptr<AlgoStopReasons<MadsStopType>> stopReasons,
                  const std::shared_ptr<RunParameters>& runParams,
                  const std::shared_ptr<PbParameters>& pbParams,
                  bool barrierInitializedFromCache = true,
                  bool useOnlyLocalFixedVariables = false )
      : Algorithm(parentStep, stopReasons, runParams, pbParams, useOnlyLocalFixedVariables),
    _hasUserSearchMethod(false),
    _hasUserPollMethod(false),
    _hasUserFreePollMethod(false)
    {
        init(barrierInitializedFromCache);
    }

    /// Helper for hot restart
    void hotRestartOnUserInterrupt() override;

    /// For suggest and observe PyNomad interface
    NOMAD::ArrayOfPoint suggest() override;
    void observe(const std::vector<NOMAD::EvalPoint>& evalPointList) override;
    
    
    /// Insert extra search methods. To be accesses
    void insertSearchMethod(size_t pos, const std::shared_ptr<SearchMethodBase>& searchMethod)
    {
        _extraSearchMethods.push_back(std::pair<size_t,std::shared_ptr<SearchMethodBase>>(pos,searchMethod));
    }
    
    std::vector<std::pair<std::size_t,std::shared_ptr<SearchMethodBase>>> & accessExtraSearchMethods()
    {
        return _extraSearchMethods;
    }

    /// \brief Set user method callback
    DLL_ALGO_API void addCallback(const CallbackType& callbackType,
                     const UserPollMethodCbFunc& userPollCbFunc);
	DLL_ALGO_API void addCallback(const CallbackType& callbackType,
                     const UserSearchMethodCbFunc& userSearchCbFunc);
    /// \brief Set user method post eval callback
	DLL_ALGO_API void addCallback(const CallbackType& callbackType,
                     const UserMethodEndCbFunc& userCbFunc) const;

    /// \brief Run user poll method callback to produce direction
    bool runCallback(const CallbackType& callbackType,
                     const Step& step,
                     std::list<Direction> & dir,
                     const size_t n) const;

    /// \brief Run user search method callback to produce trial points
    bool runCallback(const CallbackType& callbackType,
                     const Step& step,
                     EvalPointSet & trialPoints) const;

    /// \brief Run user method post eval callback to produce direction
    bool runCallback(const CallbackType& callbackType,
                     const Step& step) const;

    bool hasUserSearchMethod() const {return _hasUserSearchMethod;}
    bool hasUserPollMethod() const {return _hasUserPollMethod;}
    bool hasUserFreePollMethod() const {return _hasUserFreePollMethod;}

private:
    ///  Initialization of class, to be used by Constructor.
    /**
    \param barrierInitializedFromCache  Flag to initialized barrier from cache or not -- \b IN.
    */
    void init(bool barrierInitializedFromCache);

    /// Algorithm execution for single-objective.
    /**
     Overrides the default algorithm's run
     \return \c true if a full success was found, \c false otherwise
     */
    virtual bool runImp() override;

    /// Helper for start()
    void readInformationForHotRestart() override;
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_MADS__

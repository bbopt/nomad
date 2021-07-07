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
#ifndef __NOMAD_4_0_STOPREASON__
#define __NOMAD_4_0_STOPREASON__

#include <map>
#include "../Util/Exception.hpp"

#include "../nomad_nsbegin.hpp"


// All enums must begin with STARTED and be ended by LAST

/// Stop type that can happen any time
enum class BaseStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    MAX_TIME_REACHED        ,  ///< Max time
    INITIALIZATION_FAILED   ,  ///< Initialization failed to complete
    ERROR                   ,  ///< Error
    UNKNOWN_STOP_REASON     ,  ///< Unknown
    CTRL_C                  ,  ///< Ctrl-C
    HOT_RESTART             ,  ///< Hot restart interruption
    USER_STOPPED            ,  ///< User-stopped in a callback function
    LAST
};

/// Stop type that can happend during MADS
enum class MadsStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    MESH_PREC_REACHED       ,  ///< Mesh minimum precision stop criterion
    MIN_MESH_SIZE_REACHED   ,  ///< Min mesh size stop criterion
    MIN_FRAME_SIZE_REACHED   , ///< Min frame size stop criterion
    X0_FAIL                 ,  ///< Problem with starting point evaluation
    PONE_SEARCH_FAILED      ,  ///< Phase one search did not return a feasible point.
    LAST
};


/// Stop type that happen during the Phase One search of Mads (sub-algo)
enum class PhaseOneStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    NO_FEAS_PT              ,  ///< No feasible solution obtained during PhaseOne search
    MADS_FAIL               ,  ///< Mads fail
    LAST
};

/// Stop type that happen during SSD-Mads (super-algo)
enum class SSDMadsStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    X0_FAIL                 ,  ///< Problem with starting point evaluation
    SUBPB_MADS_FAIL         ,  ///< Subproblem Mads fail
    LAST
};

/// Stop type for Latin Hypercube
enum class LHStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    NO_POINTS_GENERATED     ,  ///< No points generated by Latin Hypercube
    ALL_POINTS_EVALUATED    ,  ///< No more points to evaluate
    LAST
};


/// Stop type for all model based searches (sgtelib or quad) or optimization
enum class ModelStopType : int
{
    STARTED                  ,  ///< Started
    ORACLE_FAIL              ,  ///< Oracle failed generating points
    MODEL_OPTIMIZATION_FAIL  ,  ///< Model optimization has failed
    INITIAL_FAIL             ,  ///< Cannot initialize model
    NOT_ENOUGH_POINTS        ,  ///< Not enough points to build model
    NO_NEW_POINTS_FOUND      ,  ///< Models optimization did not find new points
    EVAL_FAIL                ,  ///< Problem with Sgtelib evaluation
    X0_FAIL                  ,  ///< Problem with starting point evaluation
    MODEL_SINGLE_PASS_COMPLETED , ///< A single pass has been completed
    ALL_POINTS_EVALUATED     ,  ///< No more points to evaluate
    LAST
};


/// Stop type for Nelder Mead
/**
 \todo check the stop type
 */
enum class NMStopType : int
{
    STARTED                     ,  ///< Started (no stop)
    TOO_SMALL_SIMPLEX           ,  ///< Not used
    SIMPLEX_RANK_INSUFFICIENT   ,  ///< Not used
    INITIAL_FAILED              ,  ///< No valid point during initialization
    REFLECT_FAILED              ,
    EXPANSION_FAILED            ,
    OUTSIDE_CONTRACTION_FAILED  ,
    INSIDE_CONTRACTION_FAILED   ,
    SHRINK_FAILED               ,
    UNDEFINED_STEP              ,
    INSERTION_FAILED            ,
    X0_FAILED                   ,
    NM_SINGLE_COMPLETED         ,
    NM_STOP_ON_SUCCESS          ,
    NM_STOP_NO_SHRINK           ,
    LAST
};

/// Stop type for VNS
enum class VNSStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    X0_FAILED               ,  ///< Problem with starting point evaluation
    INITIAL_FAILED          ,  ///<  Pb during initialization
    SUBPB_MADS_FAILED       ,  ///< Subproblem Mads failed
    SHAKING_FAILED          ,  ///< Shaking incumbents failed
    SINGLE_PASS_COMPLETED   ,  ///< A single pass has been completed
    LAST
};


/// Stop type that can happen during evaluation (global conditions)
enum class EvalGlobalStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    MAX_BB_EVAL_REACHED     ,  ///< Max number of blackbox evaluations
    MAX_SURROGATE_EVAL_OPTIMIZATION_REACHED,///< Max number of static surrogate evaluations
    MAX_EVAL_REACHED        ,  ///< Max number of total evaluations
    MAX_BLOCK_EVAL_REACHED  ,  ///< Max number of block eval reached
    LAST
};


/// Stop type that can happen during evaluation (conditions for a main thread)
enum class EvalMainThreadStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    LAP_MAX_BB_EVAL_REACHED,   ///< Max number of blackbox evaluations for a sub algorithm run (lap run)
    SUBPROBLEM_MAX_BB_EVAL_REACHED,   ///< Max number of blackbox evaluations for a subproblem run (E.g. SSD-Mads)
    OPPORTUNISTIC_SUCCESS   ,  ///< Success found and opportunistic strategy is used
    EMPTY_LIST_OF_POINTS    ,  ///< Tried to eval an empty list
    ALL_POINTS_EVALUATED    ,  ///< No more points to evaluate
    MAX_MODEL_EVAL_REACHED  ,  ///< Max number of quad or sgtelib model evaluations
    LAST
};


/// Stop type that can happen at the end of an iteration
enum class IterStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    MAX_ITER_REACHED        ,  ///< Max number of iterations
    STOP_ON_FEAS            ,  ///< Stop because a feasible point is reached
    PHASE_ONE_COMPLETED     ,  ///< Stop because PhaseOne is done
    LAST
};

/// Template class for the stop reason of a stop type.
/**
 The possible stop types can be generic: ::IterStopType,
 ::EvalGlobalStopType, ::EvalMainThreadStopType, BaseStopType, or specific to
 an algorithm: ::MadsStopType, ::PhaseOneStopType, ::NMStopType ,....

 The default stop type is STARTED (no stop). A stop reason different than STARTED indicates what is the cause of termination. Some stop reasons indicate a normal termination that do not need propagation and others must be propagated to stop an algorithm (see StopReason::checkTerminate()).

 */
template<typename T >
class StopReason
{

private:

    T _stopReason;  ///< The stop reason stored as a stop type

    /// Dictionnary to translate a stop type into a string.
    /**
     We have template specializations of this function for each stop type.
     This function is called to display the stop reason.
     */
    std::map<T,std::string> & dict() const;

    /// Helper for constructor (check sanity)
    void testValidity() const
    {
        if ( ! std::is_enum<T>::value )
            throw Exception(__FILE__,__LINE__,"The templated stop reason is not an enum.");

        if ( dict().size() == 0 )
            throw Exception(__FILE__,__LINE__,"Dictionary not filled.");


        if ( (int) T::STARTED != 0 )
            throw Exception(__FILE__,__LINE__,"First StopType in enum must be STARTED.");

        if ( (int) T::LAST != dict().size() )
        {
            std::string s = "Not enough elements in enum dictionary (";
            s += std::to_string(dict().size()) + "), expecting " + std::to_string((int)T::LAST);
            throw Exception(__FILE__,__LINE__,s);
        }

        for (int i = (int) T::STARTED; i < (int ) T::LAST; i++)
        {
            typename std::map<T,std::string>::iterator it = dict().find( (T) i );

            if ( it== dict().end() )
                throw Exception(__FILE__,__LINE__,"All enum elements must be in dictionary.");

        }

    }

public:
    /// Constructor
    /**
     Upon construction, the validity of the stop reason is verified (sanity check). \n
     By default the stop reason is set to STARTED.
     */
    explicit StopReason ()
    {
        testValidity();

        _stopReason = T::STARTED;

    }

    /// Destructor
    virtual ~StopReason() {}

    /// The stop reason
    T get() const
    {
        return _stopReason;
    }

    /// Set the stop reason if it is listed in dictionnary
    void set (T s)
    {
        typename std::map<T,std::string>::iterator it = dict().find(s);

        if ( it==dict().end() )
            throw Exception(__FILE__,__LINE__,"Stop reason not found.");

        _stopReason = s;
    }

    /// Reset the stop reason to the default STARTED state.
    void setStarted ()
    {
        _stopReason = T::STARTED;
    }

    /// Check if it is in STARTED state.
    bool isStarted () const
    {
        return ( _stopReason == T::STARTED );
    }

    /// Translate the stop reason into a string for display
    std::string getStopReasonAsString() const
    {
        typename std::map<T,std::string>::iterator it = dict().find( _stopReason );
        return it->second;
    }

    /// Check if the stop reason requires a termination.
    /**
     This is implemented for each stop type (template specialization of the function). \n
     Except for ::EvalMainThreadStopType, a stop reason different than STARTED indicates that an algorithm or a sub-algorithm must terminate.

     \return \c true if a termination is required, \c false otherwise.
     */
    bool checkTerminate () const ;

};


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_STOPREASON__

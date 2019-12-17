/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
#ifndef __NOMAD400_STOPREASON__
#define __NOMAD400_STOPREASON__

#include <map>
#include "../Util/Exception.hpp"

#include "../nomad_nsbegin.hpp"


// All enums must begin with STARTED and be ended by LAST

/// Stop type that can happen any time
enum class BaseStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    MAX_TIME_REACHED        ,  ///< Max time
    INITIALIZATION_FAILED   ,  ///< Initilization failed to complete
    ERROR                   ,  ///< Error
    UNKNOWN_STOP_REASON     ,  ///< Unknown
    CTRL_C                  ,  ///< Ctrl-C
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

/// Stop type for Latin Hypercube
enum class LHStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    NO_POINTS_GENERATED     ,  ///< No points generated by Latin Hypercube
    ALL_POINTS_EVALUATED    ,  ///< No more points to evaluate
    LAST
};

/// Stop type for Sgtelib model optimization
enum class SgtelibModelStopType : int
{
    STARTED                 ,  ///< Started
    ORACLE_FAIL             ,  ///< Oracle failed generating points
    MODEL_OPTIMIZER_FAIL    ,  ///< Model optimizer failed
    NO_POINTS               ,  ///< No points to build model
    NO_NEW_POINTS_FOUND     ,  ///< Models optimization did not find new points
    EVAL_FAIL               ,  ///< Problem with Sgtelib evaluation
    X0_FAIL                 ,  ///< Problem with starting point evaluation
    ALL_POINTS_EVALUATED    ,  ///< No more points to evaluate
    LAST
};

/// Stop type for Nelder Mead
/**
 \todo check the stop type
 */

enum class NMStopType : int
{
    STARTED                     ,  ///< Started (no stop)
    TOO_SMALL_SIMPLEX           ,
    SIMPLEX_RANK_INSUFFICIENT   ,
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
    LAST
};

/// Stop type that can happen during evaluation
enum class EvalStopType : int
{
    STARTED                 ,  ///< Started (no stop)
    MAX_BB_EVAL_REACHED     ,  ///< Max number of blackbox evaluations
    LAP_MAX_BB_EVAL_REACHED,   ///< Max number of blackbox evaluations for a sub algorithm run (lap run)
    MAX_EVAL_REACHED        ,  ///< Max number of total evaluations
    OPPORTUNISTIC_SUCCESS   ,  ///< Success found and opportunistic strategy is used
    EMPTY_LIST_OF_POINTS    ,  ///< Tried to eval an empty list
    ALL_POINTS_EVALUATED    ,  ///< No more points to evaluate
    MAX_BLOCK_EVAL_REACHED  ,  ///< Max number of block eval reached
    MAX_SGTE_EVAL_REACHED   ,  ///< Max number of surrogate evaluations
    LAST
};

/// Stop type that can happen at the end of an iteration
enum class IterStopType : int
{
    STARTED                  ,  ///< Started (no stop)
    MAX_ITER_REACHED         ,   ///< Max number of iterations
    STOP_ON_FEAS            ,  ///< Stop because a feasible point is reached.
    LAST
};

/// Template class for the stop reason of a stop type.
/**
 The possible stop types can be generic: ::IterStopType,
 ::EvalStopType, BaseStopType, or specific to
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
     Except for ::EvalStopType, a stop reason different than STARTED indicates that an algorithm or a sub-algorithm must terminate.

     \return \c true if a termination is required, \c false otherwise.
     */
    bool checkTerminate () const ;

};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_STOPREASON__

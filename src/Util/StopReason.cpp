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
#include "../Util/StopReason.hpp"

//// Dictionary funct for BaseStopType
template<> std::map<NOMAD::BaseStopType,std::string> & NOMAD::StopReason<NOMAD::BaseStopType>::dict() const
{
    static std::map<NOMAD::BaseStopType,std::string> dictionary= {
        {NOMAD::BaseStopType::STARTED,"Started"},   // Set a the begining of a Step
        {NOMAD::BaseStopType::MAX_TIME_REACHED,"Maximum allowed time reached"},
        {NOMAD::BaseStopType::INITIALIZATION_FAILED,"Initialization failure"},
        {NOMAD::BaseStopType::ERROR,"Error"},
        {NOMAD::BaseStopType::UNKNOWN_STOP_REASON,"Unknown"},
        {NOMAD::BaseStopType::CTRL_C,"Ctrl-C"},
        {NOMAD::BaseStopType::USER_STOPPED,"User-stopped in a callback function"}

    };
    return dictionary;
}

// Returns true only to terminate an algorithm (Mads, ...)
template<> bool NOMAD::StopReason<NOMAD::BaseStopType>::checkTerminate() const
{
    switch ( _stopReason )
    {
        case NOMAD::BaseStopType::MAX_TIME_REACHED:
        case NOMAD::BaseStopType::INITIALIZATION_FAILED:
        case NOMAD::BaseStopType::ERROR:
        case NOMAD::BaseStopType::UNKNOWN_STOP_REASON:
        case NOMAD::BaseStopType::CTRL_C:
        case NOMAD::BaseStopType::USER_STOPPED:
            return true;
            break;
        case NOMAD::BaseStopType::STARTED:
            return false;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All stop types must be checked for algo terminate");
    }
    return false;
}

// Dictionary function for MadsStopType
template<> std::map<NOMAD::MadsStopType,std::string> & NOMAD::StopReason<NOMAD::MadsStopType>::dict() const
{
    static std::map<NOMAD::MadsStopType,std::string> dictionary = {
        {NOMAD::MadsStopType::STARTED,"Started"},   // Set at the begining of a Step
        {NOMAD::MadsStopType::MESH_PREC_REACHED,"Mesh minimum precision reached"},
        {NOMAD::MadsStopType::MIN_MESH_SIZE_REACHED,"Min mesh size reached"},
        {NOMAD::MadsStopType::MIN_FRAME_SIZE_REACHED,"Min frame size reached"},
        {NOMAD::MadsStopType::PONE_SEARCH_FAILED,"Phase one search did not return a feasible point."},
        {NOMAD::MadsStopType::X0_FAIL,"Problem with starting point evaluation"}
    };
    return dictionary;
}


// Returns true only to terminate an algorithm (Mads, ...)
template<> bool NOMAD::StopReason<NOMAD::MadsStopType>::checkTerminate() const
{
    switch ( _stopReason )
    {
        case NOMAD::MadsStopType::MESH_PREC_REACHED:
        case NOMAD::MadsStopType::MIN_MESH_SIZE_REACHED:
        case NOMAD::MadsStopType::MIN_FRAME_SIZE_REACHED:
        case NOMAD::MadsStopType::PONE_SEARCH_FAILED:
        case NOMAD::MadsStopType::X0_FAIL:
            return true;
            break;
        case NOMAD::MadsStopType::STARTED:
            return false;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All stop types must be checked for algo terminate");
    }
    return false;
}

// Dictionary function for PhaseOneStopType
template<> std::map<NOMAD::PhaseOneStopType,std::string> & NOMAD::StopReason<NOMAD::PhaseOneStopType>::dict() const
{
    static std::map<NOMAD::PhaseOneStopType,std::string> dictionary = {
        {NOMAD::PhaseOneStopType::STARTED,"Started"},   // Set a the begining of a Step
        {NOMAD::PhaseOneStopType::NO_FEAS_PT,"No feasible point obtained by PhaseOne search"},
        {NOMAD::PhaseOneStopType::MADS_FAIL,"Mads has terminated but no feasible point obtained"}
    };
    return dictionary;
}


// Returns true only to terminate an algorithm using PhaseOne (Mads, ...)
template<> bool NOMAD::StopReason<NOMAD::PhaseOneStopType>::checkTerminate() const
{
    switch ( _stopReason )
    {
        case NOMAD::PhaseOneStopType::NO_FEAS_PT:   //
        case NOMAD::PhaseOneStopType::MADS_FAIL:   //
            return true;
            break;
        case NOMAD::PhaseOneStopType::STARTED:
            return false;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All stop types must be checked for terminate");
    }
    return false;
}

// Dictionary function for LHStopType
template<> std::map<NOMAD::LHStopType,std::string> & NOMAD::StopReason<NOMAD::LHStopType>::dict() const
{
    static std::map<NOMAD::LHStopType,std::string> dictionary = {
        {NOMAD::LHStopType::STARTED,"Started"},   // Set a the begining of an EvaluatorControl Run
        {NOMAD::LHStopType::NO_POINTS_GENERATED, "No points generated by Latin Hypercube"},
        {NOMAD::LHStopType::ALL_POINTS_EVALUATED,"No more points to evaluate"}
    };
    return dictionary;
}


// Returns true when Latin Hypercube Sampling is complete, or no points generated
template<> bool NOMAD::StopReason<NOMAD::LHStopType>::checkTerminate() const
{
    switch ( _stopReason )
    {
        case NOMAD::LHStopType::ALL_POINTS_EVALUATED:
        case NOMAD::LHStopType::NO_POINTS_GENERATED:
            return true;
            break;
        case NOMAD::LHStopType::STARTED:
            return false;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All stop types must be checked for algo terminate");
    }
    return false;
}

// Dictionary function for NMStopType
template<> std::map<NOMAD::NMStopType,std::string> & NOMAD::StopReason<NOMAD::NMStopType>::dict() const
{
    static std::map<NOMAD::NMStopType,std::string> dictionary = {
        {NOMAD::NMStopType::STARTED,"Started"},   // Set a the begining of an EvaluatorControl Run
        {NOMAD::NMStopType::TOO_SMALL_SIMPLEX, "Simplex Y is too small"},
        {NOMAD::NMStopType::SIMPLEX_RANK_INSUFFICIENT, "Rank of the matrix DZ is too small"},
        {NOMAD::NMStopType::INITIAL_FAILED,"Initialization has failed"},
        {NOMAD::NMStopType::REFLECT_FAILED,"Reflect step has failed"}              ,
        {NOMAD::NMStopType::EXPANSION_FAILED,"Expansion step has failed"}            ,
        {NOMAD::NMStopType::OUTSIDE_CONTRACTION_FAILED,"Outside conctraction step has failed"}  ,
        {NOMAD::NMStopType::INSIDE_CONTRACTION_FAILED,"Inside contraction step failed"}   ,
        {NOMAD::NMStopType::SHRINK_FAILED,"Shrink step has failed"}               ,
        {NOMAD::NMStopType::UNDEFINED_STEP,"Unknown step"}              ,
        {NOMAD::NMStopType::INSERTION_FAILED,"Insertion of points has failed"}            ,
        {NOMAD::NMStopType::X0_FAILED,"No X0 provided or cannot evaluate X0"},
        {NOMAD::NMStopType::NM_SINGLE_COMPLETED,"NM with a single iteration is completed"},
        {NOMAD::NMStopType::NM_STOP_ON_SUCCESS,"NM iterations stopped on eval success"}
    };
    return dictionary;
}


// Returns true only when Nelder Mead algorithm is complete
template<> bool NOMAD::StopReason<NOMAD::NMStopType>::checkTerminate() const
{

    switch ( _stopReason )
    {
        case NOMAD::NMStopType::INITIAL_FAILED:
        case NOMAD::NMStopType::X0_FAILED:
        case NOMAD::NMStopType::REFLECT_FAILED:
        case NOMAD::NMStopType::EXPANSION_FAILED:
        case NOMAD::NMStopType::INSIDE_CONTRACTION_FAILED:
        case NOMAD::NMStopType::OUTSIDE_CONTRACTION_FAILED:
        case NOMAD::NMStopType::SHRINK_FAILED:
        case NOMAD::NMStopType::INSERTION_FAILED:
        case NOMAD::NMStopType::NM_SINGLE_COMPLETED:
        case NOMAD::NMStopType::NM_STOP_ON_SUCCESS:
        case NOMAD::NMStopType::UNDEFINED_STEP:
            return true;
            break;
        case NOMAD::NMStopType::STARTED:
            return false;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All NM stop types must be checked for algo terminate");
    }
    return false;
}


// Dictionary function for IterStopType
template<> std::map<NOMAD::IterStopType,std::string> & NOMAD::StopReason<NOMAD::IterStopType>::dict() const
{
    static std::map<NOMAD::IterStopType,std::string> dictionary = {
        {NOMAD::IterStopType::STARTED,"Started"},   // Set a the begining of a Step task
        {NOMAD::IterStopType::MAX_ITER_REACHED,"Maximum number of iterations reached"},
        {NOMAD::IterStopType::STOP_ON_FEAS,"A feasible point is reached"}
    };
    return dictionary;
}


// Returns true only to terminate an algorithm (Mads, ...)
template<> bool NOMAD::StopReason<NOMAD::IterStopType>::checkTerminate() const
{
    switch ( _stopReason )
    {
        case NOMAD::IterStopType::MAX_ITER_REACHED:
        case NOMAD::IterStopType::STOP_ON_FEAS:
            return true;
            break;
        case NOMAD::IterStopType::STARTED:
            return false;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All stop types must be checked for algo terminate");

    }
    return false;
}

// Dictionary function for EvalStopType
template<> std::map<NOMAD::EvalStopType,std::string> & NOMAD::StopReason<NOMAD::EvalStopType>::dict() const
{
    static std::map<NOMAD::EvalStopType,std::string> dictionary = {
        {NOMAD::EvalStopType::STARTED,                  "Started"},   // Set a the begining of an EvaluatorControl Run
        {NOMAD::EvalStopType::MAX_BB_EVAL_REACHED,      "Max number of blackbox evaluations"},
        {NOMAD::EvalStopType::LAP_MAX_BB_EVAL_REACHED,  "Max number of blackbox evaluations for a sub algorithm run (lap run)"},
        {NOMAD::EvalStopType::MAX_EVAL_REACHED,         "Max number of total evaluations"},
        {NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS,    "Success found and opportunistic strategy is used"},
        {NOMAD::EvalStopType::EMPTY_LIST_OF_POINTS,     "Tried to eval an empty list"},
        {NOMAD::EvalStopType::ALL_POINTS_EVALUATED,     "No more points to evaluate"},
        {NOMAD::EvalStopType::MAX_BLOCK_EVAL_REACHED,   "Maximum number of block eval reached"},
        {NOMAD::EvalStopType::MAX_SGTE_EVAL_REACHED,    "Max number of surrogate evaluations reached"}
    };
    return dictionary;
}

template<> bool NOMAD::StopReason<NOMAD::EvalStopType>::checkTerminate() const
{
    // Returns true only to terminate an algorithm (Mads, ...)
    switch ( _stopReason )
    {
        case NOMAD::EvalStopType::MAX_BB_EVAL_REACHED:
        case NOMAD::EvalStopType::LAP_MAX_BB_EVAL_REACHED:
        case NOMAD::EvalStopType::MAX_EVAL_REACHED:
        case NOMAD::EvalStopType::MAX_BLOCK_EVAL_REACHED:
            return true;
            break;
        case NOMAD::EvalStopType::STARTED:
        case NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS:
        case NOMAD::EvalStopType::EMPTY_LIST_OF_POINTS:
        case NOMAD::EvalStopType::ALL_POINTS_EVALUATED:
        case NOMAD::EvalStopType::MAX_SGTE_EVAL_REACHED:
            return false;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All stop types must be checked for algo terminate");
    }
    return false;
}


// Dictionary function for SgtelibModelStopType
template<> std::map<NOMAD::SgtelibModelStopType,std::string> & NOMAD::StopReason<NOMAD::SgtelibModelStopType>::dict() const
{
    static std::map<NOMAD::SgtelibModelStopType,std::string> dictionary = {
        {NOMAD::SgtelibModelStopType::STARTED, "Started"},   // Set a the begining of a Step
        {NOMAD::SgtelibModelStopType::ORACLE_FAIL, "Oracle failed generating points"},

        {NOMAD::SgtelibModelStopType::MODEL_OPTIMIZER_FAIL, "Model Optimizer has failed"},
        {NOMAD::SgtelibModelStopType::NO_POINTS, "No points to build model"},
        {NOMAD::SgtelibModelStopType::NO_NEW_POINTS_FOUND, "Models optimization did not find new points"},
        {NOMAD::SgtelibModelStopType::EVAL_FAIL, "Problem with Sgtelib Model evaluation"},
        {NOMAD::SgtelibModelStopType::X0_FAIL, "Problem with starting point evaluation"},
        {NOMAD::SgtelibModelStopType::ALL_POINTS_EVALUATED,"No more points to evaluate"}
    };
    return dictionary;
}


// Returns true only to terminate an SgtelibModel algorithm
template<> bool NOMAD::StopReason<NOMAD::SgtelibModelStopType>::checkTerminate() const
{
    switch ( _stopReason )
    {
        case NOMAD::SgtelibModelStopType::STARTED:
        case NOMAD::SgtelibModelStopType::ALL_POINTS_EVALUATED:
            return false;
            break;
        case NOMAD::SgtelibModelStopType::MODEL_OPTIMIZER_FAIL:
        case NOMAD::SgtelibModelStopType::NO_POINTS:
        case NOMAD::SgtelibModelStopType::NO_NEW_POINTS_FOUND:
        case NOMAD::SgtelibModelStopType::EVAL_FAIL:
        case NOMAD::SgtelibModelStopType::X0_FAIL:
            return true;
            break;
        default:
            NOMAD::Exception ( __FILE__, __LINE__,"All stop types must be checked for algo terminate");
    }
    return false;
}

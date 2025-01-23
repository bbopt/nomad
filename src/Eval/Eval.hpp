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
/**
 \file   Eval.hpp
 \brief  Evaluation at a point
 \author Viviane Rochon Montplaisir
 \date   March 2017
 \see    Eval.cpp
 */

#ifndef __NOMAD_4_5_EVAL__
#define __NOMAD_4_5_EVAL__

#include <functional>   // For std::function

#include "../Eval/BBOutput.hpp"
#include "../Param/EvalParameters.hpp"
#include "../Type/ComputeType.hpp"
#include "../Type/CompareType.hpp"

#include "../nomad_nsbegin.hpp"

/**
* Idea to implement:
* Eval differentiates between a BB evaluation that failed
 (failure of the black box) and an evaluation that was
 interrupted (failure cause is external to the black box). \n
* There is also states for evaluations that are not yet started,
* and evaluations that are currently running.
*/


/// Type for an evaluation status
enum class EvalStatusType
{
    EVAL_NOT_STARTED,       ///< Evaluation has not been done yet. Initial status.
    EVAL_FAILED,            ///< Evaluation failure. Do not re-submit.
    EVAL_ERROR,             ///< Evaluation did not proceed normally. May be submitted again.
    EVAL_USER_REJECTED,     ///< Evaluation was rejected by user. May be submitted again
    EVAL_USER_ACCEPTED,     ///< Evaluation was accepted by user. The status will be updated upon evaluation.
    EVAL_OK,                ///< Correct evaluation
    EVAL_IN_PROGRESS,       ///< Evaluation in progress
    EVAL_WAIT,              ///< Evaluation in progress for another instance of the same point: Wait for evaluation to be done.
    EVAL_STATUS_UNDEFINED   ///< Undefined evaluation status
};


/// Utility to convert an eval status to a string.
DLL_EVAL_API std::string enumStr(const EvalStatusType evalStatus);

/// Struct containing specific information for multiobjective optimization
struct MOInfo
{
    ArrayOfDouble fvalues;
    ArrayOfDouble intermediateVal;
    Double combineFValue;
};


/// Class for the representation of an evaluation at a point.
/**
 * \note We have a separate Eval from EvalPoint so that a Point can have multiple evaluations.
 */
class DLL_EVAL_API Eval {

private:
    EvalStatusType _evalStatus;         ///< The evaluation status.
    EvalStatusType _preEvalStatus;      ///< The pre-evaluation status. Used by user to reject or accept a point.
    BBOutput _bbOutput;                 ///<  The blackbox evaluation output.
    BBOutputTypeList _bbOutputTypeList; ///< List of output types: OBJ, PB, EB etc.
    bool _bbOutputComplete;             ///< All bbo outputs have a valid value for functions (OBJ, PB and EB).
    std::unique_ptr<MOInfo> _moInfo; ///< Multiobjective information; precomputed to have more performance
    
public:

    /*---------------*/
    /* Class Methods */
    /*---------------*/

    /// Constructor #1.
    explicit Eval();

    /// Constructor #2.
    /**
     \param params      The parameters for defining the behavior of the class -- \b IN.
     \param bbOutput    The output of the blackbox -- \b IN.
    */
    explicit Eval(const std::shared_ptr<EvalParameters>& params,
                  BBOutput bbOutput);

    /// Copy constructor.
    /**
     \param eval The copied object .
     */
    Eval(const Eval& eval);

    /// Affectation operator.
    /**
     \param eval The right-hand side object -- \b IN.
     \return           \c *this as the result of the affectation.
     */
    Eval& operator= (const Eval& eval);

    // Destructor.
    // Needed to avoid valgrind memory leak warnings.
    virtual ~Eval() {}

    /*---------*/
    /* Get/Set */
    /*---------*/

    // f and h are always recomputed.
    Double getF(const NOMAD::FHComputeTypeS& fhComputeType) const;
    const ArrayOfDouble& getFs(const NOMAD::FHComputeTypeS& fhComputeType) const;
    Double getH(const NOMAD::FHComputeTypeS& fhComputeType) const;

    EvalStatusType getEvalStatus() const { return _evalStatus; }
    EvalStatusType getPreEvalStatus() const { return _preEvalStatus; }
    void setEvalStatus(const EvalStatusType &evalStatus) { _evalStatus = evalStatus; }
    void setPreEvalStatus(const EvalStatusType &evalStatus) { _preEvalStatus = evalStatus; }
    bool isEvalOk () const { return _evalStatus == EvalStatusType::EVAL_OK; }

    bool isBBOutputComplete() const { return _bbOutputComplete; }
    BBOutput getBBOutput() const { return _bbOutput; }
    ArrayOfDouble getBBOutputByType( const BBOutputType & bboType );
    BBOutputTypeList getBBOutputTypeList() const { return _bbOutputTypeList; }
    void setBBOutputTypeList(const BBOutputTypeList& bbOutputTypeList)
    {
        _bbOutputTypeList = bbOutputTypeList;
        _bbOutputComplete = _bbOutput.isComplete(bbOutputTypeList);
    }

    std::string getBBO() const { return _bbOutput.getBBO(); }

    /// Set blackbox output
    void setBBO(const std::string &bbo,
                const BBOutputTypeList &bbOutputTypeList,
                const bool evalOk = true);

    /*---------------*/
    /* Other methods */
    /*---------------*/
    
    
    /// Assess if the point is feasible.
    /**
     \param fhComputeType How to compute f and h -- \b IN.
     \return                    \c true if no constraint is broken, \c false otherwise.
     */
    bool isFeasible(const NOMAD::FHComputeTypeS& fhComputeType) const;

    /// Can this point be re-evaluated? Based on the eval status only.
    bool canBeReEvaluated() const;

    /** Should this point be saved to cache file? Based on the eval status only.
     * These eval statuses are good: EVAL_OK, EVAL_FAILED, EVAL_USER_REJECTED,
     * EVAL_ERROR.
     * These eval statuses are not good:
     * EVAL_NOT_STARTED, EVAL_IN_PROGRESS, EVAL_WAIT, EVAL_STATUS_UNDEFINED.
    */
    bool goodForCacheFile() const;

    /*------------*/
    /* Comparison */
    /*------------*/
    /// Comparison operator \c ==.
    /**
     \param eval   The right-hand side object -- \b IN.
     \return       \c true if  \c *this \c == \c e, \c false if not.
     */
    bool operator==(const Eval &eval) const;

    /// Comparison operator \c !=.
    /**
     \param eval   The right-hand side object -- \b IN.
     \return       \c true if  \c *this \c != \c e, \c false if not.
     */
    bool operator!= (const Eval &eval) const { return !(*this == eval); }

    /// Comparison operator \c < for dominance
    /**
     \param eval    The right-hand side object -- \b IN.
     \return        \c true if  \c *this \c < \c e, \c false if not.
     \note          This operator must not be used to find/store EvalPoints in the cache or in any set.
     */
    bool operator<(const Eval& eval) const;

    /// Dominance
    /**
     Dominance as by definition 12.3 in the Book of Audet and Hare, Derivative-Free and Blackbox Optimization,
     https://doi.org/10.1007/978-3-319-68913-5
     * The feasible point x dominates the feasible point y \n
       when f(x) < f(y).
     * The infeasible point x dominates the infeasible point y \n
       when f(x) <= f(y) and h(x) <= h(y),
       with at least one strict inequality.
     * Otherwise, return false
     \param eval The right-hand side object -- \b IN.
     \param computeType How to compute f and h -- \b IN.
     \return     A boolean equal to \c true if  \c *this \c dominates \c eval.
     */
    bool dominates(const Eval& eval, const NOMAD::FHComputeTypeS& computeType) const;

    /// Comparison of multiobjective vectors
    /**
     \param eval The right-hand side object -- \b IN.
     \param onlyfvalues Flag which indicates if h-value must be taken into account. By default, set to false -- \b IN.
     \param fhComputeType How to compute f and h -- \b IN
     \return A compareType flag which can take the following value: "EQUAL", "INDIFFERENT", "DOMINATED", "DOMINATING", "UNDEFINED"
     */ 
    NOMAD::CompareType compMO(const Eval& eval,
                              const NOMAD::FHComputeTypeS& fhComputeType,
                              const bool onlyfvalues = false) const;

    /// Comparison of 2 evaluations.
    /**
     The comparison is used to find the best feasible and infeasible points in the cache.
     \note This is different than dominance.
     \param eval1   First eval -- \b IN.
     \param eval2   Second eval -- \b IN.
     \param computeType     Which compute type to look at                     -- \b IN.
     \return        If \c eval1 dominates \c eval2, return \c true. If eval1 and eval2 are both infeasible, and eval1.getH() < eval2.getH(), return \c true.
    */
    static bool compEvalFindBest(const Eval &eval1,
                                 const Eval &eval2,
                                 const FHComputeTypeS& computeType);


    /// Compute success type of one eval with respect to another one.
    /**
     \param eval1   First eval -- \b IN.
     \param eval2   Second eval -- \b IN.
     \param hMax    The max infeasibility for keeping points in barrier -- \b IN.
     \param computeType How to compute f and h -- \b IN.
     \return        The success type of the first eval with respect to the second one.
     */
    static SuccessType computeSuccessType(const Eval* eval1,
                                          const Eval* eval2,
                                          const NOMAD::FHComputeTypeS& computeType,
                                          const Double& hMax = INF);

    ///  Display of eval
    /**
     \param computeType How to compute f and h -- \b IN.
     \param prec   Display precision -- \b IN.
     \return A formatted eval as a string
     */
    std::string display(const NOMAD::FHComputeTypeS & computeType = defaultFHComputeTypeS, const int prec = DISPLAY_PRECISION_STD) const;
    
    
    /// Helper for setBBO or to be used by evaluator control
    /**
     Add default value (feasible) for revealed constraints added by some algos (like DiscoMads).
     */
    void updateForRevealedConstraints();

    /// Reset DMULTI_COMBINE_F value
    /**
     Reinitialize DMULTI_COMBINE_F to a non-defined value.
     */
    void resetDMultiCombineFValue()
    {
        _moInfo->combineFValue = NOMAD::Double();
    }
    
    /// Reset F values
    /**
     Reinitialize Fvalues to a non-defined value.
     */
    void resetFValues()
    {
        _moInfo->fvalues = NOMAD::ArrayOfDouble();
    }
    
    

private:
    /// Helpers for getF() and getH()
    Double computeHStandard( NOMAD::HNormType hNormType) const;
    Double computeFPhaseOne( NOMAD::HNormType hNormType) const;
    

    
};


/// Definition for evaluation unique pointer
typedef std::unique_ptr<Eval> EvalUPtr;

/**
 * \brief Output raw eval status
 *
 * Does not do the same as enumStr.
 */
std::ostream& operator<<(std::ostream& out, const EvalStatusType &evalStatus);


/// Input eval status
std::istream& operator>>(std::istream& is, EvalStatusType &evalStatus);


#include "../nomad_nsend.hpp"
#endif  // __NOMAD_4_5_EVAL__

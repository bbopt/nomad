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
 \file   Eval.hpp
 \brief  Evaluation at a point
 \author Viviane Rochon Montplaisir
 \date   March 2017
 \see    Eval.cpp
 */

#ifndef __NOMAD400_EVAL__
#define __NOMAD400_EVAL__

#include <functional>   // For std::function

#include "../Eval/BBOutput.hpp"
#include "../Param/EvalParameters.hpp"

#include "../nomad_nsbegin.hpp"

/**
* Idea to implement:
* Eval diffentiates between a BB evaluation that failed
 (failure of the black box) and an evaluation that was
 interrupted (failure cause is external to the black box). \n
* There is also states for evaluations that are not yet started,
* and evaluations that are currently running.
*/


/// Type for an evaluation status
/**
 * \todo These statuses should be made 100% clear.
 * Add something in the name that says if we can re-evaluate the point or not.
 * Even add some User Cases.
 */
enum class EvalStatusType
{
    EVAL_NOT_STARTED,       ///< Evaluation has not been done yet. Initial status.
    EVAL_FAILED,            ///< Evaluation failure. Do not re-submit.
    EVAL_ERROR,             ///< Evaluation did not proceed normally. May be submitted again.
    EVAL_USER_REJECTED,     ///< Evaluation was rejected by user. May be submitted again
    EVAL_CONS_H_OVER,       ///< Evaluation was rejected because constraint violation was higher than hMax. May be submitted again.
    EVAL_OK,                ///< Correct evaluation
    EVAL_IN_PROGRESS,       ///< Evaluation in progress
    EVAL_STATUS_UNDEFINED   ///< Undefined evaluation status
};


/// Utility to convert an eval status to a string.
std::string enumStr(const EvalStatusType evalStatus);


/// Class for the representation of an evaluation at a point.
/**
 * \note We have a separate Eval from EvalPoint so that a Point can have multiple evaluations.
 * \note The class stores value of f and h and manages the computation h (see Eval::_computeH and Eval::_computeHComponent ).
 */
class Eval {

private:
    bool _toBeRecomputed; ///< _f and _h are buffered values. Might have to be recomputed.
    Double _f;   ///< Value of the evaluation
    Double _h;   ///< Value of the constraint violation
    EvalStatusType _evalStatus;  ///> The evaluation status.
    BBOutput _bbOutput;  ///<  The blackbox evaluation output.

    static std::function<SuccessType(const Eval* eval1, const Eval* eval2, const Double& hMax)> _computeSuccessType;  ///< The function called to compute success type.

    /// The function that computes h from an Eval and a list of blackbox output types.
    /**
     * For a given Eval there can be several ways to compute infeasibility
     depending on the listed blackbox output types. For example, during the
     PhaseOneSearch of Mads or after the PhaseOneSearch the computation is
     different. This function sums the infeasibility measure of all constraints
     (see _computeHComponent function).
     */
    static std::function<Double(const Eval& eval, const
                                BBOutputTypeList &bbOutputTypeList)> _computeH;

  /// The function computes the infeasibility measure for a constraint.
    /**
     * The infeasibility is computed for a given blackbox output (Double) and a
     given output type (::BBOutputType). \n
     * A default function is provided in Eval::defaultComputeHComponent that return max(g_i,0)^2. A user can provide a different function using the static Eval::setComputeHComponent function. \n
     * This allows to change the computation of h or relax the constraint bounds ( hMim = 0 -> hMin>0).
     */
    static std::function<Double(const BBOutputType &bbOutputType,
                                size_t index,
                                const Double &bbo)> _computeHComponent;

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
    explicit Eval(std::shared_ptr<EvalParameters> params,
                  const BBOutput &bbOutput);

    /// Copy constructor.
    /**
     \param eval The copied object .
     */
    Eval(const Eval& eval);

    // Destructor.
    // Needed to avoid valgrind memory leak warnings.
    virtual ~Eval() {}

    /*---------*/
    /* Get/Set */
    /*---------*/

    Double getF() const;
    void setF(const Double &f);

    Double getH() const;
    void setH(const Double &h);

    EvalStatusType getEvalStatus() const { return _evalStatus; }
    void setEvalStatus(const EvalStatusType &evalStatus) { _evalStatus = evalStatus; }

    BBOutput getBBOutput() const { return _bbOutput; }
    void setBBOutput(const BBOutput &bbOutput);

    /// Set blackbox output and recompute objective and infeasibility
    void setBBOutputAndRecompute(const BBOutput &bbOutput,
                                 const BBOutputTypeList &bbOutputType);

    std::string getBBO() const { return _bbOutput.getBBO(); }

    /// Set blackbox output and recompute objective and infeasibility
    void setBBO(const std::string &bbo,
                const BBOutputTypeList &bbOutputType,
                const bool evalOk = true);

    /*---------------*/
    /* Other methods */
    /*---------------*/

    bool toBeRecomputed() const { return _toBeRecomputed; }
    void toRecompute(bool toBeRecomputed) { _toBeRecomputed = toBeRecomputed; }

    /// Compute objective function value.
    /**
     \param bbOutputTypeList    The list of types of blackbox outputs -- \b IN.
     \return                    The objective value
     */
    Double computeF(const BBOutputTypeList &bbOutputTypeList) const;

    /// Function to compute infeasibility by aggregating the contribution of each constraint into h.
    /**
     * This is the default function for Eval::_computeH.

     \param eval                The evaluation to consider -- \b IN.
     \param bbOutputTypeList    The list of types of blackbox outputs -- \b IN.
     \return                    The aggregated constraint value
     */
    static Double defaultComputeH(const Eval& eval, const BBOutputTypeList &bbOutputTypeList);

    /// Function to compute each constraint contribution.
    /**
     * This is the default function for Eval::_computeHComponent.
     * If constraint is not verified (bbo>0), the contribution to h is bbo^2 for PB constraint and Infinity for EB constraint.

     \param bbOutputType  The type of the considered constraint -- \b IN.
     \param index         The index of the blackbox output (not used here) -- \b IN.
     \param bbo           The blackbox output for the considered constraint -- \b IN.
     \return              The component of h for constraint referred by index.
     */
    static Double defaultComputeHComponent( const BBOutputType & bbOutputType ,
                                            size_t index __attribute__((unused)),
                                            const Double &bbo );


    /// Compute hPB - infeasibily h computed as if all EB were PB.
    /**
     \param eval                The evaluation -- \b IN.
     \param bbOutputTypeList    The list of types of blackbox outputs (including EB) -- \b IN.
     \return                    The aggregated constraint value
     */
    static Double computeHPB(const Eval& eval, const BBOutputTypeList &bbOutputTypeList);

    /**
     \note In order for this function to be const, Eval::_h is not be recomputed.
     If Eval::_toBeRecomputed is true, a warning is issued.
     */
    bool isFeasible() const;

    /// Can this point be re-evaluated? Based on the eval status only.
    bool canBeReEvaluated() const;

    /** Should this point be saved to cache file? Based on the eval status only.
     * These eval statuses are good: EVAL_OK, EVAL_FAILED, EVAL_USER_REJECTED,
     * EVAL_CONS_H_OVER, EVAL_ERROR.
     * These eval statuses are not good:
     * EVAL_NOT_STARTED, EVAL_IN_PROGRESS, EVAL_STATUS_UNDEFINED.
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
     Dominace as by definition 12.3 in the Book of Audet and Hare, Derivative-Free and Blackbox Optimization,
     https://doi.org/10.1007/978-3-319-68913-5
     * The feasible point x dominates the feasible point y \n
       when f(x) < f(y).
     * The infeasible point x dominates the infeasible point y \n
       when f(x) <= f(y) and h(x) <= h(y),
       with at least one strict inequality.
     * Otherwise, return false
     \param eval The right-hand side object -- \b IN.
     \return     A boolean equal to \c true if  \c *this \c dominates \c eval.
     */
    bool dominates(const Eval & eval) const;

    /// Comparison of 2 evaluations.
    /**
     The comparison is used to find the best feasible and infeasible points in the cache.
     \note This is different than dominance.

     \param eval1   First eval -- \b IN.
     \param eval2   Second eval -- \b IN.
     \return        If \c eval1 dominates \c eval2, return \c true. If eval1 and eval2 are both infeasible, and eval1.getH() < eval2.getH(), return \c true.
    */
    static bool compEvalFindBest(const Eval &eval1,
                                 const Eval &eval2);

    /// Compute success type of one eval with respect to another one.
    /**
     This is the default function for Eval::_computeSuccessType.

     \param eval1   First eval -- \b IN.
     \param eval2   Second eval -- \b IN.
     \param hMax    The max infeasibility for keeping points in barrier -- \b IN.
     \return        The success type of the first eval with respect to the second one.
     */
    static SuccessType defaultComputeSuccessType(const Eval* eval1,
                                                 const Eval* eval2,
                                                 const Double& hMax = INF);

    /// Compute success type of one eval with respect to another one.
    /**
     This is NOT the default function for Eval::_computeSuccessType.
     It is used only for PhaseOne. Requires to call ComputeSuccessType::setComputeSuccessTypeFunction

     \param eval1   First eval -- \b IN.
     \param eval2   Second eval -- \b IN.
     \param hMax    The max infeasibility for keeping points in barrier -- \b IN.
     \return        The success type of the first eval with respect to the second one.
     */
    static SuccessType computeSuccessTypePhaseOne(const Eval* eval1,
                                                  const Eval* eval2,
                                                  const Double& hMax  __attribute__((unused)));

    /// Set which function should be used to compute success type
    static void setComputeSuccessTypeFunction(const std::function<SuccessType(
                                                    const Eval* eval1,
                                                    const Eval* eval2,
                                                    const Double& hMax)>& comp)
    {
        _computeSuccessType = comp;
    }

    /// Set which function to use for infeasibility (h) computation.
    /**
     Needed during PhaseOne to set Eval::_computeH to
     Eval::computeSuccessTypePhaseOne instead of default Eval::defaultComputeH.
     */
    static void setComputeHFunction(const std::function<Double(
                                                       const Eval& eval,
                                                       const BBOutputTypeList &bbOutputTypeList)>& computeHfunc)
    {
        _computeH = computeHfunc;
    }

    /// Set which function to use for computing the components of infeasibility.
    /**
     Needed when a custom function is provided by user to replace
     Eval::defaultComputeHComponent. This can be used to relax
     some constraints for specific optimization problem. This is doable
     in library mode.
     */
    static void setComputeHComponentFunction(const std::function<Double( const BBOutputType &bbOutputType,
                                                                        size_t index ,
                                                                        const Double & bbo)>& computeHComponentfunc)
    {
        _computeHComponent = computeHComponentfunc;
    }


    /// \brief Display of eval
    /// \return A formatted eval as a string
    std::string display() const;
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
#endif  // __NOMAD400_EVAL__

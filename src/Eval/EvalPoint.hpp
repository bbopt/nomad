/**
 \file   EvalPoint.hpp
 \brief  Evaluation point
 \author Viviane Rochon Montplaisir
 \date   April 2017
 \see    EvalPoint.cpp
 */

#ifndef __NOMAD400_EVALPOINT__
#define __NOMAD400_EVALPOINT__

#ifdef USE_UNORDEREDSET
#include <unordered_set>
#else
#include <set>
#endif

#include "../Eval/Eval.hpp"
#include "../Math/Point.hpp"
#include "../Type/EvalType.hpp"

#include "../nomad_nsbegin.hpp"


/// Class for the representation of an evaluation point.
/**
 An evaluation point gathers the point coordinates \c x, and the blackbox
 outputs at these coordinates \c f(x).
*/
class EvalPoint : public Point
{
private:

    static size_t _currentTag;  ///< Value of the current tag

    EvalUPtr _eval;       ///< Value of the evaluation (truth / blackbox)

    EvalUPtr _evalSgte;   ///< Value of the surrogate evaluation

    mutable size_t  _tag; ///< Tag: Ordinal representing the order of creation

    short    _numberEval; ///< Number of times \c *this point has been evaluated (blackbox only)

    std::shared_ptr<Point> _pointFrom; ///< The frame center which generated \c *this point (blackbox only)

public:

    /*---------------*/
    /* Class Methods */
    /*---------------*/

    /// Constructor #1.
    explicit EvalPoint();

    /// Constructor #2.
    /**
     \param n Number of variables -- \b IN.
     */
    explicit EvalPoint(size_t n);

    /// Constructor #3.
    /**
      \param x Coordinates of the eval point -- \b IN.
      */
    explicit EvalPoint(const Point& x);

    /// Copy constructor.
    /**
     \param evalPoint The copied object -- \b IN.
     */
    EvalPoint(const EvalPoint& evalPoint);

private:
    /// Helper for copy constructor and others
    void copyMembers(const EvalPoint &evalPoint);

public:

    /// Affectation operator.
    /**
     \param evalPoint The right-hand side object -- \b IN.
     \return           \c *this as the result of the affectation.
     */
    EvalPoint& operator= (const EvalPoint& evalPoint);

    /// Destructor.
    virtual ~EvalPoint();

    /*---------*/
    /* Get/Set */
    /*---------*/
    /// Get Point part of this EvalPoint
    const Point* getX() const { return dynamic_cast<const Point*>(this); }

    /// Get the Eval part of this EvalPoint, using the right EvalType (BB or SGTE)
    Eval* getEval(const EvalType& evalType = EvalType::BB) const;

    /// Set the Eval part of this EvalPoint, using the right EvalType (BB or SGTE)
    void setEval(const Eval& eval, const EvalType& evalType);

    /// Clear the surrogate evaluation of \c *this
    void clearEvalSgte() { _evalSgte = nullptr; }

    /// Clear the surrogate evaluation of a point
    static void clearEvalSgte(EvalPoint& evalPoint) { evalPoint.clearEvalSgte(); }

    /// Get the objective function value of Eval of this EvalType
    Double getF(const EvalType& evalType = EvalType::BB) const;

    /// Set the objective function value of the Eval of this EvalType
    void setF(const Double f, const EvalType& evalType);

    /// Get the infeasibility measure of the Eval of this EvalType
    Double getH(const EvalType& evalType = EvalType::BB) const;

    /// Set the infeasibility measure of the Eval of this EvalType
    void setH(const Double &h, const EvalType& evalType);

    /// Get the blackbox output for the Eval of this EvalType as a \c string
    std::string getBBO(const EvalType& evalType) const;

    /// Set the blackbox output for the Eval of this EvalType from a \c string.
    /**
     \param bbo             The string containg the raw result of the blackbox evaluation -- \b IN.
     \param bboutputtypes   The list of blackbox output types -- \b IN.
     \param evalType        Blackbox or surrogate evaluation  -- \b IN.
     \param evalOk          Flag for evaluation status  -- \b IN.
    */
    void setBBO(const std::string &bbo,
                const BBOutputTypeList &bboutputtypes,
                const EvalType& evalType = EvalType::BB,
                const bool evalOk = true);

    /// Set the true or surrogate blackbox output from a \c string.
    /**
     \param bbo             The string containg the raw result of the blackbox evaluation -- \b IN.
     \param sBBOutputTypes  The blackbox output types coded as a single string -- \b IN.
     \param evalType        Blackbox or surrogate evaluation  -- \b IN.
     \param evalOk          Flag for evaluation status  -- \b IN.
     */
    void setBBO(const std::string &bbo,
                const std::string &sBBOutputTypes,
                const EvalType& evalType = EvalType::BB,
                const bool evalOk = true);

    /// Set the true or surrogate blackbox output.
    /**
     \param bbo             A blackbox evaluation output -- \b IN.
     \param evalType        Blackbox or surrogate evaluation  -- \b IN.
     \param evalOk          Flag for evaluation status  -- \b IN.
     */
    void setBBO(const BBOutput &bbo,
                const EvalType& evalType = EvalType::BB,
                const bool evalOk = true);

    /// Get evaluation status of the Eval of this EvalType
    EvalStatusType getEvalStatus(const EvalType& evalType) const;

    /// Set evaluation status of the Eval of this EvalType
    void setEvalStatus(const EvalStatusType &evalStatus, const EvalType& evalType);

    size_t getTag() const { return _tag; }
    void setTag(const size_t tag) const { _tag = tag; } ///< Sets mutable _tag
    void updateTag() const; ///< Modifies mutable _tag, and increments static _currentTag
    static void resetCurrentTag(); ///< Reset tag numbers: Use with caution. Expected to be used in unit tests  and runner only.

    short getNumberEval() const { return _numberEval; }
    void setNumberEval(const short numEval) { _numberEval = numEval; }
    void incNumberEval() { _numberEval++; }

    /// Get the Point which was the center when this point was generated
    const std::shared_ptr<Point> getPointFrom() const { return _pointFrom; }

    /// Get the Point which was the center when this point was generated
    /**
     Returns a Point in the Subspace defined by the fixedVariable
     */
    const std::shared_ptr<Point> getPointFrom(const Point& fixedVariable) const;

    /// Set the Point for which this point was generated
    void setPointFrom(const std::shared_ptr<Point> pointFrom);

    /// Set the Point for which this point was generated
    /**
     Use the fixedVariable to convert pointFrom from Subspace dimension to the full dimension.
     */
    void setPointFrom(std::shared_ptr<Point> pointFrom,
                      const Point& fixedVariable);

    /// Get evaluation feasibility flag f the Eval of this EvalType
    bool isFeasible(const EvalType& evalType) const;

    /// Recompute f and h, for a given list of blackbox output type considering that raw blackbox output is set.
    /**
     \param bbOutputType    The list of blackbox output types.
    */
    void recomputeFH(const BBOutputTypeList &bbOutputType);

    /// Comparison operator used by NM algorithm.
    /**
     \param rhs     Second eval points to compare      -- \b IN.
     \param evalType        Blackbox or surrogate evaluation  -- \b IN.
     \return        \c true if \c *this dominates x.
     */
    bool dominates(const EvalPoint& rhs, const EvalType& evalType) const;

    /// Convert a point from sub space to full space using fixed variables.
    /**
     \remark The evaluation part of \c *this is unchanged.
     */
    EvalPoint makeFullSpacePointFromFixed(const Point &fixedVariable) const;

    /// Convert a point from full space to sub space using fixed variables
    /**
     \remark The evaluation part of \c *this is unchanged.
     */
    EvalPoint makeSubSpacePointFromFixed(const Point &fixedVariable) const;

    /*----------------------*/
    /* Comparison operators */
    /*----------------------*/

    /// Comparison operator \c ==.
    /**
     \param evalPoint   The right-hand side object -- \b IN.
     \return            \c true if  \c *this \c == \c p, \c false if not.
     */
    bool operator== (const EvalPoint& evalPoint) const;

    /// Comparison operator \c !=.
    /**
     \param evalPoint   The right-hand side object -- \b IN.
     \return            \c false if  \c *this \c == \c p, \c true if not.
     */
    bool operator!= (const EvalPoint& evalPoint) const { return !(*this == evalPoint); }

    /// Comparison operator \c <, used for set ordering.
    /**
     \param x       Right-hand side object -- \b IN.
     \return        \c true if \c *this \c < \c x, \c false if not..
     */
    bool operator< (const EvalPoint& x) const;


    /*---------------*/
    /* Class methods */
    /*---------------*/
    bool isEvalOk(const EvalType& evalType) const;

    /// Display with or without format
    std::string display(const ArrayOfDouble &format = ArrayOfDouble()) const override;

    /// Display both true and surrogate evaluations.  Useful for debugging
    std::string displayAll() const;

    /// Function to test if evaluation is required.
    /**
     * Depending on the status of the Eval, should we evaluate
     * (possibly re-evaluate) this point?

     \param maxPointEval    The maximum number of point evaluations  -- \b IN.
     \param evalType        Blackbox or surrogate evaluation  -- \b IN.
     \return                \c true if evaluation is required and \c false otherwise.
     */
    bool toEval(short maxPointEval, const EvalType& evalType) const;

    /**
    \warning It is unclear if the caller wants to verify if the base point is defined,
    or if f is defined. To avoid mistakes and confusion, throw an error.
     */
    bool isDefined() const override
    {
        throw Exception(__FILE__,__LINE__,"Error: Calling EvalPoint::isDefined(). Choose ArrayOfDouble::isDefined() or Double::isDefined() instead.");
    }

    // Determine if an evalpoint has a sgte eval.
    static bool hasSgteEval(const EvalPoint& evalPoint);
    // Determine if an evalpoint has a bb (regular) eval.
    static bool hasBbEval(const EvalPoint& evalPoint);
};


/// Display useful values so that a new EvalPoint could be constructed using these values.
std::ostream& operator<<(std::ostream& os, const EvalPoint &evalPoint);

/// Get these values from stream
std::istream& operator>>(std::istream& is, EvalPoint &evalPoint);

/// Definition for eval point pointer
typedef std::shared_ptr<EvalPoint> EvalPointPtr;

/// Definition for block (vector) of EvalPointPtr
typedef std::vector<EvalPointPtr> Block;

/// Utility to find Point in EvalPoint vector
bool findInList(const Point& point, const std::vector<EvalPoint>& evalPointList,
                EvalPoint& foundEvalPoint);


/// Class for eval point compare.
/**
 * Compare the Point parts only.
 * Trying to insert an EvalPoint in a set that has already a Point defined
 * for that EvalPoint will return false.
 */
class EvalPointCompare
{
public:
    bool operator() (const EvalPoint& lhs, const EvalPoint& rhs) const
    {
        return Point::weakLess(*(lhs.getX()),*(rhs.getX()));
    }
};

/// Definition for EvalPointSet
#ifdef USE_UNORDEREDSET
    typedef std::unordered_set<EvalPoint,
                               std::hash<EvalPoint>,
                               std::equal_to<EvalPoint>> EvalPointSet;
#else
    typedef std::set<EvalPoint, EvalPointCompare> EvalPointSet;
#endif

/// Definition for EvalPointList
typedef std::vector<EvalPoint> EvalPointList;

void convertPointListToSub(EvalPointList &evalPointList,  const Point& fixedVariable);
void convertPointListToFull(EvalPointList &evalPointList, const Point& fixedVariable);


#include "../nomad_nsend.hpp"

#ifdef USE_UNORDEREDSET
/**
 * \warning If we use unordered_set, then we must define hash.
 * Template specialization for std::hash<class T> to T=EvalPoint
 */
namespace std {
    template <>
    struct hash<EvalPoint>
    {
        public:
        size_t operator()(const EvalPoint& evalPoint) const;
    };

/**
 * If we use unordered_set, then we must define function call operator to test equality.
 * Template specialization of std::equal_to<class T> to T=EvalPoint
 */
    template <>
    class equal_to<EvalPoint>
    {
        public:
        bool operator()(const EvalPoint& lhs, const EvalPoint& rhs) const;
    };
}
#endif // USE_UNORDEREDSET


#endif // __NOMAD400_EVALPOINT__

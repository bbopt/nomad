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
 \file   EvalPoint.hpp
 \brief  Evaluation point
 \author Viviane Rochon Montplaisir
 \date   April 2017
 \see    EvalPoint.cpp
 */

#ifndef __NOMAD_4_5_EVALPOINT__
#define __NOMAD_4_5_EVALPOINT__

#ifdef USE_UNORDEREDSET
#include <unordered_set>
#else
#include <set>
#endif

#include "../Eval/Eval.hpp"
#include "../Eval/MeshBase.hpp"
#include "../Math/Point.hpp"
#include "../Type/ComputeType.hpp"
#include "../Type/EvalType.hpp"
#include "../Type/StepType.hpp"

#include "../nomad_platform.hpp"
#include "../nomad_nsbegin.hpp"


/// Class for the representation of an evaluation point.
/**
 An evaluation point gathers the point coordinates \c x, and the blackbox
 outputs at these coordinates \c f(x).
*/
class DLL_EVAL_API EvalPoint : public Point
{
private:

    static int     _currentTag;  ///< Value of the current tag

    EvalUPtr _eval[3];  ///< Value of the evaluation for each eval type. The index number must correspond to EvalType enum order
                        ///< Index 0: access to BB
                        ///< Index 1: access to SURROGATE
                        ///< Index 2: access to MODEL

    mutable int                 _tag; ///< Tag: Ordinal representing the order of creation

    int                         _threadAlgo;    ///< Main thread that generated this point

    short                       _numberBBEval; ///< Number of times \c *this point has been evaluated (blackbox only)

    bool                        _evalFromCacheFile; ///<  Evaluation (BB or SURROGATE obtained from reading a cache file
    
    std::shared_ptr<EvalPoint>  _pointFrom; ///< The frame center which generated \c *this point (blackbox only). Full space.

    StepTypeList                _genSteps;   ///< Steps and algorithms that generated this point

    std::shared_ptr<Direction>  _direction; ///< True direction that generated this point. Full dimension.
    Double                      _angle;     ///< Angle of that direction with last successful dir

    std::shared_ptr<MeshBase> _mesh;
    
    int                        _isRevealing;  ///< for DiscoMads, 0 : not revealing, 1: revealing, 2: just revealed at the last eval
    
    bool                       _userFailEvalCheck; ///< for DiscoMads handling hidden constraints. True: fail eval point has been checked (not necessarily EvalOk after)
             
    
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
    /// Helper for constructors
    void initEval();

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

    /// Get the Eval part of this EvalPoint, using the right EvalType (BB or MODEL)
    Eval* getEval(EvalType evalType) const;


    /// Get the single eval with the given eval status. If no eval is defined, return nullptr. If more than one eval has the eval status, throw exception.
    NOMAD::EvalType getSingleEvalType(NOMAD::EvalStatusType evalStatusType) const;
    
    /// Set the Eval part of this EvalPoint, using the right EvalType (BB or MODEL)
    void setEval(const Eval& eval, EvalType evalType);
    
    /// Clear the eval type evaluation of \c *this
    void clearEval(EvalType evalType)
    {
        _eval[(size_t) evalType].reset();
    }

    /// Clear the model evaluation of a point
    static void clearModelEval(EvalPoint& evalPoint) { evalPoint.clearEval(EvalType::MODEL); }

    /// Get the objective function value of Eval of this EvalType,
    /// using the given ComputeType.
    Double getF(const FHComputeType& completeComputeType ) const;

    /// Get the objective function vector values of Eval of this EvalType,
    /// using the given ComputeType.
    const ArrayOfDouble& getFs(const FHComputeType& completeComputeType ) const;

    /// Get the infeasibility measure of the Eval of this EvalType
    Double getH(const FHComputeType& completeComputeType) const;

    /// Get the blackbox output for the Eval of this EvalType as a \c string
    std::string getBBO(EvalType evalType) const;

    /// Set the blackbox output for the Eval of this EvalType from a \c string.
    /**
     \param bbo                 The string containing the raw result of the blackbox evaluation -- \b IN.
     \param bbOutputTypeList    The list of blackbox output types -- \b IN.
     \param evalType            Blackbox or model evaluation  -- \b IN.
     \param evalOk              Flag for evaluation status  -- \b IN.
    */
    void setBBO(const std::string &bbo,
                const BBOutputTypeList& bbOutputTypeList,
                EvalType evalType = EvalType::LAST,
                const bool evalOk = true);

    /// Set the true or model blackbox output from a \c string.
    /**
     \param bbo             The string containing the raw result of the blackbox evaluation -- \b IN.
     \param sBBOutputTypes  The blackbox output types coded as a single string -- \b IN.
     \param evalType        Blackbox or model evaluation  -- \b IN.
     \param evalOk          Flag for evaluation status  -- \b IN.
     */
    void setBBO(const std::string &bbo,
                const std::string &sBBOutputTypes = "",
                EvalType evalType = EvalType::LAST,
                const bool evalOk = true);

    void setBBOutputType(const BBOutputTypeList& bbOutputType, const EvalType evalType) const;
    void setBBOutputType(const BBOutputTypeList& bbOutputType);
    

    /// Get evaluation status of the Eval of this EvalType
    EvalStatusType getEvalStatus(EvalType evalType) const;
    
    /// Get pre evaluation status of the Eval of this EvalType
    EvalStatusType getPreEvalStatus(EvalType evalType) const;

    /// Set evaluation status of the Eval of this EvalType
    void setEvalStatus(EvalStatusType evalStatus, EvalType evalType);
    
    /// Set pre evaluation status of the Eval of this EvalType
    void setPreEvalStatus(EvalStatusType evalStatus, EvalType evalType);

    int getTag() const { return _tag; }
    void setTag(const int tag) const { _tag = tag; } ///< Sets mutable _tag
    void updateTag() const; ///< Modifies mutable _tag, and increments static _currentTag
    static void resetCurrentTag(); ///< Reset tag numbers: Use with caution. Expected to be used in unit tests  and runner only.
    static int getCurrentTag() { return _currentTag;} ///< Access to current tag

    int getThreadAlgo() const { return _threadAlgo; }
    void setThreadAlgo(const int threadAlgo) { _threadAlgo = threadAlgo; }

    short getNumberBBEval() const { return _numberBBEval; }
    void setNumberBBEval(const short numBBEval) { _numberBBEval = numBBEval; }
    void incNumberBBEval() { _numberBBEval++; }
    
    bool getEvalIsFromCacheFile() const { return _evalFromCacheFile ;}
    void setEvalIsFromCacheFile(bool flag) { _evalFromCacheFile = flag; }

    /// Get the Point which was the center when this point was generated
    std::shared_ptr<EvalPoint> getPointFrom() const { return _pointFrom; }

    /// Push the step type
    void addGenStep(const StepType& stepType, bool inherit=true);
    /// Get the downmost step type
    StepType getGenStep() const;
    /// Get vector of all step types
    const StepTypeList& getGenSteps() const;
    /// Set vector of all step types
    void setGenSteps(const StepTypeList& genSteps);
    /// Check if this EvalPoint was generated by PhaseOne
    bool getGenByPhaseOne() const;
    /// Algo comment to be printed out in DISPLAY_STATS
    std::string getComment() const;

    /// Get Direction from which the point was generated.
    /// Value is set when setPointFrom() is called.
    const std::shared_ptr<Direction> getDirection() const { return _direction; }
    /// Get/Set angle of direction with direction of last success
    const Double& getAngle() const { return _angle; }
    void setAngle(const Double& angle) { _angle = angle; }

    
    void setMesh(const MeshBasePtr& mesh);
    MeshBasePtr getMesh() const { return _mesh; }
    
    /// Get revealing status (for DiscoMads)
    int getRevealingStatus() const {return _isRevealing;}
    
    /// Set revealing status (for DiscoMads)
    /**
     \param revealingStatus Flag for revealing status IN.
    */
    void setRevealingStatus(const int& revealingStatus){_isRevealing=revealingStatus;};
    
    /// Get user fail eval check status (for DiscoMads)
    bool getUserFailEvalCheck() const { return _userFailEvalCheck; }
    
    /// Set user fail eval check status (for DiscoMads)
    /**
     \param userFailEvalCheck Flag for user fail check IN.
     */
    void setUserFailEvalCheck(bool userFailEvalCheck){_userFailEvalCheck=userFailEvalCheck;};

    /// Get revealed constraint value (for some algorithms like DiscoMads)
    NOMAD::Double getRevealedConstraint() const;

    /// Set revealed constraint value (for some algorithms like DiscoMads) 
    void setRevealedConstraint(const NOMAD::Double &constraintValue) const;

    /// Get the Point which was the center when this point was generated
    /**
     Returns a Point in the Subspace defined by the fixedVariable
     */
    std::shared_ptr<EvalPoint> getPointFrom(const Point& fixedVariable) const;

    /// Set the Point for which this point was generated
    /**
     Use the fixedVariable to convert pointFrom from Subspace dimension to the full dimension, if needed.
     */
    void setPointFrom(const std::shared_ptr<EvalPoint>& pointFrom,
                      const Point& fixedVariable);

    /// Get evaluation feasibility flag f the Eval of this EvalType
    bool isFeasible(const FHComputeType & completeComputeType) const;

    /// Get feasibility with respect to EB constraints only
    /**
     Returns true if point satisfies all EB constraints
     */
    bool isEBOk(NOMAD::EvalType evalType) const;

    /// Comparison operator used by NM algorithm.
    /**
     \param rhs         Second eval points to compare      -- \b IN.
     \param completeComputeType   Which type of f, h computation (eval type, compute type and h norm type)  -- \b IN.
     \return        \c true if \c *this dominates x.
     */
    bool dominates(const EvalPoint& rhs,
                   const FHComputeType& completeComputeType) const;

    /// Comparison operator used in Multiobjective Optimization
    /**
     \param onlyfvalues Flag which indicates if h-value must be taken into account. By default, set to false -- \b IN.
     \param completeComputeType   Which type of f, h computation (eval type, compute type and h norm type)  -- \b IN.
     \return A compareType flag which can take the following value: "EQUAL", "INDIFFERENT", "DOMINATED", "DOMINATING", "UNDEFINED"
     */ 
    NOMAD::CompareType compMO(const EvalPoint& rhs,
                              const FHComputeType& completeComputeType,
                              bool onlyfvalues = false) const;

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


    /*---------------*/
    /* Class methods */
    /*---------------*/
    bool isEvalOk(EvalType evalType) const;

    /// Display with or without format
    std::string display(const NOMAD::FHComputeTypeS& computeType,
                        const ArrayOfDouble &pointFormat = ArrayOfDouble(),
                        const int &solFormat = NOMAD::DISPLAY_PRECISION_FULL,
                        const bool surrogateAsBB = false) const;

    std::string display(const ArrayOfDouble &pointFormat = ArrayOfDouble(),
                        const int &solFormat = NOMAD::DISPLAY_PRECISION_FULL) const;

    std::string displayForCache(const ArrayOfDouble &pointFormat) const ;

    /// Display both true and model evaluations. Useful for debugging
    std::string displayAll(const NOMAD::FHComputeTypeS& computeType = NOMAD::defaultFHComputeTypeS) const;

    /// Function to test if evaluation is required.
    /**
     * Depending on the status of the Eval, should we evaluate
     * (possibly re-evaluate) this point?

     \param maxPointEval    The maximum number of point evaluations  -- \b IN.
     \param evalType        Blackbox or model evaluation  -- \b IN.
     \return                \c true if evaluation is required and \c false otherwise.
     */
    bool toEval(short maxPointEval, EvalType evalType) const;

    /**
    \warning It is unclear if the caller wants to verify if the base point is defined, or if f is defined. To avoid mistakes and confusion, throw an error.
     */
    bool isDefined() const override
    {
        throw Exception(__FILE__,__LINE__,"Error: Calling EvalPoint::isDefined(). Choose ArrayOfDouble::isDefined() or Double::isDefined() instead.");
    }

    /// Determine if an eval point has a bb (regular) eval.
    static bool hasBbEval(const EvalPoint& evalPoint);
    
    /// Determine if an eval point has a model eval.
    static bool hasModelEval(const EvalPoint& evalPoint);
    
    /// Determine if an eval point has a static surrogate eval.
    static bool hasSurrogateEval(const EvalPoint& evalPoint);
    /// Used for phase one
    static bool isPhaseOneSolution(const EvalPoint& evalPoint, const FHComputeType & completeComputeType);
    
    ///  Get the rank of the generating directions of a vector of eval points
    /**
     \param vectEvalPoints   The eval points   -- \b IN.
     \return                 The rank of generating directions
     */
    static size_t getRank(const std::vector<NOMAD::EvalPoint> & vectEvalPoints);

    /// Reset DMULTI_COMBINE_F value of an EvalPoint
    /**
     Reinitialize DMULTI_COMBINE_F to a non-defined value.
     */
    void resetDMultiCombineFValue();
    
    /// Reset F values of an EvalPoint
    /**
     Reinitialize Fvalues to a non-defined value.
     */
    void resetFValues();
};


/// Display useful values so that a new EvalPoint could be constructed using these values.
DLL_EVAL_API std::ostream& operator<<(std::ostream& os,
                         const EvalPoint &evalPoint);

/// Get these values from stream
DLL_EVAL_API std::istream& operator>>(std::istream& is, EvalPoint &evalPoint);

/// Definition for eval point shared pointer
typedef std::shared_ptr<EvalPoint> EvalPointPtr;
typedef std::shared_ptr<const EvalPoint> EvalPointCstPtr;

/// Definition for block (vector) of EvalPointPtr
typedef std::vector<EvalPointPtr> Block;

/// Utility to find Point in EvalPoint vector
DLL_EVAL_API bool findInList(const Point& point,
                const std::vector<EvalPoint>& evalPointList,
                EvalPoint& foundEvalPoint);


/// Class for eval point compare.
/**
 * Compare the Point parts only.
 * Trying to insert an EvalPoint in a set that has already a Point defined
 * for that EvalPoint will return false.
 */
class DLL_EVAL_API EvalPointCompare
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
    static const EvalPointSet emptyEvalPointSet = {}; /// For convenience. Used when calling functions with EvalPointSet as arguments.

/// Definition for EvalPointList
typedef std::list<EvalPoint> EvalPointList;

DLL_EVAL_API void convertPointListToSub(std::vector<EvalPoint> &evalPointList,  const Point& fixedVariable);
DLL_EVAL_API void convertPointListToFull(std::vector<EvalPoint> &evalPointList, const Point& fixedVariable);


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


#endif // __NOMAD_4_5_EVALPOINT__

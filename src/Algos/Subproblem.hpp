/**
 \file   Subproblem.hpp
 \brief  Subproblem of lesser dimension than the original problem
 \author Viviane Rochon Montplaisir
 \date   February 2019
 */
#ifndef __NOMAD400_SUBPROBLEM__
#define __NOMAD400_SUBPROBLEM__

#include "../Math/Point.hpp"
#include "../Param/PbParameters.hpp"

#include "../nomad_nsbegin.hpp"

/// Class to define an optimization subproblem
/**
*  Subproblem of lesser dimension than the original problem
*
* - Sets up the new parameters
* - Keeps the necessary information to bridge the gap between subproblem and
    original problem
*/
class Subproblem
{
private:
    /**
      * The elements of this point that have defined values are fixed value
      * "variables". The elements that are undefined are for true variables.
      * This Point is always in full dimension.
     */
    const Point  _fixedVariable;
    size_t       _dimension;   ///< Dimension of the subproblem.

    /**
     Reference to the original problem's PbParameters.
     */
    const std::shared_ptr<PbParameters>  _refPbParams;

    /**
     PbParameters converted to subdimension
     */
    std::shared_ptr<PbParameters>        _subPbParams;

public:
    /// Constructor
    /**
     Pb parameters will be recomputed as dimension has changed.
     */
    explicit Subproblem(const std::shared_ptr<PbParameters> refPbParams,
                        const Point& fullFixedVariable)
      : _fixedVariable(fullFixedVariable),
        _dimension(refPbParams->getAttributeValue<size_t>("DIMENSION")),
        _refPbParams(refPbParams),
        _subPbParams(nullptr)
    {
        init();
    }

    /// Destructor
    virtual ~Subproblem();

    // Get/Set

    const Point& getFixedVariable() const { return _fixedVariable; }
    std::shared_ptr<PbParameters> getPbParams() const { return _subPbParams; }

private:
    /// Helper for constructor calls to Subproblem::setupProblemParameters
    void init();

    /// Helper for constructor
    /**
     Construct the subproblem parameters (X0, LB, UB, mesh sizes, variable groups...) based on Subproblem::_fixedVariable
     \note If a new parameter with dimension (ex. a parameter of type ArrayOfDouble, Point, or Dimension)
     is added to the class PbParameters, this method will break.
     Currently supported parameters: X0 LOWER_BOUND UPPER_BOUND BB_INPUT_TYPE INITIAL_MESH_SIZE INITIAL_FRAME_SIZE MIN_MESH_SIZE MIN_FRAME_SIZE GRANULARITY VARIABLE_GROUP
    */
    void setupProblemParameters();


    ///  Helper for setupProblemParameters()
    void resetVariableGroupsAgainstFixedVariables(ListOfVariableGroup& lvg, const Point& fixedVar) const;
};


#include "../nomad_nsend.hpp"


#endif  // __NOMAD400_SUBPROBLEM__

/**
 \file   ComputeSuccessType.hpp
 \brief  Comparison methods for EvalPoints
 \author Viviane Rochon Montplaisir
 \date   April 2017 / September 2020
 \see    ComputeSuccessType.cpp
 */

#ifndef __NOMAD400_COMPUTESUCCESSTYPE__
#define __NOMAD400_COMPUTESUCCESSTYPE__

#include "../Eval/EvalPoint.hpp"


#include "../nomad_nsbegin.hpp"
/// Definition for compute success type function.
/**
 A function of this type compares two EvalPoints, and returns the SuccessType resulting from the comparison. The function is a member of ComputeSuccessType class and set using ComputeSuccessType::setComputeSuccessTypeFunction. \n For example, computing success type is changed when doing PhaseOne, or optimizing a surrogate instead of blackbox.
*/
typedef std::function<SuccessType(const EvalPointPtr &p1,
                                  const EvalPointPtr &p2,
                                  const Double& hMax)> ComputeSuccessFunction;

class ComputeSuccessType
{
private:
    /** The function to compute success type
     */
    ComputeSuccessFunction _computeSuccessType;

public:

    /// Constructor 1
    //ComputeSuccessType(ComputeSuccessFunction computeSuccessFunction defaultComputeSuccessType)
    ComputeSuccessType(ComputeSuccessFunction computeSuccessFunction)
      : _computeSuccessType(computeSuccessFunction)
    {}

    /// Constructor 2
    ComputeSuccessType(const EvalType& evalType)
    {
        setDefaultComputeSuccessTypeFunction(evalType);
    }

    void setComputeSuccessTypeFunction(const ComputeSuccessFunction &computeSuccessFunction)
    {
        _computeSuccessType = computeSuccessFunction;
    }

    /// Set default function for comparing EvalPoints, depending if the evaluation is surrogate or blackbox
    void setDefaultComputeSuccessTypeFunction(const EvalType& evalType);

    /// Function call operator
    /**
     \param p1      First eval point -- \b IN.
     \param p2      Second eval point -- \b IN.
     \param hMax    Max acceptable infeasibility to keep point in barrier -- \b IN.
     \return        Success type of p1 over p2, considering hMax
     */
    SuccessType operator()(const EvalPointPtr& p1,
                           const EvalPointPtr& p2,
                           const Double& hMax = INF);


    /// Function for default compute
    /**
     \param evalPoint1 First eval queue point -- \b IN.
     \param evalPoint2 Second eval queue point -- \b IN.
     \param hMax       Max acceptable infeasibility to keep point in barrier   -- \b IN.
     \return           Success type.
     */
    static SuccessType defaultComputeSuccessType(const EvalPointPtr& evalPoint1,
                                                 const EvalPointPtr& evalPoint2,
                                                 const Double& hMax = INF);

    /// Function to compute success type when in PhaseOne.
    /**
     \param evalPoint   First eval queue point -- \b IN.
     \param xInf        Second eval queue point -- \b IN.
     \param hMax        Max acceptable infeasibility to keep point in barrier -- \b IN.
     \return            Success type.
     */
    static SuccessType computeSuccessTypePhaseOne(const EvalPointPtr& evalPoint,
                                                  const EvalPointPtr& xInf,
                                                  const Double& hMax __attribute__((unused)));

    /// Function to compute success type for a surrogate evaluation.
    /**
     \param evalPoint1  First eval queue point -- \b IN.
     \param evalPoint2  Second eval queue point -- \b IN.
     \param hMax        Max acceptable infeasibility to keep point in barrier   -- \b IN.
     \return            Success type.
     */
    static SuccessType computeSuccessTypeSgte(const EvalPointPtr& evalPoint1,
                                              const EvalPointPtr& evalPoint2,
                                              const Double& hMax = INF);

};
#include "../nomad_nsend.hpp"

#endif // __NOMAD400_COMPUTESUCCESSTYPE__

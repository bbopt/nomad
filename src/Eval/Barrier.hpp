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
#ifndef __NOMAD_4_3_BARRIER__
#define __NOMAD_4_3_BARRIER__

#include "../Eval/BarrierBase.hpp"
#include "../Eval/EvalPoint.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for single objective barrier following algorithm 12.2 of DFBO.
class DLL_EVAL_API Barrier : public BarrierBase
{

public:
    /// Constructor
    /**
     * hMax will be updated during optimization.
     \param hMax            The max of h to keep a point in the barrier -- \b IN.
     \param fixedVariable   The fixed variables have a fixed value -- \b IN.
     \param evalType        Type of evaluation (BB or MODEL) -- \b IN.
     \param computeType  Type of function computation (standard, phase-one or user) -- \b IN.
     \param evalPointList   Additional points to consider in building the barrier -- \b IN.
     \param barrierInitializedFromCache Flag to initialize the barrier from cache -- \b IN.
     */
    Barrier(const Double& hMax = INF,
            const Point& fixedVariable = Point(),
            EvalType evalType = EvalType::BB,
            ComputeType computeType = ComputeType::STANDARD,
            const std::vector<EvalPoint>& evalPointList = std::vector<EvalPoint>(),
            bool barrierInitializedFromCache= true)
      : BarrierBase(hMax)
    {
        init(fixedVariable, evalType, computeType, barrierInitializedFromCache);
        init(fixedVariable, evalType, evalPointList, computeType);
    }
    
    Barrier(const Barrier & b)
    {
        // Do not copy the barrier points. Do not initialize from cache.
        
        _hMax = b._hMax;
        
    }
    
    std::shared_ptr<BarrierBase> clone() const override {
      return std::make_shared<Barrier>(*this);
    }
    
    /*-----------------*/
    /* Feasible points */
    /*-----------------*/

    /// Update ref best feasible and ref best infeasible values.
    void updateRefBests() override;

   ///  Get the first feasible point in the barrier.
   /**
    * If there is no feasible point, return a \c nullptr
    \return A single feasible eval point.
    */
   EvalPointPtr getFirstXFeas() const override;


    /*-------------------*/
    /* Infeasible points */
    /*-------------------*/

   ///  Get the first infeasible point in the barrier.
   /**
    * If there is no infeasible point, return a \c nullptr
    \return A single infeasible eval point.
    */
    EvalPointPtr getFirstXInf() const override;

    /*---------------*/
    /* Other methods */
    /*---------------*/

    /// Set the hMax of the barrier
    /**
     \param hMax    The hMax -- \b IN.
    */
    void setHMax(const Double &hMax) override;

    ///  xFeas and xInf according to given points.
    /* \param evalPointList vector of EvalPoints  -- \b IN.
     * \param keepAllPoints keep all good points, or keep just one point as in NOMAD 3 -- \b IN.
     * \return true if the Barrier was updated, false otherwise
     * \note Input EvalPoints are already in subproblem dimention
     */
    SuccessType getSuccessTypeOfPoints(const EvalPointPtr xFeas,
                                       const EvalPointPtr xInf,
                                       EvalType evalType,
                                       ComputeType computeType) override;

    /// Update xFeas and xInf according to given points.
    /* \param evalPointList vector of EvalPoints  -- \b IN.
     * \param keepAllPoints keep all good points, or keep just one point as in NOMAD 3 -- \b IN.
     * \return true if the Barrier was updated, false otherwise
     * \note Input EvalPoints are already in subproblem dimention
     */
    bool updateWithPoints(const std::vector<EvalPoint>& evalPointList,
                          EvalType evalType,
                          ComputeType computeType,
                          const bool keepAllPoints = false) override;

    /// Return the barrier as a string.
    /* May be used for information, or for saving a barrier. In the former case,
     * it may be useful to set parameter max to a small value (e.g., 4). In the
     * latter case, INF_SIZE_T should be used so that all points are saved.
     * \param max Maximum number of feasible and infeasible points to display
     * \return A vector of string describing the barrier
     */
    std::vector<std::string> display(const size_t max = INF_SIZE_T) const override;

private:

    /**
     * \brief Helper function for constructor.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     \param fixedVariable   The fixed variables have a fixed value     -- \b IN.
     \param evalType        Which eval (Blackbox or Model) to use to verify feasibility  -- \b IN.
     \param computeType    Which compute type (standard, phase-one or user) must be available to find in cache  -- \b IN.
     \param barrierInitializedFromCache  Flag to initialize barrier from cache or not. -- \b IN.
     */
    void init(const Point& fixedVariable,
              EvalType evalType,
              ComputeType computeType,
              bool barrierInitializedFromCache) override;

    /**
     * \brief Helper function for constructor.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     \param fixedVariable   The fixed variables have a fixed value     -- \b IN.
     \param evalType        Which eval (Blackbox or Model) to use to verify feasibility  -- \b IN.
     \param evalPointList   Additional points to consider to construct barrier. -- \b IN.
     \param computeType    Which compute type (standard, phase-one or user) must be available to find in cache  -- \b IN.
     */
    void init(const Point& fixedVariable,
              EvalType evalType,
              const std::vector<EvalPoint>& evalPointList,
              ComputeType computeType);
    



    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkXFeas(const EvalPoint &xFeas,
                    EvalType evalType,
                    ComputeType computeType = ComputeType::STANDARD) override;

    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkXFeasIsFeas(const EvalPoint &xFeas,
                          EvalType evalType,
                          ComputeType computeType = ComputeType::STANDARD);

   /**
    * \brief Helper function for insertion.
    *
    * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
    */
   void checkXInf(const EvalPoint &xInf, EvalType evalType) override;

    
};




#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_3_BARRIER__

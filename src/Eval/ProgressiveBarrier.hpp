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
#ifndef __NOMAD_4_5_PROGRESSIVEBARRIER__
#define __NOMAD_4_5_PROGRESSIVEBARRIER__

#include "../Eval/BarrierBase.hpp"
#include "../Eval/EvalPoint.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for single objective progressive barrier (for Mads) following algorithm 12.2 of DFBO.
class DLL_EVAL_API ProgressiveBarrier : public BarrierBase
{
private:
    bool _incumbentsAndHMaxUpToDate;
    
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
    ProgressiveBarrier(const Double& hMax = INF,
            const Point& fixedVariable = Point(),
            EvalType evalType = EvalType::BB,
            FHComputeTypeS computeType = defaultFHComputeTypeS,
            const std::vector<EvalPoint>& evalPointList = std::vector<EvalPoint>(),
            bool barrierInitializedFromCache= true)
      : BarrierBase(evalType, computeType, hMax),
        _incumbentsAndHMaxUpToDate(false)
    {
        init(fixedVariable, barrierInitializedFromCache);
        init(fixedVariable,evalPointList);
    }
    
    
    
    // Copy constructor
    ProgressiveBarrier(const ProgressiveBarrier & b) : NOMAD::BarrierBase(b)
    {
        // Do not copy the barrier points. Do not initialize from cache.
    }
    
    std::shared_ptr<BarrierBase> clone() const override {
      return std::make_shared<ProgressiveBarrier>(*this);
    }
    
    /*-----------------*/
    /* Feasible points */
    /*-----------------*/

    /// Update ref best feasible and ref best infeasible values.
    void updateRefBests() override;


    /*---------------*/
    /* Other methods */
    /*---------------*/

    /// Set the hMax of the barrier
    /**
     \param hMax    The hMax -- \b IN.
    */
    void setHMax(const Double &hMax) override;

    /// SuccessType of xFeas and xInf according to given points.
    /* \param xFeas Feasible point -- \b IN.
     * \param XInf Infeasible point -- \b IN.
     * \return SuccessType of points.
     * \note Input EvalPoints are already in subproblem dimension
     */
    SuccessType getSuccessTypeOfPoints(const EvalPointPtr xFeas,
                                       const EvalPointPtr xInf) override;

    /// Update xFeas and xInf according to given points.
    /* \param evalPointList vector of EvalPoints  -- \b IN.
     * \param keepAllPoints \b IN.
     * \return true if the barrier feasible and/or infeasible incumbents are changed, false otherwise
     * \note Input EvalPoints are already in subproblem dimension
     */
    bool updateWithPoints(const std::vector<EvalPoint>& evalPointList,
                          const bool keepAllPoints = false,
                          const bool updateInfeasibleIncumbentAndHmax = false ) override;
    
    
    
    EvalPointPtr getCurrentIncumbentInf() const override;
    EvalPointPtr getCurrentIncumbentFeas() const override;

    /// Return the barrier as a string.
    /* May be used for information, or for saving a barrier. In the former case,
     * it may be useful to set parameter max to a small value (e.g., 4). In the
     * latter case, INF_SIZE_T should be used so that all points are saved.
     * \param max Maximum number of feasible and infeasible points to display
     * \return A string describing the barrier
     */
    std::vector<std::string> display(const size_t max = INF_SIZE_T) const override;

private:

    /**
     * \brief Helper function for constructor.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     \param fixedVariable   The fixed variables have a fixed value     -- \b IN.
     \param barrierInitializedFromCache  Flag to initialize barrier from cache or not. -- \b IN.
     */
    void init(const Point& fixedVariable,
              bool barrierInitializedFromCache) override;

    /**
     * \brief Helper function for constructor.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     \param fixedVariable   The fixed variables have a fixed value     -- \b IN.
     \param evalPointList   Additional points to consider to construct barrier. -- \b IN.
     */
    void init(const Point& fixedVariable,
              const std::vector<EvalPoint>& evalPointList);

    
    /** Helper for updateWithPoints
     *
     */
    bool dominates(const NOMAD::ArrayOfDouble & f1, const NOMAD::Double & h1, const NOMAD::ArrayOfDouble & f2, const NOMAD::Double & h2) const;


    /** Helper for updateWithPoints
     * Used for updating hMax
     * Get h just below hmax among all xInf
     */
    NOMAD::Double getWorstHInBarrier() const;
    
protected:

    /** Helper for updateWithPoints
        * Set the infeasible incumbent(s) from xInf
     */
    bool setInfeasibleIncumbents() ;


};




#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_5_PROGRESSIVEBARRIER__

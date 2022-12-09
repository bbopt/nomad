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
#ifndef __NOMAD_4_3_DMULTIMADSBARRIER
#define __NOMAD_4_3_DMULTIMADSBARRIER

#include "../../Eval/BarrierBase.hpp"
#include "../../Eval/EvalPoint.hpp"

#include "../../Type/BBInputType.hpp"

#include "../../nomad_nsbegin.hpp"
#include "Exception.hpp"

/// Class for multiobjective barrier following DMultiMADS algorithm
class DMultiMadsBarrier : public BarrierBase
{
private:

    EvalPointPtr _currentIncumbentFeas;  ///< current feasible of interest for DmultiMads
    EvalPointPtr _currentIncumbentInf;   ///< current infeasible of interest for DMultiMads

    ArrayOfDouble _fixedVariables; ///< The fixed variables. Use the fixed variables if the mesh and barrier points do not have the same dimension. The fixed variables are used to access the frame size for unfixed variables. This situation occurs for the EvaluatorControl barrier with fixed variables.
    
    std::vector<EvalPointPtr> _xFilterInf; ///< Stores current non dominated infeasible solutions. Can hold a copy of EvalPoint in xFeas and XInf.

    /// Dimension of the objective vectors in the barrier.
    /**
     Used for verification only.
    */
    size_t _nobj;

    /// BBInputTypes
    /*
     * Used by get mesh max frame size and  for verification (selection of frame centers by the mesh)
    */
    BBInputTypeList _bbInputsType;

    size_t _incumbentSelectionParam;


public:
    /// Constructor
    /**
     * hMax will be updated during optimization.
     \param nbObj            The number of objectives -- \b IN.
     \param hMax            The max of h to keep a point in the barrier -- \b IN.
     \param fixedVariables   The fixed variables have a fixed value -- \b IN.
     \param evalType        Type of evaluation (BB or MODEL) -- \b IN.
     \param computeType     Type of function computation (standard, phase-one or user) -- \b IN.
     \param evalPointList   Additional points to consider in building the barrier -- \b IN.
     \param barrierInitializedFromCache Flag to initialize the barrier from cache -- \b IN.
     */
    DMultiMadsBarrier(size_t nbObj,
                      const Double& hMax = INF,
              const size_t incumbentSelectionParam = 1,
              const Point& fixedVariables = Point(),
              EvalType evalType = EvalType::BB,
              ComputeType computeType = ComputeType::STANDARD,
              const std::vector<EvalPoint>& evalPointList = std::vector<EvalPoint>(),
              bool barrierInitializedFromCache= true,
              const BBInputTypeList bbInputsType= std::vector<BBInputType>())
      : BarrierBase(hMax),
        _nobj(nbObj), 
        _currentIncumbentFeas(nullptr),
        _currentIncumbentInf(nullptr),
        _fixedVariables(fixedVariables),
        _bbInputsType(bbInputsType),
        _incumbentSelectionParam(incumbentSelectionParam)
    {
        checkHMax();

        // The number of objectives is initialized via the call to the cache.
        init(fixedVariables, evalType, computeType, barrierInitializedFromCache); // Initialize with cache (if flag is true)
        init(fixedVariables, evalType, evalPointList, computeType); // Initializae with a list of points

        if (computeType == NOMAD::ComputeType::STANDARD && _nobj == 1)
        {
            std::string s = "Error: Construction of a DMultiMadsBarrier with number of objectives equal to 1. ";
            s += "In this case, use Barrier";
            throw NOMAD::Exception(__FILE__,__LINE__,s);
        }
    }

    DMultiMadsBarrier(const DMultiMadsBarrier & b)
    {
        // Do not copy the barrier points
        _nobj = b._nobj;
        _hMax = b._hMax;
        _bbInputsType = b._bbInputsType;
        _fixedVariables = b._fixedVariables;
    }
    
    std::shared_ptr<BarrierBase> clone() const override {
        return std::make_shared<DMultiMadsBarrier>(*this);
    }

    /// Update ref best feasible and ref best infeasible values.
    ///  Not used. Triggers an exception
    void updateRefBests() override ;

    void clearXFeas() override;

    ///  Get the point that is the best feasible point in the barrier
    NOMAD::EvalPointPtr getCurrentIncumbentFeas() const { return _currentIncumbentFeas ; }
    
    /// Function used to get/set the primary and secondary frame centers to generate trial points (poll and search). For DMultiMads use the barrier current incumbent.
    /**
     * If there is no feasible point, return a \c nullptr
     \return A single feasible eval point.
     */
    NOMAD::EvalPointPtr getFirstXFeas() const override {    return _currentIncumbentFeas; }
    
    ///  Get the point that is the best infeasible point in the barrier
    NOMAD::EvalPointPtr getCurrentIncumbentInf() const { return _currentIncumbentInf ; }

    // Function used to get/set the primary and secondary frame centers to generate trial points (poll and search). For DMultiMads use the barrier current incumbent.
    /**
     * If there is no infeasible point, return a \c nullptr
     \return A single infeasible eval point.
     */
    NOMAD::EvalPointPtr getFirstXInf() const override
    {
        return _currentIncumbentInf;
    }

    
    
    ///  Get the i non dominated and best feasible incumbent.
    /**
     * If there is no feasible point or index i is superior to the
     * number of feasible incumbents, raise an exception.
     \return A single feasible eval point.
     */
    EvalPoint& getXFeas(size_t i)
    {
        if (i >= _xFeas.size())
        {
            std::string s = "Error: try to get access to " + std::to_string(i) + " element of ";
            s += "DMultiMadsBarrier but number of feasible elements is " + std::to_string(_xFeas.size()) + ".";
            throw NOMAD::Exception(__FILE__,__LINE__,s);
        }
        return *_xFeas[i];
    }

    /// Remove all infeasible points (including the filter points)
    /// from the barrier.
    void clearXInf() override ;

    ///  Get the i non dominated and best infeasible incumbent.
    /**
     * If there is no infeasible point or index i is superior to the
     * number of infeasible incumbents, raise an exception.
     \return A single infeasible eval point.
     */
    EvalPoint& getXInf(size_t i)
    {
        if (i >= _xInf.size())
        {
            std::string s = "Error: try to get access to " + std::to_string(i) + " element of ";
            s += "DMultiMadsBarrier but number of infeasible elements is " + std::to_string(_xInf.size()) + ".";
            throw NOMAD::Exception(__FILE__,__LINE__,s);
        }
        return *_xInf[i];
    }

    size_t nbXFilterInf()
    {
        return _xFilterInf.size();
    }
    
    const std::vector<EvalPointPtr> & getAllXFilterInf() const
    {
        return _xFilterInf;
    }
    
    
    
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
     * \note After calling this function, all xInf elements of the MO Barrier
     * are below _hMax.
     */
    bool updateWithPoints(const std::vector<EvalPoint>& evalPointList,
                          EvalType evalType,
                          ComputeType computeType,
                          const bool keepAllPoints = false) override;
    
    
    /// Update current feas and infeas incumbents. Called by DMultiMadsUpdate and when updating the barrier
    void updateCurrentIncumbents()
    {
        updateCurrentIncumbentFeas();
        updateCurrentIncumbentInf();
    }
    

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
     \param fixedVariables   The fixed variables have a fixed value     -- \b IN.
     \param evalType        Which eval (Blackbox or Model) to use to verify feasibility  -- \b IN.
     \param computeType    Which compute type (standard, phase-one or user) must be available to find in cache  -- \b IN.
     \param barrierInitializedFromCache  Flag to initialize barrier from cache or not. -- \b IN.
     */
    void init(const Point& fixedVariables,
              EvalType evalType,
              ComputeType computeType,
              bool barrierInitializedFromCache) override;

    /**
     * \brief Helper function for constructor.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     \param fixedVariables   The fixed variables have a fixed value     -- \b IN.
     \param evalType        Which eval (Blackbox or Model) to use to verify feasibility  -- \b IN.
     \param evalPointList   Additional points to consider to construct barrier. -- \b IN.
     \param computeType    Which compute type (standard, phase-one or user) must be available to find in cache  -- \b IN.
     */
    void init(const Point& fixedVariables,
              EvalType evalType,
              const std::vector<EvalPoint>& evalPointList,
              ComputeType computeType);

    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong with the mesh parameters.
     * Will remain silent otherwise.
     */
    void checkMeshParameters(const EvalPoint &evalPoint) const;
    
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

    /**
     * \brief Helper function for infeasible point candidate search.
     */
    NOMAD::EvalPointPtr getFirstXInfNoXFeas() const;

    /**
     * \brief Helper function for infeasible point candidate search.
     */
    NOMAD::EvalPointPtr getXInfMinH() const;

    /**
     * \brief Filter infeasible incumbent solutions when hMax is set.
     *
     */
    void updateXInfAndFilterInfAfterHMaxSet();
    
    NOMAD::Double getMeshMaxFrameSize(const NOMAD::EvalPointPtr pt) const;
    
    /// Helper for updateWithPoints
    bool updateFeasWithPoint(const EvalPoint& evalPoint,
                          EvalType evalType,
                          ComputeType computeType,
                          const bool keepAllPoints );
    
    /// Helper for updateWithPoints
    bool updateInfWithPoint(const EvalPoint& evalPoint,
                          EvalType evalType,
                          ComputeType computeType,
                          const bool keepAllPoints,
                          const bool feasHasBeenUpdated);
    
    /// Helper for updateCurrentIncumbents
    void updateCurrentIncumbentFeas();
    
    /// Helper for updateCurrentIncumbents
    void updateCurrentIncumbentInf();
    
    

};

/// Display useful values so that a new Barrier could be constructed using these values.
std::ostream& operator<<(std::ostream& os, const DMultiMadsBarrier& barrier);

/// Get barrier values from stream
std::istream& operator>>(std::istream& is, DMultiMadsBarrier& barrier);

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_3_DMULTIMADSBARRIER__

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
#ifndef __NOMAD_4_4_BARRIERBASE__
#define __NOMAD_4_4_BARRIERBASE__

#include "../Eval/EvalPoint.hpp"

#include "../nomad_platform.hpp"
#include "../nomad_nsbegin.hpp"

/// Generic class for barrier following algorithm 12.2 of DFBO.
class DLL_EVAL_API BarrierBase
{
protected:

    std::vector<EvalPointPtr> _xFeas;  ///< Current feasible incumbent solutions
    std::vector<EvalPointPtr> _xInf;   ///< Current infeasible barrier points (contains infeasible incumbents) with h<=hMax
    
    std::vector<EvalPointPtr> _xIncFeas;   ///< Current feasible incumbent solutions. Can be a subset of _xFeas but for now, xIncFeas and xFeas are the same
    std::vector<EvalPointPtr> _xIncInf;   ///< Current infeasible incumbent solutions (subset of _xInf if it is defined). For now, we consider a vector (maybe DMultiMads needs it)
    
    EvalPointPtr _refBestFeas;      ///< Previous first feasible incumbent
    EvalPointPtr _refBestInf;       ///< Previous first infeasible incumbent
                                    ///< NB: can be above the barrier threshold
    
    Double _hMax;                   ///< Maximum acceptable value for h

    /// Dimension of the points in the barrier.
    /**
     * Used for verification only.
     * To be reviewed when we address category variables.
       /see _n in CacheBase.
     */
    size_t _n;

public:
    /// Constructor
    /**
     * hMax will be updated during optimization.
     \param hMax            The max of h to keep a point in the barrier -- \b IN.
     */
    BarrierBase(const Double& hMax = INF)
      : _hMax(hMax),
        _n(0)
    {}

    /// Copy Constructor
    /**
       Copy barrier some parameters. Barrier points are not copied.
     */
    BarrierBase(const BarrierBase & b)
    {
        _hMax = b._hMax;
    }
    
    // Use clone to create a barrier of the same type (for example, ProgressiveBarrier, DiscoMadsBarrier or DMultiMadsBarrier)
    virtual std::shared_ptr<BarrierBase> clone() const = 0;
    
    /*-----------------*/
    /* Feasible points */
    /*-----------------*/
    /// Get all feasible points in the barrier
    /**
     \return All the eval points that are feasible.
     */
    const std::vector<EvalPointPtr>& getAllXFeas() const { return _xFeas; }
    
    ///  Get the current incumbent feasible point in the barrier.
    /**
     * If there is no feasible point, return a \c nullptr
     \return A single feasible eval point.
     */
    virtual EvalPointPtr getCurrentIncumbentFeas() const =0;
    
    
    ///  Get all incumbent feasible points in the barrier
    /**
     \return All the eval points that are feasible incumbents (with same f and h).
     */
    const std::vector<EvalPointPtr>& getAllXIncFeas() const { return _xIncFeas; }
    
    ///  Get the point that was previously the first feasible point in the barrier.
    /**
     * If there is no feasible point, return a \c nullptr
     \return A single feasible eval point.
     */
    EvalPointPtr getRefBestFeas() const { return _refBestFeas; }
    void setRefBestFeas(const EvalPointPtr refBestFeas) { _refBestFeas = refBestFeas; }
    
    /// Update ref best feasible and ref best infeasible values.
    virtual void updateRefBests() = 0;

    /// Number of feasible points in the barrier.
    size_t nbXFeas() const {return _xFeas.size();}

    /// Remove feasible points from the barrier.
    virtual void clearXFeas();

    /*-------------------*/
    /* Infeasible points */
    /*-------------------*/
    ///  Get all infeasible points in the barrier
    /**
     \return All the eval points that are infeasible.
     */
    const std::vector<EvalPointPtr>& getAllXInf() const { return _xInf; }
    
    ///  Get all incumbent infeasible points in the barrier
    /**
     \return All the eval points that are infeasible incumbents (with same f and h).
     */
    const std::vector<EvalPointPtr>& getAllXIncInf() const { return _xIncInf; }

    
    ///  Get the currnent infeasible incumbent.
    /**
     * If there is no infeasible point, return a \c nullptr
     \return A single infeasible eval point.
     */
    virtual EvalPointPtr getCurrentIncumbentInf() const = 0;
    
    
    ///  Get the point that was previously the first infeasible point in the barrier.
    /**
     * If there is no feasible point, return a \c nullptr
     \return A single feasible eval point.
     */
    EvalPointPtr getRefBestInf() const { return _refBestInf; }
    void setRefBestInf(const EvalPointPtr refBestInf) { _refBestInf = refBestInf; }
    
    /// Number of infeasible points in the barrier.
    size_t nbXInf() const { return _xInf.size() ;}

    /// Remove infeasible points from the barrier.
    virtual void clearXInf() ;

    /*---------------*/
    /* Other methods */
    /*---------------*/
    /// Get all feasible and infeasible points ptr
    std::vector<EvalPoint> getAllPoints() const ;
    
    /// Make a copy of all feasible and infeasible points
    std::vector<EvalPointPtr> getAllPointsPtr() const ;

    /// Get first of all feasible and infeasible points.
    /** If there are feasible points, returns first feasible point.
      * else, returns first infeasible incumbent. */
    const EvalPointPtr getFirstPoint() const;
    
    /// Get the current hMax of the barrier.
    Double getHMax() const { return _hMax; }

    /// Set the hMax of the barrier
    /**
     \param hMax    The hMax -- \b IN.
     */
    virtual void setHMax(const Double &hMax) = 0;

    ///  xFeas and xInf according to given points.
    /* \param evalPointList vector of EvalPoints  -- \b IN.
     * \param keepAllPoints keep all good points, or keep just one point as in NOMAD 3 -- \b IN.
     * \return true if the Barrier was updated, false otherwise
     * \note Input EvalPoints are already in subproblem dimention
     */
    virtual SuccessType getSuccessTypeOfPoints(const EvalPointPtr xFeas,
                                               const EvalPointPtr xInf,
                                               EvalType evalType,
                                               ComputeType computeType) = 0;

    /// Update xFeas and xInf according to given points.
    /* \param evalPointList vector of EvalPoints  -- \b IN.
     * \param keepAllPoints keep all good points, or keep just one point as in NOMAD 3 -- \b IN.
     * \return true if the Barrier was updated, false otherwise
     * \note Input EvalPoints are already in subproblem dimention
     */
    virtual bool updateWithPoints(
                          const std::vector<EvalPoint>& evalPointList,
                          EvalType evalType,
                          ComputeType computeType,
                          const bool keepAllPoints = false,
                          const bool updateInfeasibleIncumbentAndHmax = false) = 0;
    
    /// Return the barrier as a string.
    /* May be used for information, or for saving a barrier. In the former case,
     * it may be useful to set parameter max to a small value (e.g., 4). In the
     * latter case, INF_SIZE_T should be used so that all points are saved.
     * \param max Maximum number of feasible and infeasible points to display
     * \return A string describing the barrier
     */
    virtual std::vector<std::string> display(const size_t max = INF_SIZE_T) const =0;
    
    
    bool findPoint(const Point & point, EvalPoint & foundEvalPoint);
    
protected:
    
    /**
     * \brief Helper function for init/constructor.
     */
    void setN();
    
    /**
     * \brief Helper function for init/setHMax.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkHMax();

    /**
     * \brief Helper function for init/constructor.
     *
     * Throw an exception if the Cache has not been instantiated yet. Will remain silent otherwise.
     */
    void checkCache();
    
    
    std::vector<EvalPointPtr>::iterator findEvalPoint(std::vector<EvalPointPtr>::iterator first, std::vector<EvalPointPtr>::iterator last, const EvalPoint & p  );
    
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
    virtual void init(const Point& fixedVariable,
                      EvalType evalType,
                      ComputeType computeType,
                      bool barrierInitializedFromCache) = 0;
    
    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkXFeas(const EvalPoint &xFeas,
                    EvalType evalType,
                    ComputeType computeType = ComputeType::STANDARD) ;
    
    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    virtual void checkXFeasIsFeas(const EvalPoint &xFeas,
                          EvalType evalType,
                          ComputeType computeType = ComputeType::STANDARD);
    
    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkXInf(const EvalPoint &xInf, EvalType evalType);
    
    
};

/// Display useful values so that a new Barrier could be constructed using these values.
DLL_EVAL_API std::ostream& operator<<(std::ostream& os, const BarrierBase& barrier);

/// Get barrier values from stream
DLL_EVAL_API std::istream& operator>>(std::istream& is, BarrierBase& barrier);

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_4_BARRIERBASE__

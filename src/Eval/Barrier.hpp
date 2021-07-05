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
#ifndef __NOMAD_4_0_BARRIER__
#define __NOMAD_4_0_BARRIER__

#include "../Eval/EvalPoint.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for barrier following algorithm 12.2 of DFBO.
class Barrier
{
private:

    std::vector<EvalPoint> _xFeas;  ///< Current feasible incumbent solutions
    std::vector<EvalPoint> _xInf;   ///< Current infeasible incumbent solutions

    EvalPointPtr _refBestFeas;      ///< Previous first feasible incumbent
    EvalPointPtr _refBestInf;       ///< Previous first infeasible incumbent

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
     \param fixedVariable   The fixed variables have a fixed value -- \b IN.
     \param evalType        Type of evaluation (BB or MODEL) -- \b IN.
     \param evalPointList   Additional points to consider in building the barrier -- \b IN.
     */
    Barrier(const Double& hMax = INF,
            const Point& fixedVariable = Point(),
            const EvalType& evalType = EvalType::BB,
            const ComputeType& computeType = ComputeType::STANDARD,
            const std::vector<EvalPoint>& evalPointList = std::vector<EvalPoint>(),
            bool barrierInitializedFromCache= true)
      : _xFeas(),
        _xInf(),
        _refBestFeas(nullptr),
        _refBestInf(nullptr),
        _hMax(hMax),
        _n(0)
    {
        init(fixedVariable, evalType, evalPointList, computeType, barrierInitializedFromCache);
    }

    /*-----------------*/
    /* Feasible points */
    /*-----------------*/
    /// Get all feasible points in the barrier
    /**
     \return All the eval points that are feasible.
     */
    const std::vector<EvalPoint>& getAllXFeas()    const { return _xFeas; }

    /// Update ref best feasible and ref best infeasible values.
    void updateRefBests();

    ///  Get the first feasible point in the barrier.
    /**
     * If there is no feasible point, return a \c nullptr
     \return A single feasible eval point.
     */
    EvalPointPtr getFirstXFeas() const;

    ///  Get the point that was previously the first feasible point in the barrier.
    /**
     * If there is no feasible point, return a \c nullptr
     \return A single feasible eval point.
     */
    EvalPointPtr getRefBestFeas() const { return _refBestFeas; }
    void setRefBestFeas(const EvalPointPtr refBestFeas) { _refBestFeas = refBestFeas; }

    /// Number of feasible points in the barrier.
    size_t nbXFeas() const { return _xFeas.size(); }

    /// Add a feasible point in the barrier.
    /**
     * If the point is feasible it is added, if not an exception is triggered.
     \param xFeas       The eval point to add -- \b IN.
     \param evalType    Which eval (Blackbox or Model) of the EvalPoint to use to verify feasibility  -- \b IN.
     */
    void addXFeas(const EvalPoint &xFeas, const EvalType& evalType, const ComputeType& computeType = ComputeType::STANDARD);

    /// Remove feasible points from the barrier.
    void clearXFeas();

    /*-------------------*/
    /* Infeasible points */
    /*-------------------*/
    ///  Get all infeasible points in the barrier
    /**
     \return All the eval points that are infeasible.
     */
    const std::vector<EvalPoint>& getAllXInf() const { return _xInf; }

    ///  Get the first infeasible point in the barrier.
    /**
     * If there is no infeasible point, return a \c nullptr
     \return A single infeasible eval point.
     */
    EvalPointPtr getFirstXInf() const;

    ///  Get the point that was previously the first infeasible point in the barrier.
    /**
     * If there is no feasible point, return a \c nullptr
     \return A single feasible eval point.
     */
    EvalPointPtr getRefBestInf() const { return _refBestInf; }
    void setRefBestInf(const EvalPointPtr refBestInf) { _refBestInf = refBestInf; }

    /// Number of infeasible points in the barrier.
    size_t nbXInf() const { return _xInf.size(); }

    /// Add an infeasible point in the barrier.
    /**
     * If the point is nullptr an exception is triggered.
     \param xInf   The eval point to add -- \b IN.
     */
    void addXInf(const EvalPoint &xInf, const EvalType& evalType);

    /// Remove infeasible points from the barrier.
    void clearXInf();

    /*---------------*/
    /* Other methods */
    /*---------------*/
    /// Get all feasible and infeasable points
    std::vector<EvalPoint> getAllPoints() const;

    /// Get first of all feasible and infeasible points.
    /** If there are feasible points, returns first feasible point.
      * else, returns first infeasible point. */
    const EvalPoint& getFirstPoint() const;

    /// Get the current hMax of the barrier.
    Double getHMax() const { return _hMax; }

    /// Set the hMax of the barrier
    /**
     \param hMax    The hMax -- \b IN.
     */
    void setHMax(const Double &hMax);

    ///  xFeas and xInf according to given points.
    /* \param evalPointList vector of EvalPoints  -- \b IN.
     * \param keepAllPoints keep all good points, or keep just one point as in NOMAD 3 -- \b IN.
     * \return true if the Barrier was updated, false otherwise
     * \note Input EvalPoints are already in subproblem dimention
     */
    SuccessType getSuccessTypeOfPoints(const std::shared_ptr<EvalPoint> & xFeas,
                                       const std::shared_ptr<EvalPoint> & xInf,
                                       const EvalType& evalType,
                                       const ComputeType& computeType);

    /// Update xFeas and xInf according to given points.
    /* \param evalPointList vector of EvalPoints  -- \b IN.
     * \param keepAllPoints keep all good points, or keep just one point as in NOMAD 3 -- \b IN.
     * \return true if the Barrier was updated, false otherwise
     * \note Input EvalPoints are already in subproblem dimention
     */
    bool updateWithPoints(const std::vector<EvalPoint>& evalPointList,
                          const EvalType& evalType,
                          const ComputeType& computeType,
                          const bool keepAllPoints = false);

    /// Return the barrier as a string.
    /* May be used for information, or for saving a barrier. In the former case,
     * it may be useful to set parameter max to a small value (e.g., 4). In the
     * latter case, INF_SIZE_T should be used so that all points are saved.
     * \param max Maximum number of feasible and infeasible points to display
     * \return A string describing the barrier
     */
    std::string display(const size_t max = INF_SIZE_T) const;

private:

    /**
     * \brief Helper function for constructor.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     \param fixedVariable   The fixed variables have a fixed value     -- \b IN.
     \param evalType        Which eval (Blackbox or Model) to use to verify feasibility  -- \b IN.
     \param evalPointList   Additional points to consider to construct barrier. -- \b IN.
     \param barrierInitializedFromCache  Flag to initialize barrier from cache or not. -- \b IN.
     */
    void init(const Point& fixedVariable,
              const EvalType& evalType,
              const std::vector<EvalPoint>& evalPointList,
              const ComputeType& computeType,
              bool barrierInitializedFromCache);

    /**
     * \brief Helper function for init/constructor.
     */
    void setN();

    /**
     * \brief Helper function for init/constructor.
     *
     * Throw an exception if the Cache has not been instantiated yet. Will remain silent otherwise.
     */
    void checkCache();

    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkXFeas(const EvalPoint &xFeas,
                    const EvalType& evalType,
                    const ComputeType& computeType = ComputeType::STANDARD);

    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkXFeasIsFeas(const EvalPoint &xFeas,
                          const EvalType& evalType,
                          const ComputeType& computeType = ComputeType::STANDARD);

    /**
     * \brief Helper function for insertion.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkXInf(const EvalPoint &xInf, const EvalType& evalType);

    /**
     * \brief Helper function for init/setHMax.
     *
     * Will throw exceptions or output error messages if something is wrong. Will remain silent otherwise.
     */
    void checkHMax();
};


/// Display useful values so that a new Barrier could be constructed using these values.
std::ostream& operator<<(std::ostream& os, const Barrier& barrier);

/// Get barrier values from stream
std::istream& operator>>(std::istream& is, Barrier& barrier);


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_BARRIER__

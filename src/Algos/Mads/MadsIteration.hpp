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
#ifndef __NOMAD_4_0_MADSITERATION__
#define __NOMAD_4_0_MADSITERATION__

#include "../../Algos/Iteration.hpp"
#include "../../Algos/MeshBase.hpp"
#include "../../Eval/EvalPoint.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for MADS iteration
/**
 A MADS iteration consists of a Search step followed by a Poll step depending on the stop reasons and successes.
 */
class MadsIteration: public Iteration
{
private:
    const std::shared_ptr<MeshBase>  _mesh;        ///< Mesh on which the points are
    SuccessType                      _success;     ///< Success type of this iteration

#ifdef TIME_STATS
    /// Time counters
    DLL_ALGO_API static double  _iterTime;          ///< Total time spent running this class
    DLL_ALGO_API static double  _searchTime;        ///< Total time spent running searches
    DLL_ALGO_API static double  _searchEvalTime;    ///< Total time spent evaluating search points
    DLL_ALGO_API static double  _pollTime;          ///< Total time spent running polls
    DLL_ALGO_API static double  _pollEvalTime;      ///< Total time spent evaluating poll points
    double                      _iterStartTime;     ///< Time at which the start method was called
#endif // TIME_STATS

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param k                  The iteration number -- \b IN.
     \param mesh               The mesh of the iteration -- \b IN.
     */
    explicit MadsIteration(const Step *parentStep,
                           const size_t k,
                           const std::shared_ptr<MeshBase> mesh)
      : Iteration(parentStep, k),
        _mesh(mesh),
        _success(SuccessType::NOT_EVALUATED)
#ifdef TIME_STATS
        ,_iterStartTime(0.0)
#endif // TIME_STATS
    {
        init();
    }

    
    NOMAD::ArrayOfPoint suggest() override;

    // Gets/Sets

    /**
     The Mads algorithm iteration possesses a mesh, unlike the base iteration that has none.
     \remark Used by Step::getIterationMesh() to pass the mesh whenever needed
     */
    const std::shared_ptr<MeshBase> getMesh() const override { return _mesh; }

    /// Return current SuccessType
    const SuccessType& getSuccessType() const { return _success; }

    /// Set SuccessType member
    void setSuccessType(const SuccessType& success) { _success = success; }

#ifdef TIME_STATS
    /// Time stats
    static double getIterTime()         { return _iterTime; }
    static double getSearchTime()       { return _searchTime; }
    static double getSearchEvalTime()   { return _searchEvalTime; }
    static double getPollTime()         { return _pollTime; }
    static double getPollEvalTime()     { return _pollEvalTime; }
#endif // TIME_STATS

    /*---------------------*/
    /* Other class methods */
    /*---------------------*/


private:
    /// Helper for constructor
    void init();

    virtual void startImp() override;

    /// Implementation of the run tasks of MADS algorithm.
    /**
     Run a MADS iteration: a Search step followed by a Poll step depending on the stop reasons and successes.
     */
    virtual bool runImp() override;

#ifdef TIME_STATS
    virtual void endImp() override;
#endif // TIME_STATS
};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_MADSITERATION__

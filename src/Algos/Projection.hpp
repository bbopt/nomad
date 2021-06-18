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

#ifndef __NOMAD_4_0_PROJECTION__
#define __NOMAD_4_0_PROJECTION__

#include "../Algos/IterationUtils.hpp"
#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"


/// Manage projection of points to mesh
/**
Generalization from projection used in sgtelib in NOMAD 3. See file: Sgtelib_Model_Search.cpp
Currently not implemented, not used, does nothing.
Should be generic for any case that needs projection.
\todo Implement and use
*/
class Projection : public Step, public IterationUtils
{
private:
    EvalPointSet                _oraclePoints;
    OutputLevel                 _displayLevel;

    // Vector of EvalPoints which have a model eval
    std::vector<EvalPoint>      _cacheModelEval;

    // Mesh and frame center to project on
    std::shared_ptr<MeshBase>   _mesh;
    std::shared_ptr<EvalPoint>  _frameCenter;

    std::set<size_t>            _indexSet;
    size_t                      _nbProjTrial;

public:
    // Constructor
    explicit Projection(const Step* parentStep,
                        const EvalPointSet &oraclePoints);

    // Destructor
    virtual ~Projection();

    // Get / Set
    EvalPointSet getOraclePoints() const { return getTrialPoints(); }

private:
    void init();
    void buildIndexSet(const size_t n);

    virtual void startImp() override;
    virtual bool runImp() override;
    virtual void endImp() override;

    void generateTrialPoints() override;

    void projectPoint(const EvalPoint& oraclePoint);

    void nonProjectedPoint(const EvalPoint& oraclePoint);

    void stdProjectedPoint(const EvalPoint& oraclePoint);

    Direction computePerturbation(const EvalPoint& oraclePoint, size_t index);

    EvalPoint buildProjectionTrialPoint(const Point& xRef, const Direction& perturbation);

    void doGreedySelection(const EvalPoint& oraclePoint,
                           const EvalPointSet& trySet,
                           std::vector<bool>& keep);

/*
    void evaluateProjectionTrialPoints(const EvalPointSet& trySet,
                           const Evaluator& ev,
                           const std::vector<bool>& keep,
                           EvalPoint& bestEvalPoint);
    */

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_PROJECTION__

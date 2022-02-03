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
#ifndef __NOMAD_4_2_ALGOSTATS__
#define __NOMAD_4_2_ALGOSTATS__

#include "../Algos/Step.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Type/EvalType.hpp"
#include "../Util/Exception.hpp"
#include "../Util/fileutils.hpp"
#include "../Util/utils.hpp"

#include "../nomad_nsbegin.hpp"


/// Class to manage the trial point stats (Search methods and Poll methods)
/**
 */
class AlgoStats
{
private:
    size_t _nbEvalsDoneBB;
    size_t _nbEvalsDoneModels;
    size_t _nbEvalsDoneSurrogate;

    size_t _nbTrialPointsGeneratedBB;
    size_t _nbTrialPointsGeneratedModels;
    size_t _nbTrialPointsGeneratedSurrogate;
    
    size_t _nbGenerateCalls;
    
#ifdef TIME_STATS
    size_t _totalRealAlgoTime;
    double _startTime;
    double _totalCPUAlgoTime;
#endif // TIME_STATS

    
public:
    /// Constructor
    /**
     \param parentStep          The Step parent -- \b IN.
     */
    explicit AlgoStats(const Step* parentStep):
    _nbEvalsDoneBB(0),
    _nbEvalsDoneModels(0),
    _nbEvalsDoneSurrogate(0),
    _nbTrialPointsGeneratedBB(0),
    _nbTrialPointsGeneratedModels(0),
    _nbTrialPointsGeneratedSurrogate(0),
    _nbGenerateCalls(0)
#ifdef TIME_STATS
    ,_totalRealAlgoTime(0.0),
    _startTime(0.0),
    _totalCPUAlgoTime(0.0)
#endif // TIME_STATS
    {
    }
    
    void incrementEvalsDone(size_t nb, EvalType  evalType);

    void incrementTrialPointsGenerated(size_t nb, EvalType evalType);
    
    void incrementGenerateCalls() { _nbGenerateCalls ++ ;}
    
    size_t getNbEvalsDone(EvalType evalType) const;
    
    size_t getNbTrialPointsGenerated(EvalType  evalType) const;
        
    void startCounting();
    void endCounting();
    
    std::string display() const;
    
};
///   Display method stats.
std::ostream& operator<<(std::ostream& os, const AlgoStats& stats);

/// Get the mesh values from stream
// std::istream& operator>>(std::istream& is, TrialPointStats& stats);


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_2_ALGOSTATS__

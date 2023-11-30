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
#ifndef __NOMAD_4_4_TEMPLATEALGOSEARCHMETHOD__
#define __NOMAD_4_4_TEMPLATEALGOSEARCHMETHOD__

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/SearchMethodAlgo.hpp"
#include "../../Algos/TemplateAlgo/TemplateAlgo.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class to perform a dummy random search method using an algorithm.
/**
 Randomly generates trial points.
 
 Can be used as a TEMPLATE for a new search method: copy and rename the file and the class name. Adapt the code to your needs. It is IMPORTANT to register the new search method in ../Algos/Mads/Search.cpp (NOMAD::Search::init()).
 
 
 */
class TemplateAlgoSearchMethod final : public SearchMethodAlgo
{
private:
    
    // TEMPLATE for a new search method: A new StopType must be defined in ../Util/StopReason.hpp and ../Util/StopReason.hpp to store the new algorithm accepted stop reasons.
    std::shared_ptr<AlgoStopReasons<RandomAlgoStopType>> _randomAlgoStopReasons;
    
    // TEMPLATE for a new search method: A new Algo must be defined. Create a directory in ../Algos and use the files provided ../Algos/TemplateAlgo to create your own Algorithm.
    std::unique_ptr<TemplateAlgo> _randomAlgo;

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit TemplateAlgoSearchMethod(const Step* parentStep )
      : SearchMethodAlgo(parentStep ),
        _randomAlgoStopReasons(nullptr),
        _randomAlgo(nullptr)
    {
        init();
    }


    /**
     Execute (start, run, end) of the template algorithm (random). Returns a \c true flag if the algorithm found better point.
     */
    virtual bool runImp() override ;


private:

    /// Helper for constructor.
    /**
     Test if search is enabled or not. Set the maximum number of trial points.
     */
    void init();

    /// Generate new points (no evaluation)
    /**
     \copydoc SearchMethodAlgo::generateTrialPointsFinal 
     Iterative random generation of trial points
     */
     void generateTrialPointsFinal() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_4_TEMPLATEALGOSEARCHMETHOD__


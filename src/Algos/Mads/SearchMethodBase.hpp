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
#ifndef __NOMAD_4_0_SEARCHMETHODBASE__
#define __NOMAD_4_0_SEARCHMETHODBASE__

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for generic search method of MADS. Run by Search.
/**
 Pure virtual class from which SearchMethodSimple and SearchMethodAlgo derive.
 */
class SearchMethodBase: public Step, public IterationUtils
{
private:

    bool _enabled; ///< Should this search method be used? Modified by parameters.

    std::string _comment; ///<  Comment shown when a search method is used

public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SearchMethodBase( const Step* parentStep )
      : Step( parentStep ),
        IterationUtils ( parentStep ),
        _enabled(true),
        _comment("")
    {
        init();
    }

    bool isEnabled() const { return _enabled; }
    void setEnabled(const bool enabled) { _enabled = enabled; }

    const std::string& getComment() const { return _comment; }
    bool hasComment() const { return (!_comment.empty()); }
    void setComment(const std::string& comment) { _comment = comment; }

    /**
     - Pure virtual function.
     - The implementation of startImp function in the derived class generates trial points  (in SearchMethodSimple) OR does nothing (in SearchMethodAlgo).
     */
    virtual void startImp() override =0 ;

    /**
     - Pure virtual function.
     - The implementation of runImp function in the derived class evaluates the trial points (in SearchMethodSimple) OR launches an algo (in SearchMethodAlgo).
     */
    virtual bool runImp() override = 0 ;

    /// Implementation of endImp (not virtual)
    /**
        Call to the postProcessing function to update the Barrier
    */
    void endImp() override ;

    /// Intermediate function (not yet implementation that can generate the trial points)
    /**
     - Display before and after generation comments.
     - Launches the implementation of the search method to generate the trial points (::generateTrialPointsImp).
     - Snap the points to bounds and mesh.
     */
    void generateTrialPoints() override;

    /**
     - Pure virtual function.
     - See derived classes (SearchMethodSimple, SearchMethodAlgo) for implementations.
     */
    virtual void generateTrialPointsImp() = 0 ;


protected:
    void init();

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_SEARCHMETHODBASE__


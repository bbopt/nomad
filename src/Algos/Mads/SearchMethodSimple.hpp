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
#ifndef __NOMAD_4_5_SEARCHMETHODSIMPLE__
#define __NOMAD_4_5_SEARCHMETHODSIMPLE__

#include "../../Algos/Mads/SearchMethodBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for simple search methods (no iterations) of MADS. Run by Search.
/**
 - Pure virtual class. Final class derived of this must implement ::generateTrialPointsImp.
  - The trial points information is completed (model or surrogate evals used for sorting) is completed before evaluation.
 - The trial points evaluation (derived classes) are performed when ::runImp is called.
 - Projection on mesh and bounds is performed after ::generateTrialPointsImp is called by the base class SearchMethodBase.
 */
class SearchMethodSimple: public SearchMethodBase
{
    
protected:
    bool _dynamicEnabled; ///< Should this simple search method generate trial points? Maybe modified dynamically by method (based on success, or alternating rule, ...).

    
public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SearchMethodSimple( const Step* parentStep )
      : SearchMethodBase( parentStep ),
    _dynamicEnabled(true)
    {
    }

    /// Intermediate function called to generate the trial points
    /**
     - Call for the intermediate base SearchMethodBase::generateTrialPoints (call generateTrialPointsImp, snap on bounds and mesh).
     - Complete trial points information (sorting is done before evaluation)
     - Sanity check on generated trial points
     - Update the points with frame center
     */
    void startImp() override;

    /// Function called to evaluate the trial points
    /**
     - Evaluate the trial points and update the barrier.
     - The projection of trial points on bounds and on mesh is performed before this function is called and after the function SearchMethodBase::generateTrialPointsImp is called.
     */
    bool runImp() override;
    
    /// Implementation of endImp (not virtual)
    /**
        Must call for SearchMethodBase endImp
        Maybe '
    */
    void endImp() override ;
    
protected:
    
    virtual void updateDynamicEnabled() {} ;
    
    virtual void updateAtStepEnd() {} ;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_SEARCHMETHODSIMPLE__


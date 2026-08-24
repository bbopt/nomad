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

#ifndef __NOMAD_4_6_ADSUPDATE__
#define __NOMAD_4_6_ADSUPDATE__

#include "../../Algos/Mads/MadsUpdate.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for Step 1. of ADS algorithm: parameter update, inherited from MadsUpdate but with punctured space check for PARTIAL_SUCCESS
/**
 The update is performed when calling the AdsUpdate::run function.
 For PARTIAL_SUCCESS, it checks if the point is in punctured space:
 - If yes, add its tag to the barrier's tag list
 - If no, change PARTIAL_SUCCESS to UNSUCCESSFUL
 */
class AdsUpdate: public MadsUpdate
{
    
public:
    // Constructor
    explicit AdsUpdate(const Step* parentStep)
      : MadsUpdate(parentStep)
    {
        init();
    }

private:

    /// Helper for constructor to check for valid ancestor.
    void init();

    /// No implementation is required for start.
    virtual void    startImp() override {}

    /// Implementation of the run tasks.
    /**
     Same as MadsUpdate::runImp(), but for PARTIAL_SUCCESS:
     - Checks if the point is in punctured space
     - If yes, adds its tag to the barrier's tag list
     - If no, changes PARTIAL_SUCCESS to UNSUCCESSFUL (treats it as unsuccessful)
     */
    virtual bool    runImp()   override;


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_6_ADSUPDATE__


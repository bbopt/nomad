/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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
#ifndef __NOMAD400_SGTELIBSEARCHMETHOD__
#define __NOMAD400_SGTELIBSEARCHMETHOD__

#include "../../Algos/Mads/SearchMethod.hpp"
#ifdef USE_SGTELIB
#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#endif

#include "../../nomad_nsbegin.hpp"

// \class SgtelibSearchMethod: Search method using library Sgtelib
class SgtelibSearchMethod: public SearchMethod
{
private:
    OutputLevel _displayLevel;
#ifdef USE_SGTELIB
    std::shared_ptr<SgtelibModel> _modelAlgo;
#endif

    // TODO finish implementing and use Projection class
    /// Get best projection
    /**
     \param  incumbent      The incumbent             -- \b IN.
     \param  deltaMeshSize  Mesh size parameter       -- \b IN.
     \param  x              The oracle point          -- \b IN/OUT.
     */
    void getBestProjection(const Point& incumbent,
                           const ArrayOfDouble& deltaMeshSize,
                           std::shared_ptr<Point> x);


/*----------------------------------------------------------------------*/


public:
    /// Constructor
    /**
     /param parentStep      The parent of this search step -- \b IN.
     */
    explicit SgtelibSearchMethod(const Step* parentStep )
      : SearchMethod(parentStep),
        _displayLevel(OutputLevel::LEVEL_NORMAL)
#ifdef USE_SGTELIB
        ,_modelAlgo(nullptr)
#endif
    {
        init();
    }

private:
    void init();

    /// Generate new points to evaluate
    void generateTrialPoints() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_SGTELIBSEARCHMETHOD__


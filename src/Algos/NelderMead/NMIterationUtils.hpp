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

#ifndef __NOMAD_4_0_NMITERATIONUTILS__
#define __NOMAD_4_0_NMITERATIONUTILS__

#include "../../Algos/IterationUtils.hpp"
#include "../../Algos/NelderMead/NMIteration.hpp"
#include "../../Algos/NelderMead/NMSimplexEvalPoint.hpp"
#include "../../Algos/Step.hpp"

#include "../../nomad_nsbegin.hpp"


/// Class of utils for NM iterations.
/**
 - Manage the simplex: update the characteristics (diameter, volume and normalized volume). The diameter is max(distance(y_i,y_j)). The volume is det(y_k-y_0)/!n (k=1,..n). The normalized volume is volume/diameter^n. \n

 - Hold a variable NMIterationUtils::_currentStepType of type ::StepType for the phase of Nelder Mead algorithm.
 - Calculate the rank of DZ=[y_i-y_0] using NMIterationUtils::_rankEps as trigger (see ::getRank function)

 */
class NMIterationUtils : public IterationUtils
{
private:

    // Simplex characteristics
    double _simplexDiam, _simplexVol, _simplexVon ;
    const EvalPoint * _simplexDiamPt1; /// First point used for simplex diameter
    const EvalPoint * _simplexDiamPt2; /// Second point used for simplex diameter

    ArrayOfDouble _Delta;  /// Delta mads frame size for scaling DZ (can be undefined)

    /// Helper for NMIterationUtils::updateYCharacteristics
    void updateYDiameter ( void );

protected:
    /// The precision for the rank calculation. Default is ::DEFAULT_EPSILON.
    Double _rankEps;

    /// The step type (REFLECT, EXPAND, INSIDE_CONTRACTION, OUTSIDE_CONTRACTION)
    StepType _currentStepType;

    std::shared_ptr<NMSimplexEvalPointSet> _nmY;  ///< The Nelder Mead simplex.

    /// Update the simplex diameter and volumes from NMIterationUtils::_nmY
    void updateYCharacteristics ( void ) ;

    /// Display all the characteristics of a simplex
    /**
     The diameter is max(distance(y_i,y_j)). The volume is det(y_k-y_0)/!n (k=1,..n). The normalized volume is volume/diameter^n.
     */
    void displayYInfo ( void ) const ;

    /**
     \return The rank of DZ=[y_i-y_0]/Delta (Delta can be ones if mesh is not available)
     */
    int getRankDZ ( ) const ;

    /// Set the stop reason according to NMIterationUtils::_currentStepType
    void setStopReason ( ) const;

public:
    /// Constructor
    /**
     The simplex is obtained from NMIteration.

     \param parentStep      The calling iteration Step.
     */
    explicit NMIterationUtils(const Step* parentStep)
      : IterationUtils(parentStep),
        _simplexDiam(0),
        _simplexVol(0),
        _simplexVon(0),
        _simplexDiamPt1(nullptr),
        _simplexDiamPt2(nullptr),
        _Delta(ArrayOfDouble()),
        _rankEps(DEFAULT_EPSILON),
        _currentStepType(StepType::NM_UNSET),
        _nmY(nullptr)
    {
        auto iter = dynamic_cast<const NMIteration*>(_iterAncestor);
        if ( nullptr != iter )
        {
            _nmY = iter->getY();

            // If a Mads mesh is available, initialize the frame size.
            auto madsMesh = iter->getMesh();
            if ( nullptr !=
                madsMesh )
            {
                _Delta = madsMesh->getDeltaFrameSize();
            }
        }
    }


};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_NMITERATIONUTILS__

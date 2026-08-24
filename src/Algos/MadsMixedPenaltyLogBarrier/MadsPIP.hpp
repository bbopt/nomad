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

#ifndef __NOMAD_4_6_MADSPIP__
#define __NOMAD_4_6_MADSPIP__

#include "../../Algos/AlgoCallback.hpp"
#include "../../Algos/Mads/Mads.hpp"

#include "../../nomad_nsbegin.hpp"


/// The (M)esh (A)daptive (D)irect (S)earch algorithm with a Mixed Penalty Logarithmic Barrier approach.
/**
 */
class MadsPIP: public Mads
{
private:
    // G^log and G^ext
    // IMPORTANT: the _Glog, _GextI and _GextE are indices based on ALL outputs (obj, pb, eq).
    std::list<size_t> _Glog, _GextI, _GextE; // indices of Glog and Gext constraints
    
    
    // Beta: for rho^beta when testing and updating rho
    const double _beta = 1.000000001;
    
    // b_cint_rho: CInt constant for testing and updating rho  --> maybe make it a parameter
    const double _b_cint_rho = 1E10;
    
    // b_rhobeta_rho: rho^beta constant for testing and updating rho
    double _b_rhobeta_rho ; // Initialized with Nomad parameter
    
    // Rho: Constraint penalization weights.
    NOMAD::Double _rho=0.1;
    // theta_rho: constant to update rho based on frame size test.
    const double _theta_rho = 0.01;
    const NOMAD::Double _b_rho_cint = 1.0;
    NOMAD::Double _b_rho_cext; // Initialized when knowing f(x0): _b_rho_cext = max(1,10^log(|f(x0)|)) or set to 1 if f(x0)=0
    
    // Settings to control algorithm
    
    // Ext to log switch
    const bool _extToLogSwitch = true;
    
    // Update best point whenever rho ext or rho log changes
    const bool _updateBestPoint = true ;
    
    double _logSwitchFeasThres ; // Initialized with Nomad parameter. Used for strengthen feasibility
    

public:
    /// Constructor
    /**
     \param parentStep          The parent of this step -- \b IN.
     \param stopReasons         The stop reasons for MADS-PIP -- \b IN.
     \param runParams           The run parameters that control MADS-PIP -- \b IN.
     \param pbParams            The problem parameters that control MADS-PIP -- \b IN.
     */
    explicit MadsPIP(const Step* parentStep,
                  std::shared_ptr<AlgoStopReasons<MadsStopType>> stopReasons,
                  const std::shared_ptr<RunParameters>& runParams,
                  const std::shared_ptr<PbParameters>& pbParams )
      : Mads(parentStep, stopReasons, runParams, pbParams, false /* barrier not initialized from cache, use X0s*/ )
    {
        init();
    }

private:
    ///  Initialization of class, to be used by Constructor.
    void init();

    virtual void startImp() override;
    
    // Private member for mega iteration callback.
    // Has access to MadsPIP member variables
    void megaIterationCallback(const NOMAD::Step& step,
                               bool &stop) ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_6_MADSPIP__

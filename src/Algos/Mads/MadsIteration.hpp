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
#ifndef __NOMAD400_MADSITERATION__
#define __NOMAD400_MADSITERATION__

#include "../../Algos/Iteration.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for MADS iteration
/**
 A MADS iteration consists of a Search step followed by a Poll step depending on the stop reasons and successes.
 */
class MadsIteration: public Iteration
{
private:
    const std::shared_ptr<EvalPoint> _frameCenter; ///< Center around which the points are generated
    const std::shared_ptr<MeshBase>  _mesh;        ///< Mesh on which the points are
    SuccessType                      _success;     ///< Success type of this iteration

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param frameCenter        Frame center of this iteration -- \b IN.
     \param k                  The iteration number -- \b IN.
     \param mesh               The mesh of the iteration -- \b IN.
     */
    explicit MadsIteration(const Step *parentStep,
                           const std::shared_ptr<EvalPoint>& frameCenter,
                           const size_t k,
                           const std::shared_ptr<MeshBase> mesh)
      : Iteration(parentStep, k),
        _frameCenter(frameCenter),
        _mesh(mesh),
        _success(SuccessType::NOT_EVALUATED)
    {
        init();
    }


    /// Destructor
    /**
     When iteration is done, Flush prints output queue.
     */
    virtual ~MadsIteration()
    {
        OutputQueue::Flush();
    }

    // Gets/Sets

    /**
     The Mads algorithm iteration possesses a frame center, unlike the base iteration that has none.
     \remark Used by Step::getIterationFrameCenter() to pass the frame center whenever needed
     */
    const std::shared_ptr<EvalPoint> getFrameCenter() const override { return _frameCenter; }

    /**
     The Mads algorithm iteration possesses a mesh, unlike the base iteration that has none.
     \remark Used by Step::getIterationMesh() to pass the mesh whenever needed
     */
    const std::shared_ptr<MeshBase> getMesh() const override { return _mesh; }

    /// Return current SuccessType
    const SuccessType& getSuccessType() const { return _success; }

    /// Set SuccessType member
    void setSuccessType(const SuccessType& success) { _success = success; }

    /*---------------------*/
    /* Other class methods */
    /*---------------------*/

    /// Is this the main iteration of the current MegaIteration?
    bool isMainIteration() const override;


private:
    /// Helper for constructor
    void init();
    
    /// No start task is required.
    virtual void startImp() override {}
    
    /// Implementation of the run tasks of MADS algorithm.
    /**
     Run a MADS iteration: a Search step followed by a Poll step depending on the stop reasons and successes.
     */
    virtual bool runImp() override ;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD400_MADSITERATION__

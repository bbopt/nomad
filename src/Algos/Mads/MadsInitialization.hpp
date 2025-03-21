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
#ifndef __NOMAD_4_5_MADSINITIALIZATION__
#define __NOMAD_4_5_MADSINITIALIZATION__

#include "../../Algos/Initialization.hpp"
#include "../../Type/BBInputType.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for Mads initialization (step 0)
/**
 The run function of this step validates and evaluates X0(s).
 Initialization of the mesh is performed at construction.
 */
class MadsInitialization : public Initialization
{
    
private:
    
    BBInputTypeList _bbInputType;
    
    Double _hMax0;  ///< Initial HMax for Mads barrier
    
protected:
    MeshBasePtr _initialMesh;
    
    bool _barrierInitializedFromCache;
    bool _isUsedForDMultiMads;
    bool _isUsedForDiscoMads;

public:
    /// Constructor
    /*
     \param parentStep                   The parent of this step -- \b IN.
     \param barrierInitializedFromCache  Flag to initialize barrier from cache or not -- \b IN.
     */
    explicit MadsInitialization(const Step* parentStep, bool barrierInitializedFromCache=true, bool isUsedForDMultiMads=false, bool isUsedForDiscoMads=false)
      : Initialization(parentStep),
        _initialMesh(nullptr),
        _barrierInitializedFromCache(barrierInitializedFromCache),
        _isUsedForDMultiMads(isUsedForDMultiMads),
        _isUsedForDiscoMads(isUsedForDiscoMads)
    {
        init();
    }

    virtual ~MadsInitialization() {}

    MeshBasePtr getMesh() const { return _initialMesh; }

private:
    void init();

    bool eval_x0s();

protected:
    
    virtual bool runImp() override;

};

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_MADSINITIALIZATION__

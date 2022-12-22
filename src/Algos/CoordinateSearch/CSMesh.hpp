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


#ifndef __NOMAD_4_3_CSMESH__
#define __NOMAD_4_3_CSMESH__

#include "../../Eval/MeshBase.hpp"

#include "../../nomad_nsbegin.hpp"

/// Class for the granular mesh of CS.
/**
 This class derives from MeshBase.
 The class manages the mesh size (delta) and the frame size (Delta) for the discretization of the variable space used by CS. Each variable has its own mesh and frame sizes which but the "cell" aspect ratios is not changed in CS (unlike in Mads) as the mesh is not anisotropic. The frame size (and mesh size) for all variables can be enlarged or decreased (see GMesh::refineDeltaFrameSize and GMesh::enlargeDeltaFrameSize). A given point can be projected on the the mesh using GMesh::scaleAndProjectOnMesh. \n

 The frame size for each variable is parameterized with one or two attributes: CSMesh::_frameSize, and CSMesh::_granularity (Delta = gran * frameSize). The first attribute is for variable having a specified minimal granularity (for example, integers have a minimal granularity of 1). This ensures that variables are always a multiple of the granularity if it is defined. The mesh size is delta = Delta/2.
 
 */
class CSMesh: public MeshBase
{

    ArrayOfDouble        _initFrameSize;  ///< The initial frame size
    ArrayOfDouble        _frameSize;  ///< The current frame size
    const ArrayOfDouble  _granularity;  ///< The fixed granularity of the mesh

public:

    /// Constructor
    /**
     \param parameters  The problem parameters attributes control the mesh mechanics -- \b IN.
     */
    explicit CSMesh(std::shared_ptr<PbParameters> parameters)
      : MeshBase(parameters),
        _initFrameSize(ArrayOfDouble()),
        _frameSize(ArrayOfDouble()),
        _granularity(parameters->getAttributeValue<ArrayOfDouble>("GRANULARITY"))
    {
        init();
    }

    // Clone a CSMesh and return a pointer to MeshBase
    std::unique_ptr<MeshBase> clone() const override {
      return std::make_unique<CSMesh>(*this);
    }
    
    /*-----------*/
    /* Get / Set */
    /*-----------*/
    const ArrayOfDouble& getInitFrameSize() const { return _initFrameSize; }
    const ArrayOfDouble& getFrameSize() const { return _frameSize; }
    const ArrayOfDouble& getGranularity() const { return _granularity; }


    /*----------------------*/
    /* Other Class Methods */
    /*----------------------*/


    void checkMeshForStopping( std::shared_ptr<AllStopReasons> algoStopReason ) const override;
    

    void updatedeltaMeshSize() override;

    /**
     \copydoc MeshBase::enlargeDeltaFrameSize
     \note This implementation relies on CSMesh::_frameSize,  and CSMesh::_granularity.
     */
    bool enlargeDeltaFrameSize(const Direction& direction) override;

    /**
     \copydoc MeshBase::refineDeltaFrameSize
     \note This implementation relies on CSMesh::_frameSize Frame size is updated . Upon failure, the frame sizes for all dimensions are refined in the same way (no need for direction and anisotropy factor).
     */
    void refineDeltaFrameSize() override;

private:

    /// Helper for refineDeltaFrameSize().
    /**
     \param frameSize     The frame size  to update -- \b IN/OUT.
     \param granularity      The granularity of the mesh -- \b IN.
     */
    void refineDeltaFrameSize(Double &frameSize,
                              const Double &granularity) const;
public:
    Double getRho(const size_t i) const override;
    ArrayOfDouble getRho() const override { return MeshBase::getRho(); }

    ArrayOfDouble getdeltaMeshSize() const override;
    Double getdeltaMeshSize(const size_t i) const override;
private:

    /// Helper function
    Double getdeltaMeshSize(const Double& frameSize,
                            const Double& granularity) const;
public:

    //
    // The documentation of overriden function is provided in the base class.
    //

    ArrayOfDouble getDeltaFrameSize() const override;
    Double getDeltaFrameSize(const size_t i) const override;
    ArrayOfDouble getDeltaFrameSizeCoarser() const override;
    Double getDeltaFrameSizeCoarser(const size_t i) const override;

private:
    /// Helper function
    Double getDeltaFrameSize(const Double& granularity, const Double& frameSize) const;

public:

    void setDeltas(const ArrayOfDouble &deltaMeshSize,
                   const ArrayOfDouble &deltaFrameSize) override;

    void setDeltas(const size_t i,
                   const Double &deltaMeshSize,
                   const Double &deltaFrameSize) override;

public:

    // Add some documentation in addition to the base class
    /**
     \copydoc MeshBase::scaleAndProjectOnMesh
     \note This implementation relies on CSMesh::_frameSize
     */
    Double scaleAndProjectOnMesh(size_t i, const Double &l) const override;

    ArrayOfDouble scaleAndProjectOnMesh(const Direction &dir) const override;

    /**
     * Project the point on the mesh centered on frameCenter. No scaling.
     */
    Point projectOnMesh(const Point& point, const Point& frameCenter) const override;


private:
    /// Helper for constructor.
    void init();

    /// Initialization of granular frame size mantissa and exponent
    /**
     \param contInitFrameSize        continuous initial frame size   -- \b IN.
    */
    void initFrameSizeGranular(const ArrayOfDouble &contInitFrameSize);

};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_3_GMESH__


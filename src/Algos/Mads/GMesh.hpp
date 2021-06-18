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
/**
 * \file   GMesh.hpp
 * \brief  Class for Granular Mesh
 * \author Viviane Rochon Montplaisir
 * \date   November 2017
 */

#ifndef __NOMAD_4_0_GMESH__
#define __NOMAD_4_0_GMESH__

#include "../../Algos/MeshBase.hpp"

#include "../../nomad_nsbegin.hpp"

// Support Mesh Index (issue (feature) #381)

/// Class for the granular mesh of Mads.
/**
 This class derives from MeshBase.
 The class manages the mesh size (delta) and the frame size (Delta) for the discretization of the variable space used by Mads. Each variable has its own mesh and frame sizes which allows to increase or decrease the anisotropy of the mesh, that is changing the "cell" aspect ratios. The frame size (and mesh size) for each variable can be enlarged or decreased (see GMesh::refineDeltaFrameSize and GMesh::enlargeDeltaFrameSize). A given point can be projected on the the mesh using GMesh::scaleAndProjectOnMesh. \n

 The frame size for each variable is parameterized with two or three attributes: GMesh::_frameSizeExp (exponent), GMesh::_frameSizeMant (mantissa) and GMesh::_granularity (Delta = gran * a * 10^b with b an integer). The last attribute is for variable having a specified minimal granularity (for example, integers have a minimal granularity of 1). The mesh size is delta = 10^(b-|b-b_0|) for variables without granularity and delta = granulariy * max(1,delta) for variables with granularity.

 */
class GMesh: public MeshBase
{

    ArrayOfDouble        _initFrameSizeExp;  ///< The initial frame size exponent.
    ArrayOfDouble        _frameSizeMant;  ///< The current frame size mantissa.
    ArrayOfDouble        _frameSizeExp;  ///< The current frame size exponent.
    const ArrayOfDouble  _granularity;  ///< The fixed granularity of the mesh
    bool                 _enforceSanityChecks;   ///< Should we enforce sanity checks?

public:

    /// Constructor
    /**
     \param parameters  The problem parameters attributes control the mesh mechanics -- \b IN.
     */
    explicit GMesh(std::shared_ptr<PbParameters> parameters)
      : MeshBase(parameters),
        _initFrameSizeExp(ArrayOfDouble()),
        _frameSizeMant(ArrayOfDouble()),
        _frameSizeExp(ArrayOfDouble()),
        _granularity(parameters->getAttributeValue<ArrayOfDouble>("GRANULARITY")),
        _enforceSanityChecks(true)
    {
        init();
    }


    /*-----------*/
    /* Get / Set */
    /*-----------*/
    const ArrayOfDouble& getInitFrameSizeExp() const { return _initFrameSizeExp; }
    const ArrayOfDouble& getFrameSizeMant() const { return _frameSizeMant; }
    const ArrayOfDouble& getFrameSizeExp() const { return _frameSizeExp; }
    const ArrayOfDouble& getGranularity() const { return _granularity; }

    /**
     Performing enforced sanity checks consists in checking GMesh::checkFrameSizeIntegrity and GMesh::checkDeltasGranularity. It is requested when delta/granularity==1.
     */
    void setEnforceSanityChecks(const bool doubleCheck) { _enforceSanityChecks = doubleCheck; }

    /*----------------------*/
    /* Other Class Methods */
    /*----------------------*/


    void checkMeshForStopping( std::shared_ptr<AllStopReasons> algoStopReason ) const override;

    void updatedeltaMeshSize() override;

    /**
     \copydoc MeshBase::enlargeDeltaFrameSize
     \note This implementation relies on GMesh::_frameSizeExp, GMesh::_frameSizeMant and GMesh::_granularity.
     */
    bool enlargeDeltaFrameSize(const Direction& direction,
                               const Double& anisotropyFactor = 0.1,
                               bool anisotropicMesh = true) override;

    /**
     \copydoc MeshBase::refineDeltaFrameSize
     \note This implementation relies on GMesh::_frameSizeExp and GMesh::_frameSizeMant. Frame size is updated and mesh size can be impacted too if _frameSizeExp is decreased. Upon failure, the frame sizes for all dimensions are refined in the same way (no need for direction and anisotropy factor).
     */
    void refineDeltaFrameSize() override;

private:

    /// Helper for refineDeltaFrameSize().
    /**
     \param frameSizeMant    The frame size mantissa to update -- \b IN/OUT.
     \param frameSizeExp     The frame size exponant to update -- \b IN/OUT.
     \param granularity      The granularity of the mesh -- \b IN.
     */
    void refineDeltaFrameSize(Double &frameSizeMant,
                              Double &frameSizeExp,
                              const Double &granularity) const;
public:
    Double getRho(const size_t i) const override;
    ArrayOfDouble getRho() const override { return MeshBase::getRho(); }

    ArrayOfDouble getdeltaMeshSize() const override;
    Double getdeltaMeshSize(const size_t i) const override;
private:

    /// Helper function
    Double getdeltaMeshSize(const Double& frameSizeExp,
                            const Double& initFrameSizeExp,
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
    Double getDeltaFrameSize(const Double& granularity, const Double& frameSizeMant, const Double& frameSizeExp) const;

public:

    void setDeltas(const ArrayOfDouble &deltaMeshSize,
                   const ArrayOfDouble &deltaFrameSize) override;

    void setDeltas(const size_t i,
                   const Double &deltaMeshSize,
                   const Double &deltaFrameSize) override;

private:
    /// helper function to verify values are correct.
    void checkDeltasGranularity(const size_t i,
                                const Double &deltaMeshSize,
                                const Double &deltaFrameSize) const;

    /// helper function to verify values are correct.
    void checkFrameSizeIntegrity(const Double &frameSizeExp,
                                const Double &frameSizeMant) const;

    /// helper function to verify values are correct.
    void checkSetDeltas(const size_t i,
                        const Double &deltaMeshSize,
                        const Double &deltaFrameSize) const;

public:

    // Add some documentation in addition to the base class
    /**
     \copydoc MeshBase::scaleAndProjectOnMesh
     \note This implementation relies on GMesh::_frameSizeExp and GMesh::_frameSizeMant.
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

    /// Round frame size exponent to an int.
    /**
     \param exp     The frame size exponent -- \b IN.
     \return        The frame size exponent as an int
     */
    int roundFrameSizeExp(const Double &exp);

    /// Round frame size mantissa to an int.
    /**
     \param mant     The frame size mantissa -- \b IN.
     \return         The frame size mantissa as an int
     */
    int roundFrameSizeMant(const Double &mant);

    /// Helper for enlargeDeltaFrameSize and getDeltaFrameSizeCoarser.
    void getLargerMantExp(Double &frameSizeMant, Double &frameSizeExp) const;
};


#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_0_GMESH__

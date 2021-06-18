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
 * \file   MeshBase.hpp
 * \brief  Base class for mesh
 * \author Viviane Rochon Montplaisir
 * \date   November 2017
 */

#ifndef __NOMAD_4_0_MESHBASE__
#define __NOMAD_4_0_MESHBASE__

#include <memory>   // for shared_ptr

#include "../Math/ArrayOfDouble.hpp"
#include "../Math/Direction.hpp"
#include "../Math/Point.hpp"
#include "../Param/PbParameters.hpp"
#include "../Util/AllStopReasons.hpp"

#include "../nomad_nsbegin.hpp"

/// The generic class for meshes (discretization of design variables space).
/**
 \note To be implemented by the derived mesh used by an algorithm.

 \note This class encompasses the mesh discretization and the frame on it.

 * In the algorithm, we have delta for mesh size parameter and
  Delta for frame size parameter. To avoid confusion in function calls,
  using deltaMeshSize and DeltaFrameSize in function names.
 * Functions for delta: MeshBase::getdeltaMeshSize, MeshBase::updatedeltaMeshSize
 * Functions for Delta: MeshBase::getDeltaFrameSize, MeshBase::enlargeDeltaFrameSize, MeshBase::refineDeltaFrameSize
 */
class MeshBase
{
protected:
    const size_t _n;     ///< Dimension
    std::shared_ptr<PbParameters>  _pbParams; ///< The parameters for a mesh are in problem parameters.
    const ArrayOfDouble  _initialMeshSize;  ///< The initial mesh size.
    const ArrayOfDouble  _minMeshSize;   ///< The minimum mesh size (stopping criterion).
    const ArrayOfDouble  _initialFrameSize; ///< The initial frame size.
    const ArrayOfDouble  _minFrameSize; ///< The minimum frame size (stopping criterion).
    const ArrayOfDouble _lowerBound;
    const ArrayOfDouble _upperBound;

public:

    /// Constructor From [mesh] parameters
    /**
     \param pbParams    The mesh parameters are taken from problem parameters -- \b IN.
     */
    explicit MeshBase(std::shared_ptr<PbParameters> pbParams);

    /**
     Virtual destructor needed to avoid compilation warning:
     destructor called on non-final 'GMesh' that has virtual functions
     but non-virtual destructor [-Wdelete-non-virtual-dtor]
     */
    virtual ~MeshBase() {}

    /*-----------------------------------------------------*/
    /* Get / Set                                           */
    /* Set removed for members because they are all const. */
    /*-----------------------------------------------------*/
    size_t getSize() const { return _n; }

    const ArrayOfDouble& getInitialMeshSize() const { return _initialMeshSize; }
    const ArrayOfDouble& getMinMeshSize() const { return _minMeshSize; }
    const ArrayOfDouble& getInitialFrameSize() const { return _initialFrameSize; }
    const ArrayOfDouble& getMinFrameSize() const { return _minFrameSize; }

    /*------------------*/
    /*   Mesh methods   */
    /*------------------*/

    /// Update mesh size.
    /**
     Update mesh size (small delta) based on frame size (big Delta)
     */
    virtual void updatedeltaMeshSize() = 0;

    /// Enlarge frame size.
    /**
     * Update frame size (big Delta) after a success.

     * The successful direction is used to ensure integrity of the mesh (no variable collapse).

     \param direction          The direction of success of the iteration       -- \b IN.
     \param anisotropyFactor   Control the development of anisotropy of the mesh (if anisotropicMesh is \c true) -- \b IN.
     \param anisotropicMesh    Flag to enable or not anisotropic mesh  -- \b IN.
     \return                   \c true if the mesh has changed, \c false otherwise.
     */
    virtual bool enlargeDeltaFrameSize(const Direction& direction,
                                       const Double& anisotropyFactor = 0.1,
                                       bool anisotropicMesh = true) = 0;

    /// Refine frame size.
    /**
     Update frame size after a failure.
     */
    virtual void refineDeltaFrameSize() = 0;


    /// Check the mesh stopping condition
    virtual void checkMeshForStopping( std::shared_ptr<AllStopReasons> algoStopReason ) const =0;


    /// Access to the ratio of frame size / mesh size parameter rho^k.
    /**
     \param i       Index of the dimension of interest -- \b IN.
     \return        The ratio frame/mesh size rho^k for index i.
     */
    virtual Double getRho(const size_t i) const = 0;

    /// Access to the ratio of frame size / mesh size parameter rho^k.
    /**
     \return The ratio frame/mesh size rho^k.
     */
    virtual ArrayOfDouble getRho() const;

    /// Access to the mesh size parameter delta^k.
    /**
     \param i       Index of the dimension of interest -- \b IN.
     \return        The mesh size parameter delta^k.
     */
    virtual Double getdeltaMeshSize(const size_t i) const = 0;

    /// Access to the mesh size parameter delta^k.
    /**
     \return The mesh size parameter delta^k.
     */
    virtual ArrayOfDouble getdeltaMeshSize() const;

    // Access to the frame size parameter Delta^k.
    /**
     \param i       Index of the dimension of interest -- \b IN.
     \return        The frame size parameter Delta^k.
     */
    virtual Double getDeltaFrameSize(const size_t i) const = 0;

    // Access to the frame size parameter Delta^k.
    /**
     \return The frame size parameter Delta^k.
     */
    virtual ArrayOfDouble getDeltaFrameSize() const;

    // Access to the frame size one shift coarser than the actual frame size.
    /**
     \param i       Index of the dimension of interest -- \b IN.
     \return        The frame size parameter Delta^k.
     */
    virtual Double getDeltaFrameSizeCoarser(const size_t i) const = 0;

    // Access to the frame size one shift coarser than the actual frame size.
    /**
     \return The frame size parameter Delta^k.
     */
    virtual ArrayOfDouble getDeltaFrameSizeCoarser() const;

    /**
     Setting deltaMeshSize and DeltaFrameSize should be done together.
     This is easier and does not seem to be a constraint for
    the general case.
     */
    virtual void setDeltas(const size_t i,
                           const Double &deltaMeshSize,
                           const Double &deltaFrameSize) = 0;

    /**
     Setting deltaMeshSize and DeltaFrameSize should be done together.
     This is easier and does not seem to be a constraint for
     the general case.
     */
    virtual void setDeltas(const ArrayOfDouble &deltaMeshSize,
                           const ArrayOfDouble &deltaFrameSize);

    /// Scale and project the ith component of a vector on the mesh
    /**
     \param i   The vector component number         -- \b IN.
     \param l   The vector component value          -- \b IN.
     \return    The ith component of a vector after mesh scaling and projection.
     */
    virtual Double scaleAndProjectOnMesh(size_t i, const Double &l) const = 0;

    /// Scale and project a direction on the mesh
    /**
     \param dir   The direction to scale and project          -- \b IN.
     \return      The ith component of a vector after mesh scaling and projection.
     */
    virtual ArrayOfDouble scaleAndProjectOnMesh(const Direction &dir) const;

    /// Project the point on the mesh centered on frameCenter. No scaling.
    virtual Point projectOnMesh(const Point& point, const Point& frameCenter) const;

    /// Verify if the point is on the mesh.
    /**
    \return     \c true if the point is on the mesh, false otherwise.
    */
    bool verifyPointIsOnMesh(const Point& point, const Point& frameCenter) const;

private:
    /// Helper for constructor.
    void init();

protected:
    /**
     Verify dimension is n. Throw an exception if it is not.
     */
    void verifyDimension(const std::string& arrayName, size_t dim);
};


///   Display useful values so that a new mesh could be constructed using these values.
std::ostream& operator<<(std::ostream& os, const MeshBase& mesh);

/// Get the mesh values from stream
std::istream& operator>>(std::istream& is, MeshBase& mesh);


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_MESHBASE__

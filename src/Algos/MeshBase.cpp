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
#include "../Algos/MeshBase.hpp"

NOMAD::MeshBase::MeshBase(const std::shared_ptr<NOMAD::PbParameters> pbParams)
      : _n(pbParams->getAttributeValue<size_t>("DIMENSION")),
        _pbParams(pbParams),
        _initialMeshSize(pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_MESH_SIZE")),
        _minMeshSize(pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("MIN_MESH_SIZE")),
        _initialFrameSize(pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE")),
        _minFrameSize(pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("MIN_FRAME_SIZE")),
        _lowerBound(pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND")),
        _upperBound(pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND"))
{
    init();
}


void NOMAD::MeshBase::init()
{
    if ( _pbParams->toBeChecked() )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameters::checkAndComply() needs to be called before constructing a mesh.");
    }
}

NOMAD::ArrayOfDouble NOMAD::MeshBase::getRho(void) const
{
    NOMAD::ArrayOfDouble rho(_n);
    for (size_t i = 0; i < _n; i++)
    {
        rho[i] = getRho(i);
    }
    return rho;
}


NOMAD::ArrayOfDouble NOMAD::MeshBase::getdeltaMeshSize(void) const
{
    NOMAD::ArrayOfDouble delta(_n);
    for (size_t i = 0; i < _n; i++)
    {
        delta[i] = getdeltaMeshSize(i);
    }
    return delta;
}


NOMAD::ArrayOfDouble NOMAD::MeshBase::getDeltaFrameSize() const
{
    NOMAD::ArrayOfDouble Delta(_n);
    for (size_t i = 0; i < _n; i++)
    {
        Delta[i] = getDeltaFrameSize(i);
    }
    return Delta;
}


NOMAD::ArrayOfDouble NOMAD::MeshBase::getDeltaFrameSizeCoarser() const
{
    NOMAD::ArrayOfDouble Delta(_n);
    for (size_t i = 0; i < _n; i++)
    {
        Delta[i] = getDeltaFrameSizeCoarser(i);
    }
    return Delta;
}


void NOMAD::MeshBase::setDeltas(const NOMAD::ArrayOfDouble &deltaMeshSize,
                                const NOMAD::ArrayOfDouble &deltaFrameSize)
{
    for (size_t i = 0; i < _n; i++)
    {
        setDeltas(i, deltaMeshSize[i], deltaFrameSize[i]);
    }
}


/*-----------------------------------------------------------*/
/*              scale and project on the mesh                */
/*-----------------------------------------------------------*/
NOMAD::Double NOMAD::MeshBase::scaleAndProjectOnMesh(size_t i, const NOMAD::Double &l) const
{
    // Not defined for MeshBase.
    throw NOMAD::Exception(__FILE__, __LINE__, "scaleAndProjectOnMesh() not defined for MeshBase.");
}


NOMAD::ArrayOfDouble NOMAD::MeshBase::scaleAndProjectOnMesh(
    const NOMAD::Direction &dir) const
{
    // Not defined for MeshBase.
    throw NOMAD::Exception(__FILE__, __LINE__, "scaleAndProjectOnMesh() not defined for MeshBase.");
}


NOMAD::Point NOMAD::MeshBase::projectOnMesh(const NOMAD::Point& point,
                                            const NOMAD::Point& frameCenter) const
{
    // Not defined for MeshBase.
    throw NOMAD::Exception(__FILE__, __LINE__, "projectOnMesh() not defined for MeshBase.");
}


bool NOMAD::MeshBase::verifyPointIsOnMesh(const NOMAD::Point& point, const NOMAD::Point& center) const
{
    bool isOnMesh = true;

    for (size_t i = 0; i < point.size(); i++)
    {
        NOMAD::Double pointRebaseI = point[i];
        NOMAD::Double centerI = center[i];
        NOMAD::Double deltaI = getdeltaMeshSize(i);

        if (   (_lowerBound[i].isDefined() && _lowerBound[i] == pointRebaseI)
            || (_upperBound[i].isDefined() && _upperBound[i] == pointRebaseI))
        {
            isOnMesh = true;
        }
        else
        {
            if (!centerI.isMultipleOf(deltaI))
            {
                // Rebase point on the mesh centered on center Point
                pointRebaseI -= centerI;
            }
            if (!pointRebaseI.isMultipleOf(deltaI))
            {
                isOnMesh = false;
                break;
            }
        }
    }

    return isOnMesh;
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::MeshBase& mesh)
{
    os << "DELTA_MESH_SIZE " << mesh.getdeltaMeshSize() << std::endl;
    os << "DELTA_FRAME_SIZE " << mesh.getDeltaFrameSize() << std::endl;

    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::MeshBase& mesh)
{
    size_t n = mesh.getSize();

    // Read line by line
    std::string name;

    NOMAD::ArrayOfDouble deltaMeshSize(n), deltaFrameSize(n);
    while (is >> name && is.good() && !is.eof())
    {
        if ("DELTA_MESH_SIZE" == name)
        {
            is >> deltaMeshSize;
        }
        else if ("DELTA_FRAME_SIZE" == name)
        {
            is >> deltaFrameSize;
        }
        else
        {
            // Put back name to istream. Maybe there is a simpler way.
            for (unsigned i = 0; i < name.size(); i++)
            {
                is.unget();
            }
            break;
        }
    }

    mesh.setDeltas(deltaMeshSize, deltaFrameSize);

    return is;

}


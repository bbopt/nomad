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
 \file   Direction.cpp
 \brief  Custom class for directions (implementation)
 \author Sebastien Le Digabel and Viviane Rochon Montplaisir
 \date   March 2017
 \see    Direction.hpp
 */
#include "../Math/Direction.hpp"
#include "../Math/RNG.hpp"

// Assignment operator
NOMAD::Direction& NOMAD::Direction::operator=(const NOMAD::Direction& dir)
{
    NOMAD::ArrayOfDouble::operator=(dir);
    return *this;
}

/*-----------------------------------------*/
/* Operators for addition and substraction */
/*-----------------------------------------*/
const NOMAD::Direction& NOMAD::Direction::operator+=(const Direction& dir1)
{
    for (size_t i = 0; i < size(); i++)
    {
        _array[i] += dir1[i];
    }
    return *this;
}


const NOMAD::Direction& NOMAD::Direction::operator-=(const Direction& dir1)
{
    for (size_t i = 0; i < size(); i++)
    {
        _array[i] -= dir1[i];
    }
    return *this;
}


/*------------------------------------------------------*/
/* Squared Norm of this direction, viewed as a vector.  */
/*------------------------------------------------------*/
const NOMAD::Double NOMAD::Direction::squaredL2Norm() const
{
    NOMAD::Double sqL2 = 0;

    for (size_t i = 0; i < size(); i++)
    {
        sqL2 += _array[i] * _array[i];
    }

    return sqL2;
}


/*------------------------------------------------*/
/* Norm of this direction X, viewed as a vector.  */
/* Norm type may be L1, L2, or LINF.              */
/*------------------------------------------------*/
const NOMAD::Double NOMAD::Direction::norm(NOMAD::NormType normType) const
{
    NOMAD::Double retNorm = 0;

    switch (normType)
    {
        case NOMAD::NormType::L1:
            for (size_t i = 0; i < size(); i++)
            {
                retNorm += _array[i].abs();
            }
            break;
        case NOMAD::NormType::L2:
        default:
            retNorm = this->squaredL2Norm();
            retNorm = sqrt(retNorm.todouble());
            break;
        case NOMAD::NormType::LINF:
            for (size_t i = 0; i < size(); i++)
            {
                retNorm = NOMAD::max(retNorm, _array[i].abs());
            }
            break;
    }

    return retNorm;
}


/// Infinite norm
const NOMAD::Double NOMAD::Direction::infiniteNorm() const
{
    return norm(NormType::LINF);
}


const NOMAD::Double NOMAD::Direction::dotProduct(const NOMAD::Direction& dir1,
                                                 const NOMAD::Direction& dir2)
{
    NOMAD::Double dot = 0.0;

    size_t size = dir1.size();
    if (size != dir2.size())
    {
        std::string err = "Dot product: vectors are not of the same size: \n";
        err += dir1.display() + "\n";
        err += dir2.display();
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (size_t i = 0; i < size; i++)
    {
        dot += dir1[i] * dir2[i];
    }

    return dot;
}


const NOMAD::Double NOMAD::Direction::cos(const NOMAD::Direction& dir1,
                                          const NOMAD::Direction& dir2)
{
    NOMAD::Double cos = 0.0;

    double norm1 = dir1.norm().todouble();
    double norm2 = dir2.norm().todouble();
    if (0.0 == norm1 || 0.0 == norm2)
    {
        std::string err = "Cosine: a vector is of size 0";
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    cos = dotProduct(dir1, dir2) / (norm1 * norm2);

    return cos;
}

/*-----------------------------------------------------------*/
/*      computation of the angle with another direction      */
/*-----------------------------------------------------------*/
const NOMAD::Double NOMAD::Direction::angle(const NOMAD::Direction& dir1,
                                            const NOMAD::Direction& dir2)
{
    if (dir1.size() != dir2.size())
    {
        return NOMAD::Double();
    }

    NOMAD::Double innerProduct = 0.0, norm1 = 0.0, norm2 = 0.0;

    for (size_t i = 0; i < dir1.size(); i++)
    {
        norm1        += dir1[i] * dir1[i];
        norm2        += dir2[i] * dir2[i];
        innerProduct += dir1[i] * dir2[i];
    }

    if (norm1 == 0.0 || norm2 == 0.0)
    {
        return NOMAD::Double();
    }

    return std::acos((innerProduct / (norm1.sqrt() * norm2.sqrt())).todouble());
}


/*--------------------------------------------------*/
/*  Compute a random direction on a unit N-Sphere   */
/*  See http://en.wikipedia.org/wiki/N-sphere       */
/*--------------------------------------------------*/
void NOMAD::Direction::computeDirOnUnitSphere(NOMAD::Direction &randomDir)
{
    size_t i;
    NOMAD::Double norm, d;
    size_t n = randomDir.size();

    for (i = 0; i < n; ++i)
    {
        randomDir[i] = NOMAD::RNG::normalRand(0,1);
    }

    norm = randomDir.norm();

    if (0 == norm)
    {
        std::string err("Cannot compute a random direction on unit sphere");
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (i = 0; i < n; ++i)
    {
        randomDir[i] /= norm;
    }

}


/*----------------------------------------------------------------*/
/*  - Householder transformation to generate _nc directions from  */
/*    a given direction.                                          */
/*  - Compute also H[i+nc] = -H[i] (completion to 2n directions). */
/*  - Private method.                                              */
/*----------------------------------------------------------------*/
void NOMAD::Direction::householder(const NOMAD::Direction &dir,
                                   bool completeTo2n,
                                   NOMAD::Direction ** H)
{
    size_t n = dir.size();

    const NOMAD::Double norm2 = dir.squaredL2Norm();
    NOMAD::Double v, h2i;

    for (size_t i = 0 ; i < n ; ++i)
    {
        h2i = 2 * dir[i];
        for (size_t j = 0 ; j < n ; ++j)
        {
            // H[i]:
            (*H[i])[j] = v = (i == j) ? norm2 - h2i * dir[j] : - h2i * dir[j];

            // -H[i]:
            if ( completeTo2n )
            {
                (*H[i+n])[j] = -v;
            }
        }
    }
}


NOMAD::Direction NOMAD::operator-(const Direction &dir)
{
    size_t size = dir.size();
    NOMAD::Direction negDir(size);
    for (size_t i = 0; i < size; i++)
    {
        negDir[i] = -dir[i];
    }
    return negDir;
}


std::ostream& NOMAD::operator<< (std::ostream& out, const NOMAD::Direction& dir)
{
    out << dir.display();
    return out;
}

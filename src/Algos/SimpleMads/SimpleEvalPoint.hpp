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
/**
 \file   SimpleEvalPoint.hpp
 \brief  Evaluation point for simple mads barrier
 \author Christophe Tribes
 \date   April 2024
 \see    SimpleEvalPoint.cpp
 */

#ifndef __NOMAD_4_5_SIMPLEEVALPOINT__
#define __NOMAD_4_5_SIMPLEEVALPOINT__

#include "../../Math/Double.hpp"
#include "../../Math/Point.hpp"



#include "../../nomad_platform.hpp"
#include "../../nomad_nsbegin.hpp"


/// Class for the representation of a simple form for evaluation point.
/**
 Simple evaluation point for the point coordinates \c x, and the blackbox
 outputs at these coordinates summarized into f(x) and h(x).
*/
class DLL_ALGO_API SimpleEvalPoint : public Point
{
private:

    Double _f, _h;
    
public:

    /*---------------*/
    /* Class Methods */
    /*---------------*/

    /// Constructor #1.
    explicit SimpleEvalPoint()
    : NOMAD::Point() {};

    /// Constructor #2.
    /**
     \param n Number of variables -- \b IN.
     */
    explicit SimpleEvalPoint(size_t n)
    : NOMAD::Point(n) {};

    /// Constructor #3.
    /**
      \param x Coordinates of the eval point -- \b IN.
      */
    explicit SimpleEvalPoint(const Point& x)
    : NOMAD::Point(x) {}; // Let F and H undefined
    

    /// Copy constructor.
    /**
     \param evalPoint The copied object -- \b IN.
     */
    SimpleEvalPoint(const SimpleEvalPoint& evalPoint)
    : NOMAD::Point(evalPoint),
      _f(evalPoint._f),
      _h(evalPoint._h) {};
    

public:

    /// Affectation operator to base class.
    /**
     \param evalPoint The right-hand side object -- \b IN.
     \return           \c *this as the result of the affectation.
     */
    SimpleEvalPoint& operator= (const SimpleEvalPoint& evalPoint);


    /// Destructor.
    virtual ~SimpleEvalPoint() {};
    
    const Double & getF() const { return _f;}
    const Double & getH() const { return _h;}
    
    void setF(const Double & f) {_f = f;}
    void setH(const Double & h) {_h = h;}
    
    bool isDefined() const ;

};

//
///// Display useful values so that a new SimpleEvalPoint could be constructed using these values.
//DLL_ALGO_API std::ostream& operator<<(std::ostream& os,
//                         const SimpleEvalPoint &evalPoint);
//
///// Get these values from stream
//DLL_ALGO_API std::istream& operator>>(std::istream& is, SimpleEvalPoint &evalPoint);
//

#include "../../nomad_nsend.hpp"

#endif // __NOMAD_4_5_SIMPLEEVALPOINT__

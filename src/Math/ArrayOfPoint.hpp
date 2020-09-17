/**
 \file   ArrayOfPoint.hpp
 \brief  Vector of Points
 \author Viviane Rochon Montplaisir
 \date   August 2019
 \see    ArrayOfPoint.cpp
 */

#ifndef __NOMAD400_ARRAY_OF_POINT__
#define __NOMAD400_ARRAY_OF_POINT__

#include <vector>
#include "../Math/Point.hpp"

#include "../nomad_nsbegin.hpp"

/// Type definition for the representation of a vector of points.
typedef std::vector<Point> ArrayOfPoint;

std::ostream& operator<<(std::ostream& out, const ArrayOfPoint& aop);
// Not implemented:
//std::istream& operator>>(std::istream& in, ArrayOfPoint& aop);



#include "../nomad_nsend.hpp"
#endif // __NOMAD400_ARRAY_OF_POINT__

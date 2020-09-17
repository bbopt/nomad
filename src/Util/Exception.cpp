/**
 \file   Exception.cpp
 \brief  custom class for exceptions (implementation)
 \author Sebastien Le Digabel
 \date   2010-03-29
 \see    Exception.hpp
 */
#include "../Util/Exception.hpp"

/*----------------------------------------------------------------*/
/*                     NOMAD::Exception::what()                   */
/*----------------------------------------------------------------*/
const char * NOMAD::Exception::what() const throw()
{
    std::ostringstream oss;
    if ( ! _file.empty() || _line > 0 )
    {
        oss << "NOMAD::Exception thrown (" << _file << ", " << _line << ")";
    }
    if ( !_what.empty() )
    {
        if (!_typeMsg.empty())
        {
            oss << " " << _typeMsg;
        }
        oss << " " << _what;
    }
    _what = oss.str();

    return _what.c_str();
}

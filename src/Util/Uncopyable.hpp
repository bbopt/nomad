/**
  \file   Uncopyable.hpp
  \brief  Base class for uncopyable classes (headers)
  \author Sebastien Le Digabel
  \date   2010-04-02
*/
#ifndef __NOMAD400_UNCOPYABLE__
#define __NOMAD400_UNCOPYABLE__

#include "../nomad_nsbegin.hpp"

  /// Base class for uncopyable classes
  /**
    See Scott Meyer's Effective C++, 3rd ed., item #6.
  */
  class Uncopyable
{

  protected:

    /// Constructor.
    explicit Uncopyable  ( void ) {}

    /// Destructor.
    virtual ~Uncopyable ( void ) {}

  private:

    /// Undefined copy constructor.
    Uncopyable ( const Uncopyable & );

    /// Undefined affectation operator.
    Uncopyable & operator = ( const Uncopyable & );
  };
#include "../nomad_nsend.hpp"

#endif // __NOMAD400_UNCOPYABLE__


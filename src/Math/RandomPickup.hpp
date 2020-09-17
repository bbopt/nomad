/**
 \file   RandomPickup.hpp
 \brief  Class for randomly picking up integers
 \author Sebastien Le Digabel
 \date   2010-04-07
 \see    RandomPickup.cpp
 */

#ifndef __NOMAD400_RANDOM_PICKUP__
#define __NOMAD400_RANDOM_PICKUP__

#include <cstdlib>
#include "../Util/Uncopyable.hpp"

using namespace std;

#include "../nomad_nsbegin.hpp"


/// Class for randomly picking up integers.
/**
   - The integers are chosen in [0;n-1] and are distinct.
   - Example displaying 5 different integers in [0;4]:
   \code
   NOMAD::RandomPickup rp(5);
   for (size_t i = 0 ; i < 5 ; ++i)
   {
       std::cout << rp.pickup() << std::endl;
   }
   \endcode
*/
class RandomPickup : private Uncopyable {

private:
    size_t _n0;     ///< Initial value of \c n.
    size_t _n;      ///< Current value of \c n.
    size_t *_elems; ///< Elements that have not been chosen yet.

public:
    /// Constructor.
    /**
       \param n -- The unsigned integer \c n defining the range
                   of values that can be picked up -- \b IN.
    */
    explicit RandomPickup(const size_t n);

    /// Destructor.
    virtual ~RandomPickup() { delete [] _elems; }

    /// Get number of remaining values
    size_t getN() const { return _n; }

    /// Reset.
    void reset();

    /// Randomly pick up an element in [0;n-1].
    /**
       \return The element.
    */
    size_t pickup();

};

#include "../nomad_nsend.hpp"


#endif // __NOMAD400_RANDOM_PICKUP__

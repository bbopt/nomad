/**
 \file   BBInput.hpp
 \brief  Input of a Blackbox evaluation
 \author Viviane Rochon Montplaisir
 \date   March 2018
 \see    BBInput.cpp
 */

// Manage input to blackbox:
// Variable types - Continuous, integer, binary
// Scaling (future work)
//
// Note: As of September 2019, this class is not used by NOMAD.

#ifndef __NOMAD400_BB_INPUT__
#define __NOMAD400_BB_INPUT__


#include "../Math/Point.hpp"
#include "../Type/BBInputType.hpp"

#include "../nomad_nsbegin.hpp"


/// Class for the representation of the input to a blackbox evaluation.
class BBInput
{
    //std::string             _rawBBI;        // Actual input string (currently not implemented).

public:

    /*---------------*/
    /* Class Methods */
    /*---------------*/
    /// Constructor
    /**
      Currently does nothing.
     \param bbInputTypeList     The list of blackbox input types  -- \b IN.
     \param point               The point -- \b IN.
     */
    explicit BBInput(const BBInputTypeList& bbInputTypeList,
                     const Point& point);

    /*---------*/
    /* Get/Set */
    /*---------*/
    // To implement as needed.
    // Maybe something like this:
    //void setBBInput(const std::string bbInputString);
    //void setBBInput(const Point point);
    //void getBBInput(Point& point) const;

    // Currently does nothing.
    void display(std::ostream &out) const;

};


inline std::ostream& operator<<(std::ostream& out, const BBInput &bbinput)
{
    bbinput.display(out);
    return out;
}


#include "../nomad_nsend.hpp"
#endif  // __NOMAD400_BB_INPUT__

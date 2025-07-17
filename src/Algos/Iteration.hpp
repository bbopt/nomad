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
#ifndef __NOMAD_4_5_ITERATION__
#define __NOMAD_4_5_ITERATION__

#include "../Algos/Step.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for iteration of an Algorithm.
/**
 This is an abstract class, each algorithm must implement its own iteration.
 */
class DLL_ALGO_API Iteration: public Step
{
    
protected:
    bool _userCallbackEnabled;
    
protected:

    size_t _k; ///< Iteration number

    void init(); ///< Utility for constructor

public:
    /// Constructor
    /**
     \param parentStep         The parent of this step -- \b IN.
     \param k                  The initial iteration counter -- \b IN.
     */
    explicit Iteration(const Step *parentStep,
                       const size_t k)
      : Step( parentStep ),
        _k(k)
    {
        init();
    }

    /// Destructor
    /**
     When iteration is done, Flush prints output queue.
     */
    virtual ~Iteration();

    // Get/Set
    /// Get name
    virtual std::string getName() const override;

    /// Get iteration number
    /**
     Iteration number is incremented when calling the default Iteration::start().
     */
    size_t getK() const { return _k; }

    /**
     \return \c nullptr for algorithms that do not use a mesh. Otherwise, this function must be reimplemented in algorithm specific iteration (for example, MadsIteration, NMIteration).
     */
    virtual const MeshBasePtr getMesh() const { return nullptr; }

protected:

    /**
     This must be implemented when an algorithm has its own iteration.
     */
    virtual void startImp()    override = 0;

    /**
     This must be implemented when an algorithm has its own iteration.
     */
    virtual bool runImp()      override = 0;

    /**
     The default implement for end function displays the stop reason, calls the customized end function if provided by the user and increment the counter. \n
     If an end implementation function specific to an algorithm is required, this function MUST be called for default tasks.
     */
     virtual void endImp()      override;

private:
    
    /// Implementation to increment the iteration counter
    virtual void incrementCounters() override { _k++ ;}

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_5_ITERATION__

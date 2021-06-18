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
#ifndef __NOMAD_4_0_OUTPUTINFO__
#define __NOMAD_4_0_OUTPUTINFO__

#include <vector>
#include "../Output/StatsInfo.hpp"

#include "../nomad_nsbegin.hpp"

/// Level of output AKA Display degree
enum class OutputLevel {
    LEVEL_NOTHING = 0,  ///< Print nothing
    LEVEL_ERROR,
    LEVEL_VERY_HIGH,    ///< Errors and results
    LEVEL_WARNING,
    LEVEL_HIGH,         ///< Print high level like Warnings
    LEVEL_STATS,        ///< For stats - parameters DISPLAY_STATS and STATS_FILE
    LEVEL_NORMAL,       ///< Print medium level like global information and useful information
    LEVEL_INFO,         ///< Lots of information
    LEVEL_LOW,          ///< Print low-level information
    LEVEL_DEBUG,        ///< Print all, Debug level
    LEVEL_DEBUGDEBUG,   ///< Print even more than you asked for
    NB_LEVEL
};


/**
 All information that may be useful for one output.
Used by OutputQueue.
 */
class OutputInfo
{
private:
    std::string          _originator;
    ArrayOfString        _msg;
    OutputLevel          _outputLevel;
    bool                 _blockStart;
    bool                 _blockEnd;
    StatsInfoUPtr        _statsInfo;

public:
    // Constructor
    explicit OutputInfo(const std::string& originator, const std::string& msg,
                        OutputLevel outputLevel = OutputLevel::LEVEL_INFO,
                        bool blockStart = false, bool blockEnd = false)
      : _originator(originator),
        _msg(),
        _outputLevel(outputLevel),
        _blockStart(blockStart),
        _blockEnd(blockEnd),
        _statsInfo()
        {
            _msg.add(msg);
        }
    // Do not specify destructor, so that default is used for destructor, move, copy.
    // virtual ~OutputInfo() {}

    //
    // Get/Set
    //
    const std::string& getOriginator() const { return _originator; }
    void setOriginator(const std::string& originator) { _originator = originator; }

    const ArrayOfString& getMsg() const { return _msg; }
    void addMsg(const std::string& msg) { _msg.add(msg); }

    const OutputLevel& getOutputLevel() const { return _outputLevel; }
    void setOutputLevel(const OutputLevel& outputLevel) { _outputLevel = outputLevel; }

    StatsInfo* getStatsInfo() const { return _statsInfo.get(); }

    void setStatsInfo(StatsInfoUPtr statsInfo)
    {
        _statsInfo = std::move(statsInfo);
    }

    //
    // Other class methods
    //
    bool isBlockStart() const { return _blockStart; }
    bool isBlockEnd()   const { return _blockEnd; }

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_OUTPUTINFO__

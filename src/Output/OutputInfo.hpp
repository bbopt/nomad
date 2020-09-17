#ifndef __NOMAD400_OUTPUTINFO__
#define __NOMAD400_OUTPUTINFO__

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

    void addMsgAndSol(const std::string& msg, const Point& point);

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

#endif // __NOMAD400_OUTPUTINFO__

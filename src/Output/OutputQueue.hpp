#ifndef __NOMAD400_OUTPUTQUEUE__
#define __NOMAD400_OUTPUTQUEUE__

#include <vector>
#ifdef _OPENMP
// Using OpenMP.
// NOTE We would do a wrapper like suggested in http://bisqwit.iki.fi/story/howto/openmp/.
// This way, the code would compile even if OpenMP is not available.
// Plus, it might be clearer.
#include <omp.h>
#endif // _OPENMP

#include "../Param/DisplayParameters.hpp"
#include "../Output/OutputInfo.hpp"
#include "../Output/StatsInfo.hpp"

#include "../nomad_nsbegin.hpp"

/// Queue of all information that was not output yet.
/**
 The output queue is a singleton. Some OutputInfo can be added to the queue. The output information is displayed when calling OutputQueue::Flush and queue is emptied. \n
 The information can be send to the standard display and to a stats file. \n
 The display is formatted with indendation of blocks of information. The display parameters (DisplayParameters) are attributes of the class provided by calling OutputQueue::initParameters. \n

 The display can be limited to a maximum block/step level. (OutputQueue::_maxStepLevel). \n

 \todo Replace calls to std::cout by something more general.

 */
class OutputQueue
{
private:
    /// Private constructor
    OutputQueue();

public:
    // Destructor
     virtual ~OutputQueue();

    /// Access to singleton
    static std::unique_ptr<OutputQueue>& getInstance();

    void initParameters(const std::shared_ptr<DisplayParameters>& params);

    void add(OutputInfo outputInfo);
    static void Add(OutputInfo outputInfo)
    {
        getInstance()->add(std::move(outputInfo));
    }

    void add(const std::string& s,
             OutputLevel outputLevel = OutputLevel::LEVEL_INFO);
    static void Add(const std::string& s,
                    OutputLevel outputLevel = OutputLevel::LEVEL_INFO)
    {
        getInstance()->add(s, outputLevel);
    }

    void add(const StatsInfo& statsInfo);
    static void Add(const StatsInfo & statsInfo)
    {
        getInstance()->add(statsInfo);
    }

    /// Print all in the queue and flush.
    /**
     OutputInfo block start and block end flags will print
     _blockStart after the msg, or _blockEnd before.
     \note Example output: \n
     Start step MADS { \n
     _______Start step SEARCH { \n
     _____________Things happening in SEARCH \n
     _____________More things happening in SEARCH \n
     _______} End step SEARCH \n
     _______Start step POLL { \n
     _____________Things happening in POLL \n
     _____________More things happening in POLL \n
     _______} End step POLL \n
     } End step MADS \n
     If there are more than one line to print, flags for block
     start and end are ignored.
     */
    static void Flush()
    {
        getInstance()->flush();
    }


    size_t getMaxStepLevel() const { return _maxStepLevel; }
    void setMaxStepLevel(const size_t maxStepLevel) { _maxStepLevel = maxStepLevel; }

    bool goodLevel(const OutputLevel& outputLevel) const;
    static bool GoodLevel(const OutputLevel& outputLevel)
    {
        return getInstance()->goodLevel(outputLevel);
    }

    // Macros for output
#define OUTPUT_STATS_START if (OutputQueue::GoodLevel(OutputLevel::LEVEL_STATS)) {
#define OUTPUT_INFO_START if (OutputQueue::GoodLevel(OutputLevel::LEVEL_INFO)) {
#define OUTPUT_DEBUG_START if (OutputQueue::GoodLevel(OutputLevel::LEVEL_DEBUG)) {
#define OUTPUT_STATS_END }
#define OUTPUT_INFO_END }
#define OUTPUT_DEBUG_END }
//#define OUTPUT_STATS_START
//#define OUTPUT_INFO_START
//#define OUTPUT_DEBUG_START
//#define OUTPUT_STATS_END
//#define OUTPUT_INFO_END
//#define OUTPUT_DEBUG_END

    int getDisplayDegree() const;
    void setDisplayDegree(const int displayDegree);

    void setStatsFileName(const std::string& statsFile) { _statsFile = statsFile; }
    void initStatsFile();
    const std::string& getStatsFileName() const { return _statsFile; }

    void setStatsFileFormat(const DisplayStatsTypeList& statsFileFormat)
    {
        _statsFileFormat = statsFileFormat;
    }
    const DisplayStatsTypeList& getStatsFileFormat() const { return _statsFileFormat; }

    // Used by OutputInfo.
    const ArrayOfDouble& getSolFormat() const
    {
        return _params->getAttributeValue<ArrayOfDouble>("SOL_FORMAT");
    }

private:
#ifdef _OPENMP
    // Acquire lock before Add or Flush.
    // NOTE It does not seem relevant for the lock to be static,
    // because OutputQueue is a singleton anyway. If staticity causes problems,
    // we could remove the static keyword.
    static omp_lock_t _s_queue_lock;
#endif // _OPENMP


    static std::unique_ptr<OutputQueue> _single; ///< The singleton


    /// Queue of all the OutputInfo we have to print.
    std::vector<OutputInfo> _queue;

    /// Display parameters
    std::shared_ptr<DisplayParameters> _params;


    std::string _statsFile;
    std::ofstream _statsStream;

    /**
     Format for stats in a file (parameter STATS_FILE).
     Might include some raw strings, do not convert to DisplayStatsType.
     */
    DisplayStatsTypeList _statsFileFormat;

    /**
     Keep track of the number of lines printed to output (DISPLAY_STATS).
     Used to print stats header regularly.
     */
    size_t _statsLineCount;

    /**
     Format width for OBJ and CONS_H. May be enlarged during the run.
     */
    size_t _objWidth;

    size_t _hWidth;

    size_t _maxStepLevel;  ///< Maximum step level we want to print out.
    OutputLevel _maxOutputLevel; ///< Output level (~display degree) we want to print out
    int _indentLevel;   ///< Internal indentation level

    const std::string _blockStart; ///< Symbol for a block start.
    const std::string _blockEnd; ///< Symbol for an end block.

    void startBlock();
    void endBlock();
    void flush();
    void flushBlock(const OutputInfo &outputInfo);
    void flushStatsToStatsFile(const StatsInfo *statsInfo);
    void flushStatsToStdout(const StatsInfo *statsInfo);
    void indent(int level);
};

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_OUTPUTQUEUE__

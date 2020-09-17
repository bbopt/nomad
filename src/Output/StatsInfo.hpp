/**
 * \file   StatsInfo.hpp
 * \brief  Class for Stats info and display
 * \author Viviane Rochon Montplaisir, Christophe Tribes
 * \date   February 2018
 */

#ifndef __NOMAD400_STATSINFO__
#define __NOMAD400_STATSINFO__

#include <memory>   // For unique_ptr
#include <vector>
#include "../Math/ArrayOfDouble.hpp"
#include "../Math/Point.hpp"
#include "../Util/ArrayOfString.hpp"

#include "../nomad_nsbegin.hpp"


/// Types for DISPLAY_STATS and STATS_FILE parameters
/**
 \warning Do not modify the order
 */
enum class DisplayStatsType
{
    DS_OBJ        ,    ///< Objective (f) value
    //   (keep in first position)
    DS_CONS_H     ,    ///< Infeasibility (h) value
    DS_H_MAX      ,    ///< Max infeasibility (h) acceptable at the time of eval
    //DS_SMOOTH_OBJ ,    ///< Smoothed objective value (f~)
    //DS_SIM_BBE    ,    ///< Number of simulated bb evaluations
    DS_BBE        ,    ///< Number of bb evaluations
    DS_ALGO_BBE   ,    ///< Number of bb evaluations for a single algo run
    DS_BLK_EVA    ,    ///< Number of block evaluation calls
    DS_BLK_SIZE   ,    ///< Number of EvalPoints in the block
    DS_LAP        ,    ///< Number of evaluations since this lap started
    DS_SGTE       ,    ///< Number of surrogate evaluations since this surrogate started
    DS_TOTAL_SGTE ,    ///< Total number of surrogate evaluations
    DS_BBO        ,    ///< All blackbox outputs
    DS_EVAL       ,    ///< Number of evaluations
    DS_CACHE_HITS ,    ///< Number of cache hits
    DS_CACHE_SIZE ,    ///< Number of points in cache
    DS_ITER_NUM   ,    ///< Iteration number
    DS_TIME       ,    ///< Wall-clock time
    DS_MESH_INDEX ,    ///< Mesh index
    DS_MESH_SIZE  ,    ///< Mesh size parameter Delta^m_k
    DS_DELTA_M    ,    ///< Same as \c DS_MESH_SIZE
    DS_FRAME_SIZE ,    ///< Frame size parameter Delta^f_k
    DS_DELTA_F    ,    ///< Same as \c DS_FRAME_SIZE
    DS_SOL        ,    ///< Solution vector
    DS_THREAD_ALGO,    ///< Thread number for the algorithm 
    DS_THREAD_NUM ,    ///< Thread number in which this evaluation was done
    DS_GEN_STEP   ,    ///< Name of the step in which this point was generated
    //DS_VAR        ,    ///< One variable
    //DS_STAT_SUM   ,    ///< Stat sum
    //DS_STAT_AVG   ,    ///< Stat avg
    DS_USER         ,  ///< User-defined string
    DS_UNDEFINED       ///< Undefined value
    //   (keep in last position)
};

typedef ArrayOfString DisplayStatsTypeList;


/// Information for stats format (parameters DISPLAY_STATS and STATS_FILE).
/**
 Also holds information about stats file.
 */
class StatsInfo
{
private:
    // Stats infos
    Double          _obj;
    Double          _consH;
    Double          _hMax;
    size_t          _bbe;
    size_t          _algoBBE;
    size_t          _blkEva;
    size_t          _blkSize;
    std::string     _bbo;
    size_t          _eval;
    size_t          _cacheHits;
    size_t          _cacheSize;
    size_t          _iterNum;
    size_t          _time;
    ArrayOfDouble   _meshIndex;
    ArrayOfDouble   _meshSize;
    ArrayOfDouble   _frameSize;
    size_t          _lap;
    size_t          _sgte;
    size_t          _totalSgte;
    Point           _sol;
    int             _threadAlgoNum;
    int             _threadNum;
    bool            _relativeSuccess;   ///> Used for priting star, or when DISPLAY_ALL_EVAL is false.
    std::string     _comment;   ///> General comment, ex. Algorithm from where this point was generated.
    std::string     _genStep;   ///> Step in which this point was generated


public:
    /*---------------*/
    /* Class Methods */
    /*---------------*/

    /// Constructor
    explicit StatsInfo();

    // Destructor is not implemented, so the compilator creates all default functions.
    //virtual ~StatsInfo() {}

public:

    // Get/Set
    void setObj(const Double& obj)                  { _obj = obj; }
    void setConsH(const Double consH)               { _consH = consH; }
    void setHMax(const Double hMax)                 { _hMax = hMax; }
    void setBBE(const size_t bbe)                   { _bbe = bbe; }
    void setAlgoBBE(const size_t bbe)               { _algoBBE = bbe; }
    void setBlkEva(const size_t blkEva)             { _blkEva = blkEva; }
    void setBlkSize(const size_t blkSize)           { _blkSize = blkSize; }
    void setBBO(const std::string& bbo)             { _bbo = bbo; }
    void setEval(const size_t eval)                 { _eval = eval; }
    void setCacheHits(const size_t cacheHits)       { _cacheHits = cacheHits; }
    void setCacheSize(const size_t cacheSize)       { _cacheSize = cacheSize; }
    void setIterNum(const size_t iterNum)           { _iterNum = iterNum; }
    void setTime(const size_t time)                 { _time = time; }
    void setMeshIndex(const ArrayOfDouble meshIndex) { _meshIndex = meshIndex; }
    void setMeshSize(const ArrayOfDouble meshSize)   { _meshSize = meshSize; }
    void setFrameSize(const ArrayOfDouble frameSize) { _frameSize = frameSize; }
    void setLap(const size_t lap)                   { _lap = lap; }
    void setSgte(const size_t sgte)                 { _sgte = sgte; }
    void setTotalSgte(const size_t totalSgte)       { _totalSgte = totalSgte; }
    void setSol(const Point sol)                    { _sol = sol; }
    void setThreadAlgo(const int threadAlgoNum)     { _threadAlgoNum = threadAlgoNum; }
    void setThreadNum(const int threadNum)          { _threadNum = threadNum; }
    void setRelativeSuccess(bool relativeSuccess)   { _relativeSuccess = relativeSuccess; }
    void setComment(const std::string& comment)     { _comment = comment; }
    void setGenStep(const std::string& genStep)     { _genStep = genStep; }

    // Should this stats be printed even if DISPLAY_ALL_EVAL is false
    bool alwaysDisplay(const bool displayInfeasible, const bool displayUnsuccessful) const;

    /// Display header
    std::string displayHeader(const DisplayStatsTypeList& format,
                              const ArrayOfDouble & solFormat = ArrayOfDouble(),
                              const size_t objWidth = 0) const;

    /**
     Display with an empty format will output BBE OBJ
     */
    std::string display(const DisplayStatsTypeList& format,
                        const ArrayOfDouble & solFormat = ArrayOfDouble(),
                        const size_t objWidth = 0,
                        const size_t hWidth = 0,
                        const bool starSuccess = false,
                        const bool appendComment = true) const;

    // Convert a string in a DisplayStatsType.
    static DisplayStatsType stringToDisplayStatsType(const std::string& inputString,
                                                     std::string& format);
    static std::string displayStatsTypeToString(const DisplayStatsType& displayStatsType);

};

typedef std::unique_ptr<StatsInfo> StatsInfoUPtr;


inline std::ostream& operator<< (std::ostream& os, const DisplayStatsType& displayStatsType)
{
    return os << StatsInfo::displayStatsTypeToString(displayStatsType);
}


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_STATSINFO__

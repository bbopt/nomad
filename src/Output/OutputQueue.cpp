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

#include <fstream>
#include "../Output/OutputQueue.hpp"
#include "../Util/Exception.hpp"

// Static members initialization
#ifdef _OPENMP
omp_lock_t NOMAD::OutputQueue::_s_queue_lock;
#endif // _OPENMP

std::unique_ptr<NOMAD::OutputQueue> NOMAD::OutputQueue::_single(nullptr);

bool NOMAD::OutputQueue::_hasBeenInitialized = false;

// Private constructor
NOMAD::OutputQueue::OutputQueue()
  : _queue(),
    _params(),
    _statsFile(""),
    _statsWritten(false),
    _totalEval(0),
    _statsFileFormat(),
    _statsLineCount(0),
    _objWidth(),
    _hWidth(),
    _maxStepLevel(20),   // outputInfos with step level > 20 won't be printed.
    _maxOutputLevel(NOMAD::OutputLevel::LEVEL_DEBUGDEBUG),  // Level of details for display
    _indentLevel(0),
    _blockStart("{"),
    _blockEnd("}")
{
}


// Destructor
NOMAD::OutputQueue::~OutputQueue()
{
    // Always flush on destruction. In fact, the queue should be
    // empty at this point.
    if (! _queue.empty())
    {
        //std::cerr << "Calling destructor on a non-empty OutputQueue." << std::endl;
        flush();
    }
#ifdef _OPENMP
    omp_destroy_lock(&_s_queue_lock);
#endif // _OPENMP
    // Close stats file
    if (!_statsFile.empty())
    {
        if (!_statsWritten)
        {
            _statsStream << "no feasible solution has been found after " << NOMAD::itos(_totalEval) << " evaluations" << std::endl;
        }
        _statsStream.close();
    }
}

void NOMAD::OutputQueue::reset()
{
    // Flush the queue
    flush();

    // Close stats file
    if (!_statsFile.empty())
    {
        if (!_statsWritten)
        {
            _statsStream << "no feasible solution has been found after " << NOMAD::itos(_totalEval) << " evaluations" << std::endl;
        }
        _statsStream.close();
    }
    _statsWritten = false;
    _totalEval = 0;
    _hasBeenInitialized = false;
}


// Access to singleton
std::unique_ptr<NOMAD::OutputQueue>& NOMAD::OutputQueue::getInstance()
{
#ifdef _OPENMP
    #pragma omp critical(initOutputQueueLock)
    {
#endif // _OPENMP
        if (nullptr == _single)
        {
#ifdef _OPENMP
            omp_init_lock(&_s_queue_lock);
#endif // _OPENMP
            _single = std::unique_ptr<OutputQueue> (new OutputQueue());
        }
#ifdef _OPENMP
    } // end of critical section
#endif // _OPENMP

    return _single;
}

// Initialize display parameters.
void NOMAD::OutputQueue::initParameters(const std::shared_ptr<NOMAD::DisplayParameters>& params)
{
    // If already called, make sure to flush and close the stats file
    if(_hasBeenInitialized)
    {
       reset();
    }

    _params = params;

    if (nullptr == _params)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "OutputQueue::initParameters: Display Parameters are NULL");
    }

    setDisplayDegree(_params->getAttributeValue<int>("DISPLAY_DEGREE"));
    _maxStepLevel = _params->getAttributeValue<size_t>("DISPLAY_MAX_STEP_LEVEL");

    _objWidth = params->getAttributeValue<size_t>("OBJ_WIDTH");
    _hWidth   = _objWidth;

    // Some treatment needed on STATS_FILE.
    // First item is the name of the file.
    // Remaining items is the format.
    auto statsFileParam = _params->getAttributeValue<NOMAD::DisplayStatsTypeList>("STATS_FILE");

    std::string statsFileName = "";
    NOMAD::DisplayStatsTypeList statsFileFormat = statsFileParam;
    if (statsFileParam.size() > 1)
    {
        statsFileName = statsFileParam[0];
        statsFileFormat.erase(0);
        if (statsFileFormat.size() == 0)
        {
            // A default stats type is provided: BBE OBJ
            statsFileFormat.add(NOMAD::StatsInfo::displayStatsTypeToString(NOMAD::DisplayStatsType::DS_BBE));
            statsFileFormat.add(NOMAD::StatsInfo::displayStatsTypeToString(NOMAD::DisplayStatsType::DS_OBJ));
        }
    }
    setStatsFileName(statsFileName);
    initStatsFile();
    setStatsFileFormat(statsFileFormat);

    _hasBeenInitialized = true;

}


bool NOMAD::OutputQueue::goodLevel(const OutputLevel& outputLevel) const
{
    if (outputLevel <= _maxOutputLevel)
    {
        return true;
    }

    return (outputLevel <= NOMAD::OutputLevel::LEVEL_STATS
            && !_statsFile.empty());
}


// Set level of details in output
// 0 -> LEVEL_NOTHING
// 1 -> LEVEL_VERY_HIGH
// 2 -> LEVEL_NORMAL
// 3 -> LEVEL_INFO
// 4 -> LEVEL_DEBUG
// 5 -> LEVEL_DEBUGDEBUG
void NOMAD::OutputQueue::setDisplayDegree(const int displayDegree)
{
    NOMAD::OutputLevel outputLevel = NOMAD::OutputLevel::LEVEL_NORMAL;
    if (displayDegree == 0)
    {
        outputLevel = NOMAD::OutputLevel::LEVEL_NOTHING;
    }
    else if (displayDegree == 1)
    {
        outputLevel = NOMAD::OutputLevel::LEVEL_VERY_HIGH;
    }
    else if (displayDegree == 2)
    {
        outputLevel = NOMAD::OutputLevel::LEVEL_NORMAL;
    }
    else if (displayDegree == 3)
    {
        outputLevel = NOMAD::OutputLevel::LEVEL_INFO;
    }
    else if (displayDegree == 4)
    {
        outputLevel = NOMAD::OutputLevel::LEVEL_DEBUG;
    }
    else if (displayDegree == 5)
    {
        outputLevel = NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;
    }
    else
    {
        std::cerr << "Unrecognized display degree to set: " << displayDegree << std::endl;
    }

    _maxOutputLevel = outputLevel;
}


// Add Output info
void NOMAD::OutputQueue::add(NOMAD::OutputInfo outputInfo)
{
    if (!goodLevel(outputInfo.getOutputLevel()))
    {
        return;
    }

#ifdef _OPENMP
    // Acquire lock before adding a new element to the queue
    omp_set_lock(&_s_queue_lock);
#endif // _OPENMP
    _queue.push_back(std::move(outputInfo));
#ifdef _OPENMP
    omp_unset_lock(&_s_queue_lock);
#endif // _OPENMP
}

void NOMAD::OutputQueue::add(const std::string & s, NOMAD::OutputLevel outputLevel)
{
    if (!goodLevel(outputLevel))
    {
        return;
    }

    // Warning: No originator in this case
    OutputInfo outputInfo("", s, outputLevel);

    Add(std::move(outputInfo));
}


// Print all in the queue and flush.
void NOMAD::OutputQueue::flush()
{
    if (_queue.empty())
    {
        return;
    }

    if (_maxOutputLevel >= NOMAD::OutputLevel::LEVEL_DEBUGDEBUG)
    {
        // hyper-debug
        std::cout << "Output all " << _queue.size() << " elements." << std::endl;
    }

#ifdef _OPENMP
    // Lock queue before flush
    omp_set_lock(&_s_queue_lock);
#endif // _OPENMP

    // Info goes to Standard output
    for (std::vector<NOMAD::OutputInfo>::iterator out_it = _queue.begin();
         out_it != _queue.end(); ++out_it)
    {
        flushBlock(std::move(*out_it));
    }
    _queue.clear();
#ifdef _OPENMP
    omp_unset_lock(&_s_queue_lock);
#endif // _OPENMP

}


void NOMAD::OutputQueue::startBlock()
{
    std::cout << " " << _blockStart;
}


void NOMAD::OutputQueue::endBlock()
{
    std::cout << _blockEnd << " ";
}


void NOMAD::OutputQueue::indent(int level)
{
    // Indent each line by level.
    for (int i = 0; i < level; i++)
    {
        std::cout << "    ";
    }
}


void NOMAD::OutputQueue::flushBlock(const NOMAD::OutputInfo &outputInfo)
{
    // Output one line for each string in the vector msg.

    NOMAD::OutputLevel outputLevel = outputInfo.getOutputLevel();
    const NOMAD::StatsInfo* statsInfo = outputInfo.getStatsInfo();

    flushStatsToStatsFile(statsInfo);

    if (outputLevel > _maxOutputLevel)
    {
        // Do not display to stdout
        return;
    }

    NOMAD::ArrayOfString msg = outputInfo.getMsg();
    if (outputLevel == NOMAD::OutputLevel::LEVEL_STATS)
    {
        flushStatsToStdout(statsInfo);
    }

    else
    {
        if (outputInfo.isBlockEnd())
        {
            _indentLevel--;
            if (_indentLevel < 0)
            {
    #ifdef _OPENMP
                omp_unset_lock(&_s_queue_lock);
    #endif // _OPENMP
                throw NOMAD::Exception(__FILE__, __LINE__, "OutputQueue has more block ends than block starts.");
            }
        }

        // Verify step level is high enough in the tree to be displayed.
        if (_indentLevel <= (int)_maxStepLevel)
        {
            for (size_t i = 0; i < msg.size(); i++)
            {
                indent(_indentLevel);
                if (outputInfo.isBlockEnd())
                {
                    endBlock();
                }

                std::cout << msg[i];

                if (outputInfo.isBlockStart())
                {
                    startBlock();
                }
                std::cout << std::endl;
            }
        }
        else
        {
            // Show some points so the user knows something is happening.
            // This display is not sophisticated.
            if (_indentLevel == (int)_maxStepLevel + 1)
            {
                indent(_indentLevel);
                std::cout << "........................................" << std::endl;
            }
        }


        if (outputInfo.isBlockStart())
        {
            _indentLevel++;
        }

    }

}


void NOMAD::OutputQueue::flushStatsToStdout(const NOMAD::StatsInfo *statsInfo)
{
    // Early out
    if (nullptr == statsInfo)
    {
        return;
    }

    if (nullptr == _params)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "OutputQueue: Display Parameters are NULL");
    }

    bool displayInfeasible      = _params->getAttributeValue<bool>("DISPLAY_INFEASIBLE");
    bool displayUnsuccessful    = _params->getAttributeValue<bool>("DISPLAY_UNSUCCESSFUL");
    bool displayAllEval         = _params->getAttributeValue<bool>("DISPLAY_ALL_EVAL");
    size_t displayHeaderFreq    = _params->getAttributeValue<size_t>("DISPLAY_HEADER");
    if ( _maxOutputLevel > NOMAD::OutputLevel::LEVEL_NORMAL)
    {
        // Do not display header if display is already heavy.
        displayHeaderFreq = NOMAD::INF_SIZE_T;
    }
    auto displayStatsFormat     = _params->getAttributeValue<NOMAD::ArrayOfString>("DISPLAY_STATS");
    bool displayInteresting     = statsInfo->alwaysDisplay(displayInfeasible, displayUnsuccessful, false);

    if (displayAllEval || displayInteresting)
    {
        // Set starSuccess and appendComment to false for intermediate computations.
        // They will be set to appropriate values just before display to
        // standard output.
        bool starSuccess = false;
        bool appendComment = false;

        // Update objWidth if OBJ went larger.
        NOMAD::DisplayStatsTypeList list;
        list.add("OBJ");
        auto solFormat = _params->getAttributeValue<NOMAD::ArrayOfDouble>("SOL_FORMAT");
        std::string objString = statsInfo->display(list, solFormat, _objWidth, _hWidth, starSuccess, appendComment);
        if (objString.size() > _objWidth)
        {
            _objWidth = objString.size();
        }
        list.clear();
        list.add("CONS_H");
        std::string hString = statsInfo->display(list, solFormat, _objWidth, _hWidth, starSuccess, appendComment);
        if (hString.size() > _hWidth)
        {
            _hWidth = hString.size();
        }
        // Also consider H_MAX
        list.clear();
        list.add("H_MAX");
        hString = statsInfo->display(list, solFormat, _objWidth, _hWidth, starSuccess, appendComment);
        if (hString.size() > _hWidth)
        {
            _hWidth = hString.size();
        }

        // Display column titles every DISPLAY_HEADER lines.
        if ( (displayHeaderFreq < NOMAD::INF_SIZE_T )
             && (0 == (_statsLineCount % displayHeaderFreq)))
        {
            if (_statsLineCount > 0)
            {
                std::cout << std::endl;
            }
            std::cout << statsInfo->displayHeader(displayStatsFormat, solFormat, _objWidth) << std::endl;
        }

        // To standard output
        // Set starSuccess: On standard output, add a star for successful evaluations,
        // when either all evaluations or all unsuccessful evaluations are shown.
        starSuccess = displayAllEval || displayUnsuccessful;
        appendComment = true;
        std::cout << statsInfo->display(displayStatsFormat, solFormat, _objWidth, _hWidth, starSuccess, appendComment) << std::endl;
        _statsLineCount++;
    }
}


void NOMAD::OutputQueue::initStatsFile()
{
    if (!_statsFile.empty())
    {
        // Open stats file and clear it (trunc)
        _statsStream.close();
        _statsStream.open(_statsFile.c_str(), std::ofstream::out | std::ios::trunc);
        if (_statsStream.fail())
        {
            std::cerr << "Warning: could not open stats file " << _statsFile << std::endl;
        }
        _statsStream.setf(std::ios::fixed);
        // Set full precision on stats file.
        _statsStream.precision(NOMAD::DISPLAY_PRECISION_FULL);
    }
}


void NOMAD::OutputQueue::flushStatsToStatsFile(const NOMAD::StatsInfo *statsInfo)
{
    // Early out
    if (_statsFile.empty() || nullptr == statsInfo)
    {
        return;
    }

    if (nullptr == _params)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "OutputQueue: Display Parameters are NULL");
    }
    // Display this statsInfo (to standard output and to stats file)
    // only if parameter DISPLAY_ALL_EVAL is true, and if it
    // is interesting to display.
    bool displayInfeasible      = _params->getAttributeValue<bool>("DISPLAY_INFEASIBLE");
    bool displayUnsuccessful    = _params->getAttributeValue<bool>("DISPLAY_UNSUCCESSFUL");
    bool displayInteresting     = statsInfo->alwaysDisplay(displayInfeasible, displayUnsuccessful, true);
    auto n = _params->getAttributeValue<NOMAD::ArrayOfDouble>("SOL_FORMAT").size();
    NOMAD::ArrayOfDouble solFormatStats(n, NOMAD::DISPLAY_PRECISION_FULL);

    if (displayInteresting)
    {
        // Add stats information.
        _statsStream << statsInfo->display(_statsFileFormat, solFormatStats, 0, 0, false, false) << std::endl;
        _statsWritten = true;
    }
}




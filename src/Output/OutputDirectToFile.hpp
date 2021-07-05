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
#ifndef __NOMAD_4_0_OUTPUTDIRECTTOFILE__
#define __NOMAD_4_0_OUTPUTDIRECTTOFILE__

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

/// Direct output of info to file (history, solution,...).
/**
 The output is a singleton. Some info can be written into files. \n
 The format of output is fixed. The parameters (DisplayParameters) are attributes of the class provided by calling OutputDirectToFile::initParameters. New files to receive output must be  registered in this function.\n
 */
class OutputDirectToFile
{
private:
    /// Private constructor
    OutputDirectToFile();

public:
    // Destructor
     virtual ~OutputDirectToFile();

    /// Access to singleton
    static std::unique_ptr<OutputDirectToFile>& getInstance();

    /// Initialization of file names using display parameters
    void init(const std::shared_ptr<DisplayParameters>& params);

    /// When history and/or solution files are active, write info in solution and history file according to the flags
    void write(const StatsInfo& outInfo, bool writeInSolutionFile, bool writeInHistoryFile=true);

    static void Write(const StatsInfo & outInfo, bool writeInSolutionFile, bool writeInHistoryFile = true)
    {
        getInstance()->write(outInfo,writeInSolutionFile,writeInHistoryFile);
    }

    /// Good to write in history and/or solution files when the file names have been defined.
    bool goodToWrite() const;
    static bool GoodToWrite()
    {
        return getInstance()->goodToWrite();
    }

    void reset() { _hasBeenInitialized = false; } // Allow to reset history file and solution file. Existing files will be overwritten

    /// Not used now
    void setOutputFileFormat(const DisplayStatsTypeList& outputFileFormat)
    {
        _outputFileFormat = outputFileFormat;
    }
    const DisplayStatsTypeList& getOutputFileFormat() const { return _outputFileFormat; }

    /// Used to enable solution file (requires a solution file name)
    void enableSolutionFile() { _enabledSolutionFile = true; }

    /// Disable solution file (temporarily, for example during PhaseOne)
    void disableSolutionFile() { _enabledSolutionFile = false; }

#define OUTPUT_DIRECTTOFILE_START if (OutputDirectToFile::GoodToWrite()) {
#define OUTPUT_DIRECTTOFILE_END }

private:
#ifdef _OPENMP
    // Acquire lock before write.
    // NOTE It does not seem relevant for the lock to be static,
    // because OutputDirectToFile is a singleton anyway. If staticity causes problems,
    // we could remove the static keyword.
    DLL_UTIL_API static omp_lock_t  _s_output_lock;
#endif // _OPENMP

    /// Helper for init
    void initHistoryFile();

    DLL_UTIL_API static bool        _hasBeenInitialized;    ///< Flag for initialization (initialization cannot be performed more than once).


    DLL_UTIL_API static std::unique_ptr<OutputDirectToFile> _single;    ///< The singleton

    size_t                          _outputSize;

    /**
     Format for output in a file.
     Might include some raw strings, do not convert to DisplayStatsType.
     */
    DisplayStatsTypeList            _outputFileFormat;

    std::string                     _solutionFile;
    std::ofstream                   _solutionStream;

    std::string                     _historyFile;
    std::ofstream                   _historyStream;

    /// Even if solution file is provided we can temporarily disable solution file (PhaseOne)
    bool                            _enabledSolutionFile;

};

#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_OUTPUTDIRECTTOFILE__

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
/**
 \file   fileutils.hpp
 \brief  Utility functions about files (headers)
 \author Viviane Rochon Montplaisir
 \date   June 2017
 \see    fileutils.cpp
 */
#ifndef __NOMAD_4_0_FILEUTILS__
#define __NOMAD_4_0_FILEUTILS__

// use of 'access' or '_access', and getpid() or _getpid():
#ifdef _MSC_VER
#include <io.h>
//#include <process.h>
#else
#include <unistd.h>
#endif

#include <fstream>

#include "../Util/defines.hpp"

#include "../nomad_nsbegin.hpp"

// Copied from NOMAD_3.
/// Check if a file exists and is executable.
/**
 \param filename  A string corresponding to a file name -- \b IN.
 \return          A boolean equal to \c true if the file is executable.
 */
bool checkExeFile(const std::string &filename);


/// Check if a file exists and is readable.
/**
 \param filename  A string corresponding to a file name -- \b IN.
 \return          A boolean equal to \c true if the file exists and is readable.
 */
bool checkReadFile(const std::string &filename);


/// Check if a file exists and is writable.
/**
 \param filename  A string corresponding to a file name -- \b IN.
 \return          A boolean equal to \c true if the file exists and is writable.
 */
bool checkWriteFile(const std::string &filename);


// Get current directory
std::string curdir();


// Extract directory from the given filename.
// If there is no directory, return ".".
std::string dirname(const std::string &filename);

// Extract file root from the given filename.
// Ex. "/path/toto.txt" returns "toto"
std::string rootname(const std::string &filename);

// Extract extension from the given filename.
// Ex. "/path/toto.txt" returns ".txt"
std::string extension(const std::string &filename);

// If filename has a path, leave it.
// If it doesn't, add dirname() to it.
std::string fullpath(const std::string &filename);


// Return true if a filename is absolute, i.e., starts with '/' (DIR_SEP).
// Return false otherwise (filename is relative, or filename is file only).
bool isAbsolute(const std::string &filename);


// Add a '/' (DIR_SEP) at the end of the dirname if it does not end with one.
void ensureDirPath(std::string &dirname);


// Input a line (from a parameters file).
// Remove comments starting with '#'.
// Replace tabs by spaces.
// Trim extra spaces.
void removeComments(std::string &line);

// Add problemDir to filename.
// Add seed if desired.
void completeFileName(std::string &filename,
                      const std::string &problemDir,
                      bool addSeed = false,
                      int seed = 0);

void addSeedToFileName(size_t nSeed,
                       const std::string& sSeed,
                       std::string& filename);


// Write to file
// Use T::operator<< (examples of T: CacheSet, Mads)
// Will break if T::operator<< is not defined.
template<typename T>
bool write(const T &info, const std::string &filename)
{
    bool writeSuccess = true;
    std::ofstream fout;

    if (filename.empty())
    {
        std::cerr << "Warning: " << typeid(T).name() << ": Cannot write to file: file name is not defined.";
        writeSuccess = false;
    }

    if (writeSuccess)
    {
        fout.open(filename.c_str(), std::ofstream::out);
        if (fout.fail())
        {
            std::cerr << "Warning: " << typeid(T).name() << ": Cannot write to file " + filename << std::endl;
            writeSuccess = false;
            fout.close();
        }
    }

    if (writeSuccess)
    {
        fout.clear();
        fout << info;
        fout.close();
    }


    return writeSuccess;
}


// Read from file
// Use T::operator>> (examples of T: CacheSet, Mads)
// Will break if T::operator>> is not defined.
template<typename T>
bool read(T &info, const std::string &filename)
{
    bool readSuccess = true;
    std::ifstream fin;

    if (filename.empty())
    {
        std::cerr << "Warning: " << typeid(T).name() << ": Cannot read file: file name is not defined.";
        readSuccess = false;
    }

    if (readSuccess)
    {
        if (!checkReadFile(filename))
        {
            std::cerr << "Warning: " << typeid(T).name() << ": File does not exist or cannot be read: " + filename << std::endl;
            readSuccess = false;
        }
    }

    if (readSuccess)
    {
        fin.open(filename.c_str(), std::ifstream::out);
        if (fin.fail())
        {
            std::cerr << "Warning: " << typeid(T).name() << ": Cannot read from file " + filename << std::endl;
            readSuccess = false;
            fin.close();
        }
    }

    if (readSuccess)
    {
        fin >> info;
    }

    fin.close();

    return readSuccess;
}


// Read all file and store it in string info.
bool readAllFile(std::string &info, const std::string &filename);


#include "../nomad_nsend.hpp"

#endif // __NOMAD_4_0_FILEUTILS__

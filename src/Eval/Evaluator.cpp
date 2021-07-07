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
#include "../Eval/Evaluator.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Util/fileutils.hpp"
#include <fstream>  // For ofstream
#include <stdio.h>  // For popen
#ifndef _WIN32
#include <unistd.h> // for getpid
#else
#include <process.h>
#define getpid _getpid
#define popen  _popen
#define pclose _pclose
#endif

// Initialize statics
std::vector<std::string> NOMAD::Evaluator::_tmpFiles = std::vector<std::string>();

namespace {
    // the cleanup of temporary files at program shutdown needs to remain with this
    // translation unit to guarantee the correct destruction order
    struct TmpFilesCleanup {
        ~TmpFilesCleanup() { NOMAD::Evaluator::removeTmpFiles(); }
    } _TmpFilesCleanup;
}

//
// Constructor
//
// NOTE: The full path to BB_EXE is added during check
NOMAD::Evaluator::Evaluator(
                    const std::shared_ptr<NOMAD::EvalParameters> &evalParams,
                    const NOMAD::EvalType evalType,
                    const NOMAD::EvalXDefined evalXDefined)
  : _evalParams(evalParams),
    _evalXDefined(evalXDefined),
    _evalType(evalType)
{
}


NOMAD::Evaluator::~Evaluator()
{
}


void NOMAD::Evaluator::initializeTmpFiles(const std::string& tmpDir)
{
    // Initialize tmp files for Evaluators
    int nbThreads = 1;
#ifdef _OPENMP
    nbThreads = omp_get_max_threads();
#endif
    std::string tmppath = tmpDir;
    NOMAD::ensureDirPath(tmppath);
    // Use the pid in the file name in case two nomad run at the same time.
    int pid = getpid();
    // Create a temporary file fo blackbox input. One for each thread number,
    // for each nomad pid. Add the file names to _tmpFiles.
    _tmpFiles.clear();
    for (auto threadNum = 0; threadNum < nbThreads; threadNum++)
    {
        std::string tmpfilestr = tmppath + "nomadtmp." + std::to_string(pid) + "." + std::to_string(threadNum);
        _tmpFiles.push_back(tmpfilestr);
    }
}


void NOMAD::Evaluator::removeTmpFiles()
{
    // Remove all temporary files, so that they do not linger around.
    auto nbThreads = _tmpFiles.size();
    for (size_t i = 0; i < nbThreads; i++)
    {
        remove(_tmpFiles[i].c_str());
    }
    _tmpFiles.clear();
}


// Default eval_x: System call to a black box that was provided
// via parameter BB_EXE and set through Evaluator::setBBExe().
bool NOMAD::Evaluator::eval_x(NOMAD::EvalPoint &x,
                              const NOMAD::Double& hMax,
                              bool &countEval) const
{
    // The user might have defined his own eval_x() for NOMAD::EvalPoint.
    // In the NOMAD code, we do not use this method.
    //
    // Implemented to be used by the Runner. In the case of the Runner,
    // eval_x is redefined. When batch mode is used (for instance for
    // Styrene), this eval_x is called. So in fact we really want to
    // use the executable defined by BB_EXE.

    _evalXDefined = NOMAD::EvalXDefined::USE_BB_EVAL;

    // Create a block of one point and evaluate it.
    NOMAD::Block block;
    std::shared_ptr<NOMAD::EvalPoint> epp = std::make_shared<NOMAD::EvalPoint>(x);
    block.push_back(epp);

    std::vector<bool> countEvalVector(1, countEval);
    std::vector<bool> evalOkVector(1, false);

    // Call eval_block
    evalOkVector = eval_block(block, hMax, countEvalVector);

    // Update x and countEval
    x = *epp;
    countEval = countEvalVector[0];

    return evalOkVector[0];
}


// Default eval_block: for block
// This is used even for blocks of 1 point.
// If we never go through this eval_block(),
// it means that eval_block was redefined by the user,
// using library mode.
std::vector<bool> NOMAD::Evaluator::eval_block(NOMAD::Block &block,
                                               const NOMAD::Double &hMax,
                                               std::vector<bool> &countEval) const
{
    std::vector<bool> evalOk(block.size(), false);
    countEval.resize(block.size(), false);

    // Verify there is at least one point to evaluate
    if (0 == block.size())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Evaluator: eval_block called with an empty block");
    }

    // Verify all points are completely defined
    for (auto it = block.begin(); it != block.end(); it++)
    {
        if (!(*it)->isComplete())
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Evaluator: Incomplete point " + (*it)->display());
        }
    }

    // Start evaluation
    for (auto it = block.begin(); it != block.end(); it++)
    {
        // Debugging. EVAL should already be IN_PROGRESS.
        if (NOMAD::EvalStatusType::EVAL_IN_PROGRESS != (*it)->getEvalStatus(_evalType))
        {
#ifdef _OPENMP
            #pragma omp critical(warningEvalX)
#endif
            {
                std::cerr << "Warning: EVAL should already be IN_PROGRESS for point " << (*it)->display() << std::endl;
            }
        }
    }

    if (NOMAD::EvalXDefined::EVAL_BLOCK_DEFINED_BY_USER == _evalXDefined)
    {
        // First time that eval_block() is called.
        // Obviously, eval_block() was not redefined by user.
        // If the blackbox is external, USE_BB_EVAL is already set by the
        // constructor. Hence, eval_x() is defined by user.
        _evalXDefined = NOMAD::EvalXDefined::EVAL_X_DEFINED_BY_USER;
    }

    if (NOMAD::EvalXDefined::USE_BB_EVAL == _evalXDefined)
    {
        evalOk = evalXBBExe(block, hMax, countEval);
    }
    else if (NOMAD::EvalXDefined::EVAL_X_DEFINED_BY_USER == _evalXDefined)
    {
        for (size_t index = 0; index < block.size(); index++)
        {
            bool countEval1 = false;
            evalOk[index] = eval_x(*block[index], hMax, countEval1);
            countEval[index] = countEval1;
        }
    }
    else
    {
        std::string s = "Error: This value of EvalXDefined is not processed: ";
        s += std::to_string((int)_evalXDefined);
        throw NOMAD::Exception(__FILE__, __LINE__, s);
    }

    return evalOk;
}


// eval_x() called in batch. Use system command and use the blackbox executable
// provided by parameter BB_EXE.
std::vector<bool> NOMAD::Evaluator::evalXBBExe(NOMAD::Block &block,
                                               const NOMAD::Double &hMax,
                                               std::vector<bool> &countEval) const
{
    std::vector<bool> evalOk(block.size(), false);

    // At this point, we are for sure in batch mode.
    // Verify blackbox executable defined by BB_EXE is available and executable.
    std::string bbExe;
    switch (_evalType)
    {
        case NOMAD::EvalType::BB:
            bbExe = _evalParams->getAttributeValue<std::string>("BB_EXE");
            break;
        case NOMAD::EvalType::SURROGATE:
            bbExe = _evalParams->getAttributeValue<std::string>("SURROGATE_EXE");
            break;
        default:
            std::string err = "Evaluator: No executable supported for EvalType ";
            err += NOMAD::evalTypeToString(_evalType);
            throw NOMAD::Exception(__FILE__,__LINE__,err);
    }

    if (bbExe.empty())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Evaluator: No blackbox executable defined.");
    }

    const int threadNum = NOMAD::getThreadNum();

    // Write a temp file for x0 and give that file as argument to bbExe.
    if ((size_t)threadNum >= _tmpFiles.size())
    {
        std::cerr << "Error: Evaluator: No temp file available." << std::endl;
        // Ugly early return
        return evalOk;
    }
    std::string tmpfile = _tmpFiles[threadNum];

    // System command
    std::ofstream xfile;
    // Open xfile and clear it (trunc)
    xfile.open(tmpfile.c_str(), std::ofstream::trunc);
    if (xfile.fail())
    {
        for (auto it = block.begin(); it != block.end(); it++)
        {
            (*it)->setEvalStatus(NOMAD::EvalStatusType::EVAL_ERROR, _evalType);
            std::cerr << "Error writing point " << (*it)->display() << " to temporary file \"" << tmpfile << "\"" << std::endl;
        }
        // Ugly early return
        return evalOk;
    }

    auto evalFormat = _evalParams->getAttributeValue<NOMAD::ArrayOfDouble>("BB_EVAL_FORMAT");
    for (auto it = block.begin(); it != block.end(); it++)
    {
        std::shared_ptr<NOMAD::EvalPoint> x = (*it);
        for (size_t i = 0; i < x->size(); i++)
        {
            if (i != 0)
            {
                xfile << " ";
            }
            xfile << (*x)[i].display(static_cast<int>(evalFormat[i].todouble()));
        }
        xfile << std::endl;
    }
    xfile.close();

    std::string cmd = bbExe + " " + tmpfile;
    std::string s;
    OUTPUT_DEBUG_START
    s = "System command: " + cmd;
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_DEBUGDEBUG);
    OUTPUT_DEBUG_END

    FILE *fresult = popen(cmd.c_str(), "r");
    if (!fresult)
    {
        // Something went wrong with the evaluation.
        // Point could be re-submitted.
        for (auto it = block.begin(); it != block.end(); it++)
        {
            (*it)->setEvalStatus(NOMAD::EvalStatusType::EVAL_ERROR, _evalType);
#ifdef _OPENMP
            #pragma omp critical(warningEvalX)
#endif
            {
                std::cerr << "Warning: Evaluation error with point " << (*it)->display() << std::endl;
            }
        }
    }

    else
    {
        for (size_t index = 0; index < block.size(); index++)
        {
            std::shared_ptr<NOMAD::EvalPoint> x = block[index];

            char buffer[1024];
            char *outputLine = nullptr;

            size_t nbTries=0;
            while (nbTries < 5)
            {
                nbTries++;

                outputLine = fgets(buffer, sizeof(buffer), fresult);

                if( feof(fresult) )
                { // c-stream eof detected. Output is empty, break the loop
                    x->setEvalStatus(NOMAD::EvalStatusType::EVAL_ERROR, _evalType);
#ifdef _OPENMP
#pragma omp critical(warningEvalX)
#endif
                    {
                        std::cerr << "Warning: Evaluation error with point " << x->display() << ": output is empty" << std::endl;
                    }
                    break;
                }

                if (NULL != outputLine)
                {
                    // Evaluation succeeded. Get and process blackbox output.
                    std::string bbo(outputLine);
                    // delete trailing '\n'
                    bbo.erase(bbo.size() - 1);

                    // Process blackbox output
                    auto bbOutputTypeList = _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
                    x->setBBO(bbo, bbOutputTypeList, _evalType);
                    auto bbOutput = x->getEval(_evalType)->getBBOutput();

                    evalOk[index] = bbOutput.getEvalOk();
                    countEval[index] = bbOutput.getCountEval(bbOutputTypeList);

                    break;
                }
            }
            // The number of tries has been reached (not eof) and still cannot read output file.
            if( ! feof(fresult) && NULL == outputLine )
            {
                // Something went wrong with the evaluation.
                // Point could be re-submitted.
                x->setEvalStatus(NOMAD::EvalStatusType::EVAL_ERROR, _evalType);
#ifdef _OPENMP
#pragma omp critical(warningEvalX)
#endif
                {
                    std::cerr << "Warning: Evaluation error with point " << x->display() << ": output is empty" << std::endl;
                }
            }
        }
        // Get exit status of the bb.exe. If it is not 0, there was an error.
        int exitStatus = pclose(fresult);

        size_t index = 0;   // used to update evalOk
        for (auto it = block.begin(); it != block.end(); it++)
        {
            std::shared_ptr<NOMAD::EvalPoint> x = (*it);
            if (exitStatus)
            {
                evalOk[index] = false;
                x->setEvalStatus(NOMAD::EvalStatusType::EVAL_ERROR, _evalType);
                s = "Warning: Evaluator returned exit status ";
                s += std::to_string(exitStatus);
                s += " for point: " + x->getX()->NOMAD::Point::display();
                NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_WARNING);
#ifdef _OPENMP
                #pragma omp critical(warningEvalX)
#endif
                {
                    std::cerr << s << std::endl;
                }
            }
            else if (!evalOk[index])
            {
                x->setEvalStatus(NOMAD::EvalStatusType::EVAL_FAILED, _evalType);
            }
            else
            {
                x->setEvalStatus(NOMAD::EvalStatusType::EVAL_OK, _evalType);
            }

            index++;
        }
    }

    return evalOk;
}

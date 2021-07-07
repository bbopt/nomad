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
 * \file   StatsInfo.cpp
 * \brief  Class for Stats info and display
 * \author Viviane Rochon Montplaisir, Christophe Tribes
 * \date   February 2018
 */

#include "../Output/StatsInfo.hpp"
#include "../Util/fileutils.hpp"

// Constructor 1
NOMAD::StatsInfo::StatsInfo()
  : _obj(),
    _consH(),
    _hMax(),
    _bbe(0),
    _feasBBE(0),
    _infBBE(0),
    _nbRelativeSuccess(0),
    _PhaseOneSuccess(0),
    _algoBBE(0),
    _blkEva(0),
    _blkSize(0),
    _bbo(),
    _eval(0),
    _cacheHits(0),
    _cacheSize(0),
    _iterNum(0),
    _time(0),
    _meshIndex(),
    _meshSize(),
    _frameSize(),
    _frameCenter(),
    _direction(),
    _lap(0),
    _modelEval(0),
    _totalModelEval(0),
    _sol(),
    _surrogateEval(0),
    _threadAlgoNum(0),
    _threadNum(0),
    _relativeSuccess(false),
    _comment(""),
    _genStep(""),
    _success(NOMAD::SuccessType::NOT_EVALUATED)
{
}


bool NOMAD::StatsInfo::alwaysDisplay(const bool displayInfeasible,
                                     const bool displayUnsuccessful,
                                     const bool forStatsFile) const
{
    bool doDisplay = false;
    if (!_obj.isDefined())
    {
        doDisplay = false;
    }
    else if (_bbe <= 1 && !forStatsFile)
    {
        // Always display X0 evaluation to standard output
        doDisplay = true;
    }
    else if (displayInfeasible || (_consH.isDefined() && _consH == 0.0))
    {
        if (displayUnsuccessful || _relativeSuccess)
        {
            doDisplay = true;
        }
    }

    return doDisplay;
}


// Convert a string like "BBE", "BBO", "OBJ", "SOL", "TIME"... to
// the corresponding NOMAD::DisplayStatsType.
// "%"
NOMAD::DisplayStatsType NOMAD::StatsInfo::stringToDisplayStatsType(const std::string& inputStr, std::string& format)
{
    NOMAD::DisplayStatsType ret;
    std::string s = inputStr;

    std::string tag;
    if (NOMAD::separateFormat(s, format, tag))
    {
        s = tag;
    }

    NOMAD::toupper(s);

    if (s == "OBJ")
    {
        ret = NOMAD::DisplayStatsType::DS_OBJ;
    }
    else if (s == "CONS_H")
    {
        ret = NOMAD::DisplayStatsType::DS_CONS_H;
    }
    else if (s == "H_MAX")
    {
        ret = NOMAD::DisplayStatsType::DS_H_MAX;
    }
    else if (s == "BBE")
    {
        ret = NOMAD::DisplayStatsType::DS_BBE;
    }
    else if (s == "FEAS_BBE")
    {
        ret = NOMAD::DisplayStatsType::DS_FEAS_BBE;
    }
    else if (s == "INF_BBE")
    {
        ret = NOMAD::DisplayStatsType::DS_INF_BBE;
    }
    else if (s == "REL_SUCC")
    {
        ret = NOMAD::DisplayStatsType::DS_REL_SUCC;
    }
    else if (s == "PHASE_ONE_SUCC")
    {
        ret = NOMAD::DisplayStatsType::DS_PHASE_ONE_SUCC;
    }
    else if (s == "ALGO_BBE")
    {
        ret = NOMAD::DisplayStatsType::DS_ALGO_BBE;
    }
    else if (s == "BLK_EVA")
    {
        ret = NOMAD::DisplayStatsType::DS_BLK_EVA;
    }
    else if (s == "BLK_SIZE")
    {
        ret = NOMAD::DisplayStatsType::DS_BLK_SIZE;
    }
    else if (s == "BBO")
    {
        ret = NOMAD::DisplayStatsType::DS_BBO;
    }
    else if (s == "EVAL")
    {
        ret = NOMAD::DisplayStatsType::DS_EVAL;
    }
    else if (s == "CACHE_HITS")
    {
        ret = NOMAD::DisplayStatsType::DS_CACHE_HITS;
    }
    else if (s == "CACHE_SIZE")
    {
        ret = NOMAD::DisplayStatsType::DS_CACHE_SIZE;
    }
    else if (s == "ITER_NUM")
    {
        ret = NOMAD::DisplayStatsType::DS_ITER_NUM;
    }
    else if (s == "TIME")
    {
        ret = NOMAD::DisplayStatsType::DS_TIME;
    }
    else if (s == "MESH_INDEX")
    {
        ret = NOMAD::DisplayStatsType::DS_MESH_INDEX;
    }
    else if (s == "MESH_SIZE" || s == "DELTA_M")
    {
        ret = NOMAD::DisplayStatsType::DS_MESH_SIZE;
    }
    // POLL_SIZE and DELTA_P are for backwards compatibility.
    else if (s == "FRAME_SIZE" || s == "DELTA_F" || s == "POLL_SIZE" || s == "DELTA_P")
    {
        ret = NOMAD::DisplayStatsType::DS_FRAME_SIZE;
    }
    else if (s == "FRAME_CENTER")
    {
        ret = NOMAD::DisplayStatsType::DS_FRAME_CENTER;
    }
    else if (s == "DIRECTION")
    {
        ret = NOMAD::DisplayStatsType::DS_DIRECTION;
    }
    else if (s == "LAP")
    {
        ret = NOMAD::DisplayStatsType::DS_LAP;
    }
    else if (s == "MODEL_EVAL")
    {
        ret = NOMAD::DisplayStatsType::DS_MODEL_EVAL;
    }
    else if (s == "SOL")
    {
        ret = NOMAD::DisplayStatsType::DS_SOL;
    }
    else if (s == "SURROGATE_EVAL")
    {
        ret = NOMAD::DisplayStatsType::DS_SURROGATE_EVAL;
    }
    else if (s == "THREAD_ALGO")
    {
        ret = NOMAD::DisplayStatsType::DS_THREAD_ALGO;
    }
    else if (s == "THREAD_NUM")
    {
        ret = NOMAD::DisplayStatsType::DS_THREAD_NUM;
    }
    else if (s == "GEN_STEP")
    {
        ret = NOMAD::DisplayStatsType::DS_GEN_STEP;
    }
    else if (s == "SUCCESS_TYPE")
    {
        ret = NOMAD::DisplayStatsType::DS_SUCCESS_TYPE;
    }
    else if (s == "TOTAL_MODEL_EVAL")
    {
        ret = NOMAD::DisplayStatsType::DS_TOTAL_MODEL_EVAL;
    }
    else
    {
        // Don't throw an exception.
        // Any string is accepted for output, for instance parenthesis
        // around the SOL, free text, etc.
        // If some string should be rejected, then the type
        // DS_UNDEFINED could be used.
        ret = NOMAD::DisplayStatsType::DS_USER;
    }

    return ret;
}


std::string NOMAD::StatsInfo::displayStatsTypeToString(const NOMAD::DisplayStatsType& displayStatsType)
{
    switch(displayStatsType)
    {
        case NOMAD::DisplayStatsType::DS_OBJ:
            return "OBJ";
        case NOMAD::DisplayStatsType::DS_CONS_H:
            return "CONS_H";
        case NOMAD::DisplayStatsType::DS_H_MAX:
            return "H_MAX";
        case NOMAD::DisplayStatsType::DS_BBE:
            return "BBE";
        case NOMAD::DisplayStatsType::DS_FEAS_BBE:
            return "FEAS_BBE";
        case NOMAD::DisplayStatsType::DS_INF_BBE:
            return "INF_BBE";
        case NOMAD::DisplayStatsType::DS_REL_SUCC:
            return "REL_SUCC";
        case NOMAD::DisplayStatsType::DS_PHASE_ONE_SUCC:
            return "PHASE_ONE_SUCC";
        case NOMAD::DisplayStatsType::DS_ALGO_BBE:
            return "ALGO_BBE";
        case NOMAD::DisplayStatsType::DS_BLK_EVA:
            return "BLK_EVA";
        case NOMAD::DisplayStatsType::DS_BLK_SIZE:
            return "BLK_SIZE";
        case NOMAD::DisplayStatsType::DS_BBO:
            return "BBO";
        case NOMAD::DisplayStatsType::DS_EVAL:
            return "EVAL";
        case NOMAD::DisplayStatsType::DS_CACHE_HITS:
            return "CACHE_HITS";
        case NOMAD::DisplayStatsType::DS_CACHE_SIZE:
            return "CACHE_SIZE";
        case NOMAD::DisplayStatsType::DS_ITER_NUM:
            return "ITER_NUM";
        case NOMAD::DisplayStatsType::DS_TIME:
            return "TIME";
        case NOMAD::DisplayStatsType::DS_MESH_INDEX:
            return "MESH_INDEX";
        case NOMAD::DisplayStatsType::DS_MESH_SIZE:
            return "MESH_SIZE";
        case NOMAD::DisplayStatsType::DS_DELTA_M:
            return "DELTA_M";
        case NOMAD::DisplayStatsType::DS_FRAME_SIZE:
            return "FRAME_SIZE";
        case NOMAD::DisplayStatsType::DS_FRAME_CENTER:
            return "FRAME_CENTER";
        case NOMAD::DisplayStatsType::DS_DIRECTION:
            return "DIRECTION";
        case NOMAD::DisplayStatsType::DS_DELTA_F:
            return "DELTA_F";
        case NOMAD::DisplayStatsType::DS_LAP:
            return "LAP";
        case NOMAD::DisplayStatsType::DS_MODEL_EVAL:
            return "MODEL_EVAL";
        case NOMAD::DisplayStatsType::DS_SOL:
            return "SOL";
        case NOMAD::DisplayStatsType::DS_THREAD_ALGO:
            return "THREAD_ALGO";
        case NOMAD::DisplayStatsType::DS_THREAD_NUM:
            return "THREAD_NUM";
        case NOMAD::DisplayStatsType::DS_GEN_STEP:
            return "GEN_STEP";
        case NOMAD::DisplayStatsType::DS_SUCCESS_TYPE:
            return "SUCCESS_TYPE";
        case NOMAD::DisplayStatsType::DS_SURROGATE_EVAL:
            return "SURROGATE_EVAL";
        case NOMAD::DisplayStatsType::DS_TOTAL_MODEL_EVAL:
            return "TOTAL_MODEL_EVAL";
        case NOMAD::DisplayStatsType::DS_USER:
            return "USER";
        case NOMAD::DisplayStatsType::DS_UNDEFINED:
        default:
            return "UNDEFINED";
    }
}


std::string NOMAD::StatsInfo::display(const NOMAD::DisplayStatsTypeList& format,
                                      const NOMAD::ArrayOfDouble & solFormat,
                                      const size_t objWidth,
                                      const size_t hWidth,
                                      const bool starSuccess,
                                      const bool appendComment) const
{
    std::string out;

    // Special case: an empty string will display BBE and OBJ.
    if (format.empty())
    {
        return display(NOMAD::DisplayStatsTypeList("BBE OBJ"), solFormat, objWidth);
    }

    // Compute precision. Precision of OBJ is the maximum precision of solFormat.
    int objPrec = 0;    // no decimals
    int hPrec = -1;     // free format
    if (solFormat.isDefined())
    {
        // We are ready to live with a wrong cast
        objPrec = static_cast<int>(solFormat.max().todouble());
        if (hWidth <= objWidth)
        {
            hPrec = objPrec;
        }
    }

    // Display selected types.
    NOMAD::DisplayStatsType statsType;
    for (size_t i = 0; i < format.size(); i++)
    {
        if (i > 0)
        {
            // Add a tab before pStart "(" or after pEnd ")".
            bool addTab = ( format[i-1] == NOMAD::ArrayOfDouble::pEnd
                           || format[i] == NOMAD::ArrayOfDouble::pStart );
            out += (addTab) ? "\t" : " ";
        }

        std::string doubleFormat;
        statsType = stringToDisplayStatsType(format[i], doubleFormat);

        if (NOMAD::DisplayStatsType::DS_OBJ == statsType)
        {
            if (doubleFormat.empty())
            {
                // Width is objWidth.
                out += _obj.display(objPrec, objWidth);
            }
            else
            {
                // doubleFormat overrides formatting
                out += _obj.display(doubleFormat);
            }
        }
        else if (NOMAD::DisplayStatsType::DS_CONS_H == statsType)
        {
            if (doubleFormat.empty())
            {
                // Width is hWidth. Use hPrec for precision.
                out += _consH.display(hPrec, hWidth);
            }
            else
            {
                // doubleFormat overrides formatting
                out += _consH.display(doubleFormat);
            }
        }
        else if (NOMAD::DisplayStatsType::DS_H_MAX == statsType)
        {
            if (doubleFormat.empty())
            {
                // Width is hWidth. Use hPrec for precision.
                out += _hMax.display(hPrec, hWidth);
            }
            else
            {
                // doubleFormat overrides formatting
                out += _hMax.display(doubleFormat);
            }
        }
        else if (NOMAD::DisplayStatsType::DS_BBE == statsType)
        {
            out += NOMAD::itos(_bbe);
        }
        else if (NOMAD::DisplayStatsType::DS_FEAS_BBE == statsType)
        {
            out += NOMAD::itos(_feasBBE);
        }
        else if (NOMAD::DisplayStatsType::DS_INF_BBE == statsType)
        {
            out += NOMAD::itos(_infBBE);
        }
        else if (NOMAD::DisplayStatsType::DS_REL_SUCC == statsType)
        {
            out += NOMAD::itos(_nbRelativeSuccess);
        }
        else if (NOMAD::DisplayStatsType::DS_PHASE_ONE_SUCC == statsType)
        {
            out += NOMAD::itos(_PhaseOneSuccess);
        }
        else if (NOMAD::DisplayStatsType::DS_ALGO_BBE == statsType)
        {
            out += NOMAD::itos(_algoBBE);
        }
        else if (NOMAD::DisplayStatsType::DS_BLK_EVA == statsType)
        {
            out += NOMAD::itos(_blkEva);
        }
        else if (NOMAD::DisplayStatsType::DS_BLK_SIZE == statsType)
        {
            out += NOMAD::itos(_blkSize);
        }
        else if (NOMAD::DisplayStatsType::DS_BBO == statsType)
        {
            out += _bbo;
        }
        else if (NOMAD::DisplayStatsType::DS_EVAL == statsType)
        {
            out += NOMAD::itos(_eval);
        }
        else if (NOMAD::DisplayStatsType::DS_CACHE_HITS == statsType)
        {
            out += NOMAD::itos(_cacheHits);
        }
        else if (NOMAD::DisplayStatsType::DS_CACHE_SIZE == statsType)
        {
            out += NOMAD::itos(_cacheSize);
        }
        else if (NOMAD::DisplayStatsType::DS_ITER_NUM == statsType)
        {
            out += NOMAD::itos(_iterNum);
        }
        else if (NOMAD::DisplayStatsType::DS_TIME == statsType)
        {
            out += NOMAD::itos(_time);
        }
        else if (NOMAD::DisplayStatsType::DS_MESH_INDEX == statsType)
        {
            out += _meshIndex.display(solFormat);
        }
        else if (NOMAD::DisplayStatsType::DS_MESH_SIZE == statsType
                 || NOMAD::DisplayStatsType::DS_DELTA_M == statsType)
        {
            out += _meshSize.display(solFormat);
        }
        else if (NOMAD::DisplayStatsType::DS_FRAME_SIZE == statsType
                 || NOMAD::DisplayStatsType::DS_DELTA_F == statsType)
        {
            out += _frameSize.display(solFormat);
        }
        else if (NOMAD::DisplayStatsType::DS_FRAME_CENTER == statsType)
        {
            out += _frameCenter.display(solFormat);
        }
        else if (NOMAD::DisplayStatsType::DS_DIRECTION == statsType)
        {
            out += _direction.display(solFormat);
        }
        else if (NOMAD::DisplayStatsType::DS_LAP == statsType)
        {
            out += NOMAD::itos(_lap);
        }
        else if (NOMAD::DisplayStatsType::DS_MODEL_EVAL == statsType)
        {
            out += NOMAD::itos(_modelEval);
        }
        else if (NOMAD::DisplayStatsType::DS_SOL == statsType)
        {
            // Here, use displayNoPar() to have the same output as NOMAD 3
            // (no additional parenthesis).
            out += _sol.displayNoPar(solFormat);
        }
        else if (NOMAD::DisplayStatsType::DS_SURROGATE_EVAL == statsType)
        {
            out += NOMAD::itos(_surrogateEval);
        }
        else if (NOMAD::DisplayStatsType::DS_THREAD_ALGO == statsType)
        {
            out += NOMAD::itos(_threadAlgoNum);
        }
        else if (NOMAD::DisplayStatsType::DS_THREAD_NUM == statsType)
        {
            out += NOMAD::itos(_threadNum);
        }
        else if (NOMAD::DisplayStatsType::DS_GEN_STEP == statsType)
        {
            out += _genStep;
        }
        else if (NOMAD::DisplayStatsType::DS_SUCCESS_TYPE == statsType)
        {
            out += NOMAD::enumStr(_success);
        }
        else if (NOMAD::DisplayStatsType::DS_TOTAL_MODEL_EVAL == statsType)
        {
            out += NOMAD::itos(_totalModelEval);
        }
        else if (NOMAD::DisplayStatsType::DS_USER == statsType)
        {
            // Display original string
            out += (format[i]);
        }
        else if (NOMAD::DisplayStatsType::DS_UNDEFINED == statsType)
        {
            // Undefined type
            out += "UNDEFINED";
        }
        else
        {
            // Exception here - we should not be there,
            // unless we added a new DisplayStatsType and forgot
            // to update this function.
            throw NOMAD::Exception(__FILE__, __LINE__, "Unrecognized stats type");
        }

    }

    if (appendComment && !_comment.empty())
    {
        // Append comment if it is non-empty.
        out += " " + _comment;
    }

    if (starSuccess)
    {
        // Add an '*' if this evaluation was a success relative to the previous
        // evaluation, or relative to the mesh center if there was no previous
        // evaluation in the same pass.
        if (_relativeSuccess)
        {
            out += " *";
        }
    }

    return out;
}


std::string NOMAD::StatsInfo::displayHeader(const NOMAD::DisplayStatsTypeList& header,
                                            const NOMAD::ArrayOfDouble & solFormat,
                                            const size_t objWidth) const
{
    std::string out;

    // Do not display formatting. Separate it from the DisplayStatsType.
    NOMAD::DisplayStatsTypeList realHeader;
    for (size_t i = 0; i < header.size(); i++)
    {
        std::string format, dst;
        NOMAD::separateFormat(header[i], format, dst);
        realHeader.add(dst);
    }

    // Display header for selected types.
    // Keep it simple for now. Ideally, the headers would align with the fields.
    out += realHeader.display();

    return out;
}



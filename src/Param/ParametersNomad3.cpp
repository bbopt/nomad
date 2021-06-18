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

#include "../Param/AllParameters.hpp"

// Implementation of set and get functions for backwards compatibility with NOMAD 3.


// Algorithm and miscellaneous parameters
std::string NOMAD::AllParameters::get_problem_dir() const
{
    return getAttributeValue<std::string>("PROBLEM_DIR");
}


void NOMAD::AllParameters::set_SEED ( int t )
{
    setAttributeValue("SEED", t);
}


int NOMAD::AllParameters::get_max_bb_eval() const
{
    auto max_bb_eval = getAttributeValue<size_t>("MAX_BB_EVAL");
    int n3_max_bb_eval =-1;
    if ( max_bb_eval < NOMAD::P_INF_INT )
        n3_max_bb_eval = static_cast<int>(max_bb_eval);
    return n3_max_bb_eval;
}


void NOMAD::AllParameters::set_MAX_BB_EVAL(int bbe)
{
    // -1 is the old way for INF
    if ( bbe >= NOMAD::P_INF_INT || bbe == -1 )
        setAttributeValue("MAX_BB_EVAL", NOMAD::INF_SIZE_T );
    else
        setAttributeValue("MAX_BB_EVAL", static_cast<size_t>(bbe) );
}


void NOMAD::AllParameters::set_MAX_EVAL(int bbe)
{
    // -1 is the old way for INF
    if ( bbe >= NOMAD::P_INF_INT || bbe == -1 )
        setAttributeValue("MAX_EVAL", NOMAD::INF_SIZE_T );
    else
        setAttributeValue("MAX_EVAL", static_cast<size_t>(bbe) );
}


int NOMAD::AllParameters::get_max_iterations() const
{
    auto max_iterations = getAttributeValue<size_t>("MAX_ITERATIONS");
    int n3_max_iterations =-1;
    if ( max_iterations < NOMAD::P_INF_INT )
        n3_max_iterations = static_cast<int>(max_iterations);
    return n3_max_iterations;

}


void NOMAD::AllParameters::set_MAX_ITERATIONS(int max_iterations)
{
    // -1 is the old way for INF
    if ( max_iterations >= NOMAD::P_INF_INT || max_iterations == -1 )
        setAttributeValue("MAX_ITERATIONS", NOMAD::INF_SIZE_T );
    else
        setAttributeValue("MAX_ITERATIONS", static_cast<size_t>(max_iterations) );
}

void NOMAD::AllParameters::set_EPSILON(const NOMAD::Double &epsilon)
{
    setAttributeValue("EPSILON", epsilon);
}

const NOMAD::Double NOMAD::AllParameters::get_epsilon() const
{
    return getAttributeValue<NOMAD::Double>("EPSILON");
}

void NOMAD::AllParameters::set_UNDEF_STR(const std::string &undefStr)
{
    setAttributeValue("UNDEF_STR", undefStr);
}

const std::string NOMAD::AllParameters::get_undef_str() const
{
    return getAttributeValue<std::string>("UNDEF_STR");
}

void NOMAD::AllParameters::set_INF_STR(const std::string &infStr)
{
    setAttributeValue("INF_STR", infStr);
}

const std::string NOMAD::AllParameters::get_inf_str() const
{
    return getAttributeValue<std::string>("INF_STR");
}



// Display parameters
// ------------------

bool NOMAD::AllParameters::set_DISPLAY_DEGREE(const int displayDegree)
{
    setAttributeValue("DISPLAY_DEGREE", displayDegree);
    return true;
}

int NOMAD::AllParameters::get_display_degree() const
{
    return getAttributeValue<int>("DISPLAY_DEGREE");
}

void NOMAD::AllParameters::set_DISPLAY_ALL_EVAL(const bool displayAllEval)
{
    setAttributeValue("DISPLAY_ALL_EVAL", displayAllEval);
}

bool NOMAD::AllParameters::get_display_all_eval() const
{
    return getAttributeValue<bool>("DISPLAY_ALL_EVAL");
}

// Mesh
// -----
const NOMAD::ArrayOfDouble& NOMAD::AllParameters::get_initial_mesh_size() const
{
    return getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_MESH_SIZE");
}


const NOMAD::ArrayOfDouble& NOMAD::AllParameters::get_initial_poll_size() const
{
// The ***_POLL_SIZE parameters have been renamed ***_FRAME_SIZE in Nomad 4
    return getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE");
}


const NOMAD::ArrayOfDouble& NOMAD::AllParameters::get_min_mesh_size() const
{
    return getAttributeValue<NOMAD::ArrayOfDouble>("MIN_MESH_SIZE");
}


const NOMAD::ArrayOfDouble& NOMAD::AllParameters::get_min_poll_size() const
{
// The ***_POLL_SIZE parameters have been renamed ***_FRAME_SIZE in Nomad 4
    return getAttributeValue<NOMAD::ArrayOfDouble>("MIN_FRAME_SIZE");
}



void NOMAD::AllParameters::set_MIN_MESH_SIZE(const NOMAD::ArrayOfDouble& delta_p_min)
{
    setAttributeValue("MIN_MESH_SIZE", delta_p_min);
}

void NOMAD::AllParameters::set_MIN_POLL_SIZE(const NOMAD::ArrayOfDouble &delta_p_min)
{
// The ***_POLL_SIZE parameters have been renamed ***_FRAME_SIZE in Nomad 4
    setAttributeValue("MIN_FRAME_SIZE", delta_p_min);
}

void NOMAD::AllParameters::set_INITIAL_MESH_SIZE(const NOMAD::ArrayOfDouble & delta_m_0)
{
    setAttributeValue("INITIAL_MESH_SIZE", delta_m_0);
}

void NOMAD::AllParameters::set_INITIAL_POLL_SIZE(const NOMAD::ArrayOfDouble & delta_m_0)
{
    // The ***_POLL_SIZE parameters have been renamed ***_FRAME_SIZE in Nomad 4
    setAttributeValue("INITIAL_FRAME_SIZE", delta_m_0);
}


void NOMAD::AllParameters::set_X0(const NOMAD::Point &x0)
{
    setAttributeValue("X0", x0);
}


const NOMAD::Point & NOMAD::AllParameters::get_x0() const
{
    return getAttributeValue<NOMAD::Point>("X0");
}


const NOMAD::ArrayOfPoint & NOMAD::AllParameters::get_x0s() const
{
    return getAttributeValue<NOMAD::ArrayOfPoint>("X0");
}


int NOMAD::AllParameters::get_dimension() const
{
    return static_cast<int>(getAttributeValue<size_t>("DIMENSION"));
}


void NOMAD::AllParameters::set_DIMENSION(size_t n)
{
    setAttributeValue("DIMENSION", n);
}


const NOMAD::ArrayOfDouble& NOMAD::AllParameters::get_lb() const
{
    return getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
}


const NOMAD::ArrayOfDouble & NOMAD::AllParameters::get_ub() const
{
    return getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
}


/// Reset the bounds.
void NOMAD::AllParameters::reset_bounds()
{
    getPbParams()->resetToDefaultValue("LOWER_BOUND");
    getPbParams()->resetToDefaultValue("UPPER_BOUND");
}


/// Set all lower bounds.
void NOMAD::AllParameters::set_LOWER_BOUND(const NOMAD::ArrayOfDouble& lb)
{
    setAttributeValue("LOWER_BOUND", lb);
}


/// Set all upper bounds.
void NOMAD::AllParameters::set_UPPER_BOUND(const NOMAD::ArrayOfDouble& ub)
{
    setAttributeValue("UPPER_BOUND", ub);
}


const NOMAD::ArrayOfDouble& NOMAD::AllParameters::get_granularity() const
{
    return getAttributeValue<NOMAD::ArrayOfDouble>("GRANULARITY");
}


void NOMAD::AllParameters::set_GRANULARITY(const NOMAD::ArrayOfDouble &granularity)
{
    setAttributeValue("GRANULARITY", granularity);
}


void NOMAD::AllParameters::set_TMP_DIR(const std::string &tmpdir)
{
    setAttributeValue("TMP_DIR", tmpdir);
}


std::string NOMAD::AllParameters::get_tmp_dir() const
{
    return getAttributeValue<std::string>("TMP_DIR");
}


// Get/Set blackbox executable
void NOMAD::AllParameters::set_BB_EXE(const std::string &bbexe)
{
    setAttributeValue("BB_EXE", bbexe);
}


std::string NOMAD::AllParameters::get_bb_exe() const
{
    return getAttributeValue<std::string>("BB_EXE");
}

void NOMAD::AllParameters::set_BB_INPUT_TYPE(const NOMAD::BBInputTypeList &bbInputType)
{
    setAttributeValue("BB_INPUT_TYPE", bbInputType);
}


const NOMAD::BBInputTypeList& NOMAD::AllParameters::get_bb_input_type() const
{
    return getAttributeValue<NOMAD::BBInputTypeList>("BB_INPUT_TYPE");
}


void NOMAD::AllParameters::set_BB_OUTPUT_TYPE(const NOMAD::BBOutputTypeList &bbOutputType)
{
    setAttributeValue("BB_OUTPUT_TYPE", bbOutputType);
}


const NOMAD::BBOutputTypeList& NOMAD::AllParameters::get_bb_output_type() const
{
    return getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
}


void NOMAD::AllParameters::set_DISPLAY_STATS(const NOMAD::ArrayOfDouble stats)
{
    setAttributeValue("DISPLAY_STATS", stats);
}


NOMAD::ArrayOfDouble NOMAD::AllParameters::get_display_stats() const
{
    return getAttributeValue<NOMAD::ArrayOfDouble>("DISPLAY_STATS");
}


void NOMAD::AllParameters::set_STATS_FILE(const NOMAD::ArrayOfDouble stats)
{
    setAttributeValue("STATS_FILE", stats);
}


NOMAD::ArrayOfDouble NOMAD::AllParameters::get_stats_file() const
{
    return getAttributeValue<NOMAD::ArrayOfDouble>("STATS_FILE");
}


void NOMAD::AllParameters::resetStatsFile()
{
    getDispParams()->resetToDefaultValue("STATS_FILE");
}


bool NOMAD::AllParameters::get_add_seed_to_file_names() const
{
    return getAttributeValue<bool>("ADD_SEED_TO_FILE_NAMES");
}


void NOMAD::AllParameters::set_ADD_SEED_TO_FILE_NAMES(bool addseed)
{
    setAttributeValue("ADD_SEED_TO_FILE_NAMES", addseed);
}


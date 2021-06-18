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
#ifndef __NOMAD_4_0_PARAMETERSNOMAD3__
#define __NOMAD_4_0_PARAMETERSNOMAD3__

// File to support backwards compatibility for older Parameters functions.
// Ex. set_DISPLAY_DEGREE(), get_display_degree().
// This file is included inside class AllParameters.

    /*--------------------------------------------------------------*/
    /*           Attributes and methods listed by categories        */
    /*--------------------------------------------------------------*/

    // Algorithm and miscellaneous parameters:
    // ---------------------------------------
public:
    /// Access to the seed (compatibility NOMAD 3).
    int get_seed() const { return getAttributeValue<int>("SEED"); }

    /// Access to the \c MAX_BB_EVAL parameter (compatibility NOMAD 3).
    int get_max_bb_eval() const;

    /// Access to the maximum number of iterations (compatibility NOMAD 3).
    int get_max_iterations() const;

    // Return the directory where the problem resides (compatibility NOMAD 3).
    /**
     \return The path to the problem, when in batch mode. An empty string, when in library mode.
    */
    std::string get_problem_dir() const;

    /// Set the seed (compatibility NOMAD 3).
    void set_SEED(int seed);

    /// Set the \c MAX_BB_EVAL parameter (compatibility NOMAD 3).
    void set_MAX_BB_EVAL(int bbe);

    /// Set the \c MAX_EVAL parameter (compatibility NOMAD 3).
    void set_MAX_EVAL(int maxEval);

    /// Set the \c MAX_ITERATIONS parameter (compatibility NOMAD 3).
    void set_MAX_ITERATIONS(int max_iterations);


    /// Set the \c EPSILON parameter (compatibility NOMAD 3).
    void set_EPSILON(const Double & eps);

    /// Access to the \c EPSILON parameter (compatibility NOMAD 3).
    const Double get_epsilon() const;

    /// Set the \c UNDEF_STR parameter (compatibility NOMAD 3).
    void set_UNDEF_STR(const std::string &undefStr);

    /// Set the \c INF_STR parameter (compatibility NOMAD 3).
    void set_INF_STR(const std::string &infStr);


    /// Access to the \c UNDEF_STR parameter (compatibility NOMAD 3).
    const std::string get_undef_str() const;

    /// Access to the \c INF_STR parameter (compatibility NOMAD 3).
    const std::string get_inf_str() const;


    // Mesh:
    // -----
public:
    /// Access to the \c INITIAL_MESH_SIZE parameter (compatibility NOMAD 3).
    const ArrayOfDouble& get_initial_mesh_size() const;

    /// Access to the \c INITIAL_FRAME_SIZE parameter (compatibility NOMAD 3).
    /**
     * The POLL_SIZE parameters have been renamed
     FRAME_SIZE in Nomad 4. For compatibility with
     Nomad 3, the old function name is used.

     \return The \c INITIAL_FRAME_SIZE parameter.
     */
    const ArrayOfDouble& get_initial_poll_size() const;

    /// Access to the \c MIN_MESH_SIZE parameter (compatibility NOMAD 3).
    const ArrayOfDouble& get_min_mesh_size() const;

    /// Access to the \c MIN_FRAME_SIZE parameter (compatibility NOMAD 3).
    /**
     * The ***_POLL_SIZE parameters have been renamed
     ***_FRAME_SIZE in Nomad 4. For compatibility with
     Nomad 3, the old function name is used.

     \return The \c MIN_FRAME_SIZE parameter.
     */
    const ArrayOfDouble & get_min_poll_size() const;

    /// Set the \c MIN_MESH_SIZE parameter (compatibility NOMAD 3).
    void set_MIN_MESH_SIZE(const ArrayOfDouble &mms);

    /// Set the \c MIN_FRAME_SIZE parameter (compatibility NOMAD 3).
    /**
     * The POLL_SIZE parameters have been renamed
     FRAME_SIZE in Nomad 4. For compatibility with
     Nomad 3, the old function name is used.

     */
    void set_MIN_POLL_SIZE(const ArrayOfDouble &mfs);

    /// Set the \c INITIAL_MESH_SIZE parameter (compatibility NOMAD 3).
    void set_INITIAL_MESH_SIZE(const ArrayOfDouble &ims);

    /// Set the \c INITIAL_FRAME_SIZE parameter (compatibility NOMAD 3).
    /**
     * The POLL_SIZE parameters have been renamed
     FRAME_SIZE in Nomad 4. For compatibility with
     Nomad 3, the old function name is used.
     */
    void set_INITIAL_POLL_SIZE(const ArrayOfDouble &ifs);


    // Starting point:
    // ------------------
public:

    /// Add a starting point (compatibility NOMAD 3).
    void set_X0(const Point & x0);

    /// Access to the starting point (compatibility NOMAD 3).
    const Point & get_x0() const;

    /// Access to the starting points (compatibility NOMAD 3).
    const ArrayOfPoint & get_x0s() const;


    // Dimension:
    // ----------
    /// Access to the dimension (compatibility NOMAD 3).
    int get_dimension() const;

    /// Set the dimension (compatibility NOMAD 3).
    virtual void set_DIMENSION (size_t  n);


    // Bounds and scaling:
    // -------------------
public:

    /// Access to the lower bounds (compatibility NOMAD 3).
    const ArrayOfDouble & get_lb() const;

    /// Access to the upper bounds (compatibility NOMAD 3).
    const ArrayOfDouble & get_ub() const;

    /// Reset the bounds (compatibility NOMAD 3).
    void reset_bounds();

    /// Set all lower bounds (compatibility NOMAD 3).
    /**
     Use undefined values for variables without bounds.
     */
    void set_LOWER_BOUND(const ArrayOfDouble& lb);

    /// Set all upper bounds (compatibility NOMAD 3).
    /**
     Use undefined values for variables without bounds.
     */
    void set_UPPER_BOUND(const ArrayOfDouble& ub);


    // Granular variables:
    // -------------------

    /// Access to the granular variables (compatibility NOMAD 3).
    const ArrayOfDouble &get_granularity() const;

    /// Set the granularity of a series of variables (compatibility NOMAD 3).
    /**
     \param granularity  The granular variables; This point is of dimension \c n;
     regular variables have a granularity of zero -- \b IN.
     */
    void set_GRANULARITY(const ArrayOfDouble &granularity);


    // Blackbox parameters for nomad batch
    // -----------------------------------
    /// Set temporary directory for blackbox execution (compatibility NOMAD 3)
    void set_TMP_DIR(const std::string &tmpdir);

    /// Access to tmp_dir (compatibility NOMAD 3)
    std::string get_tmp_dir() const;

    /// Set blackbox executable (compatibility NOMAD 3)
    void set_BB_EXE(const std::string &bbexe);

    /// Access to blackbox executable (compatibility NOMAD 3)
    std::string get_bb_exe() const;

    /// Set blackbox input type (compatibility NOMAD 3)
    void set_BB_INPUT_TYPE(const BBInputTypeList &bbInputType);

    /// Access to blackbox input type (compatibility NOMAD 3)
    const std::vector<BBInputType>& get_bb_input_type() const;

    /// Set blackbox output type (compatibility NOMAD 3)
    void set_BB_OUTPUT_TYPE(const BBOutputTypeList &bbOutputType);

    /// Access to blackbox output type (compatibility NOMAD 3)
    const std::vector<BBOutputType>& get_bb_output_type() const;

    // Display parameters
    // ------------------

    /// Set display degree (compatibility NOMAD 3)
    bool set_DISPLAY_DEGREE(const int displayDegree) ;// override;

    /// Access to display degree (compatibility NOMAD 3)
    int get_display_degree() const;

    /// Set DISPLAY_ALL_EVAL (compatibility NOMAD 3)
    void set_DISPLAY_ALL_EVAL(const bool displayAllEval);

    /// Get the flag DISPLAY_ALL_EVAL (compatibility NOMAD 3)
    bool get_display_all_eval() const;

    /// set DISPLAY_STATS (compatibility NOMAD 3)
    void set_DISPLAY_STATS(const ArrayOfDouble stats);

    /// Get DISPLAY_STATS (compatibility NOMAD 3)
    ArrayOfDouble get_display_stats() const;

    /// Reset STATS_FILE (compatibility NOMAD 3)
    void resetStatsFile();

    /// Set STATS_FILE (compatibility NOMAD 3)
    void set_STATS_FILE(const ArrayOfDouble stats);

    /// Get STATS_FILE (compatibility NOMAD 3)
    ArrayOfDouble get_stats_file() const;

    /// Set ADD_SEED_TO_FILE_NAMES (compatibility NOMAD 3)
    void set_ADD_SEED_TO_FILE_NAMES(bool addseed);

    /// Get ADD_SEED_TO_FILE_NAMES (compatibility NOMAD 3)
    bool get_add_seed_to_file_names() const;

#endif // __NOMAD_4_0_PARAMETERSNOMAD3__

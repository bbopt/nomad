/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
//////////// THIS FILE MUST BE CREATED BY EXECUTING WriteAttributeDefinitionFile ////////////
//////////// DO NOT MODIFY THIS FILE MANUALLY ///////////////////////////////////////////////

#ifndef __NOMAD_4_5_EVALATTRIBUTESDEFINITION__
#define __NOMAD_4_5_EVALATTRIBUTESDEFINITION__

_definition = {
{ "BB_EVAL_FORMAT",  "NOMAD::ArrayOfDouble",  "-",  " Format of the doubles sent to the blackbox evaluator ",  " \n  \n . BB_EVAL_FORMAT is computed from the BB_INPUT_TYPE parameter. \n  \n . Gives the format precision for doubles sent to blackbox evaluator. \n  \n . CANNOT BE MODIFIED BY USER. Internal parameter. \n  \n . No default value.\n\n",  "  internal  "  , "false" , "false" , "true" },
{ "BB_EXE",  "std::string",  "",  " Blackbox executable ",  " \n  \n . Blackbox executable name \n  \n . List of strings \n  \n . Required for batch mode \n  \n . Unused in library mode \n  \n . One executable can give several outputs \n  \n . Use \' or \", and \'$\', to specify names or commands with spaces \n  \n . When the \'$\' character is put in first position of a string, it is \n   considered as global and no path will be added \n  \n . Examples \n     . BB_EXE bb.exe \n     . BB_EXE \'$nice bb.exe\' \n     . BB_EXE \'$python bb.py\' \n  \n . Default: Empty string.\n\n",  "  basic blackbox blackboxes bb exe executable executables binary output outputs batch  "  , "false" , "false" , "true" },
{ "BB_REDIRECTION",  "bool",  "true",  " Blackbox executable redirection for outputs  ",  " \n  \n . Flag to redirect blackbox executable outputs in a stream. The redirection \n   in a stream does not require an ouptut file. NOMAD interprets the outputs from \n   the stream according to BB_OUTPUT_TYPE. If the blackbox executable \n   outputs some verbose, NOMAD cannot interpret correctly the outputs. \n  \n . If the redirection is disabled. The blackbox must output its results into a \n  file having the name of the input file (usually nomadtmp.pid.threadnum) \n  completed by \".output\". The format must follow the BB_OUTPUT_TYPE. \n  For example, for BB_OUTPUT_TYPE OBJ CSTR, we must have only the two values on \n  a single line in the output file with a end-of-line. \n   \n . Disable blackbox redirection and managing output file can be convenient when \n the blackbox outputs some verbose. All the verbose is put into a temporary \n log file nomadtmp.pid.threadnum.tmplog \n   \n . This parameter has no effect when BB_EXE is not defined like in library mode. \n  \n . Examples \n     . BB_REDIRECTION false \n  \n . Default: true\n\n",  "  basic blackbox blackboxes bb exe executable executables binary output outputs batch  "  , "false" , "false" , "true" },
{ "BB_OUTPUT_TYPE",  "NOMAD::BBOutputTypeList",  "OBJ",  " Type of outputs provided by the blackboxes ",  " \n  \n . Blackbox output types \n  \n . List of types for each blackbox output \n  \n . If BB_EXE is defined, the blackbox outputs must be returned by the executable \n on a SINGLE LINE of the standard output or in an output file \n (see BB_REDIRECTION). The order of outputs must be consistent between the blackbox \n and BB_OUTPUT_TYPE. \n  \n . Available types \n     . OBJ       : objective value to minimize (define twice for bi-objective) \n     . PB        : constraint <= 0 treated with Progressive Barrier (PB) \n     . CSTR      : same as 'PB' \n     . EB        : constraint <= 0 treated with Extreme Barrier (EB) \n     . F         : constraint <= 0 treated with Filter \n     . CNT_EVAL  : 0 or 1 output: count or not the evaluation (for batch mode and Matlab interface) \n     . NOTHING   : this output is ignored \n     . EXTRA_O   : same as 'NOTHING' \n     .  -        : same as 'NOTHING' \n     . BBO_UNDEFINED: same as 'NOTHING' \n  \n . Equality constraints are not natively supported \n  \n . Extra outputs (EXTRA_O, NOTHING, BBO_UNDEFINED, ...) are not used for \n   optimization but are available for display and custom user testing \n   (see examples). \n  \n . See parameters LOWER_BOUND and UPPER_BOUND for bound constraints \n  \n . See parameter H_NORM for the infeasibility measure computation. \n  \n . See parameter H_MIN for relaxing the feasibility criterion. \n  \n . Examples \n     . BB_EXE bb.exe                   # these two lines define \n     . BB_OUTPUT_TYPE OBJ EB EB        # that bb.exe outputs three values \n  \n . Default: OBJ\n\n",  "  basic bb exe blackbox blackboxs output outputs constraint constraints type types infeasibility norm  "  , "false" , "false" , "true" },
{ "SURROGATE_EXE",  "std::string",  "",  " Static surrogate executable ",  " \n . To indicate a static surrogate executable \n  \n . List of strings \n  \n . Surrogate executable must have the same number of outputs as blackbox  \n     executable, defined by BB_OUTPUT_TYPE. \n      \n . Static surrogate evaluations can be used for sorting trial points before \n   blackbox evaluation OR for VNS Search. \n  \n . Example \n     SURROGATE_EXE surrogate.exe     # surrogate.exe is a static surrogate executable \n                                     # for BB_EXE \n . Default: Empty string.\n\n",  "  advanced static surrogate executable  "  , "true" , "false" , "true" } };

#endif

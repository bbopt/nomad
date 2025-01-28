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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONIBEX__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONIBEX__

_definition = {
{ "USE_IBEX",  "bool",  "false",  " Boolean to determine if we want to use the functionnalities of IBEX ",  " \n  \n . Argument : bool \n  \n . Determine if you want to use the fonctionnalities of IBEX \n  \n . Default: false\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" },
{ "SYSTEM_FILE_NAME",  "string",  "-",  " File with the constraints  ",  " \n  \n . Minibex file name, describing the system (i.e constraints, variables...) of the problem. \n  \n . See the documentation here for more detail : http://www.ibex-lib.org/doc/minibex.html on how to create it. \n  \n . No default value.\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" },
{ "SET_FILE",  "bool",  "false",  " Boolean to determine if the file of the set is already created ",  " \n  \n . Argument : bool \n  \n . Determine if the Set of the problem is already created. \n  \n . Default: false\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" },
{ "SET_FILE_NAME",  "string",  "-",  " File to load with the set  ",  " \n  \n . Argument : string \n  \n . Name of the Set file. \n  \n . No need to be provided if SET_FILE = false. \n  \n . No default value.\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" } };

#endif

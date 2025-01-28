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

#ifndef __NOMAD_4_5_RUNATTRIBUTESDEFINITIONDISCO__
#define __NOMAD_4_5_RUNATTRIBUTESDEFINITIONDISCO__

_definition = {
{ "DISCO_MADS_OPTIMIZATION",  "bool",  "false",  " DiscoMads optimization ",  " \n  \n . DiscoMads optimization: reveal and escape either weak discontinuities \n in revealing outputs or hidden constraints regions (with DISCO_MADS_HID_CONST). \n  \n . If DiscoMads is used for discontinuities, they are characterized by  \n DISCO_MADS_DETECTION_RADIUS and DISCO_MADS_LIMIT_RATE.  Revealing outputs \n  should be indicated by appending -R to blackbox outputs type in parameters file. \n  \n . Prescribed distance to revealing points is controlled by  \n DISCO_MADS_EXCLUSION_RADIUS. \n  \n . Revealing poll parameters are DISCO_MADS_REVEALING_POLL_NB_POINTS and \n DISCO_MADS_REVEALING_POLL_RADIUS. \n  \n . Argument: bool \n  \n . Example: DISCO_MADS_OPTIMIZATION true \n  \n . Default: false\n\n",  "  discomads mads discontinuity "  , "true" , "false" , "true" },
{ "DISCO_MADS_DETECTION_RADIUS",  "NOMAD::Double",  "1.0",  " Radius used to reveal discontinuities in DiscoMads ",  " \n  \n . DiscoMADS algorithm attempts to reveal discontinuities between couple of \n points whose distance is strictly less than this value. \n  \n . Argument: greater than or equal to 0. If 0, the DiscoMads algorithm is  \n equivalent to the Mads algorithm with a special search (the revealing poll). \n  \n . This attribute is used only when DiscoMads optimization is active and  \n DiscoMads is used for discontinuity revealation (DISCO_MADS_HID_CONST false) \n  \n . Example: DISCO_MADS_DETECTION_RADIUS 1.0 \n  \n . This radius is called rd in the article. \n . Default: 1.0\n\n",  "  advanced discomads radius detection discontinuity "  , "true" , "false" , "true" },
{ "DISCO_MADS_LIMIT_RATE",  "NOMAD::Double",  "1",  " Limit rate of change used to reveal discontinuities in DiscoMads ",  " \n  \n . When using DiscoMADS algorithm, a discontinuity is revealed if between two \n  points at distance at most DISCO_MADS_DETECTION_RADIUS the rate of change of  \n at least one revealing output exceeds DISCO_MADS_LIMIT_RATE.  \n  \n . Argument: stricly greater than 0 \n  \n . This attribute is used only when DiscoMads optimization is active and  \n DiscoMads is used for discontinuity revealation (DISCO_MADS_HID_CONST false) \n  \n . Example: DISCO_MADS_LIMIT_RATE 1.0 \n  \n . This quantity is called tau in the article. \n . Default: 1\n\n",  "  advanced discomads rate  "  , "true" , "false" , "true" },
{ "DISCO_MADS_EXCLUSION_RADIUS",  "NOMAD::Double",  "1",  " Radius of exclusion balls around revealing points in DiscoMads",  " \n  \n . When using DiscoMADS algorithm, the points at distance at most  \n  DISCO_MADS_EXCLUSION_RADIUS of revealing points are considered to be close to \n  revealing points and thus penalized. \n  \n . Argument: stricly greater than 0 \n  \n . This attribute is used only when DiscoMads optimization is active. \n  \n . Example: DISCO_MADS_EXCLUSION_RADIUS 1.0 \n  \n . This radius is called re in the article. \n . Default: 1\n\n",  "  advanced discomads radius exclusion  "  , "true" , "false" , "true" },
{ "DISCO_MADS_REVEALING_POLL_RADIUS",  "NOMAD::Double",  "2.02",  " Revealing poll radius in DiscoMads ",  " \n  \n .  When using DiscoMADS algorithm, the revealing poll evaluates random points \n in the ball of this radius around the best feasible incumbent if it exists,  \n otherwise around the best infeasible incumbent. \n  \n . Argument: Double stricly greater than exclusion radius + detection radius, \n   e.g. 1.01*(exclusion radius + detection radius). \n  \n . This attribute is used only when DiscoMads optimization is active. \n  \n . Example: DISCO_MADS_REVEALING_POLL_RADIUS 2.5 \n  \n . This radius is called rm in the article. \n . Default: 2.02\n\n",  "  advanced discomads radius poll revealing "  , "true" , "false" , "true" },
{ "DISCO_MADS_REVEALING_POLL_NB_POINTS",  "size_t",  "1",  " Number of random points sampled by the revealing poll in DiscoMads ",  " \n  \n .  When using DiscoMADS algorithm, the revealing poll evaluates this number \n of random points in a ball around the best feasible or infeasible incumbent. \n  \n . Argument: one positive integer. A common choice is to take the problem \n dimension. \n  \n . If null, the revealing poll is disabled. This should only be used for testing \n as a strictly positive value is requiered for the convergence analysis. \n  \n . This attribute is used only when DiscoMads optimization is active. \n  \n . Example: DISCO_MADS_REVEALING_POLL_NB_POINTS 1 \n  \n . Default: 1\n\n",  "  advanced discomads radius poll revealing "  , "true" , "false" , "true" },
{ "DISCO_MADS_HID_CONST",  "bool",  "false",  " Use DiscoMADS to reveal and escape hidden constraints regions ",  " \n  \n . Special use of DiscoMads to reveal hidden constraint regions instead of  \n discontinuities. Parameters DISCO_MADS_DETECTION_RADIUS and DISCO_MADS_LIMIT_RATE \n are not considered in this casem DISCO_MADS_DETECTION_RADIUS takes the value 0. \n  \n . Argument: bool \n  \n . Example: DISCO_MADS_HID_CONST true \n  \n . Default: false\n\n",  "  discomads hidden failed "  , "true" , "false" , "true" },
{ "DISCO_MADS_HID_CONST_OUTPUT_VALUE",  "NOMAD::Double",  "1E20",  " Value attributed to objective function and PB constraints for failed evaluations \n when DiscoMads is used to reveal hidden constraints regions.",  " \n  \n . This attribute is used only when DiscoMads optimization is active and DiscoMads  \n is used to reveal hidden constraints regions (DISCO_MADS_HID_CONST) \n  \n . Argument: Double strictly positive, greater than MODEL_MAX_OUTPUT and \n less than NOMAD::INF. \n  \n . Example: DISCO_MADS_HID_CONST_OUTPUT_VALUE 1e20 \n  \n . Default: 1E20\n\n",  "  discomads hidden failed "  , "true" , "false" , "true" } };

#endif

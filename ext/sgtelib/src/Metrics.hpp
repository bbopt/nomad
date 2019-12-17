/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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
/*-------------------------------------------------------------------------------------*/
/*  sgtelib - A surrogate model library for derivative-free optimization               */
/*  Version 2.0.2                                                                      */
/*                                                                                     */
/*  Copyright (C) 2012-2017  Sebastien Le Digabel - Ecole Polytechnique, Montreal      */ 
/*                           Bastien Talgorn - McGill University, Montreal             */
/*                                                                                     */
/*  Author: Bastien Talgorn                                                            */
/*  email: bastientalgorn@fastmail.com                                                 */
/*                                                                                     */
/*  This program is free software: you can redistribute it and/or modify it under the  */
/*  terms of the GNU Lesser General Public License as published by the Free Software   */
/*  Foundation, either version 3 of the License, or (at your option) any later         */
/*  version.                                                                           */
/*                                                                                     */
/*  This program is distributed in the hope that it will be useful, but WITHOUT ANY    */
/*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    */
/*  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   */
/*                                                                                     */
/*  You should have received a copy of the GNU Lesser General Public License along     */
/*  with this program. If not, see <http://www.gnu.org/licenses/>.                     */
/*                                                                                     */
/*  You can find information on sgtelib at https://github.com/bastientalgorn/sgtelib   */
/*-------------------------------------------------------------------------------------*/

#ifndef __SGTELIB_METRICS_UTILS__
#define __SGTELIB_METRICS_UTILS__

#include "Defines.hpp"
#include "Exception.hpp"
#include "Surrogate_Utils.hpp"

#include <cstring>
#include <cctype>

namespace SGTELIB {

  // Metrics
  enum metric_t {
    METRIC_EMAX,  // Max absolute error
    METRIC_EMAXCV,// Max absolute error on cross-validation value
    METRIC_RMSE,  // Root mean square error
    METRIC_ARMSE,  // Agregate Root mean square error
    METRIC_RMSECV, // Leave-one-out cross-validation
    METRIC_ARMSECV, // Agregate Leave-one-out cross-validation
    METRIC_OE,  // Order error on the training points
    METRIC_OECV,  // Order error on the cross-validation output
    METRIC_AOE,  // Agregate Order error 
    METRIC_AOECV,  // Agregate Order error on the cross-validation output
    METRIC_EFIOE,  // Order error on the cross-validation output
    METRIC_EFIOECV,  // Agregate Order error on the cross-validation output
    METRIC_LINV   // Inverse of the likelihood
  };
  const int NB_METRIC_TYPES = 11;

  DLL_API std::string       metric_type_to_str        ( const SGTELIB::metric_t );
  DLL_API SGTELIB::norm_t   metric_type_to_norm_type  ( const SGTELIB::metric_t );
  DLL_API SGTELIB::metric_t str_to_metric_type        ( const std::string & s   );

  // Info on  metric
  // Tells if a metric returns one or multiple objectives
  // (i.e. One for all the BBO OR One per BBO)
  bool one_metric_value_per_bbo ( const SGTELIB::metric_t mt );
  bool metric_uses_cv           ( const SGTELIB::metric_t mt );

  // Convert a metric to another metric that returns only 1 obj.
  SGTELIB::metric_t metric_convert_single_obj ( const SGTELIB::metric_t mt );



}
#endif

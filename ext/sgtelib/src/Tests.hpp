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

#ifndef __SGTELIB_TESTS__
#define __SGTELIB_TESTS__

#include "sgtelib.hpp"

namespace SGTELIB {

  void sand_box (void);

  // test_quick: build the surrogate and compute the metrics
  std::string test_quick (const std::string & s , const SGTELIB::Matrix & X0 );
  // test_rmsecv: build the surrogate, then build all the CV models and verif the RMSECV
  std::string test_rmsecv (const std::string & s , const SGTELIB::Matrix & X0 );
  // test_rmse: build the surrogate, then perform prediction on the TP to verif rmse
  std::string test_rmse  (const std::string & s , const SGTELIB::Matrix & X0 );
  // test_update: build a surrogate all-at-once, or point-by-point
  std::string test_update(const std::string & s , const SGTELIB::Matrix & X0 );
  // test_pxx: build a surrogate and perform prediction on XX of various sizes
  // (especially pxx > _p)
  std::string test_pxx   (const std::string & s , const SGTELIB::Matrix & X0 );
  // test_scale: build 2 surrogates with a different scale on the data. 
  std::string test_scale (const std::string & s , const SGTELIB::Matrix & X0 );
  // test_dimension: build surrogates with various sizes of p,n,m to check for dimension errors
  std::string test_dimension  (const std::string & s );
  // test_singular_data (if there are constant inputs, or outputs or Nan outputs)
  std::string test_singular_data (const std::string & s );
  std::string test_multiple_occurrences (const std::string & s );

  // test_scale: build 2 surrogates with a different scale on the data. 
  void test_many_models ( const std::string & out_file , const SGTELIB::Matrix & X0 , const SGTELIB::Matrix & Z0 );


  void test_LOWESS_times (void);

  // analyse ensembl
  void analyse_ensemble ( const std::string & s );

  SGTELIB::Matrix test_functions (const SGTELIB::Matrix & X);
  SGTELIB::Matrix test_functions_1D (const SGTELIB::Matrix & T, const int function_index);
  double          test_functions_1D (const double            t, const int function_index);
  void check_matrix_diff (const SGTELIB::Matrix * A, const SGTELIB::Matrix * B);


  void build_test_data ( const std::string & function_name , SGTELIB::Matrix & X0 , SGTELIB::Matrix & Z0 );


}

#endif

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

#include "Surrogate_PRS_EDGE.hpp"

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS_EDGE::Surrogate_PRS_EDGE ( SGTELIB::TrainingSet & trainingset,
                                                  SGTELIB::Surrogate_Parameters param) :
  SGTELIB::Surrogate_PRS ( trainingset , param ){
  #ifdef SGTELIB_DEBUG
    std::cout << "constructor PRS_EDGE\n";
  #endif
}//


/*----------------------------*/
/*          destructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS_EDGE::~Surrogate_PRS_EDGE ( void ) {

}//


/*----------------------------*/
/*          display           */
/*----------------------------*/
void SGTELIB::Surrogate_PRS_EDGE::display_private ( std::ostream & out ) const {
  out << "q: " << _q << "\n";
}//


/*--------------------------------------*/
/*               build                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_PRS_EDGE::build_private ( void ) {
  
  const int pvar = _trainingset.get_pvar(); 
  const int nvar = _trainingset.get_nvar(); 

  // Get the number of basis functions.
  _q = Surrogate_PRS::get_nb_PRS_monomes(nvar,_param.get_degree())+nvar;

  // If _q is too big or there is not enough points, then quit
  if (_q>200) return false;
  if ( (_q>pvar-1) && (_param.get_ridge()==0) ) return false;

  // Compute the exponents of the basis functions
  _M = get_PRS_monomes(nvar,_param.get_degree());

  // DESIGN MATRIX H
  _H = compute_design_matrix ( _M , get_matrix_Xs() );

  return compute_alpha();   
}//





/*-------------------------------------------------*/
/*          Compute PRS_EDGE design matrix          */
/*-------------------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate_PRS_EDGE::compute_design_matrix ( const SGTELIB::Matrix& Monomes, 
                                                                           const SGTELIB::Matrix & Xs ) {

  // Call the standard design matrix
  const SGTELIB::Matrix H_prs = SGTELIB::Surrogate_PRS::compute_design_matrix ( Monomes, Xs );
  // Add the edge basis functions
  const int p = Xs.get_nb_rows(); 
  const int n = Xs.get_nb_cols(); 
  const int nvar = _trainingset.get_nvar();

  // Add the special basis function (xi==0)
  // Loop on the input variables
  SGTELIB::Matrix H_edge ("He",p,nvar);
  int i,j,jj;
  double v;
  double xs0; // value that xs will have if x=0 (for a given input j)
  jj=0;
  for (j=0 ; j<n ; j++){
    if (_trainingset.get_X_nbdiff(j)>1){
      xs0 = _trainingset.X_scale( 0.0 , j ); 
      for (i=0 ; i<p ; i++){
        v = (double) (Xs.get(i,j)==xs0);
        H_edge.set(i,jj,v);
      }
      jj++;
    }
  }

  SGTELIB::Matrix H(H_prs);
  H.add_cols(H_edge);
  return H;

}//
















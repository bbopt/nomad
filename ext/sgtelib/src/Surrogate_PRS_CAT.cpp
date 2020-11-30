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

#include "Surrogate_PRS_CAT.hpp"

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS_CAT::Surrogate_PRS_CAT ( SGTELIB::TrainingSet & trainingset,
                                                SGTELIB::Surrogate_Parameters param) :
  SGTELIB::Surrogate_PRS ( trainingset , param ){
  #ifdef SGTELIB_DEBUG
    std::cout << "constructor PRS_CAT\n";
  #endif
}//


/*----------------------------*/
/*          destructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS_CAT::~Surrogate_PRS_CAT ( void ) {

}//

/*----------------------------*/
/*          display           */
/*----------------------------*/
void SGTELIB::Surrogate_PRS_CAT::display_private ( std::ostream & out ) const {
  out << "q: " << _q << "\n";
  out << "nb_cat: " << _nb_cat << "\n";
}//

/*--------------------------------------*/
/*             init_private             */
/*--------------------------------------*/
bool SGTELIB::Surrogate_PRS_CAT::init_private ( void ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate_PRS_CAT : init_private\n";
  #endif
  // Compute the number of categories
  _cat.clear();
  for ( int i = 0 ; i < _p ; ++i ){
    _cat.insert ( _trainingset.get_Xs(i,0) );
  }
  _nb_cat = static_cast<int> ( _cat.size() );
  return true;
}//

/*--------------------------------------*/
/*               build                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_PRS_CAT::build_private ( void ) {

  const int pvar = _trainingset.get_pvar(); 
  const int nvar = _trainingset.get_nvar(); 

  // Get the number of basis functions.
  int nb_monomes = Surrogate_PRS::get_nb_PRS_monomes(nvar-1,_param.get_degree());
  _q = nb_monomes*_nb_cat;

  // If _q is too big or there is not enough points, then quit
  if (nb_monomes>100) return false;
  if ( (_q>pvar-1) && (_param.get_ridge()==0) ) return false;

  // Compute the exponents of the basis functions (nb there is one less variable
  // as the first one is the cat)
  _M = SGTELIB::Matrix ("M",nb_monomes,1);
  _M.fill(0.0);
  _M.add_cols(get_PRS_monomes(nvar-1,_param.get_degree()));

  // DESIGN MATRIX H
  _H = compute_design_matrix ( _M , get_matrix_Xs() );

  // Compute alpha
  bool ok = compute_alpha();   
  return ok;
}//



/*-------------------------------------------------*/
/*          Compute PRS_CAT design matrix          */
/*-------------------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate_PRS_CAT::compute_design_matrix ( const SGTELIB::Matrix& Monomes, 
                                                                          const SGTELIB::Matrix & Xs ) {

  const int p = Xs.get_nb_rows(); 
  SGTELIB::Matrix H("H",p,0);
  SGTELIB::Matrix is_cat ("is_cat",p,1);
  const SGTELIB::Matrix H_prs = SGTELIB::Surrogate_PRS::compute_design_matrix ( Monomes, Xs );
  
  std::set<double>::iterator it;
  double c;
  for (it = _cat.begin(); it != _cat.end(); ++it){
    c = *it;
    for (int i=0 ; i<p ; i++){
      is_cat.set(i,0,(double)(Xs.get(i,0)==c));
    }       
    H.add_cols( Matrix::diagA_product(is_cat,H_prs) );
  }
  return H;
}//
















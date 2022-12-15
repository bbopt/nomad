/*-------------------------------------------------------------------------------------*/
/*  sgtelib - A surrogate model library for derivative-free optimization               */
/*  Version 2.0.3                                                                      */
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

#include "Surrogate_CN.hpp"

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_CN::Surrogate_CN ( SGTELIB::TrainingSet & trainingset,
                                      SGTELIB::Surrogate_Parameters param) :
  SGTELIB::Surrogate ( trainingset , param  ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "constructor CN\n";
  #endif
}//


/*----------------------------*/
/*          destructor        */
/*----------------------------*/
SGTELIB::Surrogate_CN::~Surrogate_CN ( void ) {

}//


/*--------------------------------------*/
/*              display                 */
/*--------------------------------------*/
void SGTELIB::Surrogate_CN::display_private ( std::ostream & out ) const {
  out << "(No special members)\n";
}//


/*--------------------------------------*/
/*             build_private            */
/*--------------------------------------*/
bool SGTELIB::Surrogate_CN::build_private ( void ) {
  _ready = true;
  return true;
}//

/*--------------------------------------*/
/*       predict_private (ZZs only)      */
/*--------------------------------------*/
void SGTELIB::Surrogate_CN::predict_private ( const SGTELIB::Matrix & XXs,
                                                    SGTELIB::Matrix * ZZs) {
  
  int i,imin;
  const int pxx = XXs.get_nb_rows();

  // D : distance between points of XXs and other points of the trainingset
  SGTELIB::Matrix D = _trainingset.get_distances(XXs,get_matrix_Xs(),_param.get_distance_type());

  // Data:
  const SGTELIB::Matrix & Zs = get_matrix_Zs();

  // Loop on the points of XXs
  for (i=0 ; i<pxx ; i++){
    // imin is the index of the closest neighbor of xx in Xs
    imin = D.get_min_index_row(i); 
    // Copy the output of this point
    ZZs->set_row( Zs.get_row(imin) , i);
  }

}//


// Predict only objectives (used in Surrogate Ensemble Stat)
void SGTELIB::Surrogate_CN::predict_private_objective ( const std::vector<SGTELIB::Matrix *> & XXd,
                                                                           SGTELIB::Matrix * ZZsurr_around) {
  check_ready(__FILE__,__FUNCTION__,__LINE__);

  const size_t pxx = XXd.size();
  const int nbd = XXd[0]->get_nb_rows();

  // Data:
  const SGTELIB::Matrix & Zs = get_matrix_Zs();
  SGTELIB::Matrix Zs_obj ("Zs_obj", _p, 1);
  // Get only objectives values is Zs
  for (int j=0 ; j<_m ; j++){
    if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ){
      Zs_obj = Zs.get_col(j);
      break;
    }
  }

  // Loop on all pxx points 
  for (int i=0 ; i<static_cast<int>(pxx) ; i++){
    // XXd[i] is of dimension nbd * _n
    int d,dmin;

    // D : distance between points of XXd[i] and other points of the trainingset
    SGTELIB::Matrix D = _trainingset.get_distances(*(XXd[i]), get_matrix_Xs(), _param.get_distance_type());

    // Loop on the points of XXs
    for (d=0 ; d<nbd ; d++){
      // imin is the index of the closest neighbor of xx in Xs
      dmin = D.get_min_index_row(d); 
      // Copy the output of this point
      ZZsurr_around->set( i,d, Zs.get(dmin,0) );
    }

  } // end for i

}



/*--------------------------------------*/
/*      compute cv values               */
/*--------------------------------------*/
bool SGTELIB::Surrogate_CN::compute_cv_values (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);

  if ((_Zvs) && (_Svs)) return true;


  // Init matrices
  if ( ! _Zvs){
    _Zvs = new SGTELIB::Matrix ("Zvs",_p,_m);
    _Zvs->set_name("Zvs");
  }
    
  if ( ! _Svs){
    _Svs = new SGTELIB::Matrix ("Svs",_p,_m);
    _Svs->set_name("Svs");
  }


  int i,i2,imin=0;
  double d;
  SGTELIB::Matrix D = _trainingset.get_distances(get_matrix_Xs(),get_matrix_Xs(),_param.get_distance_type());
  const SGTELIB::Matrix & Zs = get_matrix_Zs();

  // Loop on the outputs
  for (i=0 ; i<_p ; i++){
    // Find the closest point to iv (not itself)
    double dmin = SGTELIB::INF;
    // Loop on the points of the trainingset
    for (i2=0 ; i2<_p ; i2++){
      d = D.get(i,i2);
      if ( (i!=i2) && (d<dmin) ){
        dmin = d;
        imin = i2;
      }
    }
    _Zvs->set_row( Zs.get_row(imin) , i);
    _Svs->set_row( dmin             , i);
  }

  
  return true;

}//

/*--------------------------------------*/
/*       get_matrix_Zhs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_CN::get_matrix_Zhs (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  if ( ! _Zhs){
    _Zhs = new SGTELIB::Matrix(get_matrix_Zs());
  }
  return _Zhs;
}//

/*--------------------------------------*/
/*       get_matrix_Shs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_CN::get_matrix_Shs (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  if ( ! _Shs){
    _Shs = new SGTELIB::Matrix("Shs",_p,_m);
  }
  return _Shs;
}//


/*--------------------------------------*/
/*       get_matrix_Zvs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_CN::get_matrix_Zvs (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  compute_cv_values();
  return _Zvs;
}//

/*--------------------------------------*/
/*       get_matrix_Svs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_CN::get_matrix_Svs (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  compute_cv_values();
  return _Svs;
}//
















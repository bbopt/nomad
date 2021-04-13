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

#include "Surrogate_Kriging.hpp"


/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_Kriging::Surrogate_Kriging ( SGTELIB::TrainingSet & trainingset,
                                                SGTELIB::Surrogate_Parameters param) :
  SGTELIB::Surrogate ( trainingset , param ),
  _R                 ( "R",0,0             ),  
  _Ri                ( "Ri",0,0            ),  
  _H                 ( "H",0,0             ),
  _alpha             ( "alpha",0,0         ),
  _beta              ( "beta",0,0          ),
  _var               ( "var",0,0           ){
  #ifdef SGTELIB_DEBUG
    std::cout << "constructor Kriging\n";
  #endif

}//

/*----------------------------*/
/*          destructor        */
/*----------------------------*/
SGTELIB::Surrogate_Kriging::~Surrogate_Kriging ( void ) {

}//


/*----------------------------*/
/*          display           */
/*----------------------------*/
void SGTELIB::Surrogate_Kriging::display_private ( std::ostream & out ) const {
  //_alpha.display(out);  
  _beta.display(out);
  _var.display(out);
}//


/*--------------------------------------*/
/*             init_private             */
/*--------------------------------------*/
bool SGTELIB::Surrogate_Kriging::init_private ( void ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate_Kriging : init_private\n";
  #endif
  return true;
}//


/*--------------------------------------*/
/*               build                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_Kriging::build_private ( void ) {

  // The build mainly consists of computing alpha  

  // Compute scaling distance for each training point
  const int nvar = _trainingset.get_nvar();
  const SGTELIB::Matrix & Zs = get_matrix_Zs();

  _R = compute_covariance_matrix(get_matrix_Xs());
  _H = SGTELIB::Matrix::ones(_p,1);
  _Ri = _R.lu_inverse(&_detR);

  if (_detR<=0){
    _detR = +INF;
    return false;
  }


  //std::cout << "detR = "<< _detR << "\n";
  const SGTELIB::Matrix HRi  = _H.transpose()*_Ri;
  const SGTELIB::Matrix HRiH = HRi*_H;
  _beta = HRiH.cholesky_inverse() * HRi * Zs;
  _alpha = _Ri*(Zs-_H*_beta);

  _beta.set_name("beta");
  _alpha.set_name("alpha");
  
  _var = SGTELIB::Matrix("var",1,_m);
  double v;
  SGTELIB::Matrix Zj;
  SGTELIB::Matrix Vj;
  for (int j=0 ; j<_m ; j++){
    Zj = Zs.get_col(j);
    Zj = (Zj-_H*_beta.get_col(j));
    Vj = Zj.transpose() * _Ri * Zj;
    v = Vj.get(0,0) / (_p-nvar);
    if (v<0) return false;
    _var.set(0,j,v);
    
  }

  _ready = true;
  return true;    
  
}//





/*--------------------------------------*/
/*         Compute Design matrix        */
/*--------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate_Kriging::compute_covariance_matrix ( const SGTELIB::Matrix & XXs ) {

  // Xs can be, either the training set, to build the model, or prediction points.
  const int pxx = XXs.get_nb_rows();
  const SGTELIB::Matrix Xs = get_matrix_Xs();
  const SGTELIB::Matrix coef = _param.get_covariance_coef();
  const SGTELIB::Matrix D = _trainingset.get_distances(XXs,get_matrix_Xs(),_param.get_distance_type());

  SGTELIB::Matrix R ("R",pxx,_p);

  double d, cov;// dsum;
  for (int i1=0 ; i1<pxx ; i1++){
    for (int i2=0 ; i2<_p ; i2++){
      d = D.get(i1,i2);
      cov = coef[1] * pow(d,coef[0]);
      cov = exp(-cov);
      // Add noise if the distance is 0.
      if (d==0) cov = 1.0+_param.get_ridge();
      R.set(i1,i2,cov);
    }
  }
  
  return R;
}//




/*--------------------------------------*/
/*       predict (ZZs only)             */
/*--------------------------------------*/
void SGTELIB::Surrogate_Kriging::predict_private ( const SGTELIB::Matrix & XXs,
                                                     SGTELIB::Matrix * ZZs) {
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  const int pxx = XXs.get_nb_rows();
  const SGTELIB::Matrix r = compute_covariance_matrix(XXs).transpose();
  *ZZs =  SGTELIB::Matrix::ones(pxx,1)*_beta + r.transpose() * _alpha;
}//


void SGTELIB::Surrogate_Kriging::predict_private (const SGTELIB::Matrix & XXs,
                                                SGTELIB::Matrix * ZZs,
                                                SGTELIB::Matrix * std, 
                                                SGTELIB::Matrix * ei ,
                                                SGTELIB::Matrix * cdf) {
  check_ready(__FILE__,__FUNCTION__,__LINE__);

  const int pxx = XXs.get_nb_rows();
  const double fs_min = _trainingset.get_fs_min();
  const SGTELIB::Matrix r = compute_covariance_matrix(XXs).transpose();
  int i,j;

  // Predict ZZ
  if (ZZs) predict_private(XXs,ZZs);

  // Predict std
  if (std) std->fill(-SGTELIB::INF);
  else std = new SGTELIB::Matrix ("std",pxx,_m);

  double rRr;
  const double HRH = (_H.transpose()*_Ri*_H).get(0,0);

  double v;
  SGTELIB::Matrix ri;
  for (i=0 ; i<pxx ; i++){
    ri = r.get_col(i);
    rRr = (ri.transpose()*_Ri*ri).get(0,0);
    if (fabs(rRr-1)<EPSILON){
      v = fabs(rRr-1);
    }
    else{
      v = 1-rRr+(1-rRr)*(1-rRr)/HRH;
    }
    v = fabs(v);
    for (j=0 ; j<_m ; j++){
      std->set(i,j,v*_var[j]);
    }
  }

  // Prediction of statistical data
  if ( (ei) || (cdf) ){
    double v;
    if (ei)   ei->fill(-SGTELIB::INF);
    if (cdf) cdf->fill(-SGTELIB::INF);
    for (j=0 ; j<_m ; j++){
      if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ){
        // Compute CDF
        if (cdf){
          for (i=0 ; i<pxx ; i++){
            v = normcdf( fs_min , ZZs->get(i,j) , std->get(i,j) );
            if (v<0) v=0;
            cdf->set(i,j,v);
          }
        }
        if (ei){
          for (i=0 ; i<pxx ; i++){
            v = normei( ZZs->get(i,j) , std->get(i,j) , fs_min );
            if (v<0) v=0;
            ei->set(i,j,v );
          }
        }
      }// END CASE OBJ
      else if (_trainingset.get_bbo(j)==SGTELIB::BBO_CON){
        // Compute CDF
        if (cdf){
          // Scaled Feasibility Threshold
          double cs = _trainingset.Z_scale(0.0,j);
          for (i=0 ; i<pxx ; i++){
            v = normcdf( cs , ZZs->get(i,j) , std->get(i,j) );
            if (v<0) v=0;
            cdf->set(i,j,v);
          }
        }
      }// END CASE CON
    }// End for j
  }

}//





bool SGTELIB::Surrogate_Kriging::compute_cv_values (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);

  if ((_Zvs) && (_Svs)) return true;

  const SGTELIB::Matrix & Zs = get_matrix_Zs();
  const SGTELIB::Matrix RiH = _Ri*_H;
  const SGTELIB::Matrix Q = _Ri - RiH*( _H.transpose()*_Ri*_H)*RiH.transpose();
  const SGTELIB::Matrix dQ = Q.diag_inverse();
  
  // Init matrices
  if ( !  _Zvs){
    _Zvs = new SGTELIB::Matrix;
    *_Zvs = Zs - SGTELIB::Matrix::diagA_product(dQ,Q)*Zs;
    _Zvs->replace_nan(+INF);
    _Zvs->set_name("Zvs");
  }
    
  if ( !  _Svs){
    _Svs = new SGTELIB::Matrix ("Svs",_p,_m);
    double q;
    for (int i=0 ; i<_p ; i++){
      q = dQ.get(i,i);
      for (int j=0 ; j<_m ; j++){
        _Svs->set(i,j,sqrt(_var[j]*q));
      }
    }
    _Svs->replace_nan(+INF);
    _Svs->set_name("Svs");
  }
  return true;
}//


/*--------------------------------------*/
/*       get cv values                  */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_Kriging::get_matrix_Zvs (void){
  compute_cv_values();
  return _Zvs;
}
const SGTELIB::Matrix * SGTELIB::Surrogate_Kriging::get_matrix_Svs (void){
  compute_cv_values();
  return _Svs;
}


/*--------------------------------------*/
/*       compute linv                   */
/*--------------------------------------*/
void SGTELIB::Surrogate_Kriging::compute_metric_linv (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  if ( !is_defined(SGTELIB::METRIC_LINV) ){
    #ifdef SGTELIB_DEBUG
      std::cout << "Compute metric linv\n";
    #endif
    SGTELIB::Matrix v = SGTELIB::Matrix("LINV",0,_m);
    for (int j=0 ; j<_m ; j++){
      v.set(0,j, pow(_var[j],_p)*_detR );
    }
    _metrics[SGTELIB::METRIC_LINV] = v;
  }

}//







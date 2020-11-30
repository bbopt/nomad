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

#include "Surrogate_Ensemble.hpp"

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_Ensemble::Surrogate_Ensemble ( SGTELIB::TrainingSet & trainingset,
                                                  SGTELIB::Surrogate_Parameters param ) :
  SGTELIB::Surrogate ( trainingset , param ),
  _kmax              ( 0               ),
  _kready            ( 0               ),
  _active            ( NULL            ),
  _metric            ( new double [_m] ){

  #ifdef ENSEMBLE_DEBUG
    std::cout << "constructor Ensemble 1\n";
  #endif
  // Init Model list
  model_list_preset(_param.get_preset());
  // Init the weight matrix in _param
  SGTELIB::Matrix W ("W",_kmax,_m);
  W.fill(1.0/double(_kmax));
  _param.set_weight(W);
}


/*----------------------------*/
/*          destructor        */
/*----------------------------*/
SGTELIB::Surrogate_Ensemble::~Surrogate_Ensemble ( void ) {

  delete [] _active;
  delete [] _metric;

  for (int k=0 ; k<_kmax ; k++){
    if ( _surrogates.at(k) ){
      surrogate_delete( _surrogates.at(k) );
    }
  }
  _surrogates.clear();

}//

/*--------------------------------------*/
/*              display                 */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::display_private ( std::ostream & out ) const {

  out << "kmax: " << _kmax << "\n";
  out << "kready: " << _kready << "\n";

  SGTELIB::Matrix W = _param.get_weight();
  /*
  out << "W = [ ";
  for ( int k=0 ; k<_kmax ; k++) out << W.get(k,0) << " ";
  out << " ]\n";
  */
/*
  for (int k=0 ; k<_kmax ; k++){
    out <<"model[" << k << "]: " << _surrogates.at(k)->get_string() << "\n";
  }
*/
/*
  double w;
  for (int j=0 ; j<_m ; j++){
    out << "output " << j << ":\n";
    for ( int k=0 ; k<_kmax ; k++){
      out << "  [";
      out.width(2); 
      out << k;
      out << "]: ";
      out.width(12); 
      out << _surrogates.at(k)->get_metric(_param.get_metric_type(),j) << " ; w: ";
      
      w = W.get(k,j);
      if (w==0) out << "  0 %";
      else if (w<=0.01) out << " <1 %";
      else{
        w = double(round(w*100));
        out.width(3); 
        out << w << " %";
      }
      out << " ; ";
      out << _surrogates.at(k)->get_short_string();
      if (! is_ready(k))
        out << " (Not Ready)";
      out << "\n";
    }
    // Metric of the Ensemble
    out << "  =====>";
    out.width(8); 
    out << _metric[j] ;
    out << " ; weight:       N.A. ; " << get_short_string() << "\n";
  }

*/

 double w;
  for (int j=0 ; j<_m ; j++){
    out << "output " << _p << " " << j << ":";
    for ( int k=0 ; k<_kmax ; k++){    
      w = W.get(k,j);
      if (w>EPSILON) out << " " << k ;
    }
    out << "\n";
  }



}//


/*-----------------------------------------*/
/*     display model list                  */
/* (remove all the models of a given type) */
/*-----------------------------------------*/
void SGTELIB::Surrogate_Ensemble::model_list_display ( std::ostream & out ) {
  out << "model list (_kmax=" << _kmax << "):\n";
  if (_kmax==0){
    out << "model list is empty\n";
  }
  for (int k=0 ; k<_kmax ; k++){
    out <<"  Model " << k << ": " << _surrogates.at(k)->get_string() << "\n";
  }

}//


/*-----------------------------------------*/
/*     remove all models from model list   */
/*-----------------------------------------*/
void SGTELIB::Surrogate_Ensemble::model_list_remove_all ( void ){

  std::vector<SGTELIB::Surrogate *>::iterator it = _surrogates.begin();
  while (it != _surrogates.end()){
    SGTELIB::surrogate_delete(*it);
    it = _surrogates.erase(it);
  }
  _surrogates.clear();
  _kmax = 0;
}//

/*-----------------------------------------*/
/*   add one model                         */
/*-----------------------------------------*/
void SGTELIB::Surrogate_Ensemble::model_list_add ( const std::string & definition ){
  _surrogates.push_back( SGTELIB::Surrogate_Factory(_trainingset,definition) );
  _kmax++;
}//


/*--------------------------------------*/
/*             init_private             */
/*--------------------------------------*/
bool SGTELIB::Surrogate_Ensemble::init_private ( void ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate_Ensemble : init_private\n";
  #endif

  // Need at least 2 surrogates
  if (_kmax<=1){
    #ifdef ENSEMBLE_DEBUG
      std::cout << "Surrogate_Ensemble : _kmax : " << _kmax << "\n";
    #endif
    return false;
  }

  // Build them & count the number of ready
  _kready = 0;
  int k;
  for (k=0 ; k<_kmax ; k++){
    #ifdef ENSEMBLE_DEBUG
      std::cout << "Init model " << k << "/" << _kmax << ": " << _surrogates.at(k)->get_short_string();
    #endif
    if (_surrogates.at(k)->build()){
      _kready++;
      #ifdef ENSEMBLE_DEBUG
        std::cout << " (ready)\n";
      #endif
    }
  }  
  #ifdef ENSEMBLE_DEBUG
    std::cout << "Surrogate_Ensemble : _kready/_kmax : " << _kready << "/" << _kmax << "\n";
  #endif


  // Need at least 2 ready surrogates
  if (_kready<=1){
    return false;
  }

  // Init weights with selection
  compute_W_by_select();

  return true;
}//


/*--------------------------------------*/
/*               build                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_Ensemble::build_private ( void ) {

  #ifdef ENSEMBLE_DEBUG
    std::cout << "Surrogate_Ensemble : build_private\n";
  #endif

  int k;

  // computation of the weight
  switch (_param.get_weight_type()){
    case SGTELIB::WEIGHT_SELECT:
      compute_W_by_select();
      break;
    case SGTELIB::WEIGHT_WTA1:
      compute_W_by_wta1();
      break;
    case SGTELIB::WEIGHT_WTA3:
      compute_W_by_wta3();
      break;
    case SGTELIB::WEIGHT_OPTIM:
    case SGTELIB::WEIGHT_EXTERN:
      #ifdef ENSEMBLE_DEBUG
        std::cout << "Weight corrections\n";
      #endif
      {
      SGTELIB::Matrix W = _param.get_weight();
      for (k=0 ; k<_kmax ; k++){
        if (! is_ready(k)){
          W.set_row(0.0,k);
        }
      }
      W.normalize_cols();
      _param.set_weight(W);
      }
      break;
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
         "Surrogate_Ensemble::build(): undefined aggregation method." );
  }


  _out << "BUILD...\n";

  if (check_weight_vector()){
    #ifdef ENSEMBLE_DEBUG
      std::cout << "Weights non valid\n";
    #endif
    _ready = false;
    return false;
  }
  compute_active_models();
  _ready = true;


  // Memorize the value of the metric for each output
  for (int j=0 ; j<_m ; j++){
    _metric[j] = get_metric(_param.get_metric_type(),j);
  }


  #ifdef ENSEMBLE_DEBUG
    std::cout << "Surrogate_Ensemble : end build_private\n";
  #endif


  return true;
}//


/*--------------------------------------*/
/*       compute_active_models          */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::compute_active_models ( void ) {

  // Compute the array _active
  // (_active[k] is true if the model k is ready AND the weight in k is not null for 
  // at least one output)
  SGTELIB::Matrix W = _param.get_weight();
  if (! _active){
    _active = new bool [_kmax];
  }
  int k;
  for (k=0 ; k<_kmax ; k++){
    _active[k] = false;
    if ( is_ready(k) ){
      for (int j=0 ; j<_m ; j++){
        if ( (_trainingset.get_bbo(j)!=SGTELIB::BBO_DUM) && (W.get(k,j)>EPSILON) ){
          _active[k] = true;
          break;
        }
      }
    }
  }

}//

/*--------------------------------------*/
/*       compute_W_by_select            */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::compute_W_by_select ( void ) {

  // Init Weight matrix
  SGTELIB::Matrix W ("W", _kmax , _m );
  W.fill(0.0);

  int j,k;
  int k_best = 0;
  double metric_best;
  double metric;

  // Loop on the outputs
  for (j=0 ; j<_m ; j++){
    if (_trainingset.get_bbo(j)!=SGTELIB::BBO_DUM){

      metric_best = SGTELIB::INF;
      // Find the value of the best metric
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          metric = _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
          if (! isnan(metric)) {
            metric_best = std::min(metric,metric_best);
          }
        }
      }// end loop k

      // Find the number of surrogate that have this metric value
      k_best = 0;
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          metric = _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
          // If the metric is close to metric_best
          if ( fabs(metric-metric_best)<EPSILON ){
            // Give weight to this model
            W.set(k,j,1.0);
            // Increment k_best (number of surrogates such that metric=metric_best)
            k_best++;
          }
        }
      }// end loop k

      // Normalise
      if (k_best>1){
        for (k=0 ; k<_kmax ; k++){
          if (is_ready(k)){
            if ( W.get(k,j) > EPSILON ){
              W.set(k,j,1.0/double(k_best));
            }
          }
        }// end loop k
      }//end if


    }// end DUM
  }// end loop j

  _param.set_weight(W);

}//

/*--------------------------------------*/
/*       compute_W_by_wta1              */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::compute_W_by_wta1 ( void ) {

  #ifdef ENSEMBLE_DEBUG
    std::cout << "SGTELIB::Surrogate_Ensemble::compute_W_by_wta1\n"; 
  #endif

  // Init Weight matrix
  SGTELIB::Matrix W ("W", _kmax , _m );
  W.fill(0.0);

  int k;
  double metric;
  double metric_sum;
  double weight_sum;

  // Loop on the outputs
  for (int j=0 ; j<_m ; j++){
    if (_trainingset.get_bbo(j)!=SGTELIB::BBO_DUM){

      // Compute the sum of the metric on all the surrogates ready
      metric_sum = 0;
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          metric = _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
          if (isdef(metric)) metric_sum += metric;
        }
      }

      // Affect weight:
      if (metric_sum>EPSILON){
        for (k=0 ; k<_kmax ; k++){
          if (is_ready(k)){
            metric = _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
            if (isdef(metric)) W.set(k,j,1-metric/metric_sum);
            else               W.set(k,j,0.0);
          }
        }
      }
      else{
        for (k=0 ; k<_kmax ; k++){
          if (is_ready(k)) W.set(k,j,1.0);
        }
      }

      // Normalize
      weight_sum = 0;
      for (k=0 ; k<_kmax ; k++){
        weight_sum += W.get(k,j);
      }
      W.multiply_col( 1.0/weight_sum , j );


    } // End if not DUMM
  }// End loop on outputs

  _param.set_weight(W);

}//

/*--------------------------------------*/
/*       compute_W_by_wta3              */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::compute_W_by_wta3 ( void ) {

  #ifdef ENSEMBLE_DEBUG
    std::cout << "SGTELIB::Surrogate_Ensemble::compute_W_by_wta3\n";
  #endif

  // Init Weight matrix
  SGTELIB::Matrix W ("W", _kmax , _m );
  W.fill(0.0);

  int k;
  double metric;
  double metric_avg;
  double w;
  double w_sum;

  // Loop on the outputs
  for (int j=0 ; j<_m ; j++){

    // Compute the average of the metric on all the surrogates ready
    metric_avg = 0;
    for (k=0 ; k<_kmax ; k++){
      if (is_ready(k)){
        metric_avg += _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
      }
    }
    metric_avg /= _kready;

    if (metric_avg > EPSILON){

      // Normal WA3 method.
      // Affect un-normalized weight: (which means that the sum of the weight is not 1)
      w_sum = 0;
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          metric = _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
          w = pow( metric + wta3_alpha * metric_avg , wta3_beta );
          w_sum += w;
          W.set(k,j,w);
        }
      }
      // Then, normalize  
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          W.set(k,j,W.get(k,j)/w_sum);
        }
      }

    }
    else{

      // If the metric is null for all models, then set to 1/_kready 
      w = 1.0 / double(_kready);
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          W.set(k,j,w);
        }
      }

    }
  }// End loop on outputs

  _param.set_weight(W);

}//



/*--------------------------------------*/
/*       predict (ZZ only)              */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::predict_private ( const SGTELIB::Matrix & XXs,
                                                          SGTELIB::Matrix * ZZ ) {
  #ifdef ENSEMBLE_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
  #endif

  const SGTELIB::Matrix W = _param.get_weight();
  const int pxx = XXs.get_nb_rows();

  ZZ->fill(0.0);
  // Tmp matrix for model k
  SGTELIB::Matrix * ZZk = new SGTELIB::Matrix("ZZk",pxx,_m);
 
  double w;
  for (int k=0 ; k<_kmax ; k++){
    if (_active[k]){
      // Call the output for this surrogate
      _surrogates.at(k)->predict_private(XXs,ZZk);
      for (int j=0 ; j<_m ; j++){
        w = W.get(k,j);
        for (int i=0 ; i<pxx ; i++){
          ZZ->set(i,j,   ZZ->get(i,j) + w*ZZk->get(i,j)   );
        }// end loop i
      }// end loop j
    }// end if ready
  }//end loop k

  delete ZZk;
}//


/*--------------------------------------*/
/*         predict_private              */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::predict_private ( const SGTELIB::Matrix & XXs,
                                                          SGTELIB::Matrix * ZZ ,
                                                          SGTELIB::Matrix * std, 
                                                          SGTELIB::Matrix * ei ,
                                                          SGTELIB::Matrix * cdf) {
  #ifdef ENSEMBLE_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
  #endif

  const SGTELIB::Matrix W = _param.get_weight();

  // If no statistical information is required, use the simpler prediction method
  if (! (std || ei || cdf)){
    predict_private ( XXs, ZZ );
    return;
  }

  // Else, go for the big guns...
  const int pxx = XXs.get_nb_rows();

  // Init ZZ
  bool delete_ZZ = false;
  if ( ! ZZ){
    // if ZZ is not required, we build it anyway, but delete it in the end
    ZZ = new SGTELIB::Matrix ("ZZ",pxx,_m);
    delete_ZZ = true;
  }
  ZZ->fill(0.0);


  // Fill output matricres
  if (std) std->fill(0.0);
  if (ei)   ei->fill(0.0);
  if (cdf) cdf->fill(0.0);

  // Init tmp matrices
  SGTELIB::Matrix * ZZk  = new SGTELIB::Matrix ("ZZk" ,pxx,_m);
  SGTELIB::Matrix * stdk = new SGTELIB::Matrix ("stdk",pxx,_m);

  // Tmp matrix cdfk
  SGTELIB::Matrix * cdfk;
  if (cdf) cdfk = new SGTELIB::Matrix ("cdfk",pxx,_m);
  else     cdfk = NULL;

  // Same story for ei.
  SGTELIB::Matrix * eik;
  if (ei) eik = new SGTELIB::Matrix ("eik",pxx,_m);
  else    eik = NULL;

  double w,z,s;
  

  // Loop on the models
  for (int k=0 ; k<_kmax ; k++){
    if (_active[k]){

      // Call the output for this surrogate
      _surrogates.at(k)->predict_private(XXs,ZZk,stdk,eik,cdfk);

      for (int j=0 ; j<_m ; j++){
        w = W.get(k,j);
        if (w>EPSILON/_kmax){

          // Compute ZZ
          for (int i=0 ; i<pxx ; i++){
            z = ZZk->get(i,j);
            ZZ->set( i,j, ZZ->get(i,j) + w*z );
          }

          // Compute std
          if (std){          
            for (int i=0 ; i<pxx ; i++){
              z = ZZk->get(i,j);
              s = stdk->get(i,j);
              std->set(i,j, std->get(i,j) + w*(s*s + z*z) );
              //std->set(i,j, std->get(i,j) + w*z*z );
            }// end loop i
          }
    
          // EI is linear on w
          if ( (ei) && (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ) ){
            for (int i=0 ; i<pxx ; i++){
              ei->set(i,j, ei->get(i,j) + w*eik->get(i,j) );
            }// end loop i
          }
          
          // CDF is linear on w
          if (cdf){
            for (int i=0 ; i<pxx ; i++){
              cdf->set(i,j, cdf->get(i,j) + w*cdfk->get(i,j) );
            }// end loop i
          }

        }// end if w>eps
      }// end loop j

    }// end if ready
  }//end loop k


  // Correction of std
  if (std){
    for (int j=0 ; j<_m ; j++){
      for (int i=0 ; i<pxx ; i++){
        z = ZZ->get(i,j);
        s = std->get(i,j) - z*z;
        std->set(i,j, sqrt(fabs(s)) );
      }// end loop i
    }
  }

  if (delete_ZZ) delete ZZ;
  if (ZZk ) delete ZZk;
  if (stdk) delete stdk;
  if (eik ) delete eik;
  if (cdfk) delete cdfk;

}//
 
/*--------------------------------------*/
/*       get_matrix_Zvs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_Ensemble::get_matrix_Zvs (void){
  if ( ! _Zvs){
    #ifdef ENSEMBLE_DEBUG
      check_ready(__FILE__,__FUNCTION__,__LINE__);
    #endif
    const SGTELIB::Matrix W = _param.get_weight(); 
    _Zvs = new SGTELIB::Matrix("Zv",_p,_m);
    _Zvs->fill(0.0);
    int i,j;
    double wkj;

    for (int k=0 ; k<_kmax ; k++){
      if (_active[k]){
        // Call the output for this surrogate
        const SGTELIB::Matrix * Zvs_k = _surrogates.at(k)->get_matrix_Zvs();
        for ( j=0 ; j<_m ; j++){
          wkj = W.get(k,j);
          if (wkj>0){
            for ( i=0 ; i<_p ; i++){
              _Zvs->add(i,j,  wkj*Zvs_k->get(i,j) );
            }
          }
        }// end loop j
      }// end if ready
    }//end loop k

    _Zvs->set_name("Zvs");
    _Zvs->replace_nan(+INF);

  }
  return _Zvs;
}//

/*--------------------------------------*/
/*       get_matrix_Zhs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_Ensemble::get_matrix_Zhs (void){
  if ( ! _Zhs){
    #ifdef ENSEMBLE_DEBUG
      check_ready(__FILE__,__FUNCTION__,__LINE__);
    #endif
    const SGTELIB::Matrix W = _param.get_weight();
    _Zhs = new SGTELIB::Matrix("Zv",_p,_m);
    _Zhs->fill(0.0);
    int i,j;
    double wkj;

    for (int k=0 ; k<_kmax ; k++){
      if (_active[k]){
        // Call the output for this surrogate
        const SGTELIB::Matrix * Zhs_k = _surrogates.at(k)->get_matrix_Zhs();
        for ( j=0 ; j<_m ; j++){
          wkj = W.get(k,j);
          if (wkj>0){
            for ( i=0 ; i<_p ; i++){
              _Zhs->add(i,j,  wkj*Zhs_k->get(i,j) );
            }
          }
        }// end loop j
      }// end if ready
    }//end loop k

    _Zhs->set_name("Zhs");
    _Zhs->replace_nan(+INF);

  }
  return _Zhs;
}//

/*--------------------------------------*/
/*       get_matrix_Shs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_Ensemble::get_matrix_Shs (void){
  if ( ! _Shs){
    const SGTELIB::Matrix W = _param.get_weight();
    _Shs = new SGTELIB::Matrix("Zv",_p,_m);
    _Shs->fill(0.0);
    SGTELIB::Matrix col ("col",_p,1);

    int i,j;
    double wkj;
    for (int k=0 ; k<_kmax ; k++){
      if (_active[k]){
        // Call the output for this surrogate
        const SGTELIB::Matrix * Zhs_k = _surrogates.at(k)->get_matrix_Zhs();
        const SGTELIB::Matrix * Shs_k = _surrogates.at(k)->get_matrix_Shs();

        for ( j=0 ; j<_m ; j++){
          wkj = W.get(k,j);
          if (wkj>0){
            for ( i=0 ; i<_p ; i++){
              _Shs->add(i,j,  wkj*( pow(Shs_k->get(i,j),2) + pow(Zhs_k->get(i,j),2) ) );
            }
          }
        }// end loop j
      }// end if ready
    }//end loop k

    const SGTELIB::Matrix * Zhs = get_matrix_Zhs();
    _Shs->sub( Matrix::hadamard_square( *Zhs ) );
    _Shs->hadamard_sqrt();

    _Shs->set_name("Shs");
    _Shs->replace_nan(+INF);

  }
  return _Shs;
}//





/*--------------------------------------*/
/*  to know if basic model k is ready   */
/*--------------------------------------*/
bool SGTELIB::Surrogate_Ensemble::is_ready (const int k) const{
  if ((k<0) || (k>=_kmax)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (const int k): k out of range" );
  }
  return _surrogates.at(k)->is_ready();
}


/*--------------------------------------*/
/*  external set of the weight vector   */
/*    (use model k for output j)        */
/*--------------------------------------*/
/*
void SGTELIB::Surrogate_Ensemble::set_weight_vector (const int k, const int j){
  if (_param.get_weight_type() != SGTELIB::WEIGHT_EXTERN){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): Not in EXTERN mode" );
  }
  if ((k<0) or (k>=_kmax)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): k out of range" );
  }
  if ((j<0) or (j>=_m)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): k out of range" );
  }
  if (not is_ready(k)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): Surrogate not ready" );
  }

  // Set the column j to 0
  _W.set_col( 0.0 , j );
  // Select model k for output j
  _W.set(k,j,1.0); 
  // Check and reset
  reset_metrics();
  compute_active_models();
}//
*/

/*--------------------------------------*/
/*  external set of the weight vector   */
/*   (use model k for every output)     */
/*--------------------------------------*/
/*
void SGTELIB::Surrogate_Ensemble::set_weight_vector (const int k){
  if (_param.get_weight_type() != SGTELIB::WEIGHT_EXTERN){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): Not in EXTERN mode" );
  }
  if ((k<0) or (k>=_kmax)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): k out of range" );
  }
  if (not is_ready(k)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): Surrogate not ready" );
  }
  // Put _W at 0
  _W.fill(0.0);
  // Put model k at 1.0 for every output
  _W.set_row( 1.0 , k );
  // Check and reset
  reset_metrics();
  compute_active_models();
}//
*/

/*--------------------------------------*/
/*  external set of the weight vector   */
/*          (with the whole matrix)     */
/*--------------------------------------*/
/*
void SGTELIB::Surrogate_Ensemble::set_weight_vector (const SGTELIB::Matrix & W){
  if (_param.get_weight_type() != SGTELIB::WEIGHT_EXTERN){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble::set_weight_vector (k,j): Not in EXTERN mode" );
  }
  // Set _W
  _W = W;
  // Check and reset
  reset_metrics();
  compute_active_models();
}//
*/

/*--------------------------------------*/
/*   check the weight vector            */
/*--------------------------------------*/
bool SGTELIB::Surrogate_Ensemble::check_weight_vector ( void ) const {
  const SGTELIB::Matrix W = _param.get_weight();
  double s,w;
  int j,k;
  for (j=0 ; j<_m ; j++){
    if (_trainingset.get_bbo(j)!=SGTELIB::BBO_DUM){
      for (k=0 ; k<_kmax ; k++){
        w = W.get(k,j);
        if (w<-EPSILON)    return true;   
        if (w>1+EPSILON)   return true;
        if ( isnan(w) ) return true;
      }
      s = W.get_col(j).sum();
      if (fabs(s-1.0)>_kready*EPSILON) return true;
    }
  }

  return false;

}//






/*--------------------------------------*/
/*        define model list             */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble::model_list_preset ( const std::string & preset ) {


    #ifdef ENSEMBLE_DEBUG
      std::cout << "Build model list\n";
    #endif

    model_list_remove_all();

    const std::string p = toupper(preset);
    const std::string m = " METRIC_TYPE "+_param.get_metric_type_str();
    const std::string d = " DISTANCE_TYPE "+_param.get_distance_type_str();
    const std::string dm = d+m;

    if (SGTELIB::streqi(p,"DEFAULT")) {
      model_list_add("TYPE PRS DEGREE 1 RIDGE 0");
      model_list_add("TYPE PRS DEGREE 1 RIDGE 0.001");
      model_list_add("TYPE PRS DEGREE 2 RIDGE 0");
      model_list_add("TYPE PRS DEGREE 2 RIDGE 0.001");
      model_list_add("TYPE PRS DEGREE 3 RIDGE 0.0");
      model_list_add("TYPE PRS DEGREE 6 RIDGE 0.001");
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 0.1"+dm);
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 0.3"+dm);
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 1  "+dm);
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 3  "+dm);
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 10 "+dm);
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 0.3"+dm);
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 1  "+dm);
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 3  "+dm);
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 10 "+dm);
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE I1"+dm);
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE I2"+dm);
      model_list_add("TYPE CN"+dm);
    }
    else if (SGTELIB::streqi(p,"KS")) {
      model_list_add("TYPE KS KERNEL_TYPE D1 KERNEL_COEF 0.1"+d);
      model_list_add("TYPE KS KERNEL_TYPE D1 KERNEL_COEF 0.2"+d); 
      model_list_add("TYPE KS KERNEL_TYPE D1 KERNEL_COEF 0.5"+d);
      model_list_add("TYPE KS KERNEL_TYPE D1 KERNEL_COEF 1  "+d);
      model_list_add("TYPE KS KERNEL_TYPE D1 KERNEL_COEF 2  "+d);
      model_list_add("TYPE KS KERNEL_TYPE D1 KERNEL_COEF 5  "+d);
      model_list_add("TYPE KS KERNEL_TYPE D1 KERNEL_COEF 10 "+d);
    }
    else if (SGTELIB::streqi(p,"PRS")) {
      model_list_add("TYPE PRS DEGREE 1");
      model_list_add("TYPE PRS DEGREE 2");
      model_list_add("TYPE PRS DEGREE 3");
      model_list_add("TYPE PRS DEGREE 4");
      model_list_add("TYPE PRS DEGREE 5");
      model_list_add("TYPE PRS DEGREE 6");
    }
    else if (SGTELIB::streqi(p,"IS0")) {
      model_list_add("TYPE PRS_EDGE DEGREE 2");
      model_list_add("TYPE PRS_EDGE DEGREE 3");

      model_list_add("TYPE KS            KERNEL_TYPE D1 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D1 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_IS0"); 
      model_list_add("TYPE KS            KERNEL_TYPE D1 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D1 KERNEL_COEF 1   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D1 KERNEL_COEF 2   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D1 KERNEL_COEF 5   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D1 KERNEL_COEF 10  DISTANCE_TYPE NORM2_IS0");

      model_list_add("TYPE KS            KERNEL_TYPE D2 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D2 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_IS0"); 
      model_list_add("TYPE KS            KERNEL_TYPE D2 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D2 KERNEL_COEF 1   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D2 KERNEL_COEF 2   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D2 KERNEL_COEF 5   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE KS            KERNEL_TYPE D2 KERNEL_COEF 10  DISTANCE_TYPE NORM2_IS0");

      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D1 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D1 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_IS0"); 
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D1 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D1 KERNEL_COEF 1   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D1 KERNEL_COEF 2   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D1 KERNEL_COEF 5   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D1 KERNEL_COEF 10  DISTANCE_TYPE NORM2_IS0");

      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D2 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D2 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_IS0"); 
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D2 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D2 KERNEL_COEF 1   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D2 KERNEL_COEF 2   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D2 KERNEL_COEF 5   DISTANCE_TYPE NORM2_IS0");
      model_list_add("TYPE RBF  PRESET I KERNEL_TYPE D2 KERNEL_COEF 10  DISTANCE_TYPE NORM2_IS0");
    }
    else if (SGTELIB::streqi(p,"CAT")) {
      model_list_add("TYPE PRS_CAT DEGREE 2");
      model_list_add("TYPE PRS_CAT DEGREE 3");

      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_CAT"); 
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 1   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 2   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 5   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D1 KERNEL_COEF 10  DISTANCE_TYPE NORM2_CAT");

      model_list_add("TYPE KS           KERNEL_TYPE D2 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D2 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_CAT"); 
      model_list_add("TYPE KS           KERNEL_TYPE D2 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D2 KERNEL_COEF 1   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D2 KERNEL_COEF 2   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D2 KERNEL_COEF 5   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE KS           KERNEL_TYPE D2 KERNEL_COEF 10  DISTANCE_TYPE NORM2_CAT");

      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_CAT"); 
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 1   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 2   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 5   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D1 KERNEL_COEF 10  DISTANCE_TYPE NORM2_CAT");

      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D2 KERNEL_COEF 0.1 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D2 KERNEL_COEF 0.2 DISTANCE_TYPE NORM2_CAT"); 
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D2 KERNEL_COEF 0.5 DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D2 KERNEL_COEF 1   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D2 KERNEL_COEF 2   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D2 KERNEL_COEF 5   DISTANCE_TYPE NORM2_CAT");
      model_list_add("TYPE RBF PRESET I KERNEL_TYPE D2 KERNEL_COEF 10  DISTANCE_TYPE NORM2_CAT");
    }
    else if (SGTELIB::streqi(p,"SUPER1")) {
      model_list_add("TYPE KS     KERNEL_TYPE OPTIM KERNEL_COEF OPTIM"+dm);
      model_list_add("TYPE RBF    KERNEL_TYPE OPTIM KERNEL_COEF OPTIM RIDGE 0.001 PRESET I"+dm);
      model_list_add("TYPE PRS    DEGREE OPTIM RIDGE OPTIM"+m);
      model_list_add("TYPE LOWESS DEGREE OPTIM RIDGE 0.001 KERNEL_COEF OPTIM KERNEL_TYPE D1"+dm);
    }
    else if (SGTELIB::streqi(p,"SMALL")) {
      model_list_add("TYPE PRS");
      model_list_add("TYPE KS");
      model_list_add("TYPE RBF PRESET I");
    }
    else if (SGTELIB::streqi(p,"NONE")) {
      // None
    }
    else {
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
        "Surrogate_Ensemble::model_list_preset: unrecognized preset \""+preset+"\"" );
    }  

    #ifdef ENSEMBLE_DEBUG
      std::cout << "END Build model list\n";
    #endif

}//



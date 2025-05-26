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

#include "Surrogate_Ensemble_Stat.hpp"

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_Ensemble_Stat::Surrogate_Ensemble_Stat ( SGTELIB::TrainingSet & trainingset,
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
SGTELIB::Surrogate_Ensemble_Stat::~Surrogate_Ensemble_Stat ( void ) {

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
void SGTELIB::Surrogate_Ensemble_Stat::display_private ( std::ostream & out ) const {

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
void SGTELIB::Surrogate_Ensemble_Stat::model_list_display ( std::ostream & out ) {
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
void SGTELIB::Surrogate_Ensemble_Stat::model_list_remove_all ( void ){

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
void SGTELIB::Surrogate_Ensemble_Stat::model_list_add ( const std::string & definition ){
  _surrogates.push_back( SGTELIB::Surrogate_Factory(_trainingset,definition) );
  _kmax++;
}//


/*--------------------------------------*/
/*             init_private             */
/*--------------------------------------*/
bool SGTELIB::Surrogate_Ensemble_Stat::init_private ( void ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate_Ensemble_Stat : init_private\n";
  #endif

  // Need at least 2 surrogates
  if (_kmax<=1){
    #ifdef ENSEMBLE_DEBUG
      std::cout << "Surrogate_Ensemble_Stat : _kmax : " << _kmax << "\n";
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
    std::cout << "Surrogate_Ensemble_Stat : _kready/_kmax : " << _kready << "/" << _kmax << "\n";
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
bool SGTELIB::Surrogate_Ensemble_Stat::build_private ( void ) {

  #ifdef ENSEMBLE_DEBUG
    std::cout << "Surrogate_Ensemble_Stat : build_private\n";
  #endif

  int k;

  // Some parameters
  _sigma_mult_stat = _param.get_sigma_mult();
  _lambda_p_stat = _param.get_lambda_p();
  _lambda_pi_stat = _param.get_lambda_pi();

  // compute simplex or pss
  switch (_param.get_uncertainty_type()){
    // If smooth, build a simplex around x
    case SGTELIB::UNCERTAINTY_SMOOTH:
      build_simplex_private();
      if (_nbd != _n+1){
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
         "Surrogate_Ensemble_Stat::build(): Number or directions in simplex must be n+1 = "
          + std::to_string(_n+1) + " but is " + std::to_string(_nbd) );
      }
      break;
    // If nonsmooth, build a positive spanning set around x
    case SGTELIB::UNCERTAINTY_NONSMOOTH:
      build_pss_private();
      if ( _nbd < _n+1 ){
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
         "Surrogate_Ensemble_Stat::build(): Number or directions in PSS must be superior to n+1 = "
          + std::to_string(_n+1) + " but is " + std::to_string(_nbd) );
      }
      else if ( _nbd > 2*_n ){
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
         "Surrogate_Ensemble_Stat::build(): Number or directions in PSS must be inferior to 2*n = "
          + std::to_string(2*_n) + " but is " + std::to_string(_nbd) );
      }
      break;
  }
  

  // computation of the weight
  switch (_param.get_weight_type()){
    case SGTELIB::WEIGHT_SELECT:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
         "Surrogate_Ensemble_Stat::build(): WEIGHT SELECT method incompatible with Ensemble_Stat surrogate." );
      // compute_W_by_select();
      break;
    case SGTELIB::WEIGHT_SELECT2:
      compute_W_by_select_nb(2);
      break;
    case SGTELIB::WEIGHT_SELECT3:
      compute_W_by_select_nb(3);
      break;
    case SGTELIB::WEIGHT_SELECT4:
      compute_W_by_select_nb(4);
      break;
    case SGTELIB::WEIGHT_SELECT5:
      compute_W_by_select_nb(5);
      break;
    case SGTELIB::WEIGHT_SELECT6:
      compute_W_by_select_nb(6);
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
         "Surrogate_Ensemble_Stat::build(): undefined aggregation method." );
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
    std::cout << "Surrogate_Ensemble_Stat : end build_private\n";
  #endif


  return true;
}//


/*-------------------------------------*/
/*       Build Simplex directions      */
/*-------------------------------------*/
void SGTELIB::Surrogate_Ensemble_Stat::build_simplex_private ( void ) {
  
  _SET = SGTELIB::Matrix ("Set", _n+1, _n);
  SGTELIB::Matrix ones ("ones", 1, _n);
  ones.set_row(1,0);
  double factor = - (1 + 1/sqrt(_n+1)) / _n;
  // Size parameter
  _t = _param.get_size_param();

  // Add direction _t * ( e_i - factor * (1,...,1) ) for i = 1,...,_n
  for (int i=0 ; i<_n ; i++){
    _SET.set_row( _t * factor * ones , i );
    _SET.set( i,i, _SET.get(i,i) + _t ); 
  }

  // Add last direction : _t * 1/sqrt(_n+2) * (1,...,1)
  _SET.set_row( ( _t/sqrt(2*_n+2) ) * ones , _n );

  _nbd = _SET.get_nb_rows();

}


/*---------------------------------------------------*/
/*      Build Positive Spanning Set directions       */
/*  (in this case the set of coordinate directions)  */
/*---------------------------------------------------*/
void SGTELIB::Surrogate_Ensemble_Stat::build_pss_private ( void ) {
  
  _SET = SGTELIB::Matrix ("_SET", 2*_n, _n);
  // Size parameter
  _t = _param.get_size_param();

  // Add directions _t * e_i and -_t * e_i for i = 1,...,_n
  for (int i=0 ; i<_n ; i++){

    _SET.set(   2*i,i,  _t );
    _SET.set( 2*i+1,i, -_t );
  }

  _nbd = _SET.get_nb_rows();

}



/*--------------------------------------------*/
/*   Build a set of points around point x     */
/*          and save it in XXd                */
/*--------------------------------------------*/ 
void SGTELIB::Surrogate_Ensemble_Stat::build_set_around_x( const SGTELIB::Matrix & XXs,
                                                           std::vector<SGTELIB::Matrix *> & XXd) const {
                                                        
  const int pxx = XXs.get_nb_rows();
  SGTELIB::Matrix x_plus_d;


  for (int i=0 ; i<pxx ; i++){
    for (int d=0 ; d<_nbd ; d++){
      x_plus_d = XXs.get_row(i) + _SET.get_row(d); // compute x + d
      XXd[i]->set_row(x_plus_d, d); // save x + d
    }  
  }

}


/*----------------------------------------------------*/
/*       Compute simplex gradients of models          */
/*                (objectives only)                   */
/*----------------------------------------------------*/
void SGTELIB::Surrogate_Ensemble_Stat::compute_simplex_gradient( const std::vector<SGTELIB::Matrix *> & XXd,
                                                                 const std::vector<SGTELIB::Matrix *> & ZZsurr_around,
                                                                       std::vector<SGTELIB::Matrix *> & SGradsurr) const {
  
  int pxx = ZZsurr_around[0]->get_nb_rows();
  SGTELIB::Matrix coeff ("coeff", _n+1, 1);
  SGTELIB::Matrix ones ("ones", _nbd, 1);
  ones.set_col(1,0);
  SGTELIB::Matrix A ("A", _nbd, _n);
  SGTELIB::Matrix Ai ("Ai", _nbd, _n+1);
  SGTELIB::Matrix ft;


  // Loop prediction points
  for (int i=0 ; i<pxx ; i++){

    A = *(XXd[i]);
    A.add_cols(ones);
    Ai = A.SVD_inverse();

    // Loop on models
    for (int k=0 ; k<_kmax ; k++){
      ft = ( ZZsurr_around[k]->get_row(i) ).transpose();
      coeff = Ai * ft;
      SGradsurr[k]->set_row( ( coeff.get_rows(0,_n) ).transpose() , i );
    }

  } // end for i

}


/*--------------------------------------*/
/*       compute_active_models          */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble_Stat::compute_active_models ( void ) {

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
void SGTELIB::Surrogate_Ensemble_Stat::compute_W_by_select ( void ) {

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
          if ( std::fabs(metric-metric_best)<EPSILON ){
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
/*      compute_W_by_select k bests     */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble_Stat::compute_W_by_select_nb ( const int nb_bests ) {

  // Init Weight matrix
  SGTELIB::Matrix W ("W", _kmax , _m );
  W.fill(0.0);

  int j,k;
  int k_best = 0;
  double metric;
  double metric_best;
  double metric_sum;
  double weight_sum;

  // Loop on the outputs
  for (j=0 ; j<_m ; j++){
    if (_trainingset.get_bbo(j)!=SGTELIB::BBO_DUM){
      
      // Memorize all metrics and indices
      // Find the value of the best metric
      std::vector< double > metrics; // to record metrics
      std::vector< bool > selected; // to record selected models
      metric_best = SGTELIB::INF;
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          metric = _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
          metrics.push_back(metric);
          if (! isnan(metric)) {
            metric_best = std::min(metric,metric_best);
          }
        }
        else{
          metrics.push_back(0); // just to keep tabs on indices
        }
        selected.push_back(false); // no model is selected as of now
      }// end loop k

      // Find the number of surrogate that have this metric value
      k_best = 0;
      for (k=0 ; k<_kmax ; k++){
        if (is_ready(k)){
          metric = _surrogates.at(k)->get_metric(_param.get_metric_type(),j);
          // If the metric is close to metric_best
          if ( std::fabs(metric-metric_best)<EPSILON ){
            // Give weight to this model
            W.set(k,j,1.0);
            // Increment k_best (number of surrogates such that metric=metric_best)
            k_best++;
          }
        }
      }// end loop k

      if (k_best >= nb_bests){
        // Keep all the k_best models
        // Normalise as in select
        for (k=0 ; k<_kmax ; k++){
          if (is_ready(k)){
            if ( W.get(k,j) > EPSILON ){
              W.set(k,j,1.0/double(k_best));
            }
          }
        }// end loop k
      }
      else{
        // Select nb_bests models
        int index_best=-1;
        metric_sum = 0;
        bool found_a_model;
        // Do nb_bests passes and select the best available model each time
        for (int i=0 ; i<nb_bests ; i++){
          metric_best = SGTELIB::INF;
          found_a_model = false;
          for (k=0 ; k<_kmax ; k++){
            if ( is_ready(k) && !selected[k] ){
              if ( isdef(metrics[k]) && (metrics[k] < metric_best) ) {
                metric_best = metrics[k];
                index_best = k;
                found_a_model = true;
              }
            }
          }// end for k
          if (found_a_model){
            metric_sum += metric_best;
            selected[index_best] = true;
          }
        } // end for i

        // Affect weight:
        if (metric_sum>EPSILON){
          for (k=0 ; k<_kmax ; k++){
            if (selected[k]){
              // If metric of model k is equal to metric_sum
              // set weigth to 0.1
              if ( std::fabs(metrics[k]-metric_sum)<EPSILON ){
                W.set(k,j, 0.1);
              }
              else{
                W.set(k,j, 1-metrics[k]/metric_sum);
              }
            }
            else{
              W.set(k,j,0.0);
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

        // Normalise as in wta1
      }// end if (k_best >= nb_bests) else

      
    }// end DUM
  }// end loop j

  _param.set_weight(W);

}//


/*--------------------------------------*/
/*       compute_W_by_wta1              */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble_Stat::compute_W_by_wta1 ( void ) {

  #ifdef ENSEMBLE_DEBUG
    std::cout << "SGTELIB::Surrogate_Ensemble_Stat::compute_W_by_wta1\n"; 
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
void SGTELIB::Surrogate_Ensemble_Stat::compute_W_by_wta3 ( void ) {

  #ifdef ENSEMBLE_DEBUG
    std::cout << "SGTELIB::Surrogate_Ensemble_Stat::compute_W_by_wta3\n";
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
void SGTELIB::Surrogate_Ensemble_Stat::predict_private ( const SGTELIB::Matrix & XXs,
                                                               SGTELIB::Matrix * ZZ ) {
  #ifdef ENSEMBLE_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
  #endif

  const SGTELIB::Matrix W = _param.get_weight();
  const int pxx = XXs.get_nb_rows();

  ZZ->fill(0.0);
  // Tmp matrix for model k
  SGTELIB::Matrix * ZZk = new SGTELIB::Matrix("ZZk", pxx, _m);
 
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
void SGTELIB::Surrogate_Ensemble_Stat::predict_private ( const SGTELIB::Matrix & XXs,
                                                               SGTELIB::Matrix * ZZ ,
                                                               SGTELIB::Matrix * std, 
                                                               SGTELIB::Matrix * ei ,
                                                               SGTELIB::Matrix * cdf) {
  #ifdef ENSEMBLE_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
  #endif

  const SGTELIB::Matrix W = _param.get_weight();
  const uncertainty_t uncertainty_type = _param.get_uncertainty_type();

  // If no statistical information is required, use the simpler prediction method
  if (! (std || ei || cdf)){
    predict_private ( XXs, ZZ );
    return;
  }

  //std::cout << "Computing std in ENSEMBLE_STAT\n";
  //std::cout << "Uncertainty is " << _param.get_uncertainty_type() << "\n";

  // Else, go for the big guns...
  const int pxx = XXs.get_nb_rows();
  const double fs_min = _trainingset.get_fs_min();

  // Init ZZ
  bool delete_ZZ = false;
  if ( ! ZZ){
    // if ZZ is not required, we build it anyway, but delete it in the end
    ZZ = new SGTELIB::Matrix ("ZZ", pxx, _m);
    delete_ZZ = true;
  }
  ZZ->fill(0.0);


  // Fill output matrices
  if (std) std->fill(0.0);
  if (ei)   ei->fill(0.0);
  if (cdf) cdf->fill(0.0);


  // Init value matrices

  // Points around the points of XXs
  std::vector<SGTELIB::Matrix *> XXd;
  for (int i=0 ; i<pxx ; i++){
    XXd.push_back( new SGTELIB::Matrix ("XXd" + std::to_string(i), _nbd, _n) );
    // XXd[i] is of dimension _nbd * _n
  }

  // Values of the _kmax surrogates at points of XXs
  std::vector<SGTELIB::Matrix *> ZZsurr;
  for (int k=0 ; k<_kmax ; k++){
    ZZsurr.push_back( new SGTELIB::Matrix ("ZZsurr" + std::to_string(k), pxx, _m) );
    // ZZsurr[k] is of dimension pxx * _m
  }

  // Values of the _kmax surrogates at points around those of XXs (positive spanning set or simplex) - objectives only
  std::vector<SGTELIB::Matrix *> ZZsurr_around;
  for (int k=0 ; k<_kmax ; k++){
    ZZsurr_around.push_back( new SGTELIB::Matrix ("ZZsurr_around" + std::to_string(k), pxx, _nbd) );
    // ZZsurr_around[k] is of dimension pxx * _nbd
  }

  // Simplex gradients of the _kmax surrogates at points of XXs - objectives only
  std::vector<SGTELIB::Matrix *> SGradsurr;
  for (int k=0 ; k<_kmax ; k++){
    SGradsurr.push_back( new SGTELIB::Matrix ("SGradsurr" + std::to_string(k), pxx, _n) );
    // SGradsurr[k] is of dimension pxx * _n
  }  

  double w,z;

  // Build sets of points (PSS or simplex) around points of XXs
  build_set_around_x(XXs, XXd);

  // Loop on the models
  for (int k=0 ; k<_kmax ; k++){
    if (_active[k]){

      // Call the output for this surrogate
      _surrogates.at(k)->predict_private(XXs, ZZsurr[k]);

      for (int i=0 ; i<pxx ; i++){

        // Compute ZZ
        for (int j=0 ; j<_m ; j++){
          w = W.get(k,j);
          if (w>EPSILON/_kmax){
            z = ZZsurr[k]->get(i,j);
            ZZ->set( i,j, ZZ->get(i,j) + w*z );
          }

        } // end loop j
      
      } // end for i

      // Compute model values around point i of XXs - objectives only
      _surrogates.at(k)->predict_private_objective(XXd, ZZsurr_around[k]);
      
      //
      // At this point, we have all the model values we need
      //
    
    }

  } // end for k

  // If smooth uncertainty is chosen, compute simplex gradients
  if (uncertainty_type == SGTELIB::UNCERTAINTY_SMOOTH){
  //if (_param.get_uncertainty_type() == SGTELIB::UNCERTAINTY_SMOOTH){
    compute_simplex_gradient(XXd, ZZsurr_around, SGradsurr);
  }



  // Prediction of statistical data (very similar to Kriging models)
  if ( (std) || (ei) || (cdf) ){
    double v;
    if (std) std->fill(-SGTELIB::INF);
    if (ei)   ei->fill(-SGTELIB::INF);
    if (cdf) cdf->fill(-SGTELIB::INF);
  
    for (int j=0 ; j<_m ; j++){

      // Compute STD
      if (std){
        for (int i=0 ; i<pxx ; i++){
            v = compute_sigma(i,j, ZZsurr, ZZsurr_around, SGradsurr);
            std->set(i,j, v);
          }
        }

      if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ){

        // Compute CDF
        if (cdf){
          for (int i=0 ; i<pxx ; i++){
            //v = normcdf( fs_min, ZZ->get(i,j), std->get(i,j) ); // stochastic version (no sigmoid)
            v = sigmoid( fs_min, ZZ->get(i,j), std->get(i,j) , _lambda_pi_stat ); // new version
            if (v<0) v=0;
            cdf->set(i,j,v);
          }
        }
        // Compute EI
        if (ei){
          for (int i=0 ; i<pxx ; i++){
            //v = normei( ZZ->get(i,j), std->get(i,j), fs_min ); // stochastic version (no sigmoid)
            v = newei( ZZ->get(i,j), std->get(i,j), fs_min ); // new version
            if (v<0) v=0;
            ei->set(i,j, v);
          }
        }
      }// END CASE OBJ

      else if (_trainingset.get_bbo(j)==SGTELIB::BBO_CON){

        // Compute CDF
        if (cdf){
          // Scaled Feasibility Threshold
          double cs = _trainingset.Z_scale(0.0,j);
          for (int i=0 ; i<pxx ; i++){
            //v = normcdf( cs , ZZ->get(i,j) , std->get(i,j) ); // stochastic version (no sigmoid)
            v = sigmoid( cs, ZZ->get(i,j), std->get(i,j) , _lambda_p_stat ); // new version
            if (v<0) v=0;
            cdf->set(i,j, v);
          }
        }

      }// END CASE CON
    }// End for j
  }
  

  // Delete
  for (int i=0 ; i<pxx ; i++ ){
    delete XXd[i];
  }
  for (int k=0 ; k<_kmax ; k++ ){
    delete ZZsurr[k];
    delete ZZsurr_around[k];
    delete SGradsurr[k];
  }
  
  if (delete_ZZ) delete ZZ;

}//


/*-----------------------------------*/
/*       Compute sigma (std)         */
/*-----------------------------------*/
double SGTELIB::Surrogate_Ensemble_Stat::compute_sigma ( const int i, const int j,
                                                         const std::vector<SGTELIB::Matrix *> & ZZsurr,
                                                         const std::vector<SGTELIB::Matrix *> & ZZsurr_around,
                                                         const std::vector<SGTELIB::Matrix *> & SGradsurr         ) const {
  
  double sigma = 0; // result
  double num=0, denum=0, wk=0, wl=0, sigma_kl=0;
  const SGTELIB::Matrix W = _param.get_weight();
  const uncertainty_t uncertainty_type = _param.get_uncertainty_type();

  // Loop on pairs of models
  for (int k=0 ; k<_kmax-1 ; k++){
    if (_active[k]){

      for (int l=k+1 ; l<_kmax ; l++){
        if (_active[l]){

          wk = W.get(k,j);
          wl = W.get(l,j);

          if ( (wk > EPSILON/_kmax) && (wl > EPSILON/_kmax) ){
          
            // If objective
            if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ){
              // If smooth uncertainty
              if (uncertainty_type == SGTELIB::UNCERTAINTY_SMOOTH){
                sigma_kl = compute_sigma_kl_obj_smooth(i, k, l, SGradsurr);
              }
              // If nonsmooth uncertainty
              else if (uncertainty_type == SGTELIB::UNCERTAINTY_NONSMOOTH){
                sigma_kl = compute_sigma_kl_obj_nonsmooth(i, j, k, l, ZZsurr, ZZsurr_around);
              }  
            } // end if OBJ

            // If constraint
            if (_trainingset.get_bbo(j)==SGTELIB::BBO_CON){
              // If smooth uncertainty
              if (uncertainty_type == SGTELIB::UNCERTAINTY_SMOOTH){
                sigma_kl = compute_sigma_kl_con_smooth(i, j, k, l, ZZsurr);
              }
              // If nonsmooth uncertainty
              else if (uncertainty_type == SGTELIB::UNCERTAINTY_NONSMOOTH){
                sigma_kl = compute_sigma_kl_con_nonsmooth(i, j, k, l, ZZsurr);
              }  
            } // end if CON

          num += wk*wl*sigma_kl; // w.transpose * Sigma * w  in the paper
          denum += wk*wl;        // w.transpose * T * w  in the paper

          }
        }

      } // end for l
    }
  } // end for k

  if (denum < EPSILON) {
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::compute_sigma (): division by zero" );
  }
  else if (num < EPSILON) {
    sigma = 0;
  }
  else {
    sigma = num / denum;
  }

  return _sigma_mult_stat * sigma;

}


/*-------------------------------------------*/
/*       Compute sigma_kl for objective      */
/*           and smooth uncertainty          */
/*-------------------------------------------*/
double SGTELIB::Surrogate_Ensemble_Stat::compute_sigma_kl_obj_smooth(const int i, const int k, const int l,
                                                                     const std::vector<SGTELIB::Matrix *> & SGradsurr ) const {
  double cos;

  cos = ( SGradsurr[k]->get_row(i) * SGradsurr[l]->get_row(i).transpose() ).get(0,0) / ( SGradsurr[k]->get_row(i).norm() * SGradsurr[l]->get_row(i).norm() );
  
  if (cos == cos) {
    return 0.5 * (1 - cos);
  }
  else { // avoid nan case 
    return 0.5;
  }

}


/*-------------------------------------------*/
/*       Compute sigma_kl for objective      */
/*          and nonsmooth uncertainty        */
/*-------------------------------------------*/
double SGTELIB::Surrogate_Ensemble_Stat::compute_sigma_kl_obj_nonsmooth(const int i, const int j, const int k, const int l,
                                                                        const std::vector<SGTELIB::Matrix *> & ZZsurr,
                                                                        const std::vector<SGTELIB::Matrix *> & ZZsurr_around ) const {
  
  double res = 0;

  // Loop on the _nbd directions of the positive spanning set around x
  for (int d=0 ; d<_nbd ; d++){
    res +=  ( ZZsurr_around[k]->get(i,d) < ZZsurr[k]->get(i,j) ) != ( ZZsurr_around[l]->get(i,d) < ZZsurr[l]->get(i,j) ) ; // xor
  }

  res /= _nbd;

  return res;
}


/*-------------------------------------------*/
/*       Compute sigma_kl for constraint     */
/*           and smooth uncertainty          */
/*-------------------------------------------*/
double SGTELIB::Surrogate_Ensemble_Stat::compute_sigma_kl_con_smooth(const int i, const int j, const int k, const int l,
                                                                     const std::vector<SGTELIB::Matrix *> & ZZsurr) const {  
  double cs = _trainingset.Z_scale(0.0,j);
  return sigmoid ( - ( ZZsurr[k]->get(i,j) - cs ) * ( ZZsurr[l]->get(i,j) - cs ) , 1 );
}


/*-------------------------------------------*/
/*       Compute sigma_kl for constraint     */
/*          and nonsmooth uncertainty        */
/*-------------------------------------------*/
double SGTELIB::Surrogate_Ensemble_Stat::compute_sigma_kl_con_nonsmooth(const int i, const int j, const int k, const int l,
                                                                        const std::vector<SGTELIB::Matrix *> & ZZsurr) const {
  double cs = _trainingset.Z_scale(0.0,j);
  return ( ( ZZsurr[k]->get(i,j) <= cs ) != ( ZZsurr[l]->get(i,j) <= cs ) ); // xor
}

 
/*--------------------------------------*/
/*       get_matrix_Zvs                 */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_Ensemble_Stat::get_matrix_Zvs (void){
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
const SGTELIB::Matrix * SGTELIB::Surrogate_Ensemble_Stat::get_matrix_Zhs (void){
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
const SGTELIB::Matrix * SGTELIB::Surrogate_Ensemble_Stat::get_matrix_Shs (void){
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
bool SGTELIB::Surrogate_Ensemble_Stat::is_ready (const int k) const{
  if ((k<0) || (k>=_kmax)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (const int k): k out of range" );
  }
  return _surrogates.at(k)->is_ready();
}


/*--------------------------------------*/
/*  external set of the weight vector   */
/*    (use model k for output j)        */
/*--------------------------------------*/
/*
void SGTELIB::Surrogate_Ensemble_Stat::set_weight_vector (const int k, const int j){
  if (_param.get_weight_type() != SGTELIB::WEIGHT_EXTERN){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): Not in EXTERN mode" );
  }
  if ((k<0) or (k>=_kmax)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): k out of range" );
  }
  if ((j<0) or (j>=_m)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): k out of range" );
  }
  if (not is_ready(k)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): Surrogate not ready" );
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
void SGTELIB::Surrogate_Ensemble_Stat::set_weight_vector (const int k){
  if (_param.get_weight_type() != SGTELIB::WEIGHT_EXTERN){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): Not in EXTERN mode" );
  }
  if ((k<0) or (k>=_kmax)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): k out of range" );
  }
  if (not is_ready(k)){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): Surrogate not ready" );
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
void SGTELIB::Surrogate_Ensemble_Stat::set_weight_vector (const SGTELIB::Matrix & W){
  if (_param.get_weight_type() != SGTELIB::WEIGHT_EXTERN){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
               "Surrogate_Ensemble_Stat::set_weight_vector (k,j): Not in EXTERN mode" );
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
bool SGTELIB::Surrogate_Ensemble_Stat::check_weight_vector ( void ) const {
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
      if (std::fabs(s-1.0)>_kready*EPSILON) return true;
    }
  }

  return false;

}//






/*--------------------------------------*/
/*        define model list             */
/*--------------------------------------*/
void SGTELIB::Surrogate_Ensemble_Stat::model_list_preset ( const std::string & preset ) {


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
        "Surrogate_Ensemble_Stat::model_list_preset: unrecognized preset \""+preset+"\"" );
    }  

    #ifdef ENSEMBLE_DEBUG
      std::cout << "END Build model list\n";
    #endif

}//

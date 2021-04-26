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

#include "Surrogate.hpp"

using namespace SGTELIB;

/*--------------------------------------*/
/*              constructor             */
/*--------------------------------------*/

SGTELIB::Surrogate::Surrogate ( SGTELIB::TrainingSet & trainingset,
                                const SGTELIB::Surrogate_Parameters& param) :
  // set of data used to build the model 
  _trainingset ( trainingset      ) ,
  // set of parameters used to build the model
  _param       ( param            ) ,   
  // n : dimension of the input space
  _n           (_trainingset.get_input_dim()  ) ,
  // m : dimension of the output space
  _m           (_trainingset.get_output_dim() ) ,
  // number of points in the training set
  _p_ts        (0                 ) ,
  // number of points in the training set last time the model was built
  _p_ts_old    (999999999         ) ,
  // number of points considered for the model (after possible filtering of the TS)
  _p           (0                 ) ,
  // value of _p last time the model was built
  _p_old       (999999999         ) ,
  // Boolean indicating if the model is ready
  _ready       (false             ) ,
  // Scaled value of the model output at the training points 
  _Zhs         (NULL              ) ,
  // Scaled value of the model variance at the training points
  _Shs         (NULL              ) ,
  // Scaled leave-one-out cross-validation value at the training points
  _Zvs         (NULL              ) ,
  // Scaled leave-one-out cross-validation variance at the training points
  _Svs         (NULL              ) ,
  // Set of training points used for building the model (default : all the points)
  // Note that in sgtelib 2.0.1, there is not method implemented to filter the data points.
  // So all the data points are used.
  _selected_points (1,-1          ) ,
  // Map containing all the metrics
  _metrics         (              ) ,
  // Poll size max (used during the parameter optimization with the MADS algorithm)
  // This value is returned by the parameter optimization and is used for the initialization
  // of the next parameter optimization.
  _psize_max       ( 0.5          ) ,
  // output stream
  _out             (              ) ,
  // Boolean for display
  _display         ( false        ) {;
}//


SGTELIB::Surrogate::Surrogate ( SGTELIB::TrainingSet & trainingset,
                                const SGTELIB::model_t& mt ) :
  _trainingset ( trainingset      ) ,
  _param       ( mt ) ,   
  _n     (_trainingset.get_input_dim()  ) ,
  _m     (_trainingset.get_output_dim() ) ,
  _p_ts      (0                   ) ,
  _p_ts_old  (999999999           ) ,
  _p         (0                   ) ,
  _p_old     (999999999           ) ,
  _ready (false                   ) ,
  _Zhs   (NULL                    ) ,
  _Shs   (NULL                    ) ,
  _Zvs   (NULL                    ) ,
  _Svs   (NULL                    ) ,
  _selected_points (1,-1          ) ,
  _metrics         (              ) ,
  _psize_max       ( 0.5          ) ,
  _out             (              ) ,
  _display         ( false        ) {
}//

SGTELIB::Surrogate::Surrogate ( SGTELIB::TrainingSet & trainingset,
                                const std::string & s) :
  _trainingset ( trainingset      ) ,
  _param       ( s                ) ,   
  _n     (_trainingset.get_input_dim()  ) ,
  _m     (_trainingset.get_output_dim() ) ,
  _p_ts      (0                   ) ,
  _p_ts_old  (0                   ) ,
  _p         (0                   ) ,
  _p_old     (0                   ) ,
  _ready (false                   ) ,
  _Zhs   (NULL                    ) ,
  _Shs   (NULL                    ) ,
  _Zvs   (NULL                    ) ,
  _Svs   (NULL                    ) ,
  _selected_points (1,-1          ) ,
  _metrics         (              ) ,
  _psize_max       ( 0.5          ) ,
  _out             (              ) ,
  _display         ( false        ) {
}//


/*--------------------------------------*/
/*               destructor             */
/*--------------------------------------*/
SGTELIB::Surrogate::~Surrogate ( void ) {
  reset_metrics();
}//

void SGTELIB::Surrogate::info ( void ) const {
  _trainingset.info();
}//


/*--------------------------------------*/
/*              display                 */
/*--------------------------------------*/
void SGTELIB::Surrogate::display ( std::ostream & out ) const {
  out << "Surrogate: " << get_string() << "\n";
  out << "ready: " << _ready << "\n";
  out << "n: " << _n << " (input dim)\n";
  out << "m: " << _m << " (output dim)\n";
  out << "p: " << _p << " (nb points)\n";
  out << "Metrics:\n";
  std::map< metric_t , SGTELIB::Matrix >::const_iterator it;
  for (it=_metrics.begin();it!=_metrics.end();it++){
    SGTELIB::Matrix V = it->second;
    out << "  " << SGTELIB::metric_type_to_str(it->first) << " = [ ";
    for (int j=0;j<V.get_nb_cols();j++) out << V[j] << " ";
    out << "]\n";
  }
  display_private ( out );
}//

/*--------------------------------------*/
/*       erase_data                     */
/*--------------------------------------*/
void SGTELIB::Surrogate::reset_metrics ( void ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate: reset_metrics...";
  #endif

  if (_Zhs) delete _Zhs;
  _Zhs = NULL;  

  if (_Shs) delete _Shs;
  _Shs = NULL;  

  if (_Zvs) delete _Zvs;
  _Zvs = NULL;  

  if (_Svs) delete _Svs;
  _Svs = NULL;  

  _metrics.clear();

  #ifdef SGTELIB_DEBUG
    std::cout << "OK\n";
  #endif
}//

/*--------------------------------------*/
/*               build                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate::build ( void ) {

  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate build - BEGIN\n";
  #endif

  if (streqi(_param.get_output(),"NULL")){
    _display = false;
  }
  else{
    _display = true;
  } 

  // Check the parameters of the model:
  _param.check();

  // Before building the surrogate, the trainingset must be ready
  _trainingset.build();

   if (!_trainingset.is_ready()) {
     _ready = false;
     return false;
   }

  // Number of points in the training set.
  _p_ts = _trainingset.get_nb_points();
  //std::cout << _ready << " " << _p_ts << " " << _p_ts_old << "\n";
  if ( (_ready) && (_p_ts==_p_ts_old) ){
    #ifdef SGTELIB_DEBUG
      std::cout << "Surrogate build - SKIP Build\n";
    #endif
    return true;
  }
  
  // Otherwise, the model is not ready and we need to call build_private
  _ready = false;

  // Get the number of points used in the surrogate
  if ( (_selected_points.size()==1) && (_selected_points.front()==-1) )
    _p = _p_ts;
  else  
    _p = static_cast<int>(_selected_points.size());

  // Need at least 2 point to build a surrogate.
  if (_p<2) return false;

  // Delete the intermediate data and metrics 
  // (they will have to be recomputed...)
  reset_metrics();

  // Call to the private build
  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate build - BUILD_PRIVATE\n";
  #endif

  bool ok;

  // First, the model has to be initialized.
  // This step does not involve parameter optimization. 
  // For some types of model, the initialization step does nothing.
  // The initialization step is necessary, for example, for RBF models, where the "preset"
  // has to be considered first, and the kernel have to be selected before the parameter
  // optimization.
  ok = init_private();
  if ( ! ok ) return false;

  #ifdef SGTELIB_DEBUG
    std::cout << "Number of parameters to optimize : " << _param.get_nb_parameter_optimization() << "\n";
    _param.display(std::cout);
  #endif


  // Optimize parameters
  if (_param.get_nb_parameter_optimization()>0){
    ok = optimize_parameters();
    if ( ! ok ){
      _ready = false;
      return false;
    }
  }

  // Build private
  ok = build_private();
  if ( ! ok ){
    _ready = false;
    return false;
  }

  // Memorize previous number of points
  _p_ts_old = _p_ts;
  _p_old = _p;

  #ifdef SGTELIB_DEBUG
    std::cout << "Surrogate build - END\n";
  #endif
  
  if (_display){
    _out.open(_param.get_output().c_str() , std::ios::out | std::ios::app);
    if (_out.fail()) std::cout << "Out.fail1!!!\n";
    std::cout << "Write in " << _param.get_output() << "\n";
    if (_out.fail()) std::cout << "Out.fail2!!!\n";
    display(_out);
    if (_out.fail()) std::cout << "Out.fail3!!!\n";
    _out.close();
  }

  _ready = true;
  return true;
}//

bool SGTELIB::Surrogate::init_private (void) {
  // Empty initialization function
  #ifdef SGTELIB_DEBUG
    std::cout << model_type_to_str(get_type()) << " : init_private\n";
  #endif
  return true;
}


/*--------------------------------------*/
/*               check_ready            */
/*--------------------------------------*/
void SGTELIB::Surrogate::check_ready (void) const {
    check_ready("");
}//

/*--------------------------------------*/
void SGTELIB::Surrogate::check_ready (const std::string & file,
                                      const std::string & function,
                                      const int & i        ) const {
    check_ready(file+"::"+function+"::"+itos(i));
}//
/*--------------------------------------*/
void SGTELIB::Surrogate::check_ready (const std::string & s) const {
  
  // Check the tag _ready
  if ( ! _ready){
    display(std::cout);
    std::cout << "Surrogate: NOT READY! (" << s << ")\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                 "check_ready(): Not ready!" );
  }

  // Check if the trainingset is ready
  _trainingset.check_ready("From Surrogate ()");


  // Check the new number of points in the trainingset
  if (_trainingset.get_nb_points()>_p_ts){
    display(std::cout);
    std::cout << "Surrogate: NOT READY! (" << s << ")\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                 "check_ready(): Not ready!" );
  }

}//


/*--------------------------------------*/
/*               add points             */
/*--------------------------------------*/
bool SGTELIB::Surrogate::add_points ( const SGTELIB::Matrix & Xnew ,
                                      const SGTELIB::Matrix & Znew  ){
  // It would be possible to add points via the surrogate model, but it is considered
  // bad practice. So it is not allowed. Points have to be added directly via the training set.
  throw SGTELIB::Exception ( __FILE__ , __LINE__ , "add_points: forbiden." );
  return _trainingset.add_points(Xnew,Znew);
}//
/*--------------------------------------*/
bool SGTELIB::Surrogate::add_point  ( const double * xnew ,
                                      const double * znew  ){
  throw SGTELIB::Exception ( __FILE__ , __LINE__ , "add_point: forbiden." );
  return _trainingset.add_point(xnew,znew);
}//




/*=========================================================*/
/*=========================================================*/
/*||                                                     ||*/
/*||             PREDICTION METHODS                      ||*/
/*||                                                     ||*/
/*=========================================================*/
/*=========================================================*/



/*---------------------------------------------------------------------*/
/*               predict                                               */
/* XX : set of points where a prediction must be performed             */
/* ZZ : value of the model in XX                                       */
/* std : standard deviation of the model in XX                         */
/* ei : expected improvement of the model in XX                        */
/* cdf : probability that y(x) < y_0 for each point x of XX, where     */
/*       y_0 = f_min for the objective function                        */
/*       y_0 = 0 for constraint functions                              */
/*---------------------------------------------------------------------*/
void SGTELIB::Surrogate::predict ( const SGTELIB::Matrix & XX ,
                                         SGTELIB::Matrix * ZZ ,
                                         SGTELIB::Matrix * std, 
                                         SGTELIB::Matrix * ei ,
                                         SGTELIB::Matrix * cdf) {

  // Prediction requires that the model is ready.
  check_ready(__FILE__,__FUNCTION__,__LINE__);


  // Check the number of columns in XX
  if (XX.get_nb_cols() != _n){
    
    display(std::cout);
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                 "predict(): dimension error" );
  }

  // Create the (non-scaled) output matrix
  *ZZ = SGTELIB::Matrix("ZZ",XX.get_nb_rows(),_m);

  // Scale the input (XX ---> XXs)
  SGTELIB::Matrix XXs(XX);
  XXs.set_name("XXs");
  _trainingset.X_scale(XXs);

  // Init the Expected Improvement
  if (ei){
    ei->fill(-INF);
  }

  // Call the private prediction with normalized input XXs.
  // This will return normalized values of ZZ, std and ei.
  // cdf is never normalized.
  predict_private( XXs , ZZ , std , ei , cdf );

  // If nbdiff==1, put the values to 0.0
  int pxx = XX.get_nb_rows();
  if (ZZ){
    for (int j=0 ; j<_m ; j++){
      if (_trainingset.get_Z_nbdiff(j)==1){
        for (int i=0 ; i<pxx ; i++){
          ZZ->set(i,j,0.0);
        }
      }
    }
  }


  #ifdef SGTELIB_DEBUG
    if (ZZ){
      if (ZZ->has_nan()){
        ZZ->replace_nan (+INF);
      }
    }
    if (std){
      if (std->has_nan()){
        display(std::cout); 
        throw SGTELIB::Exception ( __FILE__ , __LINE__ , "predict(): std has nan" );
      }
    }
    if (ei){
      if (ei->has_nan()){
        display(std::cout); 
        throw SGTELIB::Exception ( __FILE__ , __LINE__ , "predict(): ei has nan" );
      }
    }
    if (cdf){
      if (cdf->has_nan()){
        display(std::cout); 
        throw SGTELIB::Exception ( __FILE__ , __LINE__ , "predict(): cdf has nan" );
      }
    }
  #endif

  ZZ->replace_nan (+INF);
  std->replace_nan (+INF);
  ei->replace_nan (-INF);
  cdf->replace_nan (0);

  // UnScale the output
  // Note that ZZ is unscaled with Z_unscale:
  // ZZ_unscaled = ( ZZ_scaled - b ) / a
  // while std and ei 
  // are unscaled with ZE_unscale (without the additive constant):
  // ZZ_unscaled = ZZ_scaled / a
  if (ZZ){
    ZZ->set_name("ZZ");   
    _trainingset.Z_unscale(ZZ);
  }
  if (std){
    std->set_name("std");
    _trainingset.ZE_unscale(std);
  }
  if (ei){
    ei->set_name("ei");
    _trainingset.ZE_unscale(ei);
    // ei is only computed for the OBJ output, so the other values are dummy, 
    // So we put them all to 0.    
    for (int j=0 ; j<_m ; j++){ 
      if (_trainingset.get_bbo(j)!=SGTELIB::BBO_OBJ){
        for (int i=0 ; i<pxx ; i++){ 
          ei->set(i,j,0.0);
        }
      }
    }  
  }
  if (cdf){
    // no unscaling for cdf because this is a probability.
    cdf->set_name("cdf");
  }


}//



/*--------------------------------------*/
/*       predict (ZZs,std,ei)           */
/*--------------------------------------*/
// This function is the default method to compute std, ei and cdf.
// It can be overloaded, but models PRS, RBF and KS use the default method
// because these methods are not able to compute std, ei and cdf, so they use 
// the default proxies.
// This method relies on the private method predict_private(XXs,ZZs)
// which HAS TO be overloaded (pure virtual)

// In other words, PRS, RBF and KS use the overloaded method to compute ZZ, and the default method 
// to compute std, ei and cdf.
// Kriging models use an overloaded method for ZZ and for std, ei and cdf.

// The following method receives scaled inputs (XXs) and returns scaled outputs (ZZs, std and ei)
void SGTELIB::Surrogate::predict_private (const SGTELIB::Matrix & XXs,
                                                SGTELIB::Matrix * ZZs,
                                                SGTELIB::Matrix * std, 
                                                SGTELIB::Matrix * ei ,
                                                SGTELIB::Matrix * cdf) {
  check_ready(__FILE__,__FUNCTION__,__LINE__);


  const int pxx = XXs.get_nb_rows();
  // Scaled value of f_min.
  const double fs_min = _trainingset.get_fs_min();
  int i,j;

  // Prediction of ZZs
  if ( (ZZs) || (ei) || (cdf) ){
    predict_private(XXs,ZZs);
  }

  // Prediction of statistical data
  if ( (std) || (ei) || (cdf) ){

    if (std) std->fill(-SGTELIB::INF);
    else std = new SGTELIB::Matrix("std",pxx,_m);

    if (ei)   ei->fill(-SGTELIB::INF);
    if (cdf) cdf->fill(-SGTELIB::INF);

    // Use normalized distance to closest and rmse as std
    SGTELIB::Matrix dtc = _trainingset.get_distance_to_closest(XXs);
    dtc.set_name("dtc");

    for (j=0 ; j<_m ; j++){
      // Set std (use a proxy)
      double s = get_metric(SGTELIB::METRIC_RMSE,j); 
      std->set_col( dtc+s , j );

      if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ){
        // Compute CDF
        if (cdf){
          for (i=0 ; i<pxx ; i++){
            cdf->set(i,j, normcdf( fs_min , ZZs->get(i,j) , std->get(i,j) ) );
          }
        }
        if (ei){
          for (i=0 ; i<pxx ; i++){
            ei->set(i,j, normei( ZZs->get(i,j) , std->get(i,j) , fs_min ) );
          }
        }
      }// END CASE OBJ
      else if (_trainingset.get_bbo(j)==SGTELIB::BBO_CON){
        // Compute CDF
        if (cdf){
          // Scaled Feasibility Threshold
          double cs = _trainingset.Z_scale(0.0,j);
          for (i=0 ; i<pxx ; i++){
            cdf->set(i,j, normcdf( cs , ZZs->get(i,j) , std->get(i,j) ) );
          }
        }
      }// END CASE CON

    }// End for j

  }
}//






/*--------------------------------------*/
/*               predict                */
/*--------------------------------------*/
void SGTELIB::Surrogate::predict ( const SGTELIB::Matrix & XX ,
                                         SGTELIB::Matrix * ZZ ) {

  check_ready(__FILE__,__FUNCTION__,__LINE__);



  // Check the number of columns in XX
  if (XX.get_nb_cols() != _n){
    display(std::cout); 
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                 "predict(): dimension error" );
  }
  *ZZ = SGTELIB::Matrix("ZZ",XX.get_nb_rows(),_m);

  // Scale the input
  SGTELIB::Matrix XXs(XX);
  _trainingset.X_scale(XXs);


  // Call the private prediction with normalized input XXs
  predict_private( XXs , ZZ );
  #ifdef SGTELIB_DEBUG
    if (ZZ->has_nan()){
      display(std::cout); 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                   "predict(): ZZ has nan" );
    }
  #endif

  // UnScale the output
  _trainingset.Z_unscale(ZZ);

}//





/*=========================================================*/
/*=========================================================*/
/*||                                                     ||*/
/*||                  GET MATRICES                       ||*/
/*||                                                     ||*/
/*=========================================================*/
/*=========================================================*/


/*---------------------------------------*/
/*       compute matrix Zhs              */
/* Zhs is the prediction on the training */
/* points                                */
/*---------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate::get_matrix_Zhs (void){
  if ( ! _Zhs){
    check_ready(__FILE__,__FUNCTION__,__LINE__);

    // Init
    _Zhs = new SGTELIB::Matrix("Zhs",_p,_m);
    //call the predict function on the training points
    predict_private (get_matrix_Xs(),_Zhs);
    _Zhs->replace_nan(+INF);
    _Zhs->set_name("Zhs");
  }
  return _Zhs;
}//


/*--------------------------------------*/
/*       compute matrix Shs             */
/*  (Compute the predictive std)        */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate::get_matrix_Shs (void){
  if ( ! _Shs){
    check_ready(__FILE__,__FUNCTION__,__LINE__);

    #ifdef SGTELIB_DEBUG
      std::cout << "Compute _Shs\n";
    #endif
    // Init
    _Shs = new SGTELIB::Matrix("Shs",_p,_m);
    //call the predict function on the training points
    predict_private (get_matrix_Xs(),NULL,_Shs,NULL,NULL);
    _Shs->replace_nan(+INF);
    _Shs->set_name("Shs");
  }
  return _Shs;
}//

// If no specific method is defined, consider Svs = Shs.
const SGTELIB::Matrix * SGTELIB::Surrogate::get_matrix_Svs (void){
  if ( ! _Svs){

    _Svs = new SGTELIB::Matrix("Svs",_p,_m);
    const SGTELIB::Matrix Ds = _trainingset.get_matrix_Ds();
    for (int i=0 ; i<_p ; i++){
      double dmin = +INF;
      for (int j=0 ; j<_p ; j++){
        if (i!=j){
          dmin = std::min(dmin,Ds.get(i,j));
        }
      }
      _Svs->set_row(dmin,i);
    }
  }
  return _Svs;
}//



/*--------------------------------------*/
/*       get_Xs                         */
/* Returns the scaled input for all     */
/* the selected data points selected    */
/*--------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate::get_matrix_Xs (void){
  _trainingset.build(); 
  return _trainingset.get_matrix_Xs().get_rows(_selected_points);
}//


/*--------------------------------------*/
/*       get_Zs                         */
/* Returns the scaled output for all    */
/* the selected data points selected    */
/*--------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate::get_matrix_Zs (void){
  _trainingset.build(); 
  return _trainingset.get_matrix_Zs().get_rows(_selected_points);
}//


/*-----------------------------------------------*/
/*       get_Ds                                  */
/* Ds is provided by the training set and        */
/* contains the scaled distance between any pair */
/* of data points                                */
/*-----------------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate::get_matrix_Ds (void){
  _trainingset.build(); 
  return _trainingset.get_matrix_Ds().get( _selected_points , _selected_points );
}//





/*-----------------------------------------*/
/*       get_Zh                            */
/* Zh is the value of the model at the     */
/* data points (the "h" in "Zh" stands for */
/* "hat", which is the common notation for */
/* a model                                 */
/*-----------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate::get_matrix_Zh (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  SGTELIB::Matrix Zh (*get_matrix_Zhs()); // Get scaled matrix
  _trainingset.Z_unscale(&Zh); // Unscale
  return Zh; // Return unscaled
}//



/*-------------------------------------------------*/
/*       get_Zv                                    */
/* Zv contains the leave-one-out cross-validation  */
/* values at the training points                   */
/*-------------------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate::get_matrix_Zv (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  SGTELIB::Matrix Zv (*get_matrix_Zvs()); // Get scaled matrix
  _trainingset.Z_unscale(&Zv); // Unscale
  return Zv; // Return unscaled
}//


/*--------------------------------------*/
/*       get_Sh                         */
/*--------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate::get_matrix_Sh (void){
  // Return unscaled matrix Shs
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  SGTELIB::Matrix Sh = (*get_matrix_Shs());
  _trainingset.ZE_unscale(&Sh); // Unscale (without additive constant)
  return Sh; // Return unscaled
}//

/*--------------------------------------*/
/*       get_Sv                         */
/*--------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate::get_matrix_Sv (void){
  // Return unscaled matrix Sv
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  SGTELIB::Matrix Sv (*get_matrix_Svs()); // Get scaled matrix
  _trainingset.ZE_unscale(&Sv); // Unscale (without additive constant)
  return Sv; // Return unscaled
}//





/*=========================================================*/
/*=========================================================*/
/*||                                                     ||*/
/*||                   METRICS                           ||*/
/*||                                                     ||*/
/*=========================================================*/
/*=========================================================*/



/*---------------------------------------*/
/*  check if the metric is defined       */
/*---------------------------------------*/
bool SGTELIB::Surrogate::is_defined(const SGTELIB::metric_t mt){
  // Check if the key exists in the map
  if (_metrics.find(mt)==_metrics.end()) return false;
  // Check the size of the vector
  const int metric_vector_size = _metrics[mt].get_nb_cols();
  if (metric_vector_size<=0) return false;
  return true;
}//
/*---------------------------------------*/
bool SGTELIB::Surrogate::is_defined(const SGTELIB::metric_t mt, const int j){
  if (!is_defined(mt)) return false;
  if (   (j>=_metrics[mt].get_nb_cols()) || (j>=_m) || (j<0)    ) return false;
  return true;
}//


/*--------------------------------------*/
/*       compute metric                 */
/*--------------------------------------*/
bool SGTELIB::Surrogate::compute_metric ( const metric_t mt ){

  if (is_defined(mt)) return true;

  double m;
  int j;

  // Choose if we use the Zhs or the Zvs matrix
  // Zvs is used if we want to use cross-validation.
  const SGTELIB::Matrix Zs = get_matrix_Zs();
  const SGTELIB::Matrix * Zs_compare;
  const SGTELIB::Matrix * Ss_compare;
  if (SGTELIB::metric_uses_cv(mt)){
    Zs_compare = get_matrix_Zvs(); 
    Ss_compare = get_matrix_Svs(); 
  }
  else{
    Zs_compare = get_matrix_Zhs(); 
    Ss_compare = get_matrix_Shs(); 
  }

  // Size of the metric vector
  const int vector_size = (SGTELIB::one_metric_value_per_bbo(mt))?_m:1;
  // Init the metric vector 
  SGTELIB::Matrix v ("v",1,vector_size);
  // Norm associated to a given metric.
  norm_t associated_norm;

  switch (mt){
    case SGTELIB::METRIC_EMAX:
    case SGTELIB::METRIC_EMAXCV:
    case SGTELIB::METRIC_RMSE:
    case SGTELIB::METRIC_RMSECV:
    case SGTELIB::METRIC_ARMSE:
    case SGTELIB::METRIC_ARMSECV:
      // Get the norm associated with this metric
      associated_norm = SGTELIB::metric_type_to_norm_type(mt);
      // Compute the norm of the difference
      v = (Zs-(*Zs_compare)).col_norm( associated_norm );
      if (  (mt==SGTELIB::METRIC_ARMSE) || (mt==SGTELIB::METRIC_ARMSECV)  ){
        // For "Aggregate" metrics, compute the sum for all BBO
        v = v.sum(1);
      }
      else{
        // Otherwise, unscale
        _trainingset.ZE_unscale(v);
      }

      break;

    case SGTELIB::METRIC_OE:
    case SGTELIB::METRIC_OECV:
      // Order error. See paper: 
      // Order-based error for managing ensembles of surrogates in mesh adaptive direct search
      v = compute_order_error(*Zs_compare);
      break;

    case SGTELIB::METRIC_AOE:
    case SGTELIB::METRIC_AOECV:
      // Agregate order error. See paper: 
      //Locally weighted regression models for surrogate-assisted design optimization 
      v = SGTELIB::Matrix( compute_aggregate_order_error(*Zs_compare) );
      break;

    case SGTELIB::METRIC_EFIOE:
    case SGTELIB::METRIC_EFIOECV:
      // Agregate Order error on Expected Feasible Improvement
      v = SGTELIB::Matrix( compute_aggregate_order_error( -compute_efi(*Zs_compare,*Ss_compare) ) );
      break;

    case SGTELIB::METRIC_LINV:
       // Inverse of the likelihood (this method can be overloaded)
      compute_metric_linv();
      break;

    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Metric not recognized." );
  }

  // Check bounds.
  for (j=0; j<vector_size ; j++){
    m = v[j];
    if (isnan(m)    ){ m = SGTELIB::INF; }
    if (m < -EPSILON){ m = SGTELIB::INF; }
    if (m <= 0.0    ){ m = 0.0; }
    v.set(0,j,m);
  }

  _metrics[mt] = v;

  return true;
}//



/*--------------------------------------*/
/*       get metric (general)           */
/*--------------------------------------*/
double SGTELIB::Surrogate::get_metric (SGTELIB::metric_t mt , int j){
  // If the model is not ready, return +INF
  if (!_ready) return SGTELIB::INF; 
  // If the metric is defined, return it
  if ( is_defined(mt,j) ) return _metrics[mt][j];
  // Compute the metric, 
  if ( !compute_metric(mt) ) return SGTELIB::INF;
  // Return value
  #ifdef SGTELIB_DEBUG
    std::cout << "metric " << SGTELIB::metric_type_to_str(mt) << "[" << j << "]";
    if ( is_defined(mt,j) ) std::cout << " is def: " << _metrics[mt][j] << std::endl;
    else std::cout << " NOT defined." << std::endl;
  #endif
  if ( is_defined(mt,j) ) return _metrics[mt][j];
  // Is still not defined, return INF.
  return SGTELIB::INF;
}//

SGTELIB::Matrix SGTELIB::Surrogate::get_metric (SGTELIB::metric_t mt){
  // If the model is not ready, return +INF
  if (!_ready) return SGTELIB::Matrix(SGTELIB::INF);
  // If the metric is defined, return it
  if ( is_defined(mt) ) return _metrics[mt];
  // Compute the metric, 
  if ( !compute_metric(mt) ) return SGTELIB::Matrix(SGTELIB::INF);
  // Return value
  if ( is_defined(mt) ) return _metrics[mt];
  // Is still not defined, return INF.
  return SGTELIB::Matrix(SGTELIB::INF);
}//





/*----------------------------------------------------------*/
/*     compute EFI from the predictive mean and std         */
/*----------------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate::compute_efi( const SGTELIB::Matrix & Zs,
                                                 const SGTELIB::Matrix & Ss  ){

  if (  (Zs.get_nb_cols()!=_m) || 
        (Ss.get_nb_cols()!=_m) || 
        (Zs.get_nb_rows()!=_p) || 
        (Ss.get_nb_rows()!=_p)    ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Dimension error" );
  }

  const double fmin = _trainingset.get_f_min();
  double c0, ei;

  SGTELIB::Matrix EFI ("EFI",_p,1);
  EFI.fill(1.0);

  for (int j=0 ; j<_m ; j++){
    if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ){
      for (int i=0 ; i<_p ; i++){
        // Compute Expected Improvement for point i
        ei = SGTELIB::normei(Zs.get(i,j),Ss.get(i,j),fmin);
        // Unscale Expected Improvement
        ei = _trainingset.ZE_unscale(ei,j);
        // Multiply EFI by ei
        EFI.product(i,0,ei);
      }
    }
    else if (_trainingset.get_bbo(j)==SGTELIB::BBO_CON){
      c0 = _trainingset.Z_scale(0.0,j);
      for (int i=0 ; i<_p ; i++) EFI.product(i,0, SGTELIB::normcdf(c0,Zs.get(i,j),Ss.get(i,j)) );
    }
  }
  return EFI;
}//


/*--------------------------------------*/
/*       compute linv                   */
/*--------------------------------------*/
void SGTELIB::Surrogate::compute_metric_linv (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  if ( !is_defined(SGTELIB::METRIC_LINV) ){
    #ifdef SGTELIB_DEBUG
      std::cout << "Compute metric linv\n";
    #endif
    // Init
    SGTELIB::Matrix v = SGTELIB::Matrix("v",1,_m);

    // Compute the prediction on the training points
    const SGTELIB::Matrix * Zhs = get_matrix_Zhs();
    const SGTELIB::Matrix * Shs = get_matrix_Shs();
    // True values
    const SGTELIB::Matrix Zs = get_matrix_Zs();
    double s,dz;
    double linv;
    // TODO : improve the behavior of linv for very small s.
    for (int j=0 ; j<_m ; j++){
      if (_trainingset.get_bbo(j)!=SGTELIB::BBO_DUM){
        linv = 0;
        for (int i=0 ; i<_p ; i++){
          dz = Zhs->get(i,j)-Zs.get(i,j);
          s  = Shs->get(i,j);
          s = std::max(s ,EPSILON);   
          dz= std::max(dz,EPSILON); 
          linv += -log(s) - pow(dz/s,2)/2;
        }
        linv /= _p; // normalization by the number of points
        linv -= 0.5*log(2*3.141592654); // add the normal pdf constant
        // Add this point, we have log(prod g)/p
        linv = exp(-linv);
      }
      else{
        linv = -1;
      }
      v.set(0,j,linv);
    }
    _metrics[SGTELIB::METRIC_LINV] = v;
  }

}//


/*--------------------------------------*/
/*       compute order efficiency       */
/*--------------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate::compute_order_error ( const SGTELIB::Matrix & Zpred ){

  // Compute the order-efficiency metric by comparing the 
  // values of - _Zs (in the trainingset)
  //           - Zpred (input of this function)
  // Put the results in "OE" (output of this function)

  SGTELIB::Matrix OE = SGTELIB::Matrix("OE",1,Zpred.get_nb_cols());


  int nb_fail;
  double z1,z1h,z2,z2h;
  double c0;
  const SGTELIB::Matrix Zs = get_matrix_Zs();
  
  for (int j=0 ; j<_m ; j++){
    switch (_trainingset.get_bbo(j)){
    //===============================================//
    case SGTELIB::BBO_OBJ:
      nb_fail = 0;
      for (int i1=0 ; i1<_p ; i1++){
        z1 = Zs.get(i1,j);
        z1h = Zpred.get(i1,j);
        for (int i2=0 ; i2<_p ; i2++){
          z2  = Zs.get(i2,j);
          z2h = Zpred.get(i2,j);
          if ( (z1-z2<0)^(z1h-z2h<0) ) nb_fail++;
        }
      }
      OE.set(0,j, double(nb_fail)/double(_p*_p) );
      break;
    //===============================================//
    case SGTELIB::BBO_CON:
      nb_fail = 0;
      // Compute the feasibility threshold for scaled values
      c0 = _trainingset.Z_scale(0.0,j);
      for (int i=0 ; i<_p ; i++){
        z1  = Zs.get(i,j)    - c0;
        z1h = Zpred.get(i,j)- c0;
        if ( (z1<0)^(z1h<0) ) nb_fail++;
      }
      OE.set(0,j, double(nb_fail)/double(_p) );
      break;
    //===============================================//
    case SGTELIB::BBO_DUM:
      OE.set(0,j, -1 );
      break;
    //===============================================//
    default:
      display(std::cout); 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
    //===============================================//
    }// end switch
  }// end loop on j

  return OE;
}//


/*--------------------------------------*/
/*    compute f and h for a matrix Zs   */
/*--------------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate::compute_fh (const SGTELIB::Matrix & Zs){

  const int m = Zs.get_nb_cols();
  const int p = Zs.get_nb_rows();
  // First column: f
  // Second column: h
  SGTELIB::Matrix fh ("fh",p,2);
  fh.fill(0);

  if (m==1){
    fh.set_col(Zs,0);
    return fh;
  }
  else if (m==_m){
    int i,j;
    double d,c0;
    for (j=0 ; j<_m ; j++){
      switch (_trainingset.get_bbo(j)){
      //===============================================//
      case SGTELIB::BBO_OBJ:
        // Copy the objective in the first column of 
        fh.set_col( Zs.get_col(j) , 0 );
        break;
      //===============================================//
      case SGTELIB::BBO_CON:
        c0 = _trainingset.Z_scale(0.0,j);
        for (i=0 ; i<p ; i++){
          d = Zs.get(i,j) - c0;
          if (d>0) fh.add(i,1,d*d);
        }
        break;
      //===============================================//
      case SGTELIB::BBO_DUM:
        break;
      //===============================================//
      default:
        display(std::cout); 
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
      //===============================================//
      }// end switch
    }// end loop on j
  }
  else{
    Zs.display_short(std::cout);
    Zs.display_size(std::cout);
    std::cout << _m << " " << m << " " << _p << std::endl;
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Dimension error" );
  }
  return fh;
}//

/*--------------------------------------*/
/*       compute order efficiency       */
/*--------------------------------------*/
double SGTELIB::Surrogate::compute_aggregate_order_error (const SGTELIB::Matrix & Zpred){

  // Zpred must be a matrix with _p rows, and _m or 1 columns.
  // If there is only 1 column, then this column is considered as an aggregate of the
  // objective and constraints. For example, it can be the EFI.

  const SGTELIB::Matrix fhr = compute_fh( get_matrix_Zs() );
  const SGTELIB::Matrix fhs = compute_fh( Zpred           );

  int e = 0;
  int i1,i2;
  double hr1,hr2,hs1,hs2,fr1,fr2,fs1,fs2;
  bool inf_r,inf_s;
  // i1 and i2 are the indexes of the two points that are compared.
  // fr1 and hr1 (resp. fr2 and hr2) are the real values of f and r for these points.
  // fs1 and hs1 (resp. fs2 and hs2) are the surrogate (or CV) values.
  for (i1=0 ; i1<_p ; i1++){
    fr1 = fhr.get(i1,0);
    hr1 = fhr.get(i1,1);
    fs1 = fhs.get(i1,0);
    hs1 = fhs.get(i1,1);
    for (i2=0 ; i2<_p ; i2++){
      fr2 = fhr.get(i2,0);
      hr2 = fhr.get(i2,1);
      fs2 = fhs.get(i2,0);
      hs2 = fhs.get(i2,1);
      // Compute the order for real (r) data and for surrogate (s) model
      inf_r = ( (hr1<hr2) | ( (hr1==hr2) & (fr1<fr2) ) );
      inf_s = ( (hs1<hs2) | ( (hs1==hs2) & (fs1<fs2) ) );
      // If they don't agree, increment e. (Note that ^ is the xor operator)
      if (inf_r ^ inf_s) e++;
    }
  }
  return double(e)/double(_p*_p);

}//



/*=======================================*/
/*=======================================*/
/*||                                   ||*/
/*||          EXCLUSION AREA           ||*/
/*||                                   ||*/
/*=======================================*/
/*=======================================*/

/*--------------------------------------*/
/*       get_exclusion_area_penalty     */
/*--------------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate::get_exclusion_area_penalty ( const SGTELIB::Matrix & XX , const double tc ) const{
  // Scale the input
  SGTELIB::Matrix XXs(XX);
  XXs.set_name("XXs");
  _trainingset.X_scale(XXs);
  return _trainingset.get_exclusion_area_penalty ( XXs , tc );
}//


/*--------------------------------------*/
/*       get_distance_to_closest        */
/*--------------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate::get_distance_to_closest ( const SGTELIB::Matrix & XX ) const{
  // Scale the input
  SGTELIB::Matrix XXs(XX);
  XXs.set_name("XXs");
  _trainingset.X_scale(XXs);
  return _trainingset.get_distance_to_closest ( XXs );
}//



/*=======================================*/
/*=======================================*/
/*||                                   ||*/
/*||      PARAMETER OPTIMIZATION       ||*/
/*||                                   ||*/
/*=======================================*/
/*=======================================*/


/*--------------------------------------*/
/*  optimize model parameters           */
/*--------------------------------------*/
bool SGTELIB::Surrogate::optimize_parameters ( void ) {

  // Number of parameters to optimize
  const int N = _param.get_nb_parameter_optimization();
  // Budget
  int budget = N*_param.get_budget();

  int i,j,k;
  double d;
  const bool display = false;
  if (display){
    std::cout << "Begin parameter optimization\n";
    std::cout << "Metric: " << SGTELIB::metric_type_to_str(_param.get_metric_type()) << "\n";
  }

  //-----------------------------------------
  // Bounds, log-scale and domain
  //-----------------------------------------
  // Lower and upper bound of the parameter
  SGTELIB::Matrix lb("lb",1,N);
  SGTELIB::Matrix ub("ub",1,N);
  // Log-scale: if true, then the parameter must be positive and
  // will be optimized with a log-scale. This is equivalent to optimizing
  // the log of the parameter, instead of optimizing the parameter itself.
  // This is very interesting for parameters like the ridge coefficient,
  // which can take anywhere between 1e-16 and 1.
  bool * logscale = new bool [N];
  // The "domain" indicates if the parameter is continuous, integer, boolean,
  // categorical, or "MISC".
  // MISC parameter should not be optimized.
  SGTELIB::param_domain_t * domain = new SGTELIB::param_domain_t[N];

  // Interrogating the parameter instance.
  _param.get_x_bounds ( &lb , &ub , domain , logscale );

  //-----------------------------------------
  // Compute scaling
  // The scaling is necessary to compute the magnitude of the poll directions.
  //-----------------------------------------
  SGTELIB::Matrix scaling ("scaling",1,N);
  for (i=0 ; i<N ; i++){
    if (domain[i]==SGTELIB::PARAM_DOMAIN_CONTINUOUS){
      if (logscale[i]) d = 1;
      else d = (ub[i]-lb[i])/5;
      scaling.set(0,i,d);
      if (d<EPSILON) throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Bad scaling." );
    }
    else if (domain[i]==SGTELIB::PARAM_DOMAIN_CAT){
      scaling.set(0,i,ub[i]-lb[i]);
    }
    else{
      scaling.set(0,i,1);
    }
  }

  //-------------------------------------------------------
  // Display the information about optimized parameters
  //-------------------------------------------------------
  if (display){
    std::cout << "Model: " << get_short_string() << "\n";
    std::cout << "lb: [ ";
    for (i=0 ; i<N ; i++) std::cout << lb[i] << " ";
    std::cout << "]\n";
    std::cout << "ub: [ ";
    for (i=0 ; i<N ; i++) std::cout << ub[i] << " ";
    std::cout << "]\n";
    std::cout << "scaling: [ ";
    for (i=0 ; i<N ; i++){
      std::cout << scaling[i];
      if (logscale[i]) std::cout << "(log)";
      std::cout << " ";
    }
    std::cout << "]\n";
  }


  //----------------------------------------
  // Build set of starting points
  //----------------------------------------
  const int nx0 = budget/10;
  SGTELIB::Matrix X0 ("X0",nx0,N);

  const int use_lh = true;
  // RANDOM STARTING POINTS
  for (j=0 ; j<N ; j++){
    double lbj = lb[j];
    double ubj = ub[j];
    for (i=0 ; i<nx0 ; i++){ 
      if (use_lh) d = double(i-1)/double(nx0-2);
      else d = uniform_rand();

      if (logscale[j]) d = lb[j] * pow(ubj/lbj,d);
      else d = lbj + (ubj-lbj)*d;

      X0.set(i,j,d);
    } 
  }
  if (use_lh){
    // Shuffle columns (except first column)
    if (N>1){
      int i2;
      for (j=1; j<N; j++){
        for (i=0; i<nx0; i++){
          i2 = i + (int)std::floor(uniform_rand()*(nx0-i));
          if ( (i2<i) || (i2>=nx0) ){
            std::cout << "Error in permutation indexes\n";
            exit(0);
          }
          X0.swap(i,j,i2,j);
        }
      }
    }
  }

  // Add the default values.
  X0.add_rows(_param.get_x());





  //---------------------------------------------
  // Budget, poll size, success and objectives
  //---------------------------------------------
  SGTELIB::Matrix xtry ("xtry",1,N);
  // f contains the value of the error metric, returned by the model.
  // The smallest f, the better the model.
  // fmin is the smallest value of f found so 
  double fmin = +INF;
  // p is a penalty that allows to chose between
  // two sets of parameters that have the same f value.
  // pmin is the value of the best set of parameters so far;
  // For a given set of parameters, the value of p is returned by the class 
  // Surrogate_Parameters.
  // The penalty is particularly necessary for certain classes of error metrics
  // that are piece-wise constant (for example, all the order error metrics).
  double pmin = +INF;
  // ftry and ptry are the values of f and p for the current candidate.
  double ftry, ptry;
  // the MADS iteration is a success if a better set of parameters has been found.
  bool success;
  // Initial poll size value.
  double psize = 0.5;
  // Matrix containing the poll directions.
  SGTELIB::Matrix POLL;
  // xmin: Best set of parameters so far
  SGTELIB::Matrix xmin = X0.get_row(0);

  // Init cache of evaluated points
  SGTELIB::Matrix CACHE ("CACHE",0,N);
  bool cache_hit;

  //------------------------
  // LOOP
  //------------------------
  int iter=0;
  while (budget>0){

    success = false;

    if (display){
      std::cout << "=================================================\n";
      std::cout << "Budget: " << budget  << "\n";
      // Display best solution
      std::cout << "\nCurrent xmin:\n";
      std::cout << "X=[ " ;
      for (j=0 ; j<N ; j++) std::cout << xmin[j] << " ";
      std::cout << "] => " << fmin << " / " << pmin <<  "\n\n";
    }


    if (iter){
      // Create POLL candidates
      POLL = SGTELIB::Matrix::get_poll_directions(scaling,domain,psize);
      //POLL.display(std::cout);
      for (i=0 ; i<POLL.get_nb_rows() ; i++){
        for (j=0 ; j<N ; j++){
          // Add poll directions to poll center
          d = xmin[j];
          if (logscale[j]) d *= pow(4.0,POLL.get(i,j));  //exp(POLL.get(i,j));
          else             d += POLL.get(i,j);
          xtry.set(0,j,d);    
        }// End build candidate
        POLL.set_row(xtry,i);
      } // End Create POLL
      POLL.set_name("POLL-CANDIDATES");
      //POLL.display(std::cout);
    }
    else{
      // If iter==0, then evaluate starting points
      POLL = X0;
    }



    // Evaluate POLL
    for (i=0 ; i<POLL.get_nb_rows() ; i++){

      // Candidate
      xtry = POLL.get_row(i);
      xtry.set_name("xtry");

      // Display candidate
      if (display){
        if (iter) std::cout << "X = [ " ;
        else std::cout << "X0= [ " ;
        for (j=0 ; j<N ; j++) std::cout << xtry[j] << " ";
        std::cout << "] => ";
      }

      // Snap to bounds
      for (j=0 ; j<N ; j++){
        d = xtry[j];
        // Snap to bounds
        double lbj = lb[j];
        double ubj = ub[j];
        switch (domain[j]){
          case SGTELIB::PARAM_DOMAIN_CONTINUOUS:
            if (d<lbj) d = lbj;
            if (d>ubj) d = ubj;
            break;
          case SGTELIB::PARAM_DOMAIN_INTEGER:
            d = double(round(d));
            if (d<lbj) d=lbj;
            if (d>ubj) d=ubj;
            break;
          case SGTELIB::PARAM_DOMAIN_CAT:
            k = round(d);
            while (k>ubj) k-=int(ubj-lbj);
            while (k<lbj) k+=int(ubj-lbj);
            d = double(k);
            break;
          case SGTELIB::PARAM_DOMAIN_BOOL:
            d = (d>1/2)?1.0:0.0;
            break;
          case SGTELIB::PARAM_DOMAIN_MISC:
            throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Invalid variable domain!" );
            break;
        }
        xtry.set(0,j,d);    
      }

      // Check Cache
      cache_hit = (CACHE.find_row(xtry)!=-1);
      if (cache_hit){
        if (display) std::cout << "Cache hit\n";
      }
      else{
        // --------------------------------------
        // EVALUATION of metric and penalty
        // --------------------------------------
        // Register the xtry values in the parameter of the model
        _param.set_x(xtry);
        // Check that the parameters are consistent.
        _param.check();
        // Eval the objective (metric of the model)
        ftry = eval_objective();
        // Call the parameter class to get the penalty value.
        ptry = _param.get_x_penalty();
        // Reduce evaluation budget
        budget--;
        // Add the current point to the CACHE.
        CACHE.add_rows(xtry);

        // Display f and p
        if (display){
          if (ftry>=+INF) std::cout << "+inf" ;
          else std::cout << ftry;
          std::cout << " / " ;
          if (ptry>=+INF) std::cout << "+inf" ;
          else std::cout << ptry;
        }

        // Check for success 
        // The point xtry is a success if there is an improvement in the metric,
        // or, for an equal metric, if there is an improvement in the penalty.
        if ( (ftry<fmin) || ((ftry==fmin) && (ptry<pmin)) ){
          if (display) std::cout << "(!)";
          xmin = xtry;
          fmin = ftry;
          pmin = ptry;
          success = true;
        }
        if (display) std::cout << "\n";
      } // End Evaluation (i.e. No Cache Hit)

      // For iter==0, then we evaluate all the starting points.
      // For iter>0, if xtry is a success, we do not evaluate the other points of the POLL
      // (opportunistic evaluation of the poll)
      if ( (iter) && (success) ) break;

    }// END LOOP ON POLL (for i...)

    if (iter){
      // Update poll size
      if (success) psize*=2;
      else psize/=2;
    }
    iter++;

    // Check convergence
    if (psize<1e-6) break;
    if (budget<=0) break;

  }// End of optimization


  // Set param to optimal value
  _param.set_x(xmin);
  _param.check();

  //fmin = eval_objective();    // fmin is not used after this
  eval_objective();

  if (display){
    _param.display(std::cout);
    std::cout << "End parameter optimization\n";
    std::cout << "=================================\n";
  }

  // Check for Nan
  if (xmin.has_nan() || xmin.has_inf()) return false;

  delete [] logscale;
  delete [] domain;

  // Return success
  return true;

}//


/*--------------------------------------*/
/*    Evaluation of the error metric    */
/*       for a set of parameters        */
/*--------------------------------------*/
double SGTELIB::Surrogate::eval_objective ( void ){

  reset_metrics();

  // Build model
  bool ok = build_private();
  if ( ! ok ) return +INF;

  // Get the metric type specified in the parameter.
  const SGTELIB::metric_t mt = _param.get_metric_type();

  double metric = 0;
  // metric_multiple_obj indicate if the given metric "mt"
  // is scalar (one metric for all the blackbox outputs, like AOECV)
  // or is an array (one metric for each blackbox outputs, like RMSE)
  if (SGTELIB::one_metric_value_per_bbo(mt)){
    for (int i=0 ; i<_m ; i++) metric += get_metric(mt,i);
  }
  else{
    metric = get_metric(mt,0);
  }

  if ( isnan(metric) ) return +INF;
  if ( isinf(metric) ) return +INF;
  return metric;

}//










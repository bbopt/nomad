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

#include "Metrics.hpp"


/*----------------------------------------------------------*/
std::string SGTELIB::metric_type_to_str ( const SGTELIB::metric_t mt ) {
/*----------------------------------------------------------*/
  switch (mt){
    case SGTELIB::METRIC_EMAX   : return "EMAX"   ;
    case SGTELIB::METRIC_EMAXCV : return "EMAXCV" ;
    case SGTELIB::METRIC_RMSE   : return "RMSE"   ;
    case SGTELIB::METRIC_RMSECV : return "RMSECV" ;
    case SGTELIB::METRIC_ARMSE  : return "ARMSE"  ;
    case SGTELIB::METRIC_ARMSECV: return "ARMSECV";
    case SGTELIB::METRIC_OE     : return "OE"     ; 
    case SGTELIB::METRIC_OECV   : return "OECV"   ; 
    case SGTELIB::METRIC_AOE    : return "AOE"    ;
    case SGTELIB::METRIC_AOECV  : return "AOECV"  ;
    case SGTELIB::METRIC_EFIOE  : return "EFIOE"  ;
    case SGTELIB::METRIC_EFIOECV: return "EFIOECV";
    case SGTELIB::METRIC_LINV   : return "LINV"   ;
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined metric" );
  }
}//

/*----------------------------------------------------------*/
SGTELIB::norm_t SGTELIB::metric_type_to_norm_type ( const SGTELIB::metric_t mt ){
/*----------------------------------------------------------*/
  switch (mt){
    case SGTELIB::METRIC_EMAX   : 
    case SGTELIB::METRIC_EMAXCV : 
      return SGTELIB::NORM_INF;
    case SGTELIB::METRIC_RMSE   : 
    case SGTELIB::METRIC_RMSECV :
    case SGTELIB::METRIC_ARMSE   : 
    case SGTELIB::METRIC_ARMSECV :  
      return SGTELIB::NORM_2;
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"This metric does not have an associated norm" );
  }
}//

/*----------------------------------------------------------*/
bool SGTELIB::one_metric_value_per_bbo ( const SGTELIB::metric_t mt ) {
/*----------------------------------------------------------*/
  switch (mt){
    case SGTELIB::METRIC_EMAX   : 
    case SGTELIB::METRIC_EMAXCV : 
    case SGTELIB::METRIC_RMSE   : 
    case SGTELIB::METRIC_RMSECV : 
    case SGTELIB::METRIC_OE     : 
    case SGTELIB::METRIC_OECV   : 
    case SGTELIB::METRIC_LINV   : 
      return true;
    case SGTELIB::METRIC_ARMSE  : 
    case SGTELIB::METRIC_ARMSECV: 
    case SGTELIB::METRIC_AOE    : 
    case SGTELIB::METRIC_AOECV  : 
    case SGTELIB::METRIC_EFIOE    : 
    case SGTELIB::METRIC_EFIOECV  : 
      return false;
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined metric" );
  }
}//

/*----------------------------------------------------------*/
bool SGTELIB::metric_uses_cv ( const SGTELIB::metric_t mt ) {
/*----------------------------------------------------------*/
  switch (mt){
    case SGTELIB::METRIC_EMAXCV : 
    case SGTELIB::METRIC_RMSECV : 
    case SGTELIB::METRIC_OECV   : 
    case SGTELIB::METRIC_ARMSECV: 
    case SGTELIB::METRIC_AOECV  : 
    case SGTELIB::METRIC_EFIOECV  : 
      return true;
    case SGTELIB::METRIC_EMAX   : 
    case SGTELIB::METRIC_RMSE   : 
    case SGTELIB::METRIC_OE     : 
    case SGTELIB::METRIC_LINV   : 
    case SGTELIB::METRIC_ARMSE  : 
    case SGTELIB::METRIC_AOE    : 
    case SGTELIB::METRIC_EFIOE    : 
      return false;
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined metric" );
  }
}//


/*----------------------------------------------------------*/
SGTELIB::metric_t SGTELIB::metric_convert_single_obj ( const SGTELIB::metric_t mt ) {
/*----------------------------------------------------------*/
  switch (mt){
    // Metric that do not have a "Single obj" equivalent
    case SGTELIB::METRIC_EMAX   : 
    case SGTELIB::METRIC_EMAXCV : 
    case SGTELIB::METRIC_LINV   : 
      return SGTELIB::METRIC_AOECV;
    // Metric that have a "single obj" equivalent
    case SGTELIB::METRIC_RMSE   : 
      return SGTELIB::METRIC_ARMSE;
    case SGTELIB::METRIC_RMSECV : 
      return SGTELIB::METRIC_ARMSECV;
    case SGTELIB::METRIC_OE     : 
      return SGTELIB::METRIC_AOE;
    case SGTELIB::METRIC_OECV   : 
      return SGTELIB::METRIC_AOECV;
    // Metric that are "single obj"
    case SGTELIB::METRIC_ARMSE  : 
    case SGTELIB::METRIC_ARMSECV: 
    case SGTELIB::METRIC_AOE    : 
    case SGTELIB::METRIC_AOECV  : 
    case SGTELIB::METRIC_EFIOE  : 
    case SGTELIB::METRIC_EFIOECV: 
      return mt;
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined metric" );
  }
}//

/*----------------------------------------------------------*/
SGTELIB::metric_t SGTELIB::str_to_metric_type ( const std::string & s ) {
/*----------------------------------------------------------*/
  std::string ss = SGTELIB::toupper(s);
  if ( ss=="EMAX"   ){ return SGTELIB::METRIC_EMAX    ;}
  if ( ss=="EMAXCV" ){ return SGTELIB::METRIC_EMAXCV  ;}
  if ( ss=="RMSE"   ){ return SGTELIB::METRIC_RMSE    ;}
  if ( ss=="RMSECV" ){ return SGTELIB::METRIC_RMSECV  ;}
  if ( ss=="PRESS"  ){ return SGTELIB::METRIC_RMSECV  ;}
  if ( ss=="ARMSE"  ){ return SGTELIB::METRIC_ARMSE   ;}
  if ( ss=="ARMSECV"){ return SGTELIB::METRIC_ARMSECV ;}
  if ( ss=="OE"     ){ return SGTELIB::METRIC_OE      ;}
  if ( ss=="OECV"   ){ return SGTELIB::METRIC_OECV    ;}
  if ( ss=="AOE"    ){ return SGTELIB::METRIC_AOE     ;}
  if ( ss=="AOECV"  ){ return SGTELIB::METRIC_AOECV   ;}
  if ( ss=="EFIOE"  ){ return SGTELIB::METRIC_EFIOE   ;}
  if ( ss=="EFIOECV"){ return SGTELIB::METRIC_EFIOECV ;}
  if ( ss=="LINV"   ){ return SGTELIB::METRIC_LINV    ;}
  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unrecognised string \""+s+"\" ( "+ss+" )" );
}//


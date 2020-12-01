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

#include "Surrogate_Parameters.hpp"
using namespace SGTELIB;
/*---------------------------------------------------------------*/
/*         constructor                                           */
/* Use a model_type to define the model, then use default values */
/* associated to this specific model_type.                       */
/*---------------------------------------------------------------*/
SGTELIB::Surrogate_Parameters::Surrogate_Parameters ( const model_t mt ):
  _type  ( mt ){
  set_defaults();
  check();
}//


/*-----------------------------------------------------------------------*/
/*          constructor                                                  */
/* Use a string provided by the user. Read the model_type in the string, */
/* then parse the string to read all the other parameters.               */
/*-----------------------------------------------------------------------*/
SGTELIB::Surrogate_Parameters::Surrogate_Parameters ( const std::string & s ):
  _type  ( read_model_type(s) ){
  set_defaults();
  read_string(s);
  check();
}//

/*----------------------------*/
/*          Destructor        */
/*----------------------------*/
SGTELIB::Surrogate_Parameters::~Surrogate_Parameters ( void ) {

}//




/*----------------------------------------*/
/*   Extract the model type from a string */            
/*----------------------------------------*/
SGTELIB::model_t SGTELIB::Surrogate_Parameters::read_model_type ( const std::string & model_description){
  std::string s;
  std::istringstream in_line (model_description);	
  while ( in_line >> s ){
    s = SGTELIB::toupper(s);
    if (streqi(s,"TYPE")){
      in_line >> s;
      s = SGTELIB::toupper(s);
      return SGTELIB::str_to_model_type(s);
    }
  }
  // If nothing was found
  std::cout << "model_description: " << model_description << "\n";
  throw SGTELIB::Exception ( __FILE__ , __LINE__ , "No field \"TYPE\" found.");
}//

/*--------------------------------------------------------------------*/
/* Convert the param name provided by user into a standard param name */
/*--------------------------------------------------------------------*/
std::string SGTELIB::Surrogate_Parameters::to_standard_field_name (const std::string& field){
  
  if ( streqi(field,"TYPE") )           return "TYPE";
  if ( streqi(field,"TYPE_MODEL") )     return "TYPE";
  if ( streqi(field,"MODEL_TYPE") )     return "TYPE";

  if ( streqi(field,"DEGREE") )         return "DEGREE";

  if ( streqi(field,"WEIGHT") )         return "WEIGHT_TYPE";
  if ( streqi(field,"WEIGHT_TYPE") )    return "WEIGHT_TYPE";
  if ( streqi(field,"TYPE_WEIGHT") )    return "WEIGHT_TYPE";
  if ( streqi(field,"WEIGHT_CHOICE") )  return "WEIGHT_TYPE";
  if ( streqi(field,"CHOICE_WEIGHT") )  return "WEIGHT_TYPE";

  if ( streqi(field,"RIDGE") )                return "RIDGE";
  if ( streqi(field,"RIDGE_COEF") )           return "RIDGE";
  if ( streqi(field,"RIDGE_PARAM") )          return "RIDGE";
  if ( streqi(field,"RIDGE_COEFFICIENT") )    return "RIDGE";
  if ( streqi(field,"RIDGE_PARAMETER") )      return "RIDGE";
  if ( streqi(field,"REGULARIZATION_PARAM") ) return "RIDGE";
  if ( streqi(field,"REGULARIZATION_COEF") )  return "RIDGE";

  if ( streqi(field,"KERNEL_TYPE") )    return "KERNEL_TYPE";
  if ( streqi(field,"TYPE_KERNEL") )    return "KERNEL_TYPE";
  if ( streqi(field,"KERNEL") )         return "KERNEL_TYPE";

  if ( streqi(field,"KERNEL_COEF") )         return "KERNEL_COEF";
  if ( streqi(field,"KERNEL_COEFFICIENT") )  return "KERNEL_COEF";
  if ( streqi(field,"COEF_KERNEL") )         return "KERNEL_COEF";
  if ( streqi(field,"KERNEL_SHAPE") )        return "KERNEL_COEF";
  if ( streqi(field,"COEF_SHAPE") )          return "KERNEL_COEF";
  if ( streqi(field,"COEFFICIENT_SHAPE") )   return "KERNEL_COEF";
  if ( streqi(field,"SHAPE") )               return "KERNEL_COEF";
  if ( streqi(field,"SHAPE_COEF") )          return "KERNEL_COEF";
  if ( streqi(field,"SHAPE_COEFFICIENT") )   return "KERNEL_COEF";


  if ( streqi(field,"METRIC") )         return "METRIC_TYPE";
  if ( streqi(field,"METRIC_TYPE") )    return "METRIC_TYPE";
  if ( streqi(field,"TYPE_METRIC") )    return "METRIC_TYPE";

  if ( streqi(field,"BUDGET") )              return "BUDGET";
  if ( streqi(field,"OPTIM_BUDGET") )        return "BUDGET";
  if ( streqi(field,"BUDGET_OPTIM") )        return "BUDGET";
  if ( streqi(field,"OPTIMIZATION_BUDGET") ) return "BUDGET";
  if ( streqi(field,"BUDGET_OPTIMIZATION") ) return "BUDGET";

  if ( streqi(field,"PRESET") )         return "PRESET";

  if ( streqi(field,"OUTPUT") )         return "OUTPUT";
  if ( streqi(field,"OUTPUT_FILE") )    return "OUTPUT";

  if ( streqi(field,"DISTANCE") )       return "DISTANCE_TYPE";
  if ( streqi(field,"DISTANCE_TYPE") )  return "DISTANCE_TYPE";
  if ( streqi(field,"TYPE_DISTANCE") )  return "DISTANCE_TYPE";


  std::cout << "Field: " << field << "\n";
  throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Field not recognized: \""+field+"\"" );
  return "ERROR";
}

/*----------------------------*/
/*   Read model definition    */            
/*----------------------------*/
void SGTELIB::Surrogate_Parameters::read_string (const std::string & model_description){

  std::string field;
  std::string content;
  bool content_is_optim;
  std::istringstream in_line (model_description);	
  const bool display = false;
  if (display){
    std::cout << "Model description: " << model_description << "\n";
  }
  while ( in_line >> field ){

    if (display){
      std::cout << "FIELD: " << field ;
    }
    // Convert the field name into a std field name
    field = to_standard_field_name(field);
    if (display){
      std::cout << " (" << field << ")\n";
    }   
    // Check if this field is authorized for this type of model.
    if ( ! authorized_field(field) ){
      std::cout << "model_description: " << model_description << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unauthorized field \""+field+"\" in a model of type "+model_type_to_str(_type) );
    }
    if (display){
      std::cout << "CONTENT: " << content << "\n";
    }
    // Read the content 
    if ( !(in_line >> content) )
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Missing content for field \""+field+"\"" );

    // Detect if the content is "OPTIM".
    content_is_optim = ( streqi(content,"OPTIM") || streqi(content,"OPTIMIZATION") || streqi(content,"OPTIMIZE") );
    if (display){
      std::cout << "CONTENT IS OPTIM: " << content_is_optim << "\n";
    }

    // Check if optimization is allowed for this field.
    if ((content_is_optim) && (!authorized_optim(field))) {
      std::cout << "model_description: " << model_description << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Field \""+field+"\" cannot be optimized." );
    }

    if (streqi(field,"TYPE")){
      // If the field is "TYPE", then check that the model_type is consistent with the
      // one already detected.
      if ( str_to_model_type(content) != _type ){
        std::cout << "model_description: " << model_description << "\n";
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unconsistent model type!" );
      }
    }
    else if (streqi(field,"DEGREE")){
      if (content_is_optim)
        // _degree = default value.
        _degree_status = SGTELIB::STATUS_OPTIM;
      else{
        _degree = SGTELIB::stoi(content);
        _degree_status = SGTELIB::STATUS_FIXED;
      }
    }
    else if (streqi(field,"KERNEL_TYPE")){
      if (content_is_optim)
        // _kernel_type = default value        
        _kernel_type_status = SGTELIB::STATUS_OPTIM;
      else{
        _kernel_type = SGTELIB::str_to_kernel_type(content);
        _kernel_type_status = SGTELIB::STATUS_FIXED;
      }
    }
    else if ( streqi(field,"KERNEL_COEF") ){
      if (content_is_optim)
        // _kernel_coef = default value
        _kernel_coef_status = SGTELIB::STATUS_OPTIM;
      else{
        _kernel_coef = SGTELIB::stod(content);
        _kernel_coef_status = SGTELIB::STATUS_FIXED;
      }
    }
    else if (streqi(field,"RIDGE")){
      if (content_is_optim)
        // _ridge = default value
        _ridge_status = SGTELIB::STATUS_OPTIM;
      else{
        _ridge = SGTELIB::stod(content);
        _ridge_status = SGTELIB::STATUS_FIXED;
      }
    }
    else if ( streqi(field,"DISTANCE_TYPE") ){
      if (content_is_optim)
        // _distance_type = default value
        _distance_type_status = SGTELIB::STATUS_OPTIM;
      else{
        _distance_type = str_to_distance_type(content);
        _distance_type_status = SGTELIB::STATUS_FIXED;
      }
    }
    else if ( streqi(field,"WEIGHT_TYPE") ){
      _weight_type = str_to_weight_type(content);
      if (content_is_optim) 
        _weight_status = SGTELIB::STATUS_OPTIM;
      // Note, the weight_type can be one of the following:
      // - WTA1
      // - WTA3
      // - SELECT
      // - OPTIM
      // If the weight_type is OPTIM, then weight_status is OPTIM.
      // The variable weight_type is not optimized. It is the weights that are
      // optimized.
      // This is different from, for example, the kernel_type, which can be optimized.
    }
    else if ( streqi(field,"METRIC_TYPE") ){
      _metric_type = str_to_metric_type(content);
      // Cannot be optimized
    }
    else if ( streqi(field,"BUDGET") ){
      _budget = SGTELIB::stoi(content);
      // Cannot be optimized
    }
    else if ( streqi(field,"PRESET") ){
      _preset = content;
      // Cannot be optimized
    }
    else if ( streqi(field,"OUTPUT") ){
      _output = content;
      // Cannot be optimized
    }
    else{
      std::cout << "model_description: " << model_description << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Field not recognized: \""+field+"\"" );
    }
  }// end while read
}//


/*----------------------------------------------------------*/
/* Indicate if the field name given in input is authorized  */
/* for the current type of model                            */
/*----------------------------------------------------------*/
bool SGTELIB::Surrogate_Parameters::authorized_field ( const std::string & field ) const{

  if (streqi(field,"TYPE")) return true;
  if (streqi(field,"OUTPUT")) return true;
  if (streqi(field,"METRIC_TYPE")) return true;
  if (streqi(field,"BUDGET")) return true;

  switch (_type) {
    case SGTELIB::LINEAR:
    case SGTELIB::TGP: 
    case SGTELIB::SVN: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Not implemented yet! " );

    case SGTELIB::CN:
      if (streqi(field,"DISTANCE_TYPE")) return true;
      break;

    case SGTELIB::KRIGING: 
      if (streqi(field,"RIDGE"))         return true;
      if (streqi(field,"DISTANCE_TYPE")) return true;
      break;

    case SGTELIB::PRS: 
    case SGTELIB::PRS_EDGE: 
    case SGTELIB::PRS_CAT: 
      if (streqi(field,"DEGREE")) return true;
      if (streqi(field,"RIDGE")) return true;
      break;

    case SGTELIB::KS: 
      if (streqi(field,"KERNEL_TYPE")) return true;
      if (streqi(field,"KERNEL_COEF")) return true;
      if (streqi(field,"DISTANCE_TYPE")) return true;
      break;

    case SGTELIB::RBF: 
      if (streqi(field,"KERNEL_TYPE")) return true;
      if (streqi(field,"KERNEL_COEF")) return true;
      if (streqi(field,"DISTANCE_TYPE")) return true;
      if (streqi(field,"RIDGE")) return true;
      if (streqi(field,"PRESET")) return true;

      break;

    case SGTELIB::LOWESS:
      if (streqi(field,"DEGREE")) return true;
      if (streqi(field,"RIDGE")) return true;
      if (streqi(field,"KERNEL_TYPE")) return true;
      if (streqi(field,"KERNEL_COEF")) return true;
      if (streqi(field,"PRESET")) return true;
      if (streqi(field,"DISTANCE_TYPE")) return true;
      break; 

    case SGTELIB::ENSEMBLE: 
      if (streqi(field,"WEIGHT_TYPE")) return true;
      if (streqi(field,"PRESET")) return true;
      if (streqi(field,"DISTANCE_TYPE")) return true;
      break;

    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined model type" );
      break;
  }

  return false;
}//


/*----------------------------------------------------------*/
/* Indicate if the field given in input can be optimized    */
/*----------------------------------------------------------*/
bool SGTELIB::Surrogate_Parameters::authorized_optim ( const std::string & field ) {

  if (streqi(field,"DEGREE"))        return true;
  if (streqi(field,"RIDGE"))         return true;
  if (streqi(field,"KERNEL_TYPE"))   return true;
  if (streqi(field,"KERNEL_COEF"))   return true;
  if (streqi(field,"DISTANCE_TYPE")) return true;
  if (streqi(field,"WEIGHT_TYPE"))   return true;

  if (streqi(field,"TYPE"))          return false;
  if (streqi(field,"OUTPUT"))        return false;
  if (streqi(field,"METRIC_TYPE"))   return false;
  if (streqi(field,"PRESET"))        return false;
  if (streqi(field,"BUDGET"))        return false;

  std::cout << "Field : " << field << "\n";
  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined field" );

  return false;
}//



/*-------------------------------------------------*/
/*          Check                                  */
/* Verify the consistency of the set of parameters */
/*-------------------------------------------------*/
void SGTELIB::Surrogate_Parameters::check ( void ) {
  switch (_type) {
    case SGTELIB::LINEAR:
    case SGTELIB::TGP: 
    case SGTELIB::SVN: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Not implemented yet!" );


    case SGTELIB::PRS: 
    case SGTELIB::PRS_EDGE: 
    case SGTELIB::PRS_CAT: 
      if (_degree < 0 )
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"degree must be >= 0" );
      if (_ridge<0)
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"ridge must be >= 0" );
      break;

    case SGTELIB::KRIGING: 
      if (!SGTELIB::kernel_is_decreasing(_kernel_type))
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"kernel_type must be decreasing" );
      break;


    case SGTELIB::KS: 
      if (_kernel_coef <= 0)
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"kernel_coef must be > 0" );
      if (!SGTELIB::kernel_is_decreasing(_kernel_type))
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"kernel_type must be decreasing" );
      break;

    case SGTELIB::RBF: 
      if (_kernel_coef <= 0) 
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"kernel_coef must be > 0" );
      if (_ridge<0)
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"ridge must be >= 0" );
      if ( ( ! kernel_has_parameter(_kernel_type)) && (_kernel_type_status==SGTELIB::STATUS_FIXED) ){
        // If the kernel_type is fixed by the user, and this type of kernel does not
        // require a kernel_coefficient, then the kernel_coefficient is arbitrarily set to 1 
        // and its status is set to "fixed". See the class Kernel to see which kernel type requires
        // a coefficient or not.
        _kernel_coef = 1;
        _kernel_coef_status = SGTELIB::STATUS_FIXED;
      }
      break;

    case SGTELIB::LOWESS: 
      if ( (_degree < 0) || (_degree > 2) )
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"degree for LOWESS model must be 0, 1 or 2" );
      if (_ridge<0)
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"ridge must be >= 0" );
      // The default preset for LOWESS models is DGN. 
      // The preset defines how the weight of each data point is computed.
      // D : w_i = phi(distance_i), where distance_i is the distance between the prediction point 
      //           and the data point x_i.
      // DEN : w_i = phi(distance_i/dq_i), where dq_i is the distance between the prediction point
      //             and the q^th closest data point, and dq_i is computed with empirical method.
      // DGN : w_i = phi(distance_i/dq_i), where dq_i is computed with gamma method.
      // RE : w_i = phi(rank_i), where rank_i is the rank of x_i in terms of distance 
      //            to the prediction point, and the rank_i is computed with empirical method.
      // RG : w_i = phi(rank_i), where the rank is computed with gamma method.
      // REN : w_i = same as RE but the ranks are normalized in [0,1]
      // RGN : w_i = same as RG but the ranks are normalized in [0,1]
      if ( (_preset!="D"  ) &&
           (_preset!="DEN") &&
           (_preset!="DGN") &&
           (_preset!="RE" ) &&
           (_preset!="RG" ) &&
           (_preset!="REN") &&
           (_preset!="RGN") ){
        std::cout << "LOWESS preset : " << _preset << "\n";
        std::cout << "Possible values: D, DEN, DGN, RE, RG, REN, RGN.\n";
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"preset not recognized" );
      }
      if ( ! SGTELIB::kernel_is_decreasing( _kernel_type ) )
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"kernel_type must be decreasing" );
      break;

    case SGTELIB::ENSEMBLE: 
      break;

    case SGTELIB::CN: 
      break;

    default: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }


  // Count the number of parameters to optimize
  _nb_parameter_optimization = 0;
  if (_degree_status==SGTELIB::STATUS_OPTIM)           _nb_parameter_optimization++;
  if (_kernel_type_status==SGTELIB::STATUS_OPTIM)      _nb_parameter_optimization++;
  if (_kernel_coef_status==SGTELIB::STATUS_OPTIM)      _nb_parameter_optimization++;
  if (_ridge_status==SGTELIB::STATUS_OPTIM)            _nb_parameter_optimization++;
  if (_distance_type_status==SGTELIB::STATUS_OPTIM)    _nb_parameter_optimization++;
  if (_covariance_coef_status==SGTELIB::STATUS_OPTIM)  _nb_parameter_optimization += _covariance_coef.get_numel();
  if (_weight_status==SGTELIB::STATUS_OPTIM)           _nb_parameter_optimization += _weight.get_numel();

}//

/*----------------------------*/
/*          Display           */
/*----------------------------*/
void SGTELIB::Surrogate_Parameters::display ( std::ostream & out ) const {

  out << "Type: " << SGTELIB::model_type_to_str(_type) << std::endl;

  switch (_type) {

    case SGTELIB::SVN: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
        "Not implemented yet!" );

    case SGTELIB::CN: 
      out << "Distance_type: " << distance_type_to_str(_distance_type) << std::endl;
      break;

    case SGTELIB::KRIGING: 
      out << "Covariance coefs: " << "\n";
      out << "Ridge: " << _ridge << std::endl;
      _covariance_coef.display(out);
      break;

    case SGTELIB::PRS: 
    case SGTELIB::PRS_EDGE: 
    case SGTELIB::PRS_CAT: 
      out << "Degree: " << _degree << std::endl;
      out << "Ridge: " << _ridge << std::endl;
      out << "Distance_type: " << distance_type_to_str(_distance_type) << std::endl;
      break;

    case SGTELIB::KS: 
      out << "Kernel type: " << _kernel_type << std::endl;
      out << "Kernel coef: " << _kernel_coef << std::endl;
      out << "Distance_type: " << distance_type_to_str(_distance_type) << std::endl;
      break;

    case SGTELIB::RBF: 
      out << "Kernel type: " << _kernel_type << std::endl;
      out << "Kernel coef: " << _kernel_coef << std::endl;
      out << "Distance_type: " << distance_type_to_str(_distance_type) << std::endl;
      out << "ridge: " << _ridge << std::endl;
      break;

    case SGTELIB::LOWESS: 
      out << "Kernel coef: " << _kernel_coef << std::endl;
      out << "Distance_type: " << distance_type_to_str(_distance_type) << std::endl;
      out << "degree: " << _degree << std::endl;
      out << "ridge: " << _ridge << std::endl;
      out << "Preset: " << _preset << std::endl;
      break;

    case SGTELIB::ENSEMBLE: 
      out << "Metric type: " << metric_type_to_str(_metric_type) << std::endl;
      out << "Weight type: " << weight_type_to_str(_weight_type) << std::endl;
      out << "Preset: " << _preset << std::endl;
      break;

    default: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
}//

/*----------------------------*/
/*          Set defaults           */
/*----------------------------*/
void SGTELIB::Surrogate_Parameters::set_defaults ( void ) {

  _budget = 100;
  _metric_type = SGTELIB::METRIC_AOECV;
  _distance_type = SGTELIB::DISTANCE_NORM2;
  _distance_type_status = SGTELIB::STATUS_FIXED;
  _covariance_coef_status = SGTELIB::STATUS_FIXED;
  _weight_status = SGTELIB::STATUS_MODEL_DEFINED;
  _ridge = 0.001;
  _kernel_coef = 1;
  _kernel_type = SGTELIB::KERNEL_D1; 
  _covariance_coef = SGTELIB::Matrix("COVARIANCE_COEF",0,0);
  _weight = SGTELIB::Matrix("WEIGHT",0,0);


  switch (_type) {
    case SGTELIB::LINEAR:
    case SGTELIB::TGP: 
    case SGTELIB::SVN: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Not implemented yet!" );

    case SGTELIB::KRIGING: 
      _distance_type = SGTELIB::DISTANCE_NORM2;
      _distance_type_status = SGTELIB::STATUS_FIXED;
      _ridge = 1e-16;
      _ridge_status = SGTELIB::STATUS_OPTIM;
      _covariance_coef = SGTELIB::Matrix("COVARIANCE_COEF",1,2);
      {
        const double default_exponent = 2;
        _covariance_coef.set(0,0,default_exponent);

        const double default_factor   = 1;
        _covariance_coef.set(0,1,default_factor);
      }
      _covariance_coef_status = SGTELIB::STATUS_OPTIM;
      break;

    case SGTELIB::PRS: 
    case SGTELIB::PRS_EDGE: 
    case SGTELIB::PRS_CAT: 
      _degree = 2;
      _degree_status = SGTELIB::STATUS_FIXED;
      _ridge = 0.001;
      _ridge_status = SGTELIB::STATUS_FIXED;
      break;

    case SGTELIB::KS: 
      _kernel_type = SGTELIB::KERNEL_D1; 
      _kernel_type_status = SGTELIB::STATUS_FIXED;
      _kernel_coef = 5;
      _kernel_coef_status = SGTELIB::STATUS_OPTIM;
      _distance_type = SGTELIB::DISTANCE_NORM2;
      _distance_type_status = SGTELIB::STATUS_FIXED;
      break;

    case SGTELIB::RBF: 
      _kernel_type = SGTELIB::KERNEL_I2;
      _kernel_type_status = SGTELIB::STATUS_FIXED;
      _kernel_coef = 1;
      _kernel_coef_status = SGTELIB::STATUS_OPTIM;
      _distance_type = SGTELIB::DISTANCE_NORM2;
      _distance_type_status = SGTELIB::STATUS_FIXED;
      _ridge = 0.001;
      _ridge_status = SGTELIB::STATUS_FIXED;
      _preset = "I";
      break;

    case SGTELIB::LOWESS: 
      _kernel_coef = 1.0;
      _kernel_coef_status = SGTELIB::STATUS_OPTIM;
      _kernel_type = SGTELIB::KERNEL_D1;
      _kernel_type_status = SGTELIB::STATUS_FIXED;
      _distance_type = SGTELIB::DISTANCE_NORM2;
      _distance_type_status = SGTELIB::STATUS_FIXED;
      _degree = 2;
      _degree_status = SGTELIB::STATUS_FIXED;
      _ridge = 0.001;
      _ridge_status = SGTELIB::STATUS_FIXED;
      _preset = "DGN";
      break;

    case SGTELIB::ENSEMBLE: 
      _weight_type = WEIGHT_SELECT; 
      _weight_status = SGTELIB::STATUS_MODEL_DEFINED;
      _preset = "DEFAULT";
      break;

    case SGTELIB::CN:
      break;

    default: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }

  // Default output file
  _output = "NULL";

}//

/*----------------------------*/
/*          Display           */
/*----------------------------*/
std::string SGTELIB::Surrogate_Parameters::get_string ( void ) const {

  std::string s = "TYPE "+model_type_to_str(_type);

  switch (_type) {
    case SGTELIB::LINEAR:
    case SGTELIB::TGP: 

    case SGTELIB::SVN: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Not implemented yet!" );

    case SGTELIB::KRIGING: 
      s += " DISTANCE_TYPE " + distance_type_to_str(_distance_type);
      break;

    case SGTELIB::CN: 
      s += " DISTANCE_TYPE " + distance_type_to_str(_distance_type);
      break;

    case SGTELIB::PRS: 
    case SGTELIB::PRS_EDGE: 
    case SGTELIB::PRS_CAT: 
      s += " DEGREE " + itos(_degree);
      s += " RIDGE " + dtos(_ridge);
      break;

    case SGTELIB::KS: 
      s += " KERNEL_TYPE " + kernel_type_to_str(_kernel_type);
      s += " KERNEL_SHAPE " + dtos(_kernel_coef);
      s += " DISTANCE_TYPE " + distance_type_to_str(_distance_type);
      break;

    case SGTELIB::RBF: 
      s += " PRESET " + _preset;
      s += " KERNEL_TYPE " + kernel_type_to_str(_kernel_type);
      s += " KERNEL_SHAPE " + dtos(_kernel_coef);
      s += " DISTANCE_TYPE " + distance_type_to_str(_distance_type);
      s += " RIDGE " + dtos(_ridge);
      break;

    case SGTELIB::LOWESS: 
      s += " DEGREE " + itos(_degree);
      s += " KERNEL_TYPE " + kernel_type_to_str(_kernel_type);
      s += " KERNEL_SHAPE " + dtos(_kernel_coef);
      s += " DISTANCE_TYPE " + distance_type_to_str(_distance_type);
      s += " RIDGE " + dtos(_ridge);
      break;

    case SGTELIB::ENSEMBLE: 
      s += " METRIC_TYPE " + metric_type_to_str(_metric_type);
      s += " WEIGHT_TYPE " + weight_type_to_str(_weight_type);
      s += " PRESET " + _preset;
      s += " DISTANCE_TYPE " + distance_type_to_str(_distance_type);
      break;

    default: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
  return s;

}//

/*----------------------------*/
/*          Display           */
/*----------------------------*/
std::string SGTELIB::Surrogate_Parameters::get_short_string ( void ) const {

  std::string s = model_type_to_str(_type);

  switch (_type) {
    case SGTELIB::LINEAR:
    case SGTELIB::TGP: 
    case SGTELIB::SVN: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Not implemented yet!" );

    case SGTELIB::KRIGING: 
      if (_kernel_type != SGTELIB::KERNEL_D1){
        s += " " + kernel_type_to_str(_kernel_type);
      }
      if (_distance_type != SGTELIB::DISTANCE_NORM2){
        s += " " + distance_type_to_str(_distance_type);
      }
      break;

    case SGTELIB::PRS: 
    case SGTELIB::PRS_EDGE: 
    case SGTELIB::PRS_CAT: 
      s += " " + itos(_degree);
      if (_ridge>0) s+= " R "+dtos(_ridge);
      break;

    case SGTELIB::KS: 
      s += " " + dtos(_kernel_coef);
      if (_kernel_type != SGTELIB::KERNEL_D1){
        s += " " + kernel_type_to_str(_kernel_type);
      }
      if (_distance_type != SGTELIB::DISTANCE_NORM2){
        s += " " + distance_type_to_str(_distance_type);
      }
      break;

    case SGTELIB::CN: 
      if (_distance_type != SGTELIB::DISTANCE_NORM2){
        s += " " + distance_type_to_str(_distance_type);
      }
      break;

    case SGTELIB::RBF: 
      s += " " + dtos(_kernel_coef);
      if (_kernel_type != SGTELIB::KERNEL_D1){
        s += " " + kernel_type_to_str(_kernel_type);
      }
      if (_distance_type != SGTELIB::DISTANCE_NORM2){
        s += " " + distance_type_to_str(_distance_type);
      }
      s += " " + dtos(_ridge);
      s += " "+_preset;
      break;


    case SGTELIB::LOWESS: 
      s += " " + dtos(_kernel_coef);
      if (_distance_type != SGTELIB::DISTANCE_NORM2){
        s += " " + distance_type_to_str(_distance_type);
      }
      if (_kernel_type != SGTELIB::KERNEL_D1){
        s += " " + kernel_type_to_str(_kernel_type);
      }
      if (_degree!=2) s += " Degree " + itos(_degree);
      if (_ridge>0) s+= " R "+dtos(_ridge);
      s += " "+_preset;
      break;

    case SGTELIB::ENSEMBLE: 
      s += " " + metric_type_to_str(_metric_type) + " " + weight_type_to_str(_weight_type);
      if (!streqi(_preset,"DEFAULT")) s += " "+_preset;
      break;

    default: 
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
  return s;

}//





/*------------------------------------------------------------*/
/*  get_x                                                     */
/* Creates a vector that contains the numerical values of all */
/* all the parameters that must be optmized.                  */
/*------------------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate_Parameters::get_x ( void ){

    SGTELIB::Matrix X ("X",1,_nb_parameter_optimization);
    int k=0;
    if (_degree_status        == SGTELIB::STATUS_OPTIM) X.set(0,k++,double(_degree));
    if (_ridge_status         == SGTELIB::STATUS_OPTIM) X.set(0,k++,double(_ridge));
    if (_kernel_coef_status   == SGTELIB::STATUS_OPTIM) X.set(0,k++,double(_kernel_coef));
    if (_kernel_type_status   == SGTELIB::STATUS_OPTIM) X.set(0,k++,double(_kernel_type));
    if (_distance_type_status == SGTELIB::STATUS_OPTIM) X.set(0,k++,double(_distance_type));

    if (_covariance_coef_status == SGTELIB::STATUS_OPTIM){
      for (int j=0 ; j<_covariance_coef.get_nb_cols() ; j++){
        X.set(0,k++,_covariance_coef[j]);
      }
    }

    if (_weight_status == SGTELIB::STATUS_OPTIM){
      for (int i=0 ; i<_weight.get_nb_rows() ; i++){
        for (int j=0 ; j<_weight.get_nb_cols() ; j++){
          X.set(0,k++,_weight.get(i,j));
        }
      }
    }

    if ( k != _nb_parameter_optimization){
      std::cout << "k=" << k << "\n";
      std::cout << "_nb_parameter_optimization=" << _nb_parameter_optimization << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unconcistency in the value of k." );
    }

    return X;

}//

/*-------------------------------------------------*/
/*  Set the parameters from an external value of x */
/*-------------------------------------------------*/
void SGTELIB::Surrogate_Parameters::set_x ( const SGTELIB::Matrix& X ){

    int k=0;
    if (_degree_status        == SGTELIB::STATUS_OPTIM) _degree      = int(X.get(k++));
    if (_ridge_status         == SGTELIB::STATUS_OPTIM) _ridge       = double(X.get(k++));
    if (_kernel_coef_status   == SGTELIB::STATUS_OPTIM) _kernel_coef = double(X.get(k++));
    if (_kernel_type_status   == SGTELIB::STATUS_OPTIM) _kernel_type = SGTELIB::int_to_kernel_type(int(X.get(k++)));
    if (_distance_type_status == SGTELIB::STATUS_OPTIM) _distance_type=SGTELIB::int_to_distance_type(int(X.get(k++)));

    if (_covariance_coef_status == SGTELIB::STATUS_OPTIM){
      for (int j=0 ; j<_covariance_coef.get_nb_cols() ; j++){
        _covariance_coef.set(0,j,X.get(k++));
      }
    } 

    if (_weight_status == SGTELIB::STATUS_OPTIM){
      // For each BBO
      for (int i=0 ; i<_weight.get_nb_rows() ; i++){
        for (int j=0 ; j<_weight.get_nb_cols() ; j++){
        // For each model
          _weight.set(i,j,X.get(k++));
        }
      }
      _weight.normalize_cols();
    } 

    if ( k != _nb_parameter_optimization ){
      std::cout << "k=" << k << "\n";
      std::cout << "_nb_parameter_optimization=" << _nb_parameter_optimization << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unconcistency in the value of k." );
    }

}//

/*------------------------------------------------------------------------------------*/
/*  Parameter domains and definitions                                                 */
/* Defines the bounds, domain and log-scale for each parameter that must be optimized */
/*------------------------------------------------------------------------------------*/
void SGTELIB::Surrogate_Parameters::get_x_bounds ( SGTELIB::Matrix * LB ,
                                                   SGTELIB::Matrix * UB ,
                                                   SGTELIB::param_domain_t * domain,
                                                   bool * logscale ){

    // Check that the arrays or matrices are initialized.
    if ( (!LB) || (!UB) || (!domain) || (!logscale) ){
      std::cout << LB << " " << UB << " " << domain << " " << logscale << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Pointers are NULL." );
    }


    int i,j;
    const int N = _nb_parameter_optimization;
    for (i=0 ; i<N ; i++) logscale[i] = false;

    int k=0;
    // --------- DEGREE --------------------
    if (_degree_status == SGTELIB::STATUS_OPTIM){
      LB->set(0,k,0);
      if (_type==SGTELIB::LOWESS) {
        UB->set(0,k,2);
        domain[k] = SGTELIB::PARAM_DOMAIN_INTEGER;
      }
      else {
        UB->set(0,k,6);
        domain[k] = SGTELIB::PARAM_DOMAIN_INTEGER;
      }
      k++;
    }
    // --------- RIDGE --------------------
    if (_ridge_status == SGTELIB::STATUS_OPTIM){
      LB->set(0,k,1e-16);
      UB->set(0,k,1e-1);
      domain[k] = SGTELIB::PARAM_DOMAIN_CONTINUOUS;
      logscale[k] = true;
      k++;
    }
    // --------- KERNEL_COEF --------------------
    if (_kernel_coef_status == SGTELIB::STATUS_OPTIM){
      LB->set(0,k,1e-2);
      UB->set(0,k,100);
      domain[k] = SGTELIB::PARAM_DOMAIN_CONTINUOUS;
      logscale[k] = true;
      k++;
    }
    // --------- KERNEL_TYPE --------------------
    if (_kernel_type_status == SGTELIB::STATUS_OPTIM){
      LB->set(0,k,0);
      if (_type==SGTELIB::RBF) {
        UB->set(0,k,double(NB_KERNEL_TYPES-1));
      }
      else{
        UB->set(0,k,double(NB_DECREASING_KERNEL_TYPES-1));
      }
      domain[k] = SGTELIB::PARAM_DOMAIN_CAT;
      k++;
    }
    // --------- DISTANCE_TYPE --------------------
    if (_distance_type_status == SGTELIB::STATUS_OPTIM){
      LB->set(0,k,0);
      UB->set(0,k,double(NB_DISTANCE_TYPES-1));
      domain[k] = SGTELIB::PARAM_DOMAIN_CAT;
      k++;
    }
    // --------- COVARIANCE COEF --------------------
    if (_covariance_coef_status == SGTELIB::STATUS_OPTIM){
      const int v = _covariance_coef.get_nb_cols()/2;
      for (j=0 ; j<v ; j++){
        // Exponent parameter
        LB->set(0,k,0.5);
        UB->set(0,k,3.0);
        domain[k] = SGTELIB::PARAM_DOMAIN_CONTINUOUS;
        logscale[k] = false;
        k++;
        // Factor parameter
        LB->set(0,k,1e-3);
        UB->set(0,k,1e+0);
        domain[k] = SGTELIB::PARAM_DOMAIN_CONTINUOUS;
        logscale[k] = true;
        k++;
      }
    }
    // --------- WEIGHT --------------------
    if (_weight_status == SGTELIB::STATUS_OPTIM){
      for (i=0 ; i<_weight.get_nb_rows() ; i++){
        for (j=0 ; j<_weight.get_nb_cols() ; j++){
          LB->set(0,k,0);
          UB->set(0,k,1);
          domain[k] = SGTELIB::PARAM_DOMAIN_CONTINUOUS;
          logscale[k] = false;
          k++;
        }
      }
    }

    if ( k != N){
      std::cout << "k=" << k << "\n";
      std::cout << "N=" << N << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unconcistency in the value of k." );
    }


    // ----- CHECK CONSISTENCY -------
    bool error = false;
    for (j=0 ; j<N ; j++){
      // Check bounds order
      if (LB->get(j)>=UB->get(j)){
        error=true;
        std::cout << "Variable " << j << "\n";
        std::cout << "LB (=" << LB->get(j) << ") >= UB (=" << UB->get(j) << ")\n";
      }
      // Check that only continuous variables are using a log scale
      if ( (logscale[j]) && (domain[j]!=SGTELIB::PARAM_DOMAIN_CONTINUOUS) ){
        error=true;
        std::cout << "Variable " << j << "\n";
        std::cout << "Uses logscale and is not continuous.\n";
      }
      // Check that variables with log scale have bounds of the same sign.
      if (logscale[j]){
        if (LB->get(j)*UB->get(j)<=0){
          //error=true;
          std::cout << "Variable " << j << "\n";
          std::cout << "LB =" << LB->get(j) << "\nUB =" << UB->get(j) << "\n";
          std::cout << "The bounds are not appropriate for logscale optimization.\n";
        }
      }

      // Check domain types
      switch (domain[j]){
        case SGTELIB::PARAM_DOMAIN_CONTINUOUS:
          break;
        case SGTELIB::PARAM_DOMAIN_INTEGER:
        case SGTELIB::PARAM_DOMAIN_CAT:
          if (double(round(LB->get(j)))!=LB->get(j)){
            error=true;
            std::cout << "Variable " << j << " (Integer or Categorical)\n";
            std::cout << "LB (=" << LB->get(j) << ") is not an integer\n";
          }
          if (double(round(UB->get(j)))!=UB->get(j)){
            error=true;
            std::cout << "Variable " << j << " (Integer or Categorical)\n";
            std::cout << "UB (=" << UB->get(j) << ") is not an integer\n";
          }
          break;
        case SGTELIB::PARAM_DOMAIN_BOOL:
          if (LB->get(j)==0){
            error=true;
            std::cout << "Variable " << j << " (Boolean)\n";
            std::cout << "LB (=" << LB->get(j) << ") is not 0\n";
          }
          if (UB->get(j)==1){
            error=true;
            std::cout << "Variable " << j << " (Boolean)\n";
            std::cout << "UB (=" << UB->get(j) << ") is not 1\n";
          }
          break;
        case SGTELIB::PARAM_DOMAIN_MISC:
          error=true;
          std::cout << "Variable " << j << " is MISC\n";
          break;
      }
    }
    
    if (error){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Error in definition of LB, UB or domain!" );
    }

}//

/*------------------------------------*/
/*  Check parameter X                 */
/*------------------------------------*/
bool SGTELIB::Surrogate_Parameters::check_x ( void ){

    SGTELIB::Matrix X = get_x();
    bool error = false;
    // Check dimension of X
    if (X.get_nb_rows()!=1){
      error = true;
      std::cout << "Number of rows is not 1\n";
    }
    // Check dimension of X
    const int N = _nb_parameter_optimization;
    if (X.get_nb_cols()!=N){
      error = true;
      std::cout << "Number of cols is not consistent with _nb_parameter_optimization\n";
    }

    // Get bound info.
    SGTELIB::Matrix * LB = new SGTELIB::Matrix("LB",1,N);
    SGTELIB::Matrix * UB = new SGTELIB::Matrix("UB",1,N);
    param_domain_t * domain = new param_domain_t[N];
    bool * logscale = new bool[N];
    get_x_bounds ( LB , UB , domain , logscale);

    for (int j=0 ; j<_nb_parameter_optimization ; j++){
      // Check bounds
      if (X[j]<LB->get(j)){
        error=true;
        std::cout << "X[" << j << "] < lower bound\n";
      }
      if (X[j]>UB->get(j)){
        error=true;
        std::cout << "X[" << j << "] > upper bound\n";
      }
      // Check types
      switch (domain[j]){
        case SGTELIB::PARAM_DOMAIN_CONTINUOUS:
          break;
        case SGTELIB::PARAM_DOMAIN_INTEGER:
        case SGTELIB::PARAM_DOMAIN_CAT:
          if (double(round(X[j]))!=X[j]){
            error=true;
            std::cout << "Variable " << j << " (Integer or Categorical)\n";
            std::cout << "X[" << j << "]=" << X[j] << " is not an integer\n";
          }
          break;
        case SGTELIB::PARAM_DOMAIN_BOOL:
          if ((X[j]!=0) && (X[j]!=1)){
            error=true;
            std::cout << "Variable " << j << " (Boolean)\n";
            std::cout << "X[" << j << "]=" << X[j] << " is not a boolean\n";
          }
          break;
        case SGTELIB::PARAM_DOMAIN_MISC:
          error=true;
          std::cout << "Variable " << j << " is MISC\n";
          break;
      }
    }// End loop on j

    // Check dimension of _covariance_coef
    if (_covariance_coef.get_nb_rows()>1){
      error = true;
      std::cout << "Covariance_coef should have only one row.\n";
    }
    
    if (error){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Invalid X!" );
    }

    delete LB;
    delete UB;
    delete [] domain; 

    return ( !  error);
}//



/*------------------------------------*/
/*  Display the parameters            */
/*------------------------------------*/
void SGTELIB::Surrogate_Parameters::display_x ( std::ostream & out ){
    out << "Parameter set {\n";
    out << "    Dimension(X) " << _nb_parameter_optimization << "\n";
    if (_degree_status        == SGTELIB::STATUS_OPTIM)
      out << "    Degree " << _degree << " ";
    if (_ridge_status         == SGTELIB::STATUS_OPTIM)
      out << "    Ridge " << _ridge << " ";
    if (_kernel_coef_status   == SGTELIB::STATUS_OPTIM)
      out << "    Kernel_coef " << _kernel_coef << " ";
    if (_kernel_type_status   == SGTELIB::STATUS_OPTIM)
      out << "    Kernel_type " << SGTELIB::kernel_type_to_str(_kernel_type) << " (" << _kernel_type << ") ";
    if (_distance_type_status == SGTELIB::STATUS_OPTIM)
      out << "    Distance_type " << SGTELIB::distance_type_to_str(_distance_type) << " (" << _distance_type << ") ";
    if (_covariance_coef_status == SGTELIB::STATUS_OPTIM){
      out << "    Covariance_coef [ ";
      for (int j=0 ; j<_covariance_coef.get_nb_cols() ; j++){
        out << _covariance_coef.get(0,j) << " ";
      }
      out << "]\n";
    }
    if (_weight_status == SGTELIB::STATUS_OPTIM){
      out << "    Weight [ ";
      for (int i=0 ; i<_weight.get_nb_rows() ; i++){
        if (i>0) out << "                 ";
        for (int j=0 ; j<_weight.get_nb_cols() ; j++){
          out << _weight.get(i,j) << " ";
        }
        if (i==_weight.get_nb_rows()) out << " ]";
        out << "\n";
      }
    }

    out << "}\n";
  
}//


/*----------------------------------------------*/
/*  Smoothness penalty for a set of parameters  */
/*----------------------------------------------*/
double SGTELIB::Surrogate_Parameters::get_x_penalty ( void ){
  double pen = 0;
  
  if (_degree_status        == SGTELIB::STATUS_OPTIM) pen += _degree;
  if (_ridge_status         == SGTELIB::STATUS_OPTIM) pen += log(_ridge); 
  if (_kernel_coef_status   == SGTELIB::STATUS_OPTIM) pen += log(_kernel_coef);
  if (_distance_type_status == SGTELIB::STATUS_OPTIM){
    switch (_distance_type){
      case DISTANCE_NORM2:
        pen += 0;
        break;
      case DISTANCE_NORM1:
      case DISTANCE_NORMINF:
        pen += 1;
        break;
      case DISTANCE_NORM2_IS0:
      case DISTANCE_NORM2_CAT:
        pen += 10;
        break;
    }
  }
  if (_covariance_coef_status == SGTELIB::STATUS_OPTIM){
    int v = _covariance_coef.get_nb_cols()/2;
    int ip=0;
    for (int i=0 ; i<v ; i++){ 
      // Exponent (the larger, the smoother) 
      pen -= _covariance_coef.get(ip++);
      // Factor (the smaller, the smoother)
      pen += log( _covariance_coef.get(ip++) );
    }
  }
  if (_weight_status == SGTELIB::STATUS_OPTIM){
    double wij;
    for (int i=0 ; i<_weight.get_nb_rows() ; i++){
      for (int j=0 ; j<_weight.get_nb_cols() ; j++){
        wij = _weight.get(i,j);
        pen += wij*wij;
      }
    }
  }

  if ( isinf(pen) ) pen=+INF;
  if ( isnan(pen) ) pen=+INF;
  return pen;
}






/*-------------------------------------------------*/
/*  update the dimension of _covariance_parameter  */
/*  for Kriging models                             */
/*-------------------------------------------------*/
// The matrix containing the covariance coefficients, for Kriging models, 
// is initialized with 2 components (factor and exponent). 
// In sgtelib 2.0.1, the factor and exponent are the same for all input variables.
// However, it would be possible to use a different value of factor and exponent
// for each variable (which is the case in most Kriging implementation).
// The problem is that the Surrogate_Parameters class is built before the model itself,
// and the instance Surrogate_Parameters is not able to retrieve the dimension of
// the input space from either the training set or the model. The following method allows
// the model to tell to the Surrogate_Parameters what must be the dimension of
// the set of covariance coefficients. 
// This function is not used for now, but might be used in future versions of sgtelib.
void SGTELIB::Surrogate_Parameters::update_covariance_coef ( const int v ){
  // The input value v must be equal to 1 or to the number of variables.
  // Check the old dimension
  const int v0 = _covariance_coef.get_nb_cols()/2;
  if (v<v0) throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"v < v0" );
  if (v0 == v) return;

  // Filling values
  double factor_mean = 0;
  double exponent_mean = 0;

  //  Compute the mean value for exponent and factor.
  int k = 0;
  for (int i=0 ; i<v0 ; i++){
    exponent_mean += _covariance_coef[k++];
    factor_mean   += _covariance_coef[k++];
  }
  exponent_mean /= v0;
  factor_mean /= v0;

  // Create additional columns
  SGTELIB::Matrix Add("Add",1,2);
  Add.set(0,0,exponent_mean);
  Add.set(0,1,factor_mean);

  // Add columns to _covariance_coef
  for (int i=0 ; i<v-v0 ; i++) _covariance_coef.add_cols(Add);
}//





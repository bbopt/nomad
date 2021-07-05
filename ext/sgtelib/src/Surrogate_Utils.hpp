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

#ifndef __SGTELIB_SURROGATE_UTILS__
#define __SGTELIB_SURROGATE_UTILS__

#include "Defines.hpp"
#include "Exception.hpp"
#include "Matrix.hpp"
#include <sys/stat.h>

// Helpful for compilaton on some platforms
#ifndef _WIN32
#include <sys/time.h>
#endif

// CASE Visual Studio C++ compiler
#ifdef _MSC_VER
#pragma warning(disable:4996)
#include <io.h>
#define isnan(x) _isnan(x)
#define isdigit(x) _isdigit(x)
#define isinf(x) (!_finite(x))

typedef struct timeval {
     long tv_sec;
     long tv_usec;
} timeval;

#else
#include <unistd.h>
#endif


#include <cstring>
#include <cctype>

namespace SGTELIB {

  enum distance_t {
    DISTANCE_NORM2 ,
    DISTANCE_NORM1 ,
    DISTANCE_NORMINF ,
    DISTANCE_NORM2_IS0,
    DISTANCE_NORM2_CAT
  };
  const int NB_DISTANCE_TYPES = 5;

  // model type:
  enum model_t {
    LINEAR   ,
    TGP      ,
    DYNATREE ,
    PRS      ,
    PRS_EDGE ,
    PRS_CAT  ,
    KS       ,
    CN       ,
    KRIGING  ,
    SVN      ,
    RBF      ,
    LOWESS      ,
    ENSEMBLE 
  };
  const int NB_MODEL_TYPES = 12;

  // Aggregation methods (for the Surrogate_Ensemble)
  enum weight_t {
    WEIGHT_SELECT,// Take the model with the best metrics.
    WEIGHT_OPTIM, // Optimize the metric
    WEIGHT_WTA1,  // Goel, Ensemble of surrogates 2007
    WEIGHT_WTA3,  // Goel, Ensemble of surrogates 2007
    WEIGHT_EXTERN // Belief vector is set externaly by the user.
  };
  const int NB_WEIGHT_TYPES = 5;
  
 
  // Diff in ms
  int diff_ms(timeval t1, timeval t2);

  // Compare strings
  bool streq       ( const std::string & s1 , const std::string & s2 );
  bool streqi      ( const std::string & s1 , const std::string & s2 );
  // Check if s is a substring of S
  bool string_find ( const std::string & S  , const std::string & s );
  //bool issubstring (const std::string S , const std::string s);


  // Remove useless spaces in string
  std::string deblank ( const std::string & s_input );

  // test if a file exists
  bool exists (const std::string & file);

  // Word count
  int count_words(const std::string & s );

  // add string on a new line of an existing files
  void append_file (const std::string & s , const std::string & file);

  // wait 
  void wait (double t);

  // isdef (not nan nor inf)
  bool isdef ( const double x );

  // rounding:
  int round ( double d );
  double rceil (double d);

  // relative error:
  double rel_err ( double x , double y );

  // distance between two points:
  double dist ( const double * x , const double * y , int n );

  // same sign
  bool same_sign (const double a, const double b);

  // conversion functions (to string) :
  std::string itos ( int    );
  std::string dtos ( double );
  std::string btos ( bool   );
  double stod ( const std::string & s );
  int    stoi ( const std::string & s );
  bool   stob ( const std::string & s );

  std::string toupper ( const std::string & s   );

  DLL_API std::string model_output_to_str       ( const SGTELIB::model_output_t );
  DLL_API std::string model_type_to_str         ( const SGTELIB::model_t        );
  DLL_API std::string bbo_type_to_str           ( const SGTELIB::bbo_t          );
  DLL_API std::string weight_type_to_str        ( const SGTELIB::weight_t       );
  DLL_API std::string distance_type_to_str      ( const SGTELIB::distance_t     );


  // conversion functions (from string) :
  bool isdigit                                       ( const std::string & s );
  DLL_API SGTELIB::model_t         str_to_model_type         ( const std::string & s );
  DLL_API SGTELIB::weight_t        str_to_weight_type        ( const std::string & s );
  DLL_API SGTELIB::distance_t      str_to_distance_type      ( const std::string & s );
  DLL_API SGTELIB::distance_t      int_to_distance_type      ( const int i );


  /*
  // Find the index of the smallest value in an array v of size vsize.
  int get_min_index ( const double * v , const int vsize );
  // (optional: exclude index "i_exclude" from the search)
  int get_min_index ( const double * v , const int vsize , const int i_exclude);
  */

  // Statistics
  double normcdf ( double x );
  double normcdf ( double x , double mu , double sigma );
  double normpdf ( double x );
  double normpdf ( double x , double mu , double sigma );
  double normei  ( double fh, double sh , double f_min  );
  double gammacdf   ( double x, double a, double b);
  double gammacdfinv( double f, double a, double b);
  double lower_incomplete_gamma ( const double x, const double p );

  double uniform_rand (void);
  double quick_norm_rand (void);
}

#endif

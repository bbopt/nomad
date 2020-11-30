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

#include "Surrogate_Utils.hpp"

#include <string>
#include <sstream>


/*-------------------------------*/
/*     string comparison         */
/*-------------------------------*/
bool SGTELIB::streq  ( const std::string & s1 , const std::string & s2 ){
  return !std::strcmp(s1.c_str(),s2.c_str()); 
}//

bool SGTELIB::streqi ( const std::string & s1 , const std::string & s2 ){
  const std::string s1u = SGTELIB::toupper(s1);
  const std::string s2u = SGTELIB::toupper(s2);
  return !std::strcmp(SGTELIB::toupper(s1).c_str(),s2u.c_str()); 
}//

// Check if s is in S.
bool SGTELIB::string_find ( const std::string & S , const std::string & s ){
  const std::string Su = SGTELIB::toupper(S);
  const std::string su = SGTELIB::toupper(s);
  return ( Su.find(su) < Su.size() );
}//
/*
bool SGTELIB::issubstring (const std::string S , const std::string s){
  return result = (S.find(s) != std::string::npos);
}
*/

/*-------------------------------*/
/*      deblank                  */
/*-------------------------------*/
std::string SGTELIB::deblank ( const std::string & s_input ){
  std::string s = s_input;
  // Remove leading spaces
  while ( (s.length()) && (s.at(0)==' ') ){
    s.erase(0,1);
  }
  // Remove trailing spaces
  size_t i = s.length();
  while ( (i>0) && (s.at(i-1)==' ') ) {
    s.erase(i-1,1);
    i--;
  }  
  // Remove double spaces
  i=1;
  while (i+2<s.length()){
    if ( (s.at(i)==' ') && (s.at(i+1)==' ') ){
      s.erase(i,1);
    }
    else{
      i++;
    }
  }
  return s;
}//




/*-------------------------------*/
/*  test if a file is existing   */
/*-------------------------------*/
bool SGTELIB::exists (const std::string & file) {
  struct stat buffer;   
  return (stat (file.c_str(), &buffer) == 0); 
}//


/*-------------------------------*/
/* count_words                  */
/*-------------------------------*/
int SGTELIB::count_words( const std::string & s ) {
    int word_count( 0 );
    std::stringstream ss( s );
    std::string word;
    while( ss >> word ) ++word_count;
    return word_count;
}//


/*-------------------------------*/
/* add string on a new line of   */
/*           an existing files   */
/*-------------------------------*/
void SGTELIB::append_file (const std::string & s , const std::string & file){
  std::string dummy_str;
  std::string cmd;
  if ( ! SGTELIB::exists(file)){
    cmd = "touch "+ file;
    dummy_str = system( cmd.c_str() );
  }
  cmd = "echo "+s+" >> "+file;
  //std::cout << "cmd : " << cmd << "\n";
  dummy_str = system( cmd.c_str() );
}//


/*-------------------------------*/
/*  wait                         */
/*-------------------------------*/
void SGTELIB::wait (double t) {
  // t is a number of seconds
#ifdef _MSC_VER
//    Sleep(t*1000000.0);
#else
   usleep(t*1000000.0);
#endif
}//



/*-------------------------------*/
/*  isdef (not nan nor inf)     */
/*-------------------------------*/
bool SGTELIB::isdef ( const double x ) {
  if ( isnan(x) ) return false;
  if ( isinf(x) ) return false;
  if ( fabs(x)>=SGTELIB::INF) return false;
  if ( fabs(x)>=1e+16){
    return false;
  }
  return true;
}



/*-------------------------------*/
/*  distance between two points  */
/*-------------------------------*/
double SGTELIB::dist ( const double * x , const double * y , int n ) {
  double s = 0.0;
  double d = 0;
  for ( int i = 0 ; i < n ; ++i )
  {
    d = x[i]-y[i];
    s += d*d;
  }
  return sqrt(s);
}

/*------------------*/
/*  relative error  */
/*------------------*/
double SGTELIB::rel_err ( double x , double y ) {
  if ( x*y < 0.0 )
    return 1.0;
  if ( x*y == 0.0 )
    return fabs(x-y);
  double max = fabs(x);
  if ( fabs(y) > max )
    max = fabs(y);
  return ( fabs ( x-y ) / max );
}

/*---------------------------------*/
int SGTELIB::round ( double d ) {
/*---------------------------------*/
  return static_cast<int>(d < 0.0 ? -floor(.5-d) : floor(.5+d));
}

/*------------------------------------------*/
std::string SGTELIB::itos ( int i ) {
/*------------------------------------------*/
  std::ostringstream oss;
  oss << i;
  return oss.str();
}//

/*------------------------------------------*/
std::string SGTELIB::dtos (double d ) {
/*------------------------------------------*/
  std::ostringstream oss;
  oss << d;
  return oss.str();
}//

/*------------------------------------------*/
std::string SGTELIB::btos (bool b ) {
/*------------------------------------------*/
  if (b) return "True";
  else   return "False";
}//

/*------------------------------------------*/
double SGTELIB::stod ( const std::string & s ){
/*------------------------------------------*/
  double d = atof(s.c_str()); 
  return d;
}//

/*------------------------------------------*/
int SGTELIB::stoi ( const std::string & s ){
/*------------------------------------------*/
  int d = atoi(s.c_str()); 
  return d;
}//

/*----------------------------------------------------------*/
bool SGTELIB::stob ( const std::string & s ) {
/*----------------------------------------------------------*/
  std::string ss = toupper(s);
  if ( ss=="TRUE" ) return true;
  if ( ss=="FALSE") return false;
  if ( ss=="YES"  ) return true;
  if ( ss=="NO"   ) return false;
  if ( ss=="1"    ) return true;
  if ( ss=="0"    ) return false;
  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unrecognised string \""+s+"\" ( "+ss+" )" );
}//

/*------------------------------------------*/
bool SGTELIB::isdigit ( const std::string & s ){
/*------------------------------------------*/
  std::string::const_iterator it = s.begin();
  char c;
  while (it != s.end()){
    c = *it;
    if ( ! ( ( isdigit(std::string(1,c))) || (c=='+') || (c=='-') || (c=='.') ) ){
      return false;
    }
    it++;
  }
  return true;
}//




/*-----------------------------------------------------------------*/
/*                         NOMAD::toupper - 1/2                    */
/*-----------------------------------------------------------------*/
std::string SGTELIB::toupper ( const std::string & s )
{
  std::string s2(s);
  size_t ns = s2.size();
  for ( size_t i = 0 ; i < ns ; ++i )
    s2[i] = std::toupper(s2[i]);
  return s2;
}//


/*----------------------------------------------------------*/
std::string SGTELIB::model_output_to_str ( const SGTELIB::model_output_t ot ) {
/*----------------------------------------------------------*/
  switch ( ot ) {
  case SGTELIB::NORMAL_OUTPUT: return "normal";
  case SGTELIB::FIXED_OUTPUT : return "fixed";
  case SGTELIB::BINARY_OUTPUT: return "binary";
  default:
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
}//

/*----------------------------------------------------------*/
std::string SGTELIB::model_type_to_str ( const SGTELIB::model_t t ) {
/*----------------------------------------------------------*/
  switch ( t ) {
  case SGTELIB::LINEAR   : return "LINEAR";
  case SGTELIB::TGP      : return "TGP";
  case SGTELIB::DYNATREE : return "DYNATREE";
  case SGTELIB::KS       : return "KS";
  case SGTELIB::CN       : return "CN";
  case SGTELIB::PRS      : return "PRS";
  case SGTELIB::PRS_EDGE : return "PRS_EDGE";
  case SGTELIB::PRS_CAT  : return "PRS_CAT";
  case SGTELIB::RBF      : return "RBF";
  case SGTELIB::KRIGING  : return "KRIGING";
  case SGTELIB::SVN      : return "SVN";
  case SGTELIB::LOWESS   : return "LOWESS";
  case SGTELIB::ENSEMBLE : return "ENSEMBLE";
  default:
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
}//

/*----------------------------------------------------------*/
std::string SGTELIB::distance_type_to_str ( const SGTELIB::distance_t t ) {
/*----------------------------------------------------------*/
  switch ( t ) {
    case SGTELIB::DISTANCE_NORM2      : return "NORM2";
    case SGTELIB::DISTANCE_NORM1      : return "NORM1";
    case SGTELIB::DISTANCE_NORMINF    : return "NORMINF";
    case SGTELIB::DISTANCE_NORM2_IS0  : return "NORM2_IS0";
    case SGTELIB::DISTANCE_NORM2_CAT  : return "NORM2_CAT";
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
}//


/*----------------------------------------------------------*/
std::string SGTELIB::weight_type_to_str ( const SGTELIB::weight_t cb ) {
/*----------------------------------------------------------*/
  switch (cb){
    case SGTELIB::WEIGHT_SELECT : return "SELECT";
    case SGTELIB::WEIGHT_OPTIM  : return "OPTIM" ;
    case SGTELIB::WEIGHT_WTA1   : return "WTA1"  ;
    case SGTELIB::WEIGHT_WTA3   : return "WTA3"  ;
    case SGTELIB::WEIGHT_EXTERN : return "EXTERN";
    default :
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
}//



/*----------------------------------------------------------*/
std::string SGTELIB::bbo_type_to_str ( SGTELIB::bbo_t bbot ) {
/*----------------------------------------------------------*/
  switch ( bbot ) {
    case SGTELIB::BBO_OBJ: return "OBJ";
    case SGTELIB::BBO_CON: return "CON";
    case SGTELIB::BBO_DUM: return "DUM";
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }
}//

/*----------------------------------------------------------*/
SGTELIB::model_t SGTELIB::str_to_model_type ( const std::string & s ) {
/*----------------------------------------------------------*/
  std::string ss = SGTELIB::toupper(s);
  if ( ss=="LINEAR"         ){ return SGTELIB::LINEAR; }
  if ( ss=="TGP"            ){ return SGTELIB::TGP; }
  if ( ss=="DYNATREE"       ){ return SGTELIB::DYNATREE; }
  if ( ss=="KS"             ){ return SGTELIB::KS; }
  if ( ss=="CN"             ){ return SGTELIB::CN; }
  if ( ss=="PRS"            ){ return SGTELIB::PRS; }
  if ( ss=="PRS_EDGE"       ){ return SGTELIB::PRS_EDGE; }
  if ( ss=="PRS_CAT"        ){ return SGTELIB::PRS_CAT; }
  if ( ss=="RBF"            ){ return SGTELIB::RBF; }
  if ( ss=="KRIGING"        ){ return SGTELIB::KRIGING; }
  if ( ss=="SVN"            ){ return SGTELIB::SVN; }
  if ( ss=="LWR"            ){ return SGTELIB::LOWESS; }
  if ( ss=="LOWESS"         ){ return SGTELIB::LOWESS; }
  if ( ss=="ENSEMBLE"       ){ return SGTELIB::ENSEMBLE; }
  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unrecognised string \""+s+"\" ( "+ss+" )" );
}//

/*----------------------------------------------------------*/
SGTELIB::weight_t SGTELIB::str_to_weight_type ( const std::string & s ) {
/*----------------------------------------------------------*/
  std::string ss = SGTELIB::toupper(s);
  if ( ss=="SELECT" ){ return SGTELIB::WEIGHT_SELECT;}
  if ( ss=="OPTIM"  ){ return SGTELIB::WEIGHT_OPTIM; }
  if ( ss=="WTA1"   ){ return SGTELIB::WEIGHT_WTA1;  }
  if ( ss=="WTA3"   ){ return SGTELIB::WEIGHT_WTA3;  }
  if ( ss=="EXTERN" ){ return SGTELIB::WEIGHT_EXTERN;}
  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unrecognised string \""+s+"\" ( "+ss+" )" );
}//

/*----------------------------------------------------------*/
SGTELIB::distance_t SGTELIB::str_to_distance_type ( const std::string & s ) {
/*----------------------------------------------------------*/
  std::string ss = SGTELIB::toupper(s);
  if ( ss=="NORM2"    ){ return SGTELIB::DISTANCE_NORM2; }
  if ( ss=="NORM1"    ){ return SGTELIB::DISTANCE_NORM1; }
  if ( ss=="NORMINF"  ){ return SGTELIB::DISTANCE_NORMINF; }

  if ( ss=="ISO"      ){ return SGTELIB::DISTANCE_NORM2_IS0; }
  if ( ss=="IS0"      ){ return SGTELIB::DISTANCE_NORM2_IS0; }
  if ( ss=="NORM2_ISO"){ return SGTELIB::DISTANCE_NORM2_IS0; }
  if ( ss=="NORM2_IS0"){ return SGTELIB::DISTANCE_NORM2_IS0; }

  if ( ss=="CAT"      ){ return SGTELIB::DISTANCE_NORM2_CAT; }
  if ( ss=="NORM2_CAT"){ return SGTELIB::DISTANCE_NORM2_CAT; }
  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unrecognised string \""+s+"\" ( "+ss+" )" );
}//

/*----------------------------------------------------------*/
SGTELIB::distance_t SGTELIB::int_to_distance_type ( const int i ) {
/*----------------------------------------------------------*/
  if ( (i<0) || (i>=SGTELIB::NB_DISTANCE_TYPES) ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
      "int_to_distance_type: invalid integer "+itos(i) );
  }
  switch ( i ){
    case 0: return SGTELIB::DISTANCE_NORM2; 
    case 1: return SGTELIB::DISTANCE_NORM1; 
    case 2: return SGTELIB::DISTANCE_NORMINF; 
    case 3: return SGTELIB::DISTANCE_NORM2_IS0; 
    case 4: return SGTELIB::DISTANCE_NORM2_CAT; 
    default:
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
        "int_to_kernel_type: invalid integer "+itos(i) );
  }
}//




/*----------------------------------------------*/
/*       Same sign                              */
/*----------------------------------------------*/
bool SGTELIB::same_sign (const double a, const double b) {
  return ( (a*b>0) || ( (fabs(a)<EPSILON) && (fabs(b)<EPSILON) ) );
}//


/*----------------------------------------*/
/*  Compute CUMULATIVE Density Function   */
/*  (Centered & Normalized Gaussian law)  */
/*----------------------------------------*/
double SGTELIB::normcdf ( double x ) {
  double t , t2 , v , Phi;
  if (fabs(x)<EPSILON){
    Phi = 0.5;
  }
  else{
    t = 1.0 / ( 1.0 + 0.2316419 * fabs(x) );
    t2 = t*t;
    v = exp(-x*x/2.0)*t*(0.319381530-0.356563782*t+1.781477937*t2-1.821255978*t*t2+1.330274429*t2*t2)/2.506628274631;
    Phi = (x<0.0)?v:1.0-v;
  }
  return Phi;
}//

/*----------------------------------------*/
/*  Compute CUMULATIVE Density Function   */
/*  (Gaussian law)                        */
/*----------------------------------------*/
double SGTELIB::normcdf ( double x , double mu , double sigma ) {
  if (sigma<-EPSILON){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Surrogate_Utils::normpdf: sigma is <0" );
  } 
  // Apply lower bound to sigma
  if (APPROX_CDF){
    sigma = std::max(sigma,EPSILON);
  }
  // Compute CDF
  if (sigma<EPSILON){
    // The cdf is an Heavyside function
    return (x>mu)?1.0:0.0;
  }
  else{
    // Normal case
    return normcdf( (x-mu)/sigma );
  }
}//


/*----------------------------------------*/
/*  Compute PROBABILITY Density Function  */
/*  (Centered & Normalized Gaussian law)  */
/*----------------------------------------*/
double SGTELIB::normpdf ( double x ) {
  return 0.398942280401*exp(-0.5*x*x);
}//

/*----------------------------------------*/
/*  Compute PROBABILITY Density Function  */
/*  (Gaussian law)                        */
/*----------------------------------------*/
double SGTELIB::normpdf ( double x , double mu , double sigma ) {
  if (sigma<EPSILON){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Surrogate_Utils::normpdf: sigma is NULL" );
  } 
  return normpdf( (x-mu)/sigma )/sigma;
}//

/*----------------------------------------*/
/*  Compute EI (expected improvement)     */
/*  (Gaussian law)                        */
/*----------------------------------------*/
double SGTELIB::normei ( double fh , double sh , double f_min ) {
  if (sh<-EPSILON){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Surrogate_Utils::normei: sigma is <0" );
  } 
  // Apply lower bound to sigma
  if (APPROX_CDF){
    sh = std::max(sh,EPSILON);
  }
  if (sh<EPSILON){
    // If there is no uncertainty, then:
    //    - fh<f_min => EI = f_min-fh
    //    - fh>f_min => EI = 0
    return (fh<f_min)?(f_min-fh):0;
  }
  else{
    // Normal case
    double d = (f_min-fh)/sh;
    return (f_min-fh)*normcdf(d) + sh*normpdf(d);
  }
}//

/*----------------------------------------*/
/*  CDF of gamma distribution             */
/*----------------------------------------*/
double SGTELIB::gammacdf(double x, double a, double b){
  // a : shape coef
  // b : scale coef
  if ( (a<=0) || (b<=0) ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Surrogate_Utils::gammacdf: a or b is <0" );
  }  
  if (x<EPSILON) return 0.0;

  return lower_incomplete_gamma(x/b,a);
}//

/*----------------------------------------*/
/*  Inverse CDF of gamma distribution             */
/*----------------------------------------*/
double SGTELIB::gammacdfinv(double f, double a, double b){
  // a : shape coef
  // b : scale coef
  if ( (a<=0) || (b<=0) ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Surrogate_Utils::gammacdfinv: a or b is <=0" );
  }  
  if ( (f<0) || (f>1) ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Surrogate_Utils::gammacdfinv: f<0 or f>1" );
  }  
  if (f==1.0) return INF;
  if (f==0.0) return 0;

  //std::cout << "f,a,b : " << f << " " << a << " " << b << "\n";

  double xmin = 0;
  double xmax = 1;
  double xtry;
  // Extend upper bound
int k = 0;
  while (SGTELIB::gammacdf(xmax,a,b)<f){
    xmin = xmax;
    xmax*=2.0;
    //std::cout << "up " << xmax << " " << a << " " << b << " " <<  SGTELIB::gammacdf(xmax,a,b) << "\n";
    k++;
    if (k>10) break;
  }

  while (xmax-xmin>10000*EPSILON){
    xtry = (xmin+xmax)/2.0;
    if (SGTELIB::gammacdf(xtry,a,b)>f) xmax = xtry;
    else xmin = xtry;
    //std::cout << "dichotomie : " << xtry << "\n";
  }
  return (xmin+xmax)/2.0;
}//

/*----------------------------------------*/
/*  lower incomplete gamma function       */
/*----------------------------------------*/
//  See:
//  Milton Abramowitz, Irene Stegun,Handbook of Mathematical Functions, National Bureau of Standards, 1964.
//  Stephen Wolfram, The Mathematica Book,Fourth Edition,Cambridge University Press, 1999.
double SGTELIB::lower_incomplete_gamma ( const double x, double p ){
  // Special cases
  if ( ( x < EPSILON ) || ( p < EPSILON ) ) return 0;

#ifdef _MSC_VER
#if ( _MSC_VER <= 1600 )
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                                "Surrogate_Utils:: lgamma function not supported with VisualStudio 2010 or lower " );

#endif
#endif
  double f = exp( p * log ( x ) - lgamma ( p + 1.0 ) - x );

  double dv = 1.0, v = 1.0;
  while (dv > v/1e+9) {
    dv *= x / (++p);
    v += dv;
  }
  return v*f;
    
}//

/*----------------------------------------*/
/*  difference between two timeval, in ms */
/*----------------------------------------*/
int SGTELIB::diff_ms(timeval t1, timeval t2){
  return static_cast<int>((((t1.tv_sec - t2.tv_sec) * 1000000) + (t1.tv_usec - t2.tv_usec +500))/1000);
}//

/*----------------------------------------*/
/*  uniform rand generator               */
/*----------------------------------------*/
double SGTELIB::uniform_rand (void){
  return double(rand() / double(INT_MAX));
}//

/*----------------------------------------*/
/*  quick gaussian random generator       */
/*----------------------------------------*/
double SGTELIB::quick_norm_rand (void){
  const int N = 24;
  double d = 0;
  for (int i=1 ; i<N ; i++) d+= uniform_rand();
  d -= double(N)/2.0;
  d *= sqrt(12.0/double(N));
  return d;
}//

/*----------------------------------------*/
/*  relative ceil                         */
/*----------------------------------------*/
double SGTELIB::rceil (double d){
  if (d>0) return std::ceil(d);
  else if (d<0) return std::floor(d);
  else return 0.0;
}//






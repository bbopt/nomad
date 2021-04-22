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

#include "TrainingSet.hpp"
using namespace SGTELIB;

/*--------------------------------------*/
/*              constructor             */
/*--------------------------------------*/
SGTELIB::TrainingSet::TrainingSet ( const Matrix & X ,
                                    const Matrix & Z  ) :
  _p            ( X.get_nb_rows()   ) , // Nb of points
  _n            ( X.get_nb_cols()   ) , // Nb of input
  _m            ( Z.get_nb_cols()   ) , // Nb of output
  _ready        ( false             ) , 
  _bbo          ( new bbo_t [_m]    ) , // Types of output (OBJ, CON or DUM)
  _bbo_is_def   ( false             ) , // Allows to know if _bbo has been def
  _j_obj        ( 0                 ) , // Index of the output that represents the objective
  _f_min        ( INF               ) ,   
  _fs_min       ( INF               ) , 
  _i_min        ( 0                 ) , // Index of the point where f_min is reached
  _X            ( X                 ) , // Input Data
  _Z            ( Z                 ) , // Output Data
  _Xs           ( "TrainingSet._Xs" , _p , _n ) , // Normalized Input Data
  _Zs           ( "TrainingSet._Zs" , _p , _m ) , // Normalized Output Data
  _Ds           ( "TrainingSet._Ds" , _p , _p ) , // Nb of varying input
  _nvar         ( -1                ) , // Nb of varying output
  _mvar         ( -1                ) ,
  _pvar         ( -1                ) ,
  _X_lb         ( new double   [_n] ) ,
  _X_ub         ( new double   [_n] ) ,
  _X_scaling_a  ( new double   [_n] ) ,
  _X_scaling_b  ( new double   [_n] ) ,
  _X_mean       ( new double   [_n] ) ,
  _X_std        ( new double   [_n] ) ,
  _X_nbdiff     ( new int      [_n] ) ,
  _X_nbdiff1    ( 0                 ) ,
  _X_nbdiff2    ( 0                 ) ,
  _Z_lb         ( new double   [_m] ) ,
  _Z_ub         ( new double   [_m] ) ,
  _Z_replace    ( new double   [_m] ) ,
  _Z_scaling_a  ( new double   [_m] ) ,
  _Z_scaling_b  ( new double   [_m] ) ,
  _Z_mean       ( new double   [_m] ) , // Mean of each output
  _Z_std        ( new double   [_m] ) ,
  _Zs_mean      ( new double   [_m] ) , // Mean of each normalized output
  _Z_nbdiff     ( new int      [_m] ) ,
  _Ds_mean      ( 0.0               ) {

  // Init bounds
  for (int i=0 ; i<_n ; i++){
    _X_lb[i] = 0;
    _X_ub[i] = 0;
  }  
  for (int j=1 ; j<_m ; j++){
    _Z_lb[j] = 0;
    _Z_ub[j] = 0;
  }

  // Init the _bbo with standard values:
  // First is the objective,
  // Then constraints
  _bbo[0] = BBO_OBJ;
  for (int j=1 ; j<_m ; j++){
    _bbo[j] = BBO_CON;
    _Z_lb[j] = 0;
    _Z_ub[j] = 0;
  }

}//





/*---------------------------*/
/*      copy constructor     */
/*---------------------------*/

SGTELIB::TrainingSet::TrainingSet ( const TrainingSet & C ) : 
              _p    ( C._p ) ,
              _n    ( C._n ) ,
              _m    ( C._m ) {
  info();
  throw Exception ( __FILE__ , __LINE__ ,
       "TrainingSet: copy constructor forbiden." );

}//




/*--------------------------------------*/
/*               info (debug)           */
/*--------------------------------------*/
void SGTELIB::TrainingSet::info (void) const{
  std::cout << "   ## ## TrainingSet::info  " << this << " " << _ready << " " << _p << "\n";
}

/*--------------------------------------*/
/*               destructor             */
/*--------------------------------------*/
SGTELIB::TrainingSet::~TrainingSet ( void ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "   ## ## Delete TrainingSet " << this << "\n";
  #endif
  delete [] _bbo;
  delete [] _X_lb;
  delete [] _X_ub;
  delete [] _X_scaling_a;
  delete [] _X_scaling_b;
  delete [] _X_mean;
  delete [] _X_std;
  delete [] _X_nbdiff;
  delete [] _Z_lb;
  delete [] _Z_ub;
  delete [] _Z_replace;
  delete [] _Z_scaling_a;
  delete [] _Z_scaling_b;
  delete [] _Z_mean;
  delete [] _Z_std;
  delete [] _Zs_mean;
  delete [] _Z_nbdiff;
}//

/*--------------------------------------*/
/*               operator =             */
/*--------------------------------------*/
SGTELIB::TrainingSet & SGTELIB::TrainingSet::operator = ( const SGTELIB::TrainingSet & A ) {
  A.info();
  throw Exception ( __FILE__ , __LINE__ ,
       "TrainingSet: operator \"=\" forbiden." );
  return *this;
}


/*--------------------------------------*/
/*          Set bbo type                */
/*--------------------------------------*/
void SGTELIB::TrainingSet::set_bbo_type (const std::string & line){
  // BBOT must be separated by space

  if (_bbo_is_def){
      throw Exception ( __FILE__ , __LINE__ ,
           "TrainingSet::set_bbo_type: _bbo must be defined before the first build." );
  }
  #ifdef SGTELIB_DEBUG
    std::cout << "SGTELIB::TrainingSet::set_bbo_type\n";
    std::cout << "Input string: \"" << line << "\"\n";
  #endif 

  std::string s;
  std::istringstream in_line (line);		
  
  int j = 0;
  while (in_line.good()){
  	in_line >> s;
    if (j>=_m){
      throw Exception ( __FILE__ , __LINE__ ,
           "TrainingSet::set_bbo_type: string error (j>_m)" );
    }
    if ( (streqi(s,"OBJ")) || (streqi(s,"O")) ){
      _bbo[j] = BBO_OBJ;
      _j_obj = j;
    }
    else if ( (streqi(s,"CON")) || (streqi(s,"C")) ){
      _bbo[j] = BBO_CON;
    }
    else if ( (streqi(s,"DUM")) || (streqi(s,"D")) ){
      _bbo[j] = BBO_DUM;
    }
    else{
      throw Exception ( __FILE__ , __LINE__ ,
           "TrainingSet::set_bbo_type: string error (string "+s+" not recognized)" );
    }
    j++;
  }

  // Check the number of OBJ
  double n = 0;
  for (j=0 ; j<_m ; j++){
    if (_bbo[j] == BBO_OBJ){
      n++;
    }
  }
  if (n>1){
    throw Exception ( __FILE__ , __LINE__ ,
           "TrainingSet::set_bbo_type: Only one obj is allowed" );
  }

  // Check the number of OBJ+CON
  for (j=0 ; j<_m ; j++){
    if (_bbo[j] == BBO_CON){
      n++;
    }
  }
  if (n==0){
    throw Exception ( __FILE__ , __LINE__ ,
           "TrainingSet::set_bbo_type: all outputs are \"DUM\"" );
  }

  #ifdef SGTELIB_DEBUG
    std::cout << "Output types:\n";
    for (j=0 ; j<_m ; j++){
      std::cout << j << ": " << bbo_type_to_str(_bbo[j]) << "\n";
    }
  #endif

  // nb: this method is supposed to be used only once, in the begining.
  // Plus, it may have big changes on the results, so this changes the
  // trainingset to "not ready".
  _fs_min = INF;
  _f_min  = INF;
  _bbo_is_def = true;
  _ready = false;
}//

/*--------------------------------------*/
/*          Construct                   */
/*--------------------------------------*/
void SGTELIB::TrainingSet::build ( void ){



  // check the dimensions:
  if ( _X.get_nb_rows() != _Z.get_nb_rows() )
    throw Exception ( __FILE__ , __LINE__ ,
             "TrainingSet::build(): dimension error" );

  // Check number of points
  if ( _p < 1 )
    throw Exception ( __FILE__ , __LINE__ ,
             "TrainingSet::build(): empty training set");

  if ( ! _ready){
    #ifdef SGTELIB_DEBUG
      std::cout << "TrainingSet::build BEGIN, X:(" << _p << "," << _n << ") Z:(" << _p << "," << _m << ")\n";
    #endif

    // Compute the number of varying input and output
    compute_nbdiff(_X,_X_nbdiff,_nvar);
    compute_nbdiff(_Z,_Z_nbdiff,_mvar);

    // Compute the number of input dimension for which 
    // nbdiff is greater than 1 (resp. greater than 2).
    _X_nbdiff1 = 0;
    _X_nbdiff2 = 0;
    for (int j=0 ; j<_n ; j++){
      if (_X_nbdiff[j]>1) _X_nbdiff1++;
      if (_X_nbdiff[j]>2) _X_nbdiff2++;
    }

    // Check singular data (inf and void)
    if (check_singular_data())
    {
    #ifdef SGTELIB_DEBUG
        std::cout << "SGTELIB::TrainingSet::check_singular_data(): incorrect data. Some _Z has no defined value !\n";
    #endif
        return;
    }

    // Compute bounds over columns of X and Z
    compute_bounds();

    // Compute scaling values
    compute_scaling();

    // Compute scaled matrices
    compute_scaled_matrices();

    // Build matrix of distances between each pair of points
    compute_Ds();

    // Compute fs_min
    compute_f_min();

    // The training set is now ready!
    _ready = true;

    #ifdef SGTELIB_DEBUG
      std::cout << "TrainingSet::build END\n";
    #endif

  }

  // _bbo is considered as defined. It can not be modified anymore.
  _bbo_is_def = true;

}//


/*--------------------------------------*/
/*  Check if the training set is ready  */
/*--------------------------------------*/
void SGTELIB::TrainingSet::check_ready (void) const{
  if ( ! _ready){
    std::cout << "TrainingSet: NOT READY!\n";
    throw Exception ( __FILE__ , __LINE__ , "TrainingSet::check_ready(): TrainingSet not ready. Use method TrainingSet::build()" );
  }
}//
/*--------------------------------------*/
void SGTELIB::TrainingSet::check_ready (const std::string & file,
                                                  const std::string & function,
                                                  const int & i        ) const {
    check_ready(file+"::"+function+":"+itos(i));
}//
/*--------------------------------------*/
void SGTELIB::TrainingSet::check_ready (const std::string & s) const{
  if ( ! _ready){
    std::cout << "TrainingSet: NOT READY! (" << s << ")\n";
    throw Exception ( __FILE__ , __LINE__ , "TrainingSet::check_ready(): TrainingSet not ready. Use method TrainingSet::build()" );
  }
}//


/*--------------------------------------*/
/*                 add_points           */
/*--------------------------------------*/
bool SGTELIB::TrainingSet::add_points ( const Matrix & Xnew ,
                                        const Matrix & Znew ) {

  // Check dim
  if ( Xnew.get_nb_rows() != Znew.get_nb_rows() || Xnew.get_nb_cols() != _n || Znew.get_nb_cols() != _m ){
    throw Exception ( __FILE__ , __LINE__ , "TrainingSet::add_points(): dimension error" );
  }

  // Check for nan
  if (Xnew.has_nan()){
    throw Exception ( __FILE__ , __LINE__ , "TrainingSet::add_points(): Xnew is nan" );
  }
  if (Znew.has_nan()){
    throw Exception ( __FILE__ , __LINE__ , "TrainingSet::add_points(): Znew is nan" );
  }

  // Add the points in the trainingset
  _X.add_rows(Xnew);
  _Z.add_rows(Znew);

  // Add empty rows
  int pnew = Xnew.get_nb_rows();
  _Xs.add_rows(pnew);
  _Zs.add_rows(pnew);
  _Ds.add_rows(pnew);
  _Ds.add_cols(pnew);
  // Update p
  _p += pnew;
  // Note that the trainingset needs to be updated.
  _ready = false;

  return true;
}//

/*--------------------------------------*/
/*                 add_points           */
/*--------------------------------------*/
bool SGTELIB::TrainingSet::partial_reset_and_add_points ( const Matrix & Xnew ,
                                        const Matrix & Znew )
{
    
    if ( Xnew.get_nb_cols()!= _n || Znew.get_nb_cols() != _m)
        return false;
    
    _X = Xnew ;
    _Z = Znew;
    
    _p = Xnew.get_nb_rows();
    
    _Xs = SGTELIB::Matrix( "TrainingSet._Xs" , _p , _n );
    _Zs = SGTELIB::Matrix( "TrainingSet._Zs" , _p , _m );
    _Ds = SGTELIB::Matrix( "TrainingSet._Ds" , _p , _p );

    _ready = false;
    return true;
}
    
/*--------------------------------------*/
/*                 add_point            */
/*--------------------------------------*/
bool SGTELIB::TrainingSet::add_point ( const double * xnew ,
                                       const double * znew ) {

  return add_points ( Matrix::row_vector ( xnew , _n ),
                      Matrix::row_vector ( znew , _m ) );
}//

/*---------------------------------------------------*/
/*  compute the mean and std over                    */
/*  the columns of a matrix                          */
/*---------------------------------------------------*/
bool SGTELIB::TrainingSet::check_singular_data ( void ){
    
    int i,j;
    bool e = false;
    // Check that all the _X data are defined
    for ( j = 0 ; j < _n ; j++ ) {
        for ( i = 0 ; i < _p ; i++ ) {
            if ( ! isdef(_X.get(i,j))){
                std::cout << "_X(" << i << "," << j << ") = " << _X.get(i,j) << "\n";
                e = true;
            }
        }
    }
    
    // Check that, for each output index, SOME data are defined
    bool isdef_Zj; // True if at least one value is defined for output j.
    // Loop on the output indexes
    for ( j = 0 ; j < _m ; j++ ) {
        // no def value so far
        isdef_Zj = false;
        for ( i = 0 ; i < _p ; i++ ) {
            if (isdef(_Z.get(i,j))){
                isdef_Zj = true;
                break;
            }
        }
        // if there is more than 10 points and no correct value was found, return flag error.
        if ( (_p>10) && ( ! isdef_Zj) ){
            e = true;
        }
    }
    
    return e;
    
}//

/*---------------------------------------------------*/
/*  compute the mean and std over                    */
/*  the columns of a matrix                          */
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::compute_mean_std ( void ){

  int i,j;
  double v, mu, var;
  // Loop on the inputs
  for ( j=0 ; j<_n ; j++ ) {
    // Loop on lines for MEAN computation
    mu = 0;
    for ( i=0 ; i<_p ; i++ ) {
      mu += _X.get(i,j);
    }
    mu /= _p; 
    _X_mean[j] = mu;
    // Loop on lines for VAR computation
    var = 0;
    for ( i=0 ; i<_p ; i++ ) {
      v = _X.get(i,j);
      var += (v-mu)*(v-mu);
    }
    var /= (_p-1);
    _X_std[j] = sqrt(var);
  }

  // Loop on the outputs
  for ( j=0 ; j<_m ; j++ ) {
    // Loop on lines for MEAN computation
    mu = 0;
    for ( i=0 ; i<_p ; i++ ) {
      v = _Z.get(i,j);
      if ( ! isdef(v)) v = _Z_replace[j];
      mu += v;
    }
    mu /= _p; 
    _Z_mean[j] = mu;
    // Loop on lines for VAR computation
    var = 0;
    for ( i=0 ; i<_p ; i++ ) {
      v = _Z.get(i,j);
      if ( ! isdef(v)) v = _Z_replace[j];
      var += (v-mu)*(v-mu);
    }
    var /= (_p-1);
    _Z_std[j] = sqrt(var);
  }

}//

/*---------------------------------------------------*/
/*  compute the bounds over                          */
/*  the columns of a matrix                          */
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::compute_bounds ( void ){

  int i,j;
  double v = 0;

  // Bound of X
  for ( j=0 ; j<_n ; j++ ) {
    _X_lb[j] = +INF;
    _X_ub[j] = -INF;
    // Loop on points
    for ( i=0 ; i<_p ; i++ ) {
      v = _X.get(i,j);
      _X_lb[j] = std::min(v,_X_lb[j]);
      _X_ub[j] = std::max(v,_X_ub[j]);
    }
  }

  // Bound of Z
  for ( j=0 ; j<_m ; j++ ) {
    _Z_lb[j] = +INF;
    _Z_ub[j] = -INF;
    // Loop on points
    for ( i=0 ; i<_p ; i++ ) {
      v = _Z.get(i,j);
      if ( isdef(v) ){
        _Z_lb[j] = std::min(v,_Z_lb[j]);
        _Z_ub[j] = std::max(v,_Z_ub[j]);
      }
    }

    // Compute replacement value for undef Z
    // If there are no correct bounds defined yet
    if ( ( ! isdef(_Z_lb[j])) || ( ! isdef(_Z_ub[j])) ){
      _Z_replace[j] = 1.0;
    }
    else{
      _Z_replace[j] = std::max(_Z_ub[j],0.0) + 0.1*std::max(_Z_ub[j]-_Z_lb[j],1.0);
    }

  }

}//


/*---------------------------------------------------*/
/*  compute the number of different values over      */
/*  the columns of a matrix                          */
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::compute_nbdiff ( const Matrix & MAT , 
                                            int * nbdiff,
                                            int & njvar ){
  
  int nj = MAT.get_nb_cols(); // nb of columns
  njvar = 0; // nb of columns that are not constant
  for ( int j = 0 ; j < nj ; j++ ){
    nbdiff[j] = MAT.get_nb_diff_values(j);  
    if (nbdiff[j]>1) njvar++;
  }
}//







/*---------------------------------------------------*/
/*  compute the bounds over                          */
/*  the columns of a matrix                          */
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::compute_nvar_mvar ( void ){

  // Compute _nvar
  if (_nvar!=_n){
    _nvar = 0;
    for ( int j = 0 ; j < _n ; j++ )
      if (_X_nbdiff[j] > 1) _nvar++;
  }

  // Compute _mvar
  if (_mvar!=_m){
    _mvar = 0;
    for ( int j = 0 ; j < _m ; j++ )
      if (_Z_nbdiff[j] > 1) _mvar++;
  }
}//



/*---------------------------------------------------*/
/*  compute scaling parameters                       */
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::compute_scaling ( void ){
  int j=0;

  // Neutral values
  for ( j = 0 ; j < _n ; j++ ) {
    _X_scaling_a[j] = 1;
    _X_scaling_b[j] = 0;
  }
  for ( j = 0 ; j < _m ; j++ ) {
    _Z_scaling_a[j] = 1;
    _Z_scaling_b[j] = 0;
  }

  switch (scaling_method){
  case SCALING_NONE:
    //Nothing to do!
    break;
  case SCALING_MEANSTD:
    // Compute mean and std over columns of X and Z
    compute_mean_std();
    // Compute scaling constants
    for ( j = 0 ; j < _n ; j++ ) {
      if (_X_nbdiff[j]>1) _X_scaling_a[j] = 1/_X_std[j];
      _X_scaling_b[j] = -_X_mean[j]*_X_scaling_a[j];
    }
    for ( j = 0 ; j < _m ; j++ ) {
      if (_Z_nbdiff[j]>1) _Z_scaling_a[j] = 1/_Z_std[j];
      _Z_scaling_b[j] = -_Z_mean[j]*_Z_scaling_a[j];
    }
    break;
  case SCALING_BOUNDS:
    // Compute scaling constants
    for ( j = 0 ; j < _n ; j++ ) {
      if (_X_nbdiff[j]>1) _X_scaling_a[j] = 1/(_X_ub[j]-_X_lb[j]);
      _X_scaling_b[j] = -_X_lb[j]*_X_scaling_a[j];
    }
    for ( j = 0 ; j < _m ; j++ ) {
      if (_Z_nbdiff[j]>1)  _Z_scaling_a[j] = 1/(_Z_ub[j]-_Z_lb[j]);
      _Z_scaling_b[j] = -_Z_lb[j]*_Z_scaling_a[j];
    }
    break;
  }// end switch
}//


/*---------------------------------------------------*/
/*  compute scale matrices _Xs and _Zs               */
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::compute_scaled_matrices ( void ){

  double v, mu;
  int i,j;

  // Compute _Xs
  for ( j = 0 ; j < _n ; j++ ){
    for ( i = 0 ; i < _p ; i++ ){
      v = _X.get(i,j)*_X_scaling_a[j]+_X_scaling_b[j];
      _Xs.set(i,j,v);
    }
  }

  // Compute _Zs and Mean_Zs
  for ( j = 0 ; j < _m ; j++ ){
    mu = 0;
    for ( i = 0 ; i < _p ; i++ ){
      v = _Z.get(i,j);
      if ( ! isdef(v)){
        v = _Z_replace[j];
      }
      v = v*_Z_scaling_a[j]+_Z_scaling_b[j];
      mu +=v;
      _Zs.set(i,j,v);
    }
    _Zs_mean[j] = mu/_p;
  }

}//



/*---------------------------------------------------*/
/*  compute distance matrix                          */
/*  the columns of a matrix                          */
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::compute_Ds ( void ){
  double d;
  double di1i2;
  _pvar = _p;
  _Ds_mean = 0.0;
  bool unique;
  for ( int i1 = 0 ; i1 < _p-1 ; i1++ ){
    _Ds.set(i1,i1,0.0);
    unique = true;
    for ( int i2 = i1+1 ; i2 < _p ; i2++ ){
      d = 0;
      for ( int j = 0 ; j < _n ; j++ ){
        di1i2 = _Xs.get(i1,j)-_Xs.get(i2,j);
        d += di1i2*di1i2;
      }
      d = sqrt(d);
      _Ds.set(i1,i2,d);
      _Ds.set(i2,i1,d);
      // Compute the mean distance between the points
      _Ds_mean += d;
      // If d==0, then the point i2 is not unique. 
      if (fabs(d)<EPSILON){
        unique = false;
      } 
    }
    // If there are some points equal to the point of index i2, 
    // then reduce the number of different points.
    if ( ! unique) _pvar--;
  }
  _Ds_mean /= double(_pvar*(_pvar-1)/2);

}//



/*---------------------------------------------------*/
/*  compute fs_min (scaled value of f_min)             */
/*---------------------------------------------------*/
// the lazy way....
void SGTELIB::TrainingSet::compute_f_min ( void ){

  double f;
  bool feasible;  
  // Go through all points
  for ( int i=0 ; i<_p ; i++ ){
    // Get the unscaled objective
    f = _Z.get(i,_j_obj);
    // If objective is good
    if (f<_f_min){
      // check the constraints
      feasible = true;
      for ( int j=0 ; j<_m ; j++ )
        if (_bbo[j]==BBO_CON)
          if (_Z.get(i,j)>0.0){ feasible = false; break; }
      // If the point is feasible, save the value.
      if (feasible){
        _f_min = f;
        _i_min = i;
      }
    }
  }
  // Compute the scaled objective.
  _fs_min = Z_scale( _f_min, _j_obj );

}//


/*---------------------------------------------------*/
/*  get                                              */
/*---------------------------------------------------*/
double SGTELIB::TrainingSet::get_Xs ( const int i , const int j ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(); 
    // Check index
    if ( (i<0) || (i>=_p) || (j<0) || (j>=_n) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::TrainingSet(): dimension error" );
    }
  #endif
  // Return value
  return _Xs.get(i,j);
}//
/*---------------------------------------------------*/
double SGTELIB::TrainingSet::get_Zs ( const int i , const int j ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(); 
    // Check index
    if ( (i<0) || (i>=_p) || (j<0) || (j>=_m) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::TrainingSet(): dimension error" );
    }
  #endif
  // Return value
  return _Zs.get(i,j);
}//
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::get_Xs ( const int i , double * x ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
    // Check index
    if ( (i<0) || (i>=_p) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::TrainingSet(): dimension error" );
    }
  #endif
  // Check/initialize pointer
  if ( ! x){
    x = new double [_n];
  }
  // Fill pointer
  for ( int j = 0 ; j < _n ; j++ ){
    x[j] = _Xs.get(i,j);
  }
}//
/*---------------------------------------------------*/
void SGTELIB::TrainingSet::get_Zs ( const int i , double * z ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
    // Check index
    if ( (i<0) || (i>=_p) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::get_Zs(): dimension error" );
    }
  #endif
  // Check/initialize pointer
  if ( ! z){
    z = new double [_m];
  }
  // Fill pointer
  for ( int j = 0 ; j < _m ; j++ ){
    z[j] = _Zs.get(i,j);
  }
}//
/*---------------------------------------------------*/
double SGTELIB::TrainingSet::get_Zs_mean ( const int j ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
    // Check index
    if ( (j<0) || (j>=_m) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::get_Zs_mean(): dimension error" );
    }
  #endif
  return _Zs_mean[j];
}//


/*---------------------------------------------------*/
int SGTELIB::TrainingSet::get_X_nbdiff ( const int i ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
    // Check index
    if ( (i<0) || (i>=_n) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::get_X_nbdiff(): dimension error" );
    }
  #endif
  return _X_nbdiff[i];
}//

/*---------------------------------------------------*/
const Matrix SGTELIB::TrainingSet::get_X_nbdiff ( void ) const {
  Matrix V ("NbDiff",1,_n);
  for (int j=0 ; j<_n ; j++){
    V.set(0,j,(double)_X_nbdiff[j]);
  }
  return V;
}//


/*---------------------------------------------------*/
int SGTELIB::TrainingSet::get_Z_nbdiff ( const int j ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(__FILE__,__FUNCTION__,__LINE__);
    // Check index
    if ( (j<0) || (j>=_m) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::get_Z_nbdiff(): dimension error" );
    }
  #endif
  return _Z_nbdiff[j];
}//
/*---------------------------------------------------*/
// Return the normalized distance between points i1 an i2
double SGTELIB::TrainingSet::get_Ds ( const int i1 , const int i2 ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(); 
    // Check index
    if ( (i1<0) || (i1>=_p) || (i2<0) || (i2>=_p) ){
      throw Exception ( __FILE__ , __LINE__ ,
               "TrainingSet::get_Ds(): dimension error" );
    }
  #endif
  return _Ds.get(i1,i2);
}//

/*--------------------------------------------------*/
/* compute the distances between two sets of points */
/*--------------------------------------------------*/
Matrix SGTELIB::TrainingSet::get_distances ( const Matrix & A , 
                                                      const Matrix & B , 
                                                      const distance_t dt       ) const{


  switch (dt){

    case DISTANCE_NORM1:
      return Matrix::get_distances_norm1(A,B);

    case DISTANCE_NORM2:
      return Matrix::get_distances_norm2(A,B);

    case DISTANCE_NORMINF:
      return Matrix::get_distances_norminf(A,B);

    case DISTANCE_NORM2_IS0:
      // Two points x and y are in the same "IS0-class" if (x_j==0 <=> y_j==0 for each j).
      // The distance "IS0" between two points of the same IS0-class is the norm 2 distance.
      // The distance "IS0" between two points of different IS0-class is INF.
      {
        const int n = A.get_nb_cols();
        const int pa = A.get_nb_rows();
        const int pb = B.get_nb_rows();
        double v,d;
        int ia, ib, j;
        Matrix D = Matrix::get_distances_norm2(A,B);
        double * x0 = new double [n];
        for (j=0 ; j < n ; j++){
          x0[j] = X_scale( 0.0 , j ); 
        }
        for (ia=0 ; ia < pa ; ia++){
          for (ib=0 ; ib < pb ; ib++){
            // For each value of D
            d = D.get(ia,ib);
            v = d*d;
            for (j=0 ; j < n ; j++){
              // If they are not in the same 0-class
              if (  (fabs(A.get(ia,j)-x0[j])<EPSILON) ^ (fabs(B.get(ib,j)-x0[j])<EPSILON)  ){
                v+=10000;
              }
            }
            v = sqrt(v);
            D.set(ia,ib,v);
          }
        }
        delete [] x0;
        return D;
      }

    case DISTANCE_NORM2_CAT:
      // Two points x and y are in the same "X0-class" if x_0==y_0.
      // The distance "IS0" between two points of the same X0-class is the norm 2 distance.
      // The distance "IS0" between two points of different X0-class is INF.
      {
        const int pa = A.get_nb_rows();
        const int pb = B.get_nb_rows();
        double v,d;
        int ia, ib, j;
        Matrix D = Matrix::get_distances_norm2(A,B);
        j = 0;
        for (ib=0 ; ib < pb ; ib++){
          for (ia=0 ; ia < pa ; ia++){
            // For each value of D
            d = D.get(ia,ib);
            v = d*d;
            // If they are not in the same 0-class
            if (  fabs(A.get(ia,j)-B.get(ib,j))>EPSILON  ) {
              v+=10000;
            }
            v = sqrt(v);
            D.set(ia,ib,v);
          }
        }
        return D;
      }

    default:
      throw Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }

  

}//


/*--------------------------------------*/
/*    X scale: x->y: y = a.x + b        */
/*--------------------------------------*/
void SGTELIB::TrainingSet::X_scale ( double * x ) const {
  for ( int j = 0 ; j < _n ; j++ )
    x[j] = _X_scaling_a[j] * x[j] + _X_scaling_b[j];
}//

double SGTELIB::TrainingSet::X_scale ( double x , int var_index ) const {
  return _X_scaling_a[var_index] * x + _X_scaling_b[var_index];
}//

/*--------------------------------------*/
/*    X unscale: y->x: x = (y-b)/a      */
/*--------------------------------------*/
void SGTELIB::TrainingSet::X_unscale ( double * y ) const {
  for ( int j = 0 ; j < _n ; j++ )
    y[j] = ( y[j] - _X_scaling_b[j] ) / _X_scaling_a[j];
}

double SGTELIB::TrainingSet::X_unscale ( double y , int var_index ) const {
  return ( y - _X_scaling_b[var_index] ) / _X_scaling_a[var_index];
}

/*--------------------------------------*/
/*    Z scale: z->w: w = a.z + b        */
/*--------------------------------------*/
void SGTELIB::TrainingSet::Z_scale ( double * z ) const {
  for ( int j = 0 ; j < _m ; j++ )
    z[j] = _Z_scaling_a[j] * z[j] + _Z_scaling_b[j];
}//

double SGTELIB::TrainingSet::Z_scale ( double z , int j ) const {
  return _Z_scaling_a[j] * z + _Z_scaling_b[j];
}//

/*--------------------------------------*/
/*    Z unscale: w->z: z = (w-b)/a      */
/*--------------------------------------*/
void SGTELIB::TrainingSet::Z_unscale ( double * w ) const {
  for ( int j = 0 ; j < _m ; j++ )
    w[j] = Z_unscale(w[j],j);
}//

double SGTELIB::TrainingSet::Z_unscale ( double w , int j ) const {
  if ( (boolean_rounding) && (_Z_nbdiff[j]==2) ){
    double Zs_middle = 0.0;
    if (boolean_rounding==1){
      // Threshold is midway between the biggest and smallest value of Z;
      Zs_middle = Z_scale ( (_Z_ub[j]+_Z_lb[j])/2.0 , j );
    }
    else if (boolean_rounding==2){
      // Threshold is the mean of Z;  
      Zs_middle = _Zs_mean[j];
    }
    return (w>Zs_middle)?_Z_ub[j]:_Z_lb[j];
  }
  else{
    return ( w - _Z_scaling_b[j] ) / _Z_scaling_a[j];
  }
}//

/*------------------------------------------*/
/*    ZE unscale: w->z: z = (w)/a           */
/* Used to unscale errors, std and EI */
/*------------------------------------------*/
double SGTELIB::TrainingSet::ZE_unscale ( double w , int j ) const {
  return w / _Z_scaling_a[j];
}//


/*--------------------------------------*/
/*    X scale: x->y: y = a.x + b        */
/*--------------------------------------*/
void SGTELIB::TrainingSet::X_scale ( Matrix & X ) {
  int p = X.get_nb_rows();
  int n = X.get_nb_cols();
  if (n!=_n){
    throw Exception ( __FILE__ , __LINE__ ,
                 "TrainingSet::TrainingSet(): dimension error" );
  }
  double v;
  // UnScale the output
  for (int i=0 ; i<p ; i++){
    for (int j=0 ; j<n ; j++){
      // Z 
      v = X.get(i,j);
      v = X_scale ( v , j );
      X.set(i,j,v);
    }
  }
}//

/*--------------------------------------*/
/*    Z unscale: w->z: z = (w-b)/a      */
/*--------------------------------------*/
void SGTELIB::TrainingSet::Z_unscale ( Matrix * Z ) {
  int p = Z->get_nb_rows();
  int m = Z->get_nb_cols();
  if (m!=_m){
    throw Exception ( __FILE__ , __LINE__ ,
                 "TrainingSet::TrainingSet(): dimension error" );
  }
  double v;
  // UnScale the output
  for (int i=0 ; i<p ; i++){
    for (int j=0 ; j<m ; j++){
      // Z 
      v = Z->get(i,j);
      v = Z_unscale ( v , j );
      Z->set(i,j,v);
    }
  }
}//
Matrix SGTELIB::TrainingSet::Z_unscale ( const Matrix & Z ) {
  Matrix Z2 (Z);
  Z_unscale(&Z2);
  return Z2;
}//

/*--------------------------------------*/
/*    ZE unscale: w->z: z = w/a      */
/*--------------------------------------*/
void SGTELIB::TrainingSet::ZE_unscale ( Matrix * ZE ) {
  int p = ZE->get_nb_rows();
  int m = ZE->get_nb_cols();
  if (m!=_m){
    throw Exception ( __FILE__ , __LINE__ ,
                 "TrainingSet::TrainingSet(): dimension error" );
  }
  double v;
  // UnScale the output
  for (int i=0 ; i<p ; i++){
    for (int j=0 ; j<m ; j++){
      // Z 
      v = ZE->get(i,j);
      v = ZE_unscale ( v , j );
      ZE->set(i,j,v);
    }
  }
}//
Matrix SGTELIB::TrainingSet::ZE_unscale ( const Matrix & ZE ) {
  Matrix ZE2 (ZE);
  ZE_unscale(&ZE2);
  return ZE2;
}//


/*--------------------------------------*/
/*    get d1 over d2                    */
/*--------------------------------------*/
double SGTELIB::TrainingSet::get_d1_over_d2 ( const Matrix & XXs ) const {
  if (XXs.get_nb_rows()>1){
    throw Exception ( __FILE__ , __LINE__ ,
         "TrainingSet::get_d1_over_d2: XXs must have only one line." );
  } 
  double d1 = +INF;
  double d2 = +INF;
  double d;
  double dxj;
  int i,i1,j;
  i1 = 0; // Index of the closest point 

  // If only 1 point, it is not possible to compute d2, 
  // so we use a dummy value.
  if (_p==1){
    return 1.0;
  }

  // Parcours des points
  for ( i=0 ; i<_p ; i++ ){

    // Calcul de d
    d = 0.0;
    for ( j=0 ; j<_n ; j++){
      dxj = XXs.get(0,j)-_Xs.get(i,j);
      d += dxj*dxj;
    }
    if (d==0){
      return 0.0;
    }
    if (d<d1){ 
      d2=d1;
      d1=d;
      i1=i;// Memorize index of closest point
    }
    else if ((d<d2) && (_Ds.get(i,i1)>0)){
      // nb: the point i can be kept as 2nd closest point only if it is different from the point 
      // i1, which means that the distance to this point must be non null.
      d2=d;
    }

  }
  return sqrt(d1/d2);
}//


/*--------------------------------------*/
/*    get d1 over d2                    */
/*--------------------------------------*/
double SGTELIB::TrainingSet::get_d1 ( const Matrix & XXs ) const {
  if (XXs.get_nb_rows()>1){
    throw Exception ( __FILE__ , __LINE__ ,
         "TrainingSet::get_d1: XXs must have only one line." );
  } 
  double d;
  double d1 = +INF;
  int i,j;
  double dxj;

  // Parcours des points
  for ( i=0 ; i<_p ; i++ ){

    // Calcul de d
    d = 0.0;
    for ( j=0 ; j<_n ; j++){
      dxj = XXs.get(0,j)-_Xs.get(i,j);
      d += dxj*dxj;
    }
    if (d==0){
      return 0.0;
    }
    if (d<d1){ 
      d1=d;
    }

  }
  return sqrt(d1);
}//


/*--------------------------------------*/
/*       get_exclusion_area_penalty     */
/*--------------------------------------*/
Matrix SGTELIB::TrainingSet::get_exclusion_area_penalty ( const Matrix & XXs , const double tc ) const {
  const int pxx = XXs.get_nb_rows();
  double r12,p;
  //double logtc = log(tc);

  // tc = 0 => no penalty
  // tc > 0 => infinite penalty for points of the cache
  // Small value of tc (close to 0) => penalty is null nearly everywhere
  // Large value of tc (close to 1) => penalty is non null nearly everywhere
  
  Matrix P ("P",pxx,1);
  for (int i=0 ; i<pxx ; i++){
    r12 = get_d1_over_d2( XXs.get_row(i) );
    if ( r12<tc )
      //p = std::max(0.0,-1+log(r12)/logtc);
      p = 1e+9 - r12;
    else
      p = 0.0;
    P.set(i,0,p);
  }
  return P;
}//

/*--------------------------------------*/
/*       get_distance_to_closest        */
/*--------------------------------------*/
Matrix SGTELIB::TrainingSet::get_distance_to_closest ( const Matrix & XXs ) const {
  #ifdef SGTELIB_DEBUG
    check_ready(); 
  #endif
  const int pxx = XXs.get_nb_rows();
  double d;
  Matrix P ("P",pxx,1);
  for (int i=0 ; i<pxx ; i++){
    d = get_d1 ( XXs.get_row(i) );
    P.set(i,0,d);
  }
  return P;
}//



/*--------------------------------------*/
/*       select points                  */
/*--------------------------------------*/
std::list<int> SGTELIB::TrainingSet::select_greedy ( const Matrix & X,
                                                     const int imin,
                                                     const int pS,
                                                     const double lambda0, 
                                                     const distance_t dt ){

  const int p = X.get_nb_rows();
  const int n = X.get_nb_cols();

  if ( pS<3 || pS>=p ){
    std::cout << "pS = " << pS << "\n";
    throw Exception ( __FILE__ , __LINE__ ,"TrainingSet::TrainingSet(): wrong value of pS" );
  }

  std::list<int> S;
  S.clear();
  
  int inew;
  Matrix xnew("xnew",1,n);
  Matrix x   ("x"   ,1,n);

  // Select the best point (here, the set B is only the point i_min)
  xnew = X.get_row(imin);
  // Distance vector between the set B and the cache
  Matrix dB = get_distances(X,xnew,dt);
  dB.set_name("dB");
  // Add to S
  S.push_back(imin);
  #ifdef SGTELIB_DEBUG
    std::cout << "First point : " << imin << "\n";
  #endif

  // Select the further point from B 
  // (nb : selecting one point randomly works as well)
  inew = dB.get_max_index();
  xnew = X.get_row(inew);
  // Distance vector between the set S and the cache
  Matrix dS = get_distances(X,xnew,dt);
  dS.set_name("dS");
  // Add to S
  S.push_back(inew);
  #ifdef SGTELIB_DEBUG
    std::cout << "Second point : " << inew << "\n";
  #endif
  
  // As B is in S, we can take the min of both distances
  dS = Matrix::min(dS,dB);

  // Compute lambda init :
  #ifdef SGTELIB_DEBUG
    std::cout << "Compute lambda init\n";
  #endif
  double lambda = 0;
  if (lambda0!=0){
    for (int i=0 ; i<p ; i++){
      if (dB.get(i)>0){
        lambda = std::max ( lambda , dS.get(i)/dB.get(i) );
      }
    }
    lambda *= lambda0;
  }


  // Iterative selection
  #ifdef SGTELIB_DEBUG
    std::cout << "Start greedy selection (S.size / pS = " << S.size() << " / " << pS << ")\n";
  #endif
  while ((int) S.size() < pS){
    #ifdef SGTELIB_DEBUG
      std::cout << "New iteration with lambda = " << lambda << "\n";
      (dS-lambda*dB).display(std::cout);
    #endif
    inew = (dS-lambda*dB).get_max_index();
    #ifdef SGTELIB_DEBUG
      std::cout << "inew : " << inew << "\n";
    #endif
    if (dS.get(inew)==0){
      #ifdef SGTELIB_DEBUG
        std::cout << "dS(inew) == 0 !\n";
      #endif
      // Update lambda
      lambda *= 0.99;
      if (lambda<1e-8) break;
    }
    else{
      #ifdef SGTELIB_DEBUG
        std::cout << "Add point " << inew << " to set\n";
      #endif
      // Add index in S
      S.push_back(inew);
      // Get coordinates of new points
      xnew = X.get_row(inew);
      // Update dS
      dS = Matrix::min( dS , get_distances(X,xnew,dt) );
      dS.set_name("dS");
    }
  }// End while

  return S;
}//


/*--------------------------------------*/
/*               display                */
/*--------------------------------------*/
void SGTELIB::TrainingSet::display ( std::ostream & out ) const {
  check_ready();
  int j;

  // dimensions:
  out << "Number of points, p=";
  out.width(4);
  out << _p << "  (" << _pvar << ")\n";
  out << "Input dimension,  n=";
  out.width(4);
  out << _n << "  (" << _nvar << ")\n";
  out << "Output dimension, m=";
  out.width(4);
  out << _m << "  (" << _mvar << ")\n";


  if (_ready){
    out << "X (Input matrix):\n";
    out << "___________________________________________________________________________________\n";
    out << "Dim|type|nbdiff|       mean        std|         lb         ub|         a          b|\n";
    out << "---|----|------|----------------------|----------------------|---------------------|\n";
    for ( j = 0 ; j < _n ; j++ ){
      out.width(3);
      out << j            <<"| ";
      out << " NA| ";
      out.width(5);
      out << _X_nbdiff[j] <<"| ";
      out.width(10);
      out << _X_mean[j]   <<" ";
      out.width(10);
      out << _X_std[j]    <<"| ";
      out.width(10); 
      out << _X_lb[j]     <<" ";
      out.width(10); 
      out << _X_ub[j]     <<"|";
      out.width(10); 
      out << _X_scaling_a[j]     <<" ";
      out.width(10); 
      out << _X_scaling_b[j]     <<"|\n";
    }
    out << "------------------------------------------------------------------------------------\n";

    out << "\n";

    out << "Z (Input matrix):\n";
    out << "___________________________________________________________________________________\n";
    out << "Dim|type|nbdiff|       mean        std|         lb         ub|         a          b|\n";
    out << "---|----|------|----------------------|----------------------|---------------------|\n";
    for ( j = 0 ; j < _m ; j++ ){
      out.width(3);
      out << j            <<"| ";
      out << bbo_type_to_str(_bbo[j]) << "| ";
      out.width(5);
      out << _Z_nbdiff[j] <<"| ";
      out.width(10);
      out << _Z_mean[j]   <<" ";
      out.width(10);
      out << _Z_std[j]    <<"| ";
      out.width(10); 
      out << _Z_lb[j]     <<" ";
      out.width(10); 
      out << _Z_ub[j]     <<"|";
      out.width(10); 
      out << _Z_scaling_a[j]     <<" ";
      out.width(10); 
      out << _Z_scaling_b[j]     <<"|\n";
    }
    out << "------------------------------------------------------------------------------------\n";
    std::cout << "fs_min: " << _fs_min << "\n";
    std::cout << "f_min:  " << _f_min << "\n";
  }


  out << std::endl;

}//

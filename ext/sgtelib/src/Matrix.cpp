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

#include "Matrix.hpp"

/*---------------------------*/
/*        constructor 1      */
/*---------------------------*/
SGTELIB::Matrix::Matrix ( const std::string & name ,
                          int                 nbRows    ,
                          int                 nbCols    ) : 
               _name   ( name ) ,
               _nbRows ( nbRows    ) ,
               _nbCols ( nbCols    )   {
#ifdef SGTELIB_DEBUG
  if ( _nbRows < 0 || _nbCols < 0 )
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::constructor 1: bad dimensions" );
#endif

  int i , j;

  _X = new double * [_nbRows];
  for ( i = 0 ; i < _nbRows ; ++i ) {
    _X[i] = new double [_nbCols];
    for ( j = 0 ; j < _nbCols ; ++j )
      _X[i][j] = 0.0;
  }
}//

/*---------------------------*/
/*        constructor 2      */
/*---------------------------*/
SGTELIB::Matrix::Matrix ( const std::string & name ,
                          int                 nbRows    ,
                          int                 nbCols    ,
                          double           ** A      ) :
               _name ( name ) ,
               _nbRows    ( nbRows    ) ,
               _nbCols    ( nbCols    )   {
#ifdef SGTELIB_DEBUG
  if ( _nbRows < 0 || _nbCols < 0 )
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , 
              "Matrix::constructor 2: bad dimensions" );
#endif

  int i , j;

  _X = new double * [_nbRows];
  for ( i = 0 ; i < _nbRows ; ++i ) {
    _X[i] = new double [_nbCols];
    for ( j = 0 ; j < _nbCols ; ++j )
      _X[i][j] = A[i][j];
  }
}//

/*---------------------------*/
/*        constructor 3      */
/*---------------------------*/
SGTELIB::Matrix::Matrix ( const std::string & file_name ) : 
                  _name ( "no_name" ) ,
                  _nbRows    ( 0         ) ,
                  _nbCols    ( 0         ) ,
                  _X    ( NULL      )   {
  *this = import_data(file_name);
}//


/*---------------------------*/
/*        constructor 4      */
/*---------------------------*/
SGTELIB::Matrix::Matrix (void) : 
               _name ( "" ) ,
               _nbRows    ( 0   ) ,
               _nbCols    ( 0   ) {
  _X = new double * [0];
}//

/*---------------------------*/
/*        constructor 5      */
/*---------------------------*/
SGTELIB::Matrix::Matrix (double v) : 
               _name ( "double" ) ,
               _nbRows    ( 1   ) ,
               _nbCols    ( 1   ) {
  #ifdef SGTELIB_DEBUG
    std::cout << "Matrix Constructor 5\n";
  #endif
  _X = new double * [1];
  _X[0] = new double [1];
  _X[0][0] = v;
}//


/*---------------------------*/
/*      copy constructor     */
/*---------------------------*/
SGTELIB::Matrix::Matrix ( const SGTELIB::Matrix & A ) : 
                          _name ( A._name ) ,
                          _nbRows    ( A._nbRows    ) ,
                          _nbCols    ( A._nbCols    ) {
  int i , j;
  _X = new double * [_nbRows];
  for ( i = 0 ; i < _nbRows ; ++i ) {
    _X[i] = new double [_nbCols];
    for ( j = 0 ; j < _nbCols ; ++j )
      _X[i][j] = A._X[i][j];
  }
}//




/*---------------------------*/
/*    affectation operator   */
/*---------------------------*/
SGTELIB::Matrix & SGTELIB::Matrix::operator = ( const SGTELIB::Matrix & A ) {
  
  if ( this == &A )
    return *this;

  int i , j;

  if ( _nbRows != A._nbRows || _nbCols != A._nbCols ) {

    for ( i = 0 ; i < _nbRows ; ++i )
      delete [] _X[i];
    delete [] _X;

    _nbRows = A._nbRows;
    _nbCols = A._nbCols;

    _X = new double * [_nbRows];
    for ( i = 0 ; i < _nbRows ; ++i ) {
      _X[i] = new double [_nbCols];
      for ( j = 0 ; j < _nbCols ; ++j )
        _X[i][j] = A._X[i][j];
    }
  }
  else {
    for ( i = 0 ; i < _nbRows ; ++i )
      for ( j = 0 ; j < _nbCols ; ++j )
        _X[i][j] = A._X[i][j];
  }
    
  _name = A._name;

  return *this;
}//



/*---------------------------*/
/*     import data           */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::import_data  ( const std::string & file_name ){

  std::ifstream in ( file_name.c_str() );

  if ( in.fail() ) {
    in.close();
    std::ostringstream oss;
    oss << "SGTELIB::Matrix::import_data: cannot open file " << file_name;
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , oss.str() );
  }

  std::string s;
  std::string line;
  while (std::getline(in, line)) s += line+";";

  return string_to_matrix(s);

}//


/*---------------------------------------*/
/*    affectation operator from string   */
/*---------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::string_to_matrix ( std::string s ) {

  // Replace tabs, newline etc.
  std::replace( s.begin(), s.end(), '\t', ' ');
  std::replace( s.begin(), s.end(), '\n', ';');
  std::replace( s.begin(), s.end(), '\r', ';');
  std::replace( s.begin(), s.end(), ',' , ' ');

  // Remove extra spaces
  s = SGTELIB::deblank(s);
  size_t i;
  std::string curline;

  // Find name
  std::string name = "MAT"; // default name
  i = std::min(s.find("="),s.find("["));
  if (i!=std::string::npos){
    curline = SGTELIB::deblank(s.substr(0,i));
    if (curline.size()){
      name = curline;
    }
    s = s.substr(i+1);
  }

  // Replace closing brakets by semi-colon.
  std::replace( s.begin(), s.end(), '=', ' ');
  std::replace( s.begin(), s.end(), '[', ' ');
  std::replace( s.begin(), s.end(), ']', ' ');

  // Read data
  int nbCols=-1;
  SGTELIB::Matrix M;
  while (true){

    i = s.find(";");
    if (i==std::string::npos) break;
    curline = SGTELIB::deblank(s.substr(0,i));
    s = s.substr(i+1);
    if (curline.size()){
      if (nbCols==-1){
        nbCols= SGTELIB::count_words(curline);
        M = SGTELIB::Matrix(name,0,nbCols);
      }
      M.add_rows(SGTELIB::Matrix::string_to_row(curline,nbCols));
    }

  }
  return M;

}//

/*------------------------------*/
/* convert string to row vector */
/*------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::string_to_row  ( const std::string & s , int nbCols ){
  if (nbCols<=0){
    nbCols = count_words(s);
  }
  SGTELIB::Matrix row("r",1,nbCols);
  double v;
  std::stringstream ss( s );
  int i=0;
  while( ss >> v ) row._X[0][i++] = v;
  if (i++!=nbCols){
    std::cout << "In line \"" << s << "\"\n";
    std::cout << "Found " << i << " components\n";
    std::cout << "Expected " << nbCols << " components\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
       "Matrix::string_to_row : cannot read line "+s );
  }
  return row;
}//

/*---------------------------*/
/*         operators         */
/*---------------------------*/
SGTELIB::Matrix operator * (const SGTELIB::Matrix & A , const SGTELIB::Matrix & B) {
  return SGTELIB::Matrix::product(A,B);
}
SGTELIB::Matrix operator * (const SGTELIB::Matrix & A , const double v) {
  int nbCols = A.get_nb_rows();
  int nbRows = A.get_nb_cols();
  SGTELIB::Matrix B(SGTELIB::dtos(v)+"*"+A.get_name(),nbCols,nbRows);
  int i,j;
  for ( i = 0 ; i < nbCols ; ++i ) {
    for ( j = 0 ; j < nbRows ; ++j ) {
      B.set(i,j,v*A.get(i,j));
    }
  }
  return B;
}
SGTELIB::Matrix operator * (const double v , const SGTELIB::Matrix & A) {
  return A*v;
}
SGTELIB::Matrix operator + (const SGTELIB::Matrix & A , const SGTELIB::Matrix & B) {
  return SGTELIB::Matrix::add(A,B);
}//
SGTELIB::Matrix operator + (const SGTELIB::Matrix & A , const double v) {
  const int nbCols = A.get_nb_rows();
  const int nbRows = A.get_nb_cols();
  SGTELIB::Matrix B(SGTELIB::dtos(v)+"+"+A.get_name(),nbCols,nbRows);
  int i,j;
  for ( i = 0 ; i < nbCols ; ++i ) {
    for ( j = 0 ; j < nbRows ; ++j ) {
      B.set(i,j,v+A.get(i,j));
    }
  }
  return B;
}//

SGTELIB::Matrix operator + (const double v , const SGTELIB::Matrix & A) {
  return A+v;
}//

SGTELIB::Matrix operator - (const SGTELIB::Matrix & A , const SGTELIB::Matrix & B) {
  return SGTELIB::Matrix::sub(A,B);
}//

SGTELIB::Matrix operator - (const SGTELIB::Matrix & A) {
  SGTELIB::Matrix M = A*(-1.0);
  M.set_name("(-"+A.get_name()+")");
  return M;
}//

SGTELIB::Matrix operator - (const double v , const SGTELIB::Matrix & A) {
  return v+(-A);
}//

SGTELIB::Matrix operator - (const SGTELIB::Matrix & A , const double v) {
  return A+(-v);
}//

SGTELIB::Matrix operator / (const SGTELIB::Matrix & A , const double v) {
  if (v==0){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
       "Matrix::operator /: divide by 0" );
  }
  return A * (1.0/v);
}//

/*---------------------------*/
/*         destructor        */
/*---------------------------*/
SGTELIB::Matrix::~Matrix ( void ) {
  /*
  #ifdef SGTELIB_DEBUG
    std::cout << "Delete " << _name << " of size (" << _nbCols << "," << _nbRows << ")...";
    display(std::cout);
    std::cout.flush();
  #endif
  */
  for ( int i = 0 ; i < _nbRows ; ++i )
    delete [] _X[i];
  delete [] _X; 
  /*
  #ifdef SGTELIB_DEBUG
    std::cout << "Deleted.\n";
  #endif
  */
}//

/*---------------------------*/
/*        add one row        */
/*---------------------------*/
void SGTELIB::Matrix::add_row  ( const double * row ) {

  double ** new_X = new double * [_nbRows+1];

  for ( int i = 0 ; i < _nbRows ; ++i )
    new_X[i] = _X[i];

  new_X[_nbRows] = new double [_nbCols];
  for ( int j = 0 ; j < _nbCols ; ++j )
    new_X[_nbRows][j] = row[j];
  
  delete [] _X;
  _X = new_X;
  ++_nbRows;
}//

/*-------------------------------*/
/*          add rows         */
/*---------------------------*/
void SGTELIB::Matrix::add_rows ( const Matrix & A ) {

  if ( A._nbCols != _nbCols )
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::add_rows(): bad dimensions" );

  int i , j;
  int new_nbRows = _nbRows + A._nbRows;

  double ** new_X = new double * [new_nbRows];

  for ( i = 0 ; i < _nbRows ; ++i )
    new_X[i] = _X[i];

  for ( i = _nbRows ; i < new_nbRows ; ++i ) {
    new_X[i] = new double [_nbCols];
    for ( j = 0 ; j < _nbCols ; ++j )
      new_X[i][j] = A._X[i-_nbRows][j];
  }

  delete [] _X;
  _X = new_X;
  _nbRows = new_nbRows;
}//

/*---------------------------*/
/*          add rows         */
/*---------------------------*/
void SGTELIB::Matrix::add_cols ( const Matrix & A ) {

  if ( A._nbRows != _nbRows )
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::add_cols(): bad dimensions" );

  int i , j;
  int new_nbCols = _nbCols + A._nbCols;
  double * x;

  for ( i = 0 ; i < _nbRows ; ++i ) {
    // Complete line
    x = new double [new_nbCols];
    // Original columns
    for ( j = 0 ; j < _nbCols ; ++j )
      x[j] = _X[i][j];
    // Additional columns
    for ( j = _nbCols ; j < new_nbCols ; ++j )
      x[j] = A._X[i][j-_nbCols];
    // Remove original line
    delete [] _X[i];
    // Put new line
    _X[i] = x;
  }

  _nbCols = new_nbCols;
}//

/*---------------------------------*/
/*          add empty rows         */
/*---------------------------------*/
void SGTELIB::Matrix::add_rows ( const int p ) {

  int i , j;
  int new_nbRows = _nbRows + p;

  double ** new_X = new double * [new_nbRows];

  for ( i = 0 ; i < _nbRows ; ++i )
    new_X[i] = _X[i];

  for ( i = _nbRows ; i < new_nbRows ; ++i ) {
    new_X[i] = new double [_nbCols];
    for ( j = 0 ; j < _nbCols ; ++j )
      new_X[i][j] = 0.0;
  }

  delete [] _X;

  _X = new_X;
  _nbRows = new_nbRows;
}//

/*---------------------------------*/
/*         remove last rows        */
/*---------------------------------*/
void SGTELIB::Matrix::remove_rows ( const int p ) {

  int i;
  int new_nbRows = _nbRows - p;

  double ** new_X = new double * [new_nbRows];

  for ( i = 0 ; i < new_nbRows ; ++i )
    new_X[i] = _X[i];

  for ( i = new_nbRows ; i < _nbRows ; ++i )
    delete [] _X[i];

  delete [] _X;
  _X = new_X;
  _nbRows = new_nbRows;
}//

/*---------------------------------*/
/*          add empty cols         */
/*---------------------------------*/
void SGTELIB::Matrix::add_cols ( const int p ) {

  int i , j;
  int new_nbCols = _nbCols + p;

  double * x;

  for ( i = 0 ; i < _nbRows ; ++i ) {
    x = new double [new_nbCols];
    for ( j = 0 ; j < _nbCols ; ++j )
      x[j] = _X[i][j];
    for ( j = _nbCols ; j < new_nbCols ; ++j )
      x[j] = 0.0;
    delete [] _X[i];
    _X[i] = x;
  }

  _nbCols = new_nbCols;
}//

/*-----------------------------------------*/
/*  get the fix columns (fixed variables)  */
/*-----------------------------------------*/
void SGTELIB::Matrix::get_fix_columns ( std::list<int> & fix_col ) const {
  fix_col.clear();
  for ( int j = 0 ; j < _nbCols ; ++j )
    if ( get_nb_diff_values(j) == 1 )
      fix_col.push_back(j);
}//


/*-----------------------------------------*/
/*  is the matrix symmetric                */
/*-----------------------------------------*/
bool SGTELIB::Matrix::is_sym ( void ) const {
  if (_nbCols!=_nbRows) return false;
  int i,j;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = i+1 ; j < _nbCols ; ++j ) {
      if (_X[i][j]!=_X[j][i]){
        return false;
      }
    }
  }
  return true;
}//


/*-------------------------------*/
/*  set the matrix randomly  */
/*---------------------------*/
void SGTELIB::Matrix::set_random ( double l , double u , bool round ) {
  int i , j;
  for ( i = 0 ; i < _nbRows ; ++i ){
    for ( j = 0 ; j < _nbCols ; ++j ){
      _X[i][j] = l + (u-l) * SGTELIB::uniform_rand();
      if ( round )
        _X[i][j] = SGTELIB::round ( _X[i][j] );
    }
  }
}//

/*---------------------------*/
/*  fill with value v        */
/*---------------------------*/
void SGTELIB::Matrix::fill ( double v ) {
  int i , j;
  for ( i = 0 ; i < _nbRows ; ++i ){
    for ( j = 0 ; j < _nbCols ; ++j ){
      _X[i][j] = v;
    }
  }
}//

/*---------------------------*/
/*     set element (i,j)     */
/*---------------------------*/
void SGTELIB::Matrix::set ( const int i , const int j , const double d ) {
  #ifdef SGTELIB_DEBUG
    if ( i < 0 || i >= _nbRows || j < 0 || j >= _nbCols ){
      display(std::cout);
      std::cout << "Error: try to set (" << i << "," << j << ") while dim is [" << _nbRows << "," << _nbCols << "]\n";
      std::cout.flush();
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::set(i,j): bad index" );
    }
  #endif
  _X[i][j] = d;
}//

void SGTELIB::Matrix::set_row (const SGTELIB::Matrix & T , const int i){
  #ifdef SGTELIB_DEBUG
    if ( i < 0 || i >= _nbRows || T.get_nb_rows()!=1 || T.get_nb_cols()!=_nbCols ){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::set_row: bad index" );
    }
  #endif
  for (int j=0 ; j<_nbCols ; j++){
    _X[i][j] = T.get(0,j);
  }
}//

void SGTELIB::Matrix::set_col (const SGTELIB::Matrix & T , const int j){
  #ifdef SGTELIB_DEBUG
    if ( j < 0 || j >= _nbCols || T.get_nb_rows()!=_nbRows || T.get_nb_cols()!=1 ){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::set_col: bad index" );
    }
  #endif
  for (int i=0 ; i<_nbRows ; i++){
    _X[i][j] = T.get(i,0);
  }
}//

void SGTELIB::Matrix::set_row (const double v , const int i){
  #ifdef SGTELIB_DEBUG
    if ( i < 0 || i >= _nbRows ){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::set_row: bad index" );
    }
  #endif
  for (int j=0 ; j<_nbCols ; j++){
    _X[i][j] = v;
  }
}//

void SGTELIB::Matrix::set_col (const double v , const int j){
  #ifdef SGTELIB_DEBUG
    if ( j < 0 || j >= _nbCols ){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::set_col: bad index" );
    }
  #endif
  for (int i=0 ; i<_nbRows ; i++){
    _X[i][j] = v;
  }
}//

void SGTELIB::Matrix::swap (const int i1 , const int j1 , const int i2 , const int j2 ){
  #ifdef SGTELIB_DEBUG
    if ( i1 < 0 || i1 >= _nbRows || j1 < 0 || j1 >= _nbCols || i2 < 0 || i2 >= _nbRows || j2 < 0 || j2 >= _nbCols ){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::permut: bad index" );
    }
  #endif
  double buffer = _X[i1][j1];
  _X[i1][j1] = _X[i2][j2];
  _X[i2][j2] = buffer;
}//

void SGTELIB::Matrix::multiply_row (const double v , const int i){
  #ifdef SGTELIB_DEBUG
    if ( i < 0 || i >= _nbRows ){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::set_row: bad index" );
    }
  #endif
  for (int j=0 ; j<_nbCols ; j++){
    _X[i][j] *= v;
  }
}//

void SGTELIB::Matrix::multiply_col (const double v , const int j){
  #ifdef SGTELIB_DEBUG
    if ( j < 0 || j >= _nbCols ){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::set_col: bad index" );
    }
  #endif
  for (int i=0 ; i<_nbRows ; i++){
    _X[i][j] *= v;
  }
}//

/*---------------------------*/
/*  access to element (k)    */
/*---------------------------*/
double SGTELIB::Matrix::get ( const int k ) const {
  return (*this)[k];
}//


const double & SGTELIB::Matrix::operator [] ( int k ) const {
  int i = 0 , j = 0;
  if (_nbRows==1) j=k; 
  else if (_nbCols==1) i=k;
  else throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::[k]: the matrix is not a vector" );
  return _X[i][j];
}//
double & SGTELIB::Matrix::operator [] ( int k ){
  int i = 0 , j = 0;
  if (_nbRows==1) j=k; 
  else if (_nbCols==1) i=k;
  else throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::[k]: the matrix is not a vector" );
  return _X[i][j];
}//

/*------------------------------------------*/
/*  get (access to a subpart of the matrix  */
/*------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get (const std::list<int> & list_cols , const std::list<int> & list_rows) const {
  return get_rows(list_rows).get_cols(list_cols);
}//

/*---------------------------*/
/*      get row              */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_row (const int i) const {
  #ifdef SGTELIB_DEBUG
    if ( i < 0 || i >= _nbRows )
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::get_row(i): bad index" );
  #endif
  SGTELIB::Matrix A (_name+"(i,:)",1,_nbCols);
  for (int j=0 ; j<_nbCols ; j++){
    A._X[0][j] = _X[i][j];
  }
  return A;
}//

/*---------------------------*/
/*      get col              */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_col (const int j) const {
  #ifdef SGTELIB_DEBUG
    if ( j < 0 || j >= _nbCols )
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::get_row(i): bad index" );
  #endif
  SGTELIB::Matrix A (_name+"(:,j)",_nbRows,1);
  for (int i=0 ; i<_nbRows ; i++){
    A._X[i][0] = _X[i][j];
  }
  return A;
}//

/*---------------------------*/
/*      get_rows             */
/* get rows from i1 to i2-1  */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_rows  (const int i1, const int i2) const {

  if ( (i1<0) || (i1>_nbRows) || (i2<0) || (i2>_nbRows) || (i1>=i2) ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::get_rows: bad index" );
  }

  // Otherwise, select the rows
  const int nbRows = i2-i1;
  const int nbCols = _nbCols;
  SGTELIB::Matrix A (_name+"(i1:i2-1,:)",nbRows,nbCols);

  int i,k=0;
  for ( i=i1 ; i<i2 ; i++ )
    A.set_row(get_row(i),k++);
  return A;
}//

/*---------------------------*/
/*      get_cols             */
/* get rows from i1 to i2-1  */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_cols  (const int i1, const int i2) const {

  if ( (i1<0) || (i1>_nbCols) || (i2<0) || (i2>_nbCols) || (i1>=i2) ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::get_cols: bad index" );
  }

  // Otherwise, select the rows
  const int nbCols = i2-i1;
  const int nbRows = _nbRows;
  SGTELIB::Matrix A (_name+"(:,i1:i2-1)",nbRows,nbCols);

  int i,k=0;
  for ( i=i1 ; i<i2 ; i++ )
    A.set_col(get_col(i),k++);
  return A;
}//

/*---------------------------*/
/*      get_rows             */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_rows (const std::list<int> & list_rows) const {

  // If the list has only one element and this element is -1, then all rows are returned.
  if ( (list_rows.size()==1) && (list_rows.front()==-1) ){
    return *this;
  }

  // Otherwise, select the rows
  const int nbRows = static_cast<int>(list_rows.size());
  const int nbCols = _nbCols;
  SGTELIB::Matrix A (_name+"_get_rows",nbRows,nbCols);
  
  std::list<int>::const_iterator it;

  int k=0;
  for ( it = list_rows.begin() ; it != list_rows.end() ; ++it ) {
    if ( *it < 0 || *it >= _nbRows ) {
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::get_rows: bad index" );
    }
    A.set_row(get_row(*it),k++);
  }
  return A;
}//

/*---------------------------*/
/*      get_cols             */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_cols (const std::list<int> & list_cols) const {

  // If the list has only one element and this element is -1, then all rows are returned.
  if ( (list_cols.size()==1) && (list_cols.front()==-1) ){
    return *this;
  }

  // Otherwise, select the rows
  const int nbRows = _nbRows;
  const int nbCols = static_cast<int>(list_cols.size());
  SGTELIB::Matrix A (_name+"_get_cols",nbRows,nbCols);
  
  std::list<int>::const_iterator it;

  int k=0;
  for ( it = list_cols.begin() ; it != list_cols.end() ; ++it ) {
    if ( *it < 0 || *it >= _nbCols ) {
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::get_rows: bad index" );
    }
    A.set_col(get_col(*it),k++);
  }

  return A;
}//

/*----------------------------------------------------*/
/*  count the number of different values in column j  */
/*----------------------------------------------------*/
int SGTELIB::Matrix::get_nb_diff_values ( int j ) const {
  std::set<double> s;
  for ( int i = 0 ; i < _nbRows ; ++i ){
    s.insert ( _X[i][j] );
  }
  return static_cast<int> ( s.size() );
}//

/*---------------------------*/
/*          display          */
/*---------------------------*/
void SGTELIB::Matrix::display ( std::ostream & out ) const {
  int i , j;
  out << std::endl << _name << "=[\n";
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j )
      out << "\t" << std::setw(10) << _X[i][j] << " ";
    out << ";" << std::endl;
  }
  out << "];" << std::endl;
}//

/*---------------------------*/
/*          display          */
/*---------------------------*/
void SGTELIB::Matrix::display_short ( std::ostream & out ) const {
  if (get_numel()<5) display(out);
  else{
    out << std::endl << _name << " ( " << _nbRows << " x " << _nbCols << " ) =\n[";
    out << "\t" << std::setw(10) << _X[0][0] << " ";
    if (_nbCols>2) out << "... ";
    out << "\t" << std::setw(10) << _X[0][_nbCols] << "\n";
    if (_nbRows>2) out << "\t       ...";
    if (_nbCols>2) out << "    ";
    if (_nbRows>2) out << "\t       ...\n";
    out << "\t" << std::setw(10) << _X[_nbRows-1][0] << " ";
    if (_nbCols>2) out << "... ";
    out << "\t" << std::setw(10) << _X[_nbRows-1][_nbCols] << "]\n";
  }
}//

/*---------------------------*/
/*          write            */
/*---------------------------*/
void SGTELIB::Matrix::write ( const std::string & file_name ) const {
  std::ofstream output_file;
  output_file.open (file_name.c_str());
  display(output_file);
  output_file.close();  
}//

/*-------------------------------*/
/*          display_size         */
/*-------------------------------*/
void SGTELIB::Matrix::display_size ( std::ostream & out ) const {
  out << "Matrix " << _name << " : " << _nbRows << " , " << _nbCols << "\n";
}//

/*---------------------------*/
/* double * to column vector */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::col_vector ( const double * v,
                                              const int n     )  {
  if ( ! v){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::column_vector: v is null" );
  }  
  SGTELIB::Matrix V("V",n,1);
  for (int i=0 ; i<n ; i++){
    V._X[i][0] = v[i];
  }
  return V;
}//

/*---------------------------*/
/* double * to row vector    */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::row_vector ( const double * v,
                                              const int n     )  {
  if ( ! v){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::column_vector: v is null" );
  }  
  SGTELIB::Matrix V("V",1,n);
  for (int i=0 ; i<n ; i++){
    V._X[0][i] = v[i];
  }
  return V;
}//

/*---------------------------*/
/*          product          */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::product ( const SGTELIB::Matrix & A,
                                           const SGTELIB::Matrix & B )  {

  #ifdef SGTELIB_DEBUG
    if (A.get_nb_cols()!=B.get_nb_rows()){
      std::cout << "A (" << A.get_name() << ") : " << A.get_nb_rows() << " , " << A.get_nb_cols() << "\n";
      std::cout << "B (" << B.get_name() << ") : " << B.get_nb_rows() << " , " << B.get_nb_cols() << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::product(A,B): dimension error" );
    }  
    // Init matrix
    SGTELIB::Matrix C(A.get_name()+"*"+B.get_name(),A.get_nb_rows(),B.get_nb_cols());
  #else
    SGTELIB::Matrix C("C",A.get_nb_rows(),B.get_nb_cols());
  #endif

  // Compute
  int i,j,k;
  const int nb_rows = C.get_nb_rows();
  const int nb_cols = C.get_nb_cols();
  const int nb_inter= A.get_nb_cols();
  //double v;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C._X[i][j] = 0;
    }
    for ( k = 0 ; k < nb_inter; ++k ){
      for ( j = 0 ; j < nb_cols ; ++j ){
        C._X[i][j] += A._X[i][k]*B._X[k][j];
      }
    }
  }
  return C;
}//

/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::product ( const SGTELIB::Matrix & A,
                                           const SGTELIB::Matrix & B, 
                                           const SGTELIB::Matrix & C){
  return product(A,product(B,C));
}//
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::product ( const SGTELIB::Matrix & A,
                                           const SGTELIB::Matrix & B,
                                           const SGTELIB::Matrix & C, 
                                           const SGTELIB::Matrix & D){
  return product(product(A,B),product(C,D));
}//
/*---------------------------*/

/*---------------------------------------------------*/
/* Subset product                                    */
/* multiply                                          */
/* the p first rows and q first columns of A         */
/* with the q first rows and r first columns of B.   */
/* Result is a matrix of size p/r.                   */
/*---------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::subset_product (const SGTELIB::Matrix & A,
                                                 const SGTELIB::Matrix & B,
                                                 int p,
                                                 int q,
                                                 int r){

  // Default p value
  if (p==-1){
    p=A.get_nb_rows();
  }
  // Otherwise, need to check the number of rows.
  else if (A.get_nb_rows()<p){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::subset_product: dimension error" );
  }  

  // Default q value
  if ( (q==-1) & (A.get_nb_cols()==B.get_nb_rows()) ){
    q = A.get_nb_cols();
  }
  else{
    // Check for q
    if (A.get_nb_cols()<q){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::subset_product: dimension error" );
    }  
    if (B.get_nb_rows()<q){
      throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::subset_product: dimension error" );
    }  
  }

  // Default r value
  if (r==-1){
    r = B.get_nb_cols();
  }
  else if (B.get_nb_cols()<r){
    // Check for r
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::subset_product: dimension error" );
  }  

  SGTELIB::Matrix C("A*B",p,r);
  for (int i=0 ; i<p ; i++){
    for (int j=0 ; j<r ; j++){
      for (int k=0 ; k<q ; k++){
        C._X[i][j] += A._X[i][k]*B._X[k][j];
      } 
    }
  }
  return C;

}//

/*---------------------------*/
/* Hadamard product          */
/* (Term to term product)    */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::hadamard_product ( const SGTELIB::Matrix & A,
                                                    const SGTELIB::Matrix & B )  {
  const int nb_rows = A.get_nb_rows();
  const int nb_cols = A.get_nb_cols();

  if (B.get_nb_rows()!=nb_rows){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::hadamard_product(A,B): dimension error" );
  }  
  if (B.get_nb_cols()!=nb_cols){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::hadamard_product(A,B): dimension error" );
  }  

  // Init matrix
  SGTELIB::Matrix C(A.get_name()+".*"+B.get_name(),nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C.set(i,j,A.get(i,j)*B.get(i,j));
    }
  }
  return C;
}//

/*---------------------------*/
/* Hadamard square           */
/* (Term to term square)     */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::hadamard_square ( const SGTELIB::Matrix & A )  {
  const int nb_rows = A.get_nb_rows();
  const int nb_cols = A.get_nb_cols();

  // Init matrix
  SGTELIB::Matrix C("("+A.get_name()+").^2",nb_rows,nb_cols);

  // Compute
  int i,j;
  double a;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      a = A._X[i][j];
      C._X[i][j] = a*a;
    }
  }
  return C;
}//

/*---------------------------*/
/* Hadamard square           */
/* (Term to term square)     */
/*---------------------------*/
void SGTELIB::Matrix::hadamard_square ( void )  {
  // change name
  _name = "("+_name+").^2";
  // Compute
  int i,j;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      _X[i][j] *= _X[i][j];
    }
  }
}//

/*---------------------------*/
/* Hadamard inverse          */
/* (Term to term inverse)    */
/*---------------------------*/
void SGTELIB::Matrix::hadamard_inverse ( void )  {
  // change name
  _name = "("+_name+").^-1";
  // Compute
  int i,j;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      _X[i][j] = 1/_X[i][j];
    }
  }
}//

/*---------------------------*/
/* Hadamard sqrt             */
/* (Term to term sqrt)       */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::hadamard_sqrt ( const SGTELIB::Matrix & A )  {
  const int nb_rows = A.get_nb_rows();
  const int nb_cols = A.get_nb_cols();

  // Init matrix
  SGTELIB::Matrix C("sqrt("+A.get_name()+")",nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C._X[i][j] = sqrt(fabs(A._X[i][j]));
    }
  }
  return C;
}//

/*---------------------------*/
/* Hadamard power            */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::hadamard_power ( const SGTELIB::Matrix & A , const double e )  {

  if (e==1.0) return A;

  const int nb_rows = A.get_nb_rows();
  const int nb_cols = A.get_nb_cols();

  // Init matrix
  SGTELIB::Matrix C("pow("+A.get_name()+","+dtos(e)+")",nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C._X[i][j] = pow(A._X[i][j],e);
    }
  }
  return C;
}//

/*---------------------------*/
/* Hadamard sqrt             */
/* (Term to term sqrt)       */
/*---------------------------*/
void SGTELIB::Matrix::hadamard_sqrt ( void )  {
  // change name
  _name = "sqrt("+_name+")";
  // Compute
  int i,j;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      _X[i][j] = sqrt(fabs(_X[i][j]));
    }
  }
}//

/*---------------------------------------*/
/*             diagA_product             */
/*        Product A*B, considering       */
/*    that A is diag (no verif on this)  */
/*---------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::diagA_product ( const SGTELIB::Matrix & A,
                                                 const SGTELIB::Matrix & B )  {

  // Init matrix
  const int na = A.get_nb_rows();
  const int ma = A.get_nb_cols();
  const int n = B.get_nb_rows();
  const int m = B.get_nb_cols();
  SGTELIB::Matrix C(A.get_name()+"*"+B.get_name(),n,m);

  int i,j;
  double Aii;

  if ( (na==ma) || (ma==n) ){
    // A is square, use the diag terms
    for ( i = 0 ; i < n ; i++ ) {
      Aii = A._X[i][i];
      for ( j = 0 ; j < m ; j++ ) {
        C._X[i][j] = Aii*B._X[i][j];
      }
    }
    return C;
  }
  else if ( (na==1) && (ma==n) ){
    // A is a line vector
    for ( i = 0 ; i < n ; i++ ) {
      Aii = A._X[0][i];
      for ( j = 0 ; j < m ; j++ ) {
        C._X[i][j] = Aii*B._X[i][j];
      }
    }
    return C;
  }
  else if ( (na==n) && (ma==1) ){
    // A is a col vector
    for ( i = 0 ; i < n ; i++ ) {
      Aii = A._X[i][0];
      for ( j = 0 ; j < m ; j++ ) {
        C._X[i][j] = Aii*B._X[i][j];
      }
    }
    return C;
  }
  else {
    std::cout << "A (" << A.get_name() << ") : " << A.get_nb_rows() << " , " << A.get_nb_cols() << "\n";
    std::cout << "B (" << B.get_name() << ") : " << B.get_nb_rows() << " , " << B.get_nb_cols() << "\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::diagA_product(A,B): dimension error" );
  }  

}//

/*---------------------------------------*/
/*             diagB_product             */
/*        Product A*B, considering       */
/*    that B is diag (no verif on this)  */
/*---------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::diagB_product ( const SGTELIB::Matrix & A,
                                                 const SGTELIB::Matrix & B )  {

  // Init matrix
  const int n = A.get_nb_rows();
  const int m = A.get_nb_cols();
  const int nb = B.get_nb_rows();
  const int mb = B.get_nb_cols();
  SGTELIB::Matrix C(A.get_name()+"*"+B.get_name(),n,m);

  int i,j;
  double Bjj;

  if ( (nb==mb) && (mb==n) ){
    // B is square, use the diag terms
    for ( j = 0 ; j < m ; j++ ) {
      Bjj = B._X[j][j];
      for ( i = 0 ; i < n ; i++ ) {
        C._X[i][j] = A._X[i][j]*Bjj;
      }
    }
    return C;
  }
  else if ( (nb==1) && (mb==m) ){
    // B is a line vector
    for ( j = 0 ; j < m ; j++ ) {
      Bjj = B._X[0][j];
      for ( i = 0 ; i < n ; i++ ) {
        C._X[i][j] = A._X[i][j]*Bjj;
      }
    }
    return C;
  }
  else if ( (nb==m) && (mb==1) ){
    // B is a col vector
    for ( j = 0 ; j < m ; j++ ) {
      Bjj = B._X[j][0];
      for ( i = 0 ; i < n ; i++ ) {
        C._X[i][j] = A._X[i][j]*Bjj;
      }
    }
    return C;
  }
  else {
    std::cout << "A (" << A.get_name() << ") : " << A.get_nb_rows() << " , " << A.get_nb_cols() << "\n";
    std::cout << "B (" << B.get_name() << ") : " << B.get_nb_rows() << " , " << B.get_nb_cols() << "\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::diagB_product(A,B): dimension error" );
  }  
}//


/*-------------------------------------*/
/*      transposeA_product             */
/*        Product A'*B                 */
/*-------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::transposeA_product ( const SGTELIB::Matrix & A,
                                                      const SGTELIB::Matrix & B )  {

  if (A.get_nb_rows()!=B.get_nb_rows()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::transposeA_product(A,B): dimension error" );
  }  

  // Init matrix
  SGTELIB::Matrix C(A.get_name()+"'*"+B.get_name(),A.get_nb_cols(),B.get_nb_cols());

  // Compute
  int i,j,k;
  const int nb_rows = C.get_nb_rows();
  const int nb_cols = C.get_nb_cols();
  const int nb_inter= A.get_nb_rows();
  //double v;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C._X[i][j] = 0;
      for ( k = 0 ; k < nb_inter; ++k ){
        C._X[i][j] += A._X[k][i]*B._X[k][j];
      }
    }
  }
  return C;

}//



/*---------------------------*/
/*        get matrix P       */
/* P = I - H*Ai*H'           */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_matrix_P ( const SGTELIB::Matrix & Ai,
                                                const SGTELIB::Matrix & H ){
  const int p = H.get_nb_rows();
  std::cout << "Function get_matrix_P should be avoided !!\n";
  return identity(p) - (H * Ai * H.transpose());
}//

/*---------------------------*/
/*  compute P*Zs             */
/* where P = I - H*Ai*H'     */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_matrix_PZs ( const SGTELIB::Matrix & Ai,
                                                  const SGTELIB::Matrix & H ,
                                                  const SGTELIB::Matrix & Zs){
  //return Zs - H * Ai * H' * Zs;
  return Zs - (H * Ai) * transposeA_product(H,Zs);
}//

/*---------------------------*/
/*      get matrix dPi       */
/*      dPi = diag(P)^-1     */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_matrix_dPi ( const SGTELIB::Matrix & Ai,
                                                  const SGTELIB::Matrix & H ){
  const int p = H.get_nb_rows();
  SGTELIB::Matrix dPi ("dPi",p,p);
  SGTELIB::Matrix h;
  double v;

  for ( int i = 0 ; i < p ; ++i ) {
    h = H.get_row(i);
    v = (h*Ai*h.transpose()).get(0,0);    
    v = 1.0/( 1.0 - v );
    dPi.set(i,i,v);
  }
  return dPi;
}//

/*---------------------------*/
/*  compute dPi*P*Zs         */
/* where P = I - H*Ai*H'     */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_matrix_dPiPZs ( const SGTELIB::Matrix & Ai,
                                                     const SGTELIB::Matrix & H ,
                                                     const SGTELIB::Matrix & Zs){

  const SGTELIB::Matrix HAi = H*Ai;
  //SGTELIB::Matrix dPiPZs = Zs - HAi * (H.transpose()*Zs);
  SGTELIB::Matrix dPiPZs = Zs - HAi * transposeA_product(H,Zs);

  // Take dPi into account
  const int p = H.get_nb_rows();
  const int q = H.get_nb_cols();
  double v;
  int i,j;
  for ( i = 0 ; i < p ; i++ ) {
    v = 0;
    for ( j = 0 ; j < q ; j++) v += HAi._X[i][j]*H._X[i][j];
    v = 1.0/( 1.0 - v );
    dPiPZs.multiply_row ( v , i );
  }

  return dPiPZs;
}//

/*---------------------------*/
/*  compute dPi*P*Zs         */
/* where P = I - H*Ai*H'     */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_matrix_dPiPZs ( const SGTELIB::Matrix & Ai,
                                                     const SGTELIB::Matrix & H ,
                                                     const SGTELIB::Matrix & Zs,
                                                     const SGTELIB::Matrix & ALPHA){

  const SGTELIB::Matrix HAi = H*Ai;
  //SGTELIB::Matrix dPiPZs = Zs - HAi * (H.transpose()*Zs);
  SGTELIB::Matrix dPiPZs = Zs - H*ALPHA;

  // Take dPi into account
  const int p = H.get_nb_rows();
  const int q = H.get_nb_cols();
  double v;
  int i,j;
  for ( i = 0 ; i < p ; i++ ) {
    v = 0;
    for ( j = 0 ; j < q ; j++) v += HAi._X[i][j]*H._X[i][j];
    v = 1.0/( 1.0 - v );
    dPiPZs.multiply_row ( v , i );
  }

  return dPiPZs;
}//




/*---------------------------*/
/*      get trace P          */
/*---------------------------*/
double SGTELIB::Matrix::get_trace_P ( const SGTELIB::Matrix & Ai,
                                      const SGTELIB::Matrix & H ){
  const int p = H.get_nb_rows();
  SGTELIB::Matrix h;
  double v;
  double trace = 0;

  for ( int i = 0 ; i < p ; ++i ) {
    h = H.get_row(i);
    v = (h*Ai*h.transpose()).get(0,0);    
    trace += v;
  }
  return trace;
}//

/*---------------------------*/
/*          addition         */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::add ( const SGTELIB::Matrix & A,
                                       const SGTELIB::Matrix & B )  {

  if (A.get_nb_cols()!=B.get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::add(A,B): dimension error" );
  }  
  if (A.get_nb_rows()!=B.get_nb_rows()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::add(A,B): dimension error" );
  }  

  const int nb_rows = A.get_nb_rows();
  const int nb_cols = A.get_nb_cols();

  // Init matrix
  SGTELIB::Matrix C(A.get_name()+"+"+B.get_name(),nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C._X[i][j] = A._X[i][j]+B._X[i][j];
    }
  }
  return C;
}//

/*---------------------------*/
/*          addition         */
/*---------------------------*/
void SGTELIB::Matrix::add ( const SGTELIB::Matrix & B ) {

  if ( _nbCols != B.get_nb_cols() ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Matrix::add(B): dimension error" );
  }  
  if ( _nbRows != B.get_nb_rows() ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Matrix::add(B): dimension error" );
  }  

  // Compute
  int i,j;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      _X[i][j] += B.get(i,j);
    }
  }
}//

/*---------------------------*/
/*     add and fill with 0   */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::add_fill ( const SGTELIB::Matrix & A,
                                            const SGTELIB::Matrix & B )  {

  const int nb_rows = std::max(A.get_nb_rows(),B.get_nb_rows());
  const int nb_cols = std::max(A.get_nb_cols(),B.get_nb_cols());

  // Init matrix
  SGTELIB::Matrix C(A.get_name()+"+"+B.get_name(),nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < A.get_nb_rows() ; ++i ) {
    for ( j = 0 ; j < A.get_nb_cols() ; ++j ){
      C._X[i][j] = A._X[i][j];
    }
  }
  for ( i = 0 ; i < B.get_nb_rows() ; ++i ) {
    for ( j = 0 ; j < B.get_nb_cols() ; ++j ){
      C._X[i][j] += B._X[i][j];
    }
  }
  return C;
}//

/*---------------------------*/
/*          substraction     */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::sub ( const SGTELIB::Matrix & A,
                                       const SGTELIB::Matrix & B )  {

  if (A.get_nb_cols()!=B.get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::sub(A,B): dimension error" );
  }  
  if (A.get_nb_rows()!=B.get_nb_rows()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::sub(A,B): dimension error" );
  }  

  int nb_rows = A.get_nb_rows();
  int nb_cols = A.get_nb_cols();

  // Init matrix
  SGTELIB::Matrix C(A.get_name()+"-"+B.get_name(),nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C.set(i,j,A.get(i,j)-B.get(i,j));
    }
  }
  return C;
}//

/*---------------------------*/
/*          substraction     */
/*---------------------------*/
void SGTELIB::Matrix::sub ( const SGTELIB::Matrix & B ) {

  if ( _nbCols != B.get_nb_cols() ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Matrix::sub(B): dimension error" );
  }  
  if ( _nbRows != B.get_nb_rows() ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Matrix::sub(B): dimension error" );
  }  

  // Compute
  int i,j;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      _X[i][j] -= B._X[i][j];
    }
  }
}//

/*---------------------------*/
/*          identity         */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::identity ( const int n )  {

  // Init matrix
  SGTELIB::Matrix I("I",n,n);
  // Fill with 0.0
  I.fill(0.0);
  // Put one on the diag
  for ( int i = 0 ; i < n ; ++i ) {
    I.set(i,i,1.0);
  }
  return I;
}//

/*---------------------------*/
/*          ones             */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::ones ( const int nbRows , const int nbCols ) {
  // Init matrix
  SGTELIB::Matrix matrixOnes("Ones",nbRows,nbCols);
  // Fill with 1.0
  matrixOnes.fill(1.0);
  return matrixOnes;
}//

/*---------------------------*/
/* rank                      */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::rank ( void ) const {

  // TODO: use faster method...
  if ((_nbRows>1) && (_nbCols>1)) 
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Matrix::rank: dimension error" );

  SGTELIB::Matrix R;
  if (_nbRows>1){
    R = this->transpose().rank();
    R = R.transpose();
  }
  else{
    const int m = _nbCols;
    SGTELIB::Matrix D = *this;
    R = SGTELIB::Matrix("R",1,m);
    double dmin;
    int i,j,jmin=0;
    for (i=0 ; i<m ; i++){
      dmin = +INF;
      for (j=0 ; j<m ; j++){
        if (D._X[0][j]<dmin){
          jmin = j;
          dmin = D._X[0][j];
        }
      }
      R.set(0,jmin,double(i));
      D.set(0,jmin,INF);
    }
  }
  return R;
}//

/*---------------------------*/
/*        Trace              */
/*---------------------------*/
double SGTELIB::Matrix::trace ( void ) const{
  int min_nm = std::min(_nbCols,_nbRows);
  double v = 0;
  for (int i=0 ; i<min_nm ; i++){
    // get diagonal term
    v += get(i,i);
  }
  return v;
}//

/*---------------------------*/
/*        Rmse               */
/*---------------------------*/
double SGTELIB::Matrix::rmse ( void ) const{
  double v = 0;

  // Compute
  int i,j;
  double xij;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      xij = _X[i][j];
      v += xij*xij;
    }
  }
  v /= _nbRows * _nbCols;
  v = sqrt(v);
  return v;
}//

/*---------------------------*/
/*        norm               */
/*---------------------------*/
double SGTELIB::Matrix::normsquare ( void ) const{
  double v = 0;
  // Compute
  int i,j;
  double xij;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      xij = _X[i][j];
      v += xij*xij;
    }
  }
  return v;
}//

/*----------------------------------------*/
/*        normalize_cols                  */
/* Normalizes each column so that the sum */
/* of the terms on this column is 1       */
/*----------------------------------------*/
void SGTELIB::Matrix::normalize_cols ( void ){
  int i,j;
  double d;
  for ( j = 0 ; j < _nbCols ; ++j ){
    d = 0;

    for ( i = 0 ; i < _nbRows ; ++i )
      d += _X[i][j];

    if (d==0){
      for ( i = 0 ; i < _nbRows ; ++i )
        _X[i][j] = 1/_nbRows;
    }
    else{
      for ( i = 0 ; i < _nbRows ; ++i )
        _X[i][j] /= d;
    }
  }
}//



/*------------------------------------------*/
/*  compute the norm of each column         */
/*------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::col_norm ( const norm_t nt ) const {
  
  SGTELIB::Matrix N = SGTELIB::Matrix("Norm",1,_nbCols);

  int i,j;
  for ( j = 0 ; j < _nbCols ; ++j ){
    double v = 0;
    switch (nt){
      case SGTELIB::NORM_0:
        for (i=0;i<_nbRows;++i) v += double(fabs(_X[i][j])<EPSILON);
        v /= double(_nbCols);
        break;
      case SGTELIB::NORM_1:
        for (i=0;i<_nbRows;++i) v += fabs(_X[i][j]);
        v /= double(_nbCols);
        break;
      case SGTELIB::NORM_2:
        for (i=0;i<_nbRows;++i) v += _X[i][j]*_X[i][j];
        v = sqrt(v/double(_nbCols));
        break;
      case SGTELIB::NORM_INF:
        for (i=0;i<_nbRows;++i) v = std::max(v,fabs(_X[i][j]));
        break;
    }
    N.set(0,j,v);
  }
  return N;
}//



/*---------------------------*/
/*        square norm        */
/*---------------------------*/
double SGTELIB::Matrix::norm ( void ) const{
  return sqrt(normsquare());
}//

/*---------------------------*/
/*        sum                */
/*---------------------------*/
double SGTELIB::Matrix::sum ( void ) const{
  double v = 0;

  // Compute
  int i,j;
  for ( i = 0 ; i < _nbRows ; ++i ) {
    for ( j = 0 ; j < _nbCols ; ++j ){
      v += _X[i][j];
    }
  }
  return v;
}//

/*-------------------------------------------------*/
/*        sum                                      */
/* Sum the element along one of the two directions */
/* direction == 1 => return the sum for each row   */
/* direction == 2 => return the sum for each col   */
/*-------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::sum ( const int direction ) const{

  double v;
  int i,j;
  
  if (direction == 1){
    SGTELIB::Matrix S ("S",1,_nbCols);
    for ( j = 0 ; j < _nbCols ; ++j ){
      v = 0;
      for ( i = 0 ; i < _nbRows ; ++i ) {
        v += _X[i][j];
      }
      S._X[0][j] = v;
    }
    return S;
  }
  else if (direction == 2){
    SGTELIB::Matrix S ("S",_nbRows,1);
    for ( i = 0 ; i < _nbRows ; ++i ) {
      v = 0;
      for ( j = 0 ; j < _nbCols ; ++j ){
        v += _X[i][j];
      }
      S._X[i][0] = v;
    }
    return S;
  }
  else{
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::sum(direction): direction must be 1 or 2" );
  }
}//

/*---------------------------*/
/*        mean               */
/*---------------------------*/
double SGTELIB::Matrix::mean ( void ) const{
  return sum()/(_nbRows*_nbCols);
}//

/*---------------------------*/
/*        count              */
/* Number of non null values */
/*---------------------------*/
int SGTELIB::Matrix::count ( void ) const{
  int v = 0;

  // Compute
  int i,j;
  const int nb_rows = get_nb_rows();
  const int nb_cols = get_nb_cols();
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      v += (fabs(_X[i][j])>EPSILON)? 1:0 ;
    }
  }
  return v;
}//

/*---------------------------*/
/*        Diag inverse       */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::diag_inverse ( void ) const{
  // Return a new matrix such that:
  // Non diagonal terms are null
  // Diagonal terms are the inverse of the original ones.
  // Works for square and non-square matrices.
  //
  // Does not work if a diag term is 0.

  // New matrix
  SGTELIB::Matrix DI("diag("+_name+")^-1",_nbCols,_nbRows);
  // nb: this constructor initializes the matrix to 0.0

  const int min_nm = std::min(_nbCols,_nbRows);

  double v = 0;
  for (int i=0 ; i<min_nm ; i++){
    // get diagonal term
    v = get(i,i);
    DI.set(i,i,1/v);
  }
  return DI;
}//

/*------------------------------------------------*/
/*  solve linear system with conjugate gradient   */
/*------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::conjugate_solve ( const SGTELIB::Matrix & A ,
                                                  const SGTELIB::Matrix & b ,
                                                  const SGTELIB::Matrix & x0 ,
                                                  const double tol) {

  const int n = x0.get_nb_rows();
  SGTELIB::Matrix x = x0;
  SGTELIB::Matrix r = b-A*x;
  double rr = r.normsquare();
  SGTELIB::Matrix p = r;
  SGTELIB::Matrix Ap;
  double rr_old,alpha,pAp;
  int iter = 0;
  while (iter < 100){
    Ap = A*p;
    pAp = 0;
    for (int i=0 ; i<n ; i++) pAp += p._X[i][0]*Ap._X[i][0];
    alpha = rr/pAp;
    x = x+alpha*p;
    rr_old = rr;
    r = r-alpha*Ap;
    rr = r.normsquare();
    if (rr < tol) break;
    p = r + (rr/rr_old)*p;

    Ap.set_name("Ap");
    x.set_name("x");
    r.set_name("r");
    p.set_name("p");
  }
  return x;
}//

/*---------------------------*/
/*  cholesky decomposition   */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::cholesky ( void ) const {

  if (get_nb_rows()!=get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::cholesky(): dimension error" );
  }

  const int n = get_nb_rows();
  SGTELIB::Matrix L ("L",n,n);

  double s;
  int i,j,k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < (i+1); j++) {
      s = 0;
      for (k = 0; k < j; k++){
        s += L._X[i][k] * L._X[j][k];
      }
      L._X[i][j] = (i == j) ?
         sqrt(_X[i][i] - s) :
         (1.0 / L._X[j][j] * (_X[i][j] - s));
    }
  }
  return L;
}//

/*---------------------------*/
/*  cholesky inverse         */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::cholesky_inverse ( void ) const {
  return cholesky_inverse(NULL);
}


SGTELIB::Matrix SGTELIB::Matrix::cholesky_inverse ( double * det ) const {
  SGTELIB::Matrix L  = cholesky();
  SGTELIB::Matrix Li = tril_inverse(L);

  const int n = _nbRows;

  // Compute A = Li'*Li
  // Note: by taking into account the fact that Li is tri inf,
  // It is possible to divide the cost of the computation
  // of Li'*Li by 3.
  SGTELIB::Matrix A ("A",n,n);
  int i,j,k,kmin;
  for (i=0 ; i<n ; i++){
    for (j=0 ; j<n ; j++){
      A._X[i][j] = 0;
      kmin = std::max(i,j);
      for (k=kmin ; k<n ; k++){
          A._X[i][j] += Li._X[k][i]*Li._X[k][j];
      }
    }
  }

  if (det){
    double v = 1;
    for (i=0 ; i<n ; i++) v *= L._X[i][i];
    v *= v;
    if ( isnan(v)) v=+INF;
    *det = v;
  }

  return A;
}//

/*-----------------------------------------*/
/*  Solve Upper Triangular Linear system   */
/*-----------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::triu_solve( const SGTELIB::Matrix & U , 
                                            const SGTELIB::Matrix & b ){
  const int n = U.get_nb_rows();
  if (n!=U.get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::triu_solve(): dimension error" );
  }
  if (n!=b.get_nb_rows()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::triu_solve(): dimension error" );
  }
  if (1!=b.get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::triu_solve(): dimension error" );
  }
  
  SGTELIB::Matrix x = b;

  for (int i=n-1 ; i>=0 ; i--){
    for (int j=i+1 ; j<n ; j++){
      x._X[i][0] -= U._X[i][j]*x._X[j][0];
    }
    x._X[i][0] /= U._X[i][i];
  }

  return x;
}//

/*-----------------------------------------*/
/*    Inverse Lower Triangular Matrix      */
/*-----------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::tril_inverse( const SGTELIB::Matrix & L ){
  const int n = L.get_nb_rows(); 
  SGTELIB::Matrix Li = L;
  SGTELIB::Matrix b ("b",n,1);
  
  for (int i=0 ; i<n ; i++){
    b.set(i,0,1.0);
    Li.set_col( SGTELIB::Matrix::tril_solve(L,b) , i);
    b.set(i,0,0.0);
  }

  return Li;
}//

/*-----------------------------------------*/
/*  Solve Lower Triangular Linear system   */
/*-----------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::tril_solve( const SGTELIB::Matrix & L , 
                                             const SGTELIB::Matrix & b ){
  const int n = L.get_nb_rows();
  if (n!=L.get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::tril_solve(): dimension error" );
  }
  if (n!=b.get_nb_rows()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::tril_solve(): dimension error" );
  }
  if (1!=b.get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::tril_solve(): dimension error" );
  }
  
  SGTELIB::Matrix x = b;

  for (int i=0 ; i<n ; i++){
    for (int j=0 ; j<i ; j++){
      x._X[i][0] -= L._X[i][j]*x._X[j][0];
    }
    x._X[i][0] /= L._X[i][i];
  }

  return x;
}//

/*-----------------------------------------*/
/*  Solve System with Cholesky             */
/*-----------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::cholesky_solve ( const SGTELIB::Matrix & A ,
                                                  const SGTELIB::Matrix & b ) {
  SGTELIB::Matrix L = A.cholesky();
  SGTELIB::Matrix y = tril_solve(L,b);
  SGTELIB::Matrix x = triu_solve(L.transpose(),y);
  return x;
}//

/*---------------------------*/
/*        SVD inverse        */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::SVD_inverse ( void ) const {
  
  if (get_nb_rows()!=get_nb_cols()){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
             "Matrix::SVD_inverse(): dimension error" );
  }

  // Init SVD matrices
  SGTELIB::Matrix * U;
  SGTELIB::Matrix * W;
  SGTELIB::Matrix * V;

  // Perform SVD
  std::string error_msg;
  SVD_decomposition ( error_msg , U, W, V, 1000000000 );

  // Inverse diag terms of W.
  for (int i=0 ; i<W->get_nb_rows() ; i++){
    W->set(i,i,1/W->get(i,i));
  }
  
  *U = U->transpose();
  SGTELIB::Matrix INVERSE (product ( *V , *W , *U ));
  INVERSE.set_name("inv("+_name+")");
  delete V;
  delete W;
  delete U;
  return INVERSE;
}//

/*---------------------------*/
/*        transpose          */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::transpose ( void ) const{
  SGTELIB::Matrix A (_name+"'",_nbCols,_nbRows); 
  for (int i=0 ; i<_nbCols ; i++){
    for (int j=0 ; j<_nbRows ; j++){
      A.set(i,j,_X[j][i]);
    }
  }
  return A;
}//

/*---------------------------*/
/*        diag               */
/*---------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::diag ( void ) const{

  SGTELIB::Matrix A;
  if (_nbCols==_nbRows){
    A = SGTELIB::Matrix("A",_nbRows,1);   
    for (int i=0 ; i<_nbCols ; i++) A.set(i,0,_X[i][i]);
  }
  else if ( (_nbCols==1) || (_nbRows==1) ){
    const int n=std::max(_nbCols,_nbRows);
    A = SGTELIB::Matrix("A",_nbRows,1);   
    for (int i=0 ; i<n ; i++) A.set(i,i,get(i));
  }
  else{
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Matrix::diag(): dimension error" );
  }

  A.set_name("diag("+_name+")");
  return A;
}//

/*--------------------------------------------------------------*/
/*                        SVD decomposition                     */
/*  inspired and recoded from an old numerical recipes version  */
/*--------------------------------------------------------------*/
/*                                                              */
/*           M = U . W . V'   (M: current matrix object)        */
/*                                                              */
/*           M ( nbRows x nbCols )                              */
/*           U ( nbRows x nbCols )                              */
/*           W ( nbCols x nbCols )                              */
/*           V ( nbCols x nbCols )                              */
/*                                                              */
/*           U.U' = U'.U = I if nbRows = nbCols                 */
/*           U'.U = I        if nbRows > nbCols                 */
/*                                                              */
/*           V.V' = V'.V = I                                    */
/*                                                              */
/*           W diagonal, given as a size-nbCols vector          */
/*                                                              */
/*           V is given, not V'                                 */
/*                                                              */
/*--------------------------------------------------------------*/

bool SGTELIB::Matrix::SVD_decomposition ( std::string & error_msg ,
                                          SGTELIB::Matrix *& MAT_U,  // OUT, nbRows x nbCols
                                          SGTELIB::Matrix *& MAT_W,  // OUT, nbCols x nbCols, diagonal
                                          SGTELIB::Matrix *& MAT_V,  // OUT, nbCols x nbCols
                                          int           max_mpn  ) const {

  // Dimension
  const int nbRows = _nbRows;
  const int nbCols = _nbCols;

  // init matrices for SVD
  double ** U = new double *[nbRows];
  double  * W = new double  [nbCols];
  double ** V = new double *[nbCols];
  for (int i = 0 ; i < nbCols ; ++i ) {
    U[i] = new double[nbCols];
    V[i] = new double[nbCols];
  }

  // call SVD
  bool result;
  result = this->SVD_decomposition ( error_msg , U , W , V , max_mpn );

  // Init matrix for result
  MAT_U = new SGTELIB::Matrix ("MAT_U",nbRows,nbCols);
  MAT_W = new SGTELIB::Matrix ("MAT_W",nbCols,nbCols);
  MAT_V = new SGTELIB::Matrix ("MAT_V",nbCols,nbCols);

  // Fill matrices
  for (int i=0 ; i<nbRows ; i++){
    for (int j=0 ; j<nbCols ; j++){
      MAT_U->set(i,j,U[i][j]); 
    }
  }
  for (int i=0 ; i<nbCols ; i++){
    for (int j=0 ; j<nbCols ; j++){
      MAT_V->set(i,j,V[i][j]); 
      MAT_W->set(i,j,0.0); 
    }
    MAT_W->set(i,i,W[i]); 
  }

  // Delete U, V, W
  for (int i=0 ; i<nbRows ; i++){
    delete [] U[i];
  }
  delete [] U;
  for (int j=0 ; j<nbCols ; j++){
    delete [] V[j];
  }
  delete [] V;
  delete [] W;

  return result;
}//

bool SGTELIB::Matrix::SVD_decomposition ( std::string & error_msg ,
            double     ** U         ,  // OUT, nbRows x nbCols
            double      * W         ,  // OUT, nbCols x nbCols, diagonal
            double     ** V         ,  // OUT, nbCols x nbCols
            int           max_mpn     ) const {
  const int nbRows = _nbRows;
  const int nbCols = _nbCols;
  
  error_msg.clear();

  if ( max_mpn > 0 && nbRows+nbCols > max_mpn ) {
    error_msg = "SVD_decomposition() error: nbRows+nbCols > " + SGTELIB::itos ( max_mpn );
    return false;
  }

  double * rv1   = new double[nbCols];
  double   scale = 0.0;
  double   g     = 0.0;
  double   norm  = 0.0;

  int      nm1   = nbCols - 1;

  bool   flag;
  int    i , j , k , l = 0 , its , jj , nm = 0;
  double s , f , h , tmp , c , x , y , z , absf , absg , absh;

  const int NITER = 30;

  // copy the current matrix into U:
  for ( i = 0 ; i < nbRows ; ++i )
    for ( j = 0 ; j < nbCols ; ++j )
      U[i][j] = _X[i][j];

  // Householder reduction to bidiagonal form:
  for ( i = 0 ; i < nbCols ; ++i ) {
    l      = i + 1;
    rv1[i] = scale * g;
    g      = s = scale = 0.0;
    if ( i < nbRows ) {
      for ( k = i ; k < nbRows ; ++k )
   scale += fabs ( U[k][i] );
      if ( scale != 0.0 ) {
   for ( k = i ; k < nbRows ; ++k ) {
     U[k][i] /= scale;
     s += U[k][i] * U[k][i];
  }
   f       = U[i][i];
  g       = ( f >= 0.0 ) ? -fabs(sqrt(s)) : fabs(sqrt(s));
   h       = f * g - s;
  U[i][i] = f - g;
   for ( j = l ; j < nbCols ; ++j ) {
     for ( s = 0.0 , k = i ; k < nbRows ; ++k )
       s += U[k][i] * U[k][j];
     f = s / h;
     for ( k = i ; k < nbRows ; ++k )
       U[k][j] += f * U[k][i];
   }
   for ( k = i ; k < nbRows ; ++k )
     U[k][i] *= scale;
      }
    }
    W[i] = scale * g;
    g    = s = scale = 0.0;
    if ( i < nbRows && i != nm1 ) {
      for ( k = l ; k < nbCols ; ++k )
   scale += fabs ( U[i][k] );
      if ( scale != 0.0 ) {
   for ( k = l ; k < nbCols ; ++k ) {
     U[i][k] /= scale;
     s       += U[i][k] * U[i][k];
   }
   f       = U[i][l];
  g       = ( f >= 0.0 ) ? -fabs(sqrt(s)) : fabs(sqrt(s));
   h       = f * g - s;
   U[i][l] = f - g;
   for ( k = l ; k < nbCols ; ++k )
     rv1[k] = U[i][k] / h;
   for ( j = l ; j < nbRows ; ++j ) {
     for ( s=0.0,k=l ; k < nbCols ; ++k )
       s += U[j][k] * U[i][k];
     for ( k=l ; k < nbCols ; ++k )
       U[j][k] += s * rv1[k];
   }
   for ( k = l ; k < nbCols ; ++k )
     U[i][k] *= scale;
      }
    }
    tmp  = fabs ( W[i] ) + fabs ( rv1[i] );
    norm = ( norm > tmp ) ? norm : tmp;
  }

  // accumulation of right-hand transformations:
  for ( i = nm1 ; i >= 0 ; --i ) {
    if ( i < nm1 ) {
      if ( g != 0.0 ) {
   for ( j = l ; j < nbCols ; ++j )
     V[j][i] = ( U[i][j] / U[i][l] ) / g;
   for ( j = l ; j < nbCols ; ++j ) {
     for ( s = 0.0 , k = l ; k < nbCols ; ++k )
       s += U[i][k] * V[k][j];
     for ( k = l ; k < nbCols ; ++k )
       V[k][j] += s * V[k][i];
   }
      }
      for ( j = l ; j < nbCols ; ++j )
   V[i][j] = V[j][i] = 0.0;
    }
    V[i][i] = 1.0;
    g       = rv1[i];
    l       = i;
  }

  // accumulation of left-hand transformations:
  for ( i = ( ( nbRows < nbCols ) ? nbRows : nbCols ) - 1 ; i >= 0 ; --i ) {
    l = i + 1;
    g = W[i];
    for ( j = l ; j < nbCols ; ++j )
      U[i][j] = 0.0;
    if ( g != 0.0 ) {
      g = 1.0 / g;
      for ( j = l ; j < nbCols ; ++j ) {
  for ( s = 0.0 , k = l ; k < nbRows ; ++k )
    s += U[k][i] * U[k][j];
  f = ( s / U[i][i] ) * g;
  for ( k = i ; k < nbRows ; ++k )
    U[k][j] += f * U[k][i];
      }
      for ( j = i ; j < nbRows ; ++j )
  U[j][i] *= g;
    }
    else
      for ( j = i ; j < nbRows ; ++j )
  U[j][i] = 0.0;
    ++U[i][i];
  }

  // diagonalization of the bidiagonal form:
  for ( k = nm1 ; k >= 0 ; --k ) {
    for ( its = 1 ; its <= NITER ; its++ ) {
      flag = true;
      for ( l = k ; l >= 0 ; l-- ) {
   nm = l - 1;
   if ( nm < 0 || fabs ( rv1[l]) + norm == norm ) {
     flag = false;
     break;
  }
   if ( fabs ( W[nm] ) + norm == norm )
    break;
      }
      if ( flag ) {
   c = 0.0;
   s = 1.0;
   for ( i = l ; i <= k ; i++ ) {
     f      = s * rv1[i];
     rv1[i] = c * rv1[i];
     if ( fabs(f) + norm == norm )
       break;
     g = W[i];

    absf = fabs(f);
    absg = fabs(g);
    h    = ( absf > absg ) ?
      absf * sqrt ( 1.0 + pow ( absg/absf , 2.0 ) ) :
      ( ( absg==0 ) ? 0.0 : absg * sqrt ( 1.0 + pow ( absf/absg , 2.0 ) ) );

     W[i] =  h;
     h    =  1.0 / h;
     c    =  g * h;
     s    = -f * h;
     for ( j = 0 ; j < nbRows ; ++j ) {
       y = U[j][nm];
       z = U[j][ i];
       U[j][nm] = y * c + z * s;
       U[j][ i] = z * c - y * s;
     }
   }
      }
      z = W[k];
      if ( l == k) {
  if ( z < 0.0 ) {
     W[k] = -z;
     for ( j = 0 ; j < nbCols ; j++ )
       V[j][k] = -V[j][k];
   }
   break;  // this 'break' is always active if k==0
      }
      if ( its == NITER ) {
  error_msg = "SVD_decomposition() error: no convergence in " +
              SGTELIB::itos ( NITER ) + " iterations";
  delete [] rv1;
   return false;
      }
      x  = W[l];
      nm = k - 1;
      y  = W[nm];
      g  = rv1[nm];
      h  = rv1[k];
      f  = ( (y-z) * (y+z) + (g-h) * (g+h) ) / ( 2.0 * h * y );
      
      absf = fabs(f);
      g    = ( absf > 1.0 ) ?
  absf * sqrt ( 1.0 + pow ( 1.0/absf , 2.0 ) ) :
  sqrt ( 1.0 + pow ( absf , 2.0 ) );

      f = ( (x-z) * (x+z) +
      h * ( ( y / ( f + ( (f >= 0)? fabs(g) : -fabs(g) ) ) ) - h ) ) / x;
      c = s = 1.0;

      for ( j = l ; j <= nm ; ++j ) {
   i = j + 1;
   g = rv1[i];
   y = W[i];
   h = s * g;
   g = c * g;

  absf = fabs(f);
  absh = fabs(h);
  z    = ( absf > absh ) ?
    absf * sqrt ( 1.0 + pow ( absh/absf , 2.0 ) ) :
    ( ( absh==0 ) ? 0.0 : absh * sqrt ( 1.0 + pow ( absf/absh , 2.0 ) ) );

   rv1[j] = z;
   c      = f / z;
   s      = h / z;
   f      = x * c + g * s;
   g      = g * c - x * s;
   h      = y * s;
   y     *= c;
   for ( jj = 0 ; jj < nbCols ; ++jj ) {
     x = V[jj][j];
     z = V[jj][i];
     V[jj][j] = x * c + z * s;
     V[jj][i] = z * c - x * s;
   }

  absf = fabs(f);
  absh = fabs(h);
  z    = ( absf > absh ) ?
    absf * sqrt ( 1.0 + pow ( absh/absf , 2.0 ) ) :
    ( ( absh==0 ) ? 0.0 : absh * sqrt ( 1.0 + pow ( absf/absh , 2.0 ) ) );

   W[j] = z;

   if ( z ) {
     z = 1.0 / z;
     c = f * z;
     s = h * z;
   }
   f = c * g + s * y;
   x = c * y - s * g;
   for ( jj = 0 ; jj < nbRows ; ++jj ) {
     y = U[jj][j];
     z = U[jj][i];
     U[jj][j] = y * c + z * s;
     U[jj][i] = z * c - y * s;
   }
       }
       rv1[l] = 0.0;
       rv1[k] = f;
       W  [k] = x;
     }
  }

  delete [] rv1;
  return true;
}//

/*--------------------------------*/
/* is there any NaN in the matrix */
/*--------------------------------*/
bool SGTELIB::Matrix::has_nan ( void ) const {
  int i , j;
  for ( i = 0 ; i < _nbRows ; ++i ){
    for ( j = 0 ; j < _nbCols ; ++j ){
      if ( isnan(_X[i][j])){
        return true;
      }
    }
  }
  return false;
}//

/*--------------------------------*/
/* is there any NaN in the matrix */
/*--------------------------------*/
bool SGTELIB::Matrix::has_inf ( void ) const {
  int i , j;
  for ( i = 0 ; i < _nbRows ; ++i ){
    for ( j = 0 ; j < _nbCols ; ++j ){
      if ( isinf(_X[i][j])){
        return true;
      }
    }
  }
  return false;
}//

/*--------------------------------*/
/* is there any NaN in the matrix */
/*--------------------------------*/
void SGTELIB::Matrix::replace_nan ( double d ) {
  int i , j;
  for ( i = 0 ; i < _nbRows ; ++i ){
    for ( j = 0 ; j < _nbCols ; ++j ){
      if ( isnan(_X[i][j])){
        _X[i][j] = d;
      }
    }
  }
}//

/*-------------------------------------------------*/
/* return the index of the largest value           */
/*-------------------------------------------------*/
int SGTELIB::Matrix::get_max_index (void ) {
  int i,j,k=0,kmax=0;
  double vmax = -SGTELIB::INF;
  // We use the same mono-indexation as in matlab:
  //   1 4 
  //   2 5
  //   3 6
  for ( j = 0 ; j < _nbCols ; ++j ){
    for ( i = 0 ; i < _nbRows ; ++i ){
      if (_X[i][j] > vmax){
        vmax = _X[i][j];
        kmax = k;
      }
      k++;
    }
  }
  return kmax;
}//

/*-------------------------------------------------*/
/* min & max of a matrix                           */
/*-------------------------------------------------*/
double SGTELIB::Matrix::min (void) {
  double d = +INF;
  int i,j;
  for ( j = 0 ; j < _nbCols ; ++j ){
    for ( i = 0 ; i < _nbRows ; ++i ){
      d = std::min(d,_X[i][j]);
    }
  }
  return d;
}//

double SGTELIB::Matrix::max (void) {
  double d = -INF;
  int i,j;
  for ( j = 0 ; j < _nbCols ; ++j ){
    for ( i = 0 ; i < _nbRows ; ++i ){
      d = std::max(d,_X[i][j]);
    }
  }
  return d;
}//

/*-------------------------------------------------*/
/* min and max of two matrices                     */
/*-------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::max ( const SGTELIB::Matrix & A , 
                                       const SGTELIB::Matrix & B ){
  const int nb_rows = A.get_nb_rows();
  const int nb_cols = A.get_nb_cols();

  if (B.get_nb_rows()!=nb_rows){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::max(A,B): dimension error" );
  }  
  if (B.get_nb_cols()!=nb_cols){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::max(A,B): dimension error" );
  }  

  // Init matrix
  SGTELIB::Matrix C("max("+A.get_name()+";"+B.get_name()+")",nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C._X[i][j] = std::max( A._X[i][j] , B._X[i][j] );
    }
  }
  return C;
}//

SGTELIB::Matrix SGTELIB::Matrix::min ( const SGTELIB::Matrix & A , 
                                       const SGTELIB::Matrix & B ){
  const int nb_rows = A.get_nb_rows();
  const int nb_cols = A.get_nb_cols();

  if (B.get_nb_rows()!=nb_rows){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::min(A,B): dimension error" );
  }  
  if (B.get_nb_cols()!=nb_cols){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "Matrix::min(A,B): dimension error" );
  }  

  // Init matrix
  SGTELIB::Matrix C("min("+A.get_name()+";"+B.get_name()+")",nb_rows,nb_cols);

  // Compute
  int i,j;
  for ( i = 0 ; i < nb_rows ; ++i ) {
    for ( j = 0 ; j < nb_cols ; ++j ){
      C._X[i][j] = std::min( A._X[i][j] , B._X[i][j] );
    }
  }
  return C;
}//

/*-------------------------------------------------*/
/* return the index of the smallest value          */
/*-------------------------------------------------*/
int SGTELIB::Matrix::get_min_index (void ) {
  int i,j,k=0,kmin=0;
  double vmin = +SGTELIB::INF;
  // We use the same mono-indexation as in matlab:
  //   1 4 
  //   2 5
  //   3 6
  for ( j = 0 ; j < _nbCols ; ++j ){
    for ( i = 0 ; i < _nbRows ; ++i ){
      if (_X[i][j] < vmin){
        vmin = _X[i][j];
        kmin = k;
      }
      k++;
    }
  }
  return kmin;
}//

/*-------------------------------------------------*/
/* return the index of the smallest value on row i */
/*-------------------------------------------------*/
int SGTELIB::Matrix::get_min_index_row ( const int i ) {
  int j, jmin=0;
  double vmin = +SGTELIB::INF;
  for ( j = 0 ; j < _nbCols ; ++j ){
    if (_X[i][j] < vmin){
      vmin = _X[i][j];
      jmin = j;
    }
  }
  return jmin;
}//

/*-------------------------------------------------*/
/* return the index of the smallest value on row i */
/*-------------------------------------------------*/
int SGTELIB::Matrix::get_min_index_col ( const int j ) {
  int i, imin=0;
  double vmin = +SGTELIB::INF;
  for ( i = 0 ; i < _nbRows ; ++i ){
    if (_X[i][j] < vmin){
      vmin = _X[i][j];
      imin = i;
    }
  }
  return imin;
}//

/*-------------------------------------------------*/
/* Return NORM2 distance                           */
/*-------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_distances_norm2 ( const SGTELIB::Matrix & A , 
                                                       const SGTELIB::Matrix & B ){
  const int n = A.get_nb_cols();
  if ( B.get_nb_cols()!=n ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "get_distances_norm2: dimension error" );
  }

  const int pa = A.get_nb_rows();
  const int pb = B.get_nb_rows();
  SGTELIB::Matrix D = SGTELIB::Matrix("D",pa,pb);
  double v,d;
  int ia, ib, j;

  for (ia=0 ; ia < pa ; ia++){
    for (ib=0 ; ib < pb ; ib++){
      // Distance between the point ia of the cache and the point ib of the matrix XXs
      v = 0;
      for (j=0 ; j < n ; j++){
        d = A._X[ia][j]-B._X[ib][j];
        v += d*d;
      }
      D._X[ia][ib] = sqrt(v);
    }
  }
  return D;
}//

/*-------------------------------------------------*/
/* Return NORM1 distance                           */
/*-------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_distances_norm1 ( const SGTELIB::Matrix & A , 
                                                       const SGTELIB::Matrix & B ){
  const int n = A.get_nb_cols();
  if ( B.get_nb_cols()!=n ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "get_distances_norm2: dimension error" );
  }

  const int pa = A.get_nb_rows();
  const int pb = B.get_nb_rows();
  SGTELIB::Matrix D = SGTELIB::Matrix("D",pa,pb);
  double v;
  int ia, ib, j;

  for (ia=0 ; ia < pa ; ia++){
    for (ib=0 ; ib < pb ; ib++){
      // Distance between the point ia of the cache and the point ib of the matrix XXs
      v = 0;
      for (j=0 ; j < n ; j++){
        v += fabs(A._X[ia][j]-B._X[ib][j]);
      }
      D._X[ia][ib] = v;
    }
  }
  return D;
}//

/*-------------------------------------------------*/
/* Return NORMINF distance                         */
/*-------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_distances_norminf ( const SGTELIB::Matrix & A , 
                                                         const SGTELIB::Matrix & B ){
  const int n = A.get_nb_cols();
  if ( B.get_nb_cols()!=n ){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "get_distances_norm2: dimension error" );
  }

  const int pa = A.get_nb_rows();
  const int pb = B.get_nb_rows();
  SGTELIB::Matrix D = SGTELIB::Matrix("D",pa,pb);
  double v;
  int ia, ib, j;

  for (ia=0 ; ia < pa ; ia++){
    for (ib=0 ; ib < pb ; ib++){
      // Distance between the point ia of the cache and the point ib of the matrix XXs
      v = 0;
      for (j=0 ; j < n ; j++){
        v = std::max( v , fabs(A._X[ia][j]-B._X[ib][j]) );
      }
      D._X[ia][ib] = v;
    }
  }
  return D;
}//


/*-------------------------------------------------*/
/* find_row                                        */
/* Check if the matrix has a row identical to R    */
/*-------------------------------------------------*/
int SGTELIB::Matrix::find_row (SGTELIB::Matrix & R){

  // If the matrix is empty, return false.
  if (_nbRows==0) return -1;

  // Then, check dimensions.
  if (R.get_nb_rows()!=1) 
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "find_row: dimension error" );
  if (R.get_nb_cols()!=_nbCols)
    throw SGTELIB::Exception ( __FILE__ , __LINE__ , "find_row: dimension error" );

  // Look for the row. 
  int i,j;
  bool diff;
  for (i=0 ; i<_nbRows ; i++){
    diff = false;    
    for (j=0 ; j<_nbCols ; j++){
      if (_X[i][j]!=R._X[0][j]){
        diff = true;
        break;
      }
    }
    if ( ! diff) return i;
  }
  return -1;




}//


/*-------------------------------------------------*/
/* Generate poll directions                        */
/*-------------------------------------------------*/
SGTELIB::Matrix SGTELIB::Matrix::get_poll_directions ( const SGTELIB::Matrix& scaling,
                                                       const SGTELIB::param_domain_t * domain,
                                                       double psize ) {

  int i,j,k;
  // i : index of the poll direction (rows of D and POLL)
  // (Each line of D and POLL is a poll direction)
  // j : index of the optimization variable (columns of D and POLL)
  // k : integer buffer, when needed.
  double d;

  // Number of variables
  const int N = scaling.get_nb_cols();
  SGTELIB::Matrix D("D",N,N);

  // Number of continuous variables
  int Ncont = 0;
  for (j=0 ; j<N ; j++){
    if (domain[j]==SGTELIB::PARAM_DOMAIN_CONTINUOUS) Ncont++;
  }

  // Generate directions for continuous variables
  if (Ncont>0){
    // Generate one random direction
    SGTELIB::Matrix v("v",1,N);
    for (j=0 ; j<N ; j++){
      if (domain[j]==SGTELIB::PARAM_DOMAIN_CONTINUOUS){
        v._X[0][j] = SGTELIB::quick_norm_rand();
      }
    }

    // Normalize v (Euclidian Norm)
    v = v/v.norm();

    // Build D (Householder matrix)
    for (i=0 ; i<N ; i++){
      if (domain[i]==SGTELIB::PARAM_DOMAIN_CONTINUOUS){
        for (j=0 ; j<N ; j++){
          D._X[i][j] = double(i==j)-2*v[i]*v[j];
        }
      }
    }
  } // END if (Ncont>0)



  double msize = std::min(psize*psize,psize);
  double rho = psize/msize;
  // Normalize directions
  for (i=0 ; i<N ; i++){

    // Fill continous dimensions with rand if necessary
    if (domain[i]!=SGTELIB::PARAM_DOMAIN_CONTINUOUS){
      for (j=0 ; j<N ; j++){
        if (domain[j]==SGTELIB::PARAM_DOMAIN_CONTINUOUS){
          D._X[i][j] = 2*uniform_rand()-1;
        }
      }
    }

    // Find max asb
    d = 0;
    for (j=0 ; j<N ; j++) d = std::max( d , fabs(D._X[i][j]) );

    // Scale continuous dimensions
    for (j=0 ; j<N ; j++){
      if (domain[j]==SGTELIB::PARAM_DOMAIN_CONTINUOUS){
        D._X[i][j] = scaling[j]*msize*SGTELIB::rceil(rho*D._X[i][j]/d);
      }
    }

    // Add extended POLL for discrete values
    if ( (domain[i]==SGTELIB::PARAM_DOMAIN_INTEGER) ||
         (domain[i]==SGTELIB::PARAM_DOMAIN_BOOL) ){
      D._X[i][i] = (i%2==0)?-1:+1;
    }
    else if (domain[i]==SGTELIB::PARAM_DOMAIN_CAT){
      D._X[i][i] = SGTELIB::rceil(uniform_rand()*scaling[i]);
    }

  }


  // Add opposite directions and sort
  SGTELIB::Matrix POLL("POLL-DIR",2*N,N);
  k = 0;
  for (i=0 ; i<N ; i++){
    if (domain[i]==SGTELIB::PARAM_DOMAIN_CONTINUOUS){
        POLL.set_row(D.get_row(i),k++);
        POLL.set_row(-D.get_row(i),k++);
    }
  }
  for (i=0 ; i<N ; i++){
    if (domain[i]!=SGTELIB::PARAM_DOMAIN_CONTINUOUS){
        POLL.set_row(D.get_row(i),k++);
        POLL.set_row(-D.get_row(i),k++);
    }
  }

  if (k!=2*N){
    std::cout << "k,N : " << k << " " << N << "\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Unconcistency in the value of k." );
  }

  return POLL;
}//

/*-------------------------------------------------*/
/* Swap two rows                                   */
/*-------------------------------------------------*/
void SGTELIB::Matrix::swap_rows(const int i1 , const int i2){
  double buffer;
  for (int j=0 ; j<_nbCols ; j++){
    buffer    = _X[i1][j];
    _X[i1][j] = _X[i2][j];
    _X[i2][j] = buffer;
  }
}//

/*-------------------------------------------------*/
/* LU inverse                                      */
/*-------------------------------------------------*/

SGTELIB::Matrix SGTELIB::Matrix::lu_inverse ( void ) const {
  return lu_inverse(NULL);
}

SGTELIB::Matrix SGTELIB::Matrix::lu_inverse ( double * det ) const{

  const int N = _nbRows;
  SGTELIB::Matrix A (*this);

  int i,j,k,ip=0;
  double pivot,pivot_max;

  // Permuation vector
  int * P = new int [N];
  for (i=0 ; i<N ; i++) P[i]=i;

  // LU factorization (in-place)
  for (k=0 ; k<N-1 ; k++){

    // Find pivot
    pivot_max = -1;
    for (i=k ; i<N ; i++){
      pivot = A._X[k][i];
      if (pivot<0) pivot*=-1;
      if (pivot>pivot_max){
        ip = i;
        pivot_max = pivot;
      }
    }

    // Swap rows of A and P
    if (ip!=k){
      A.swap_rows(ip,k);
      i=P[ip]; P[ip]=P[k]; P[k]=i;
    }

    // Gaussian elimination
    for (j=k+1 ; j<N ; j++){
      pivot = A._X[j][k]/A._X[k][k];
      A._X[j][k] = pivot;
      for (i=k+1 ; i<N ; i++) A._X[j][i] -= pivot*A._X[k][i];
    }
  }

  // Construct the whole matrix P (under the name Ai)
  SGTELIB::Matrix Ai ("Ai",N,N);
  for (i=0 ; i<N; i++) Ai._X[i][P[i]] = 1;



  if (det){
    // Compute the determinant of the matrix that is inverted
    double v = 1;
    // Compute the determinant of U
    for (i=0 ; i<N ; i++) v *= A._X[i][i];
    // Comput ethe determinant of P
    i = 0;
    while (i<N){
      if (P[i]!=i){
        j = P[i];
        P[i] = P[j];
        P[j] = j;
        v *= -1;
      }
      else i++;
    }
    *det = v;
  }
  
  // Triangular inversion for each column of Ai.
  SGTELIB::Matrix y;
  for (k=0 ; k<N; k++){
    y = Ai.get_col(k);
    
    // Tri-L solve
    for (i=0 ; i<N ; i++){
      for (j=0 ; j<i ; j++){
        y._X[i][0] -= A._X[i][j]*y._X[j][0];
      }
    }

    // Tri-U solve
    for (int i=N-1 ; i>=0 ; i--){
      for (int j=i+1 ; j<N ; j++){
        y._X[i][0] -= A._X[i][j]*y._X[j][0];
      }
      y._X[i][0] /= A._X[i][i];
    }

    Ai.set_col(y,k);
  }

  delete [] P;
  
  Ai.set_name(_name+"^-1");
  return Ai;

}//




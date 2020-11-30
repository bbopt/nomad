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

#include "Surrogate_PRS.hpp"

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS::Surrogate_PRS ( SGTELIB::TrainingSet & trainingset,
                                        SGTELIB::Surrogate_Parameters param) :
  SGTELIB::Surrogate ( trainingset , param ),
  _q                 ( 0           ),
  _M                 ( "M",0,0     ),
  _H                 ( "H",0,0     ),
  _Ai                ( "Ai",0,0    ),
  _alpha             ( "alpha",0,0 ){
  #ifdef SGTELIB_DEBUG
    std::cout << "constructor PRS\n";
  #endif
}//


/*----------------------------*/
/*          destructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS::~Surrogate_PRS ( void ) {

}//


/*----------------------------*/
/*          display           */
/*----------------------------*/
void SGTELIB::Surrogate_PRS::display_private ( std::ostream & out ) const {
  out << "q: " << _q << "\n";
}//


/*--------------------------------------*/
/*               build                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_PRS::build_private ( void ) {
  
  const int pvar = _trainingset.get_pvar(); 
  const int nvar = _trainingset.get_nvar(); 

  // Get the number of basis functions.
  _q = Surrogate_PRS::get_nb_PRS_monomes(nvar,_param.get_degree());

  // If _q is too big or there is not enough points, then quit
  if (_q>200) return false;
  if ( (_q>pvar-1) && (_param.get_ridge()==0) ) return false;

  // Compute the exponents of the basis functions
  _M = get_PRS_monomes(nvar,_param.get_degree());

  // DESIGN MATRIX H
  _H = compute_design_matrix ( _M , get_matrix_Xs() );

  // Compute alpha
  if ( !  compute_alpha()) return false;

  _ready = true; 
  return true;
}//

/*--------------------------------------*/
/*          Compute PRS matrix          */
/*--------------------------------------*/
const SGTELIB::Matrix SGTELIB::Surrogate_PRS::compute_design_matrix ( const SGTELIB::Matrix& Monomes, 
                                                                      const SGTELIB::Matrix & Xs ) {

  const int n = Xs.get_nb_cols(); // Nb of points in the matrix X given in argument
  const int p = Xs.get_nb_rows(); // Nb of points in the matrix X given in argument
  double v;
  int i,j,jj,k,exponent;

  const int nbMonomes = Monomes.get_nb_rows();

  // Init the design matrix  
  SGTELIB::Matrix H("H",p,nbMonomes);
  // Current basis function (vector column to construct 1 basis function)
  SGTELIB::Matrix h("h",p,1);

  // j is the corresponding index among all input (j in [0;n-1])
  // jj is the index of the input variabe amongst the varying input (jj in [0;nvar-1])
  // k is the index of the monome (ie: the basis function) (k in [0;q-1])
  // i is the index of a point (i in [0;p-1])
  // Loop on the monomes
  for (k=0 ; k<nbMonomes ; k++){
    h.fill(1.0);
    jj=0;
    // Loop on the input variables
    for (j=0 ; j<n ; j++){
      if (_trainingset.get_X_nbdiff(j)>1){
        exponent = int(Monomes.get(k,jj)); 
        if (exponent>0){
          for (i=0 ; i<p ; i++){
            v = h.get(i,0);
            v *= pow(Xs.get(i,jj),exponent);
            h.set(i,0,v);
          }
        }
        jj++;
      }
    }
    H.set_col(h,k);
  }
  return H;
}//


/*--------------------------------------*/
/*       compute alpha                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_PRS::compute_alpha ( void ){

  const SGTELIB::Matrix   Ht = _H.transpose();
  const SGTELIB::Matrix & Zs = get_matrix_Zs();

  // Ridge
  double r = _param.get_ridge();

  // COMPUTE COEFS
  if (r>0){
    //_Ai = (Ht*_H+r*SGTELIB::Matrix::identity(_q)).SVD_inverse();
    _Ai = (Ht*_H+r*SGTELIB::Matrix::identity(_q)).cholesky_inverse();
  }
  else{
    //_Ai = (Ht*_H).SVD_inverse();
    _Ai = (Ht*_H).cholesky_inverse();
  }

  _alpha = _Ai * (Ht * Zs);
  _alpha.set_name("alpha");
  if (_alpha.has_nan()){
    return false;
  }
  return true;
}

/*--------------------------------------*/
/*       predict (ZZs only)             */
/*--------------------------------------*/
void SGTELIB::Surrogate_PRS::predict_private ( const SGTELIB::Matrix & XXs,
                                                     SGTELIB::Matrix * ZZs ) {
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  *ZZs = compute_design_matrix(_M,XXs) * _alpha;
}//

/*--------------------------------------*/
/*       compute Zvs                    */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_PRS::get_matrix_Zvs (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  // Not necessary. Zv is computed in "build".
  if ( !  _Zvs){
    _Zvs = new SGTELIB::Matrix;
    // Projection matrix
    const SGTELIB::Matrix & Zs = get_matrix_Zs();
    SGTELIB::Matrix dPiPZs = SGTELIB::Matrix::get_matrix_dPiPZs(_Ai,_H,Zs);

    // dPi is the inverse of the diag of P 
    // Compute _Zv = Zs - dPi*P*Zs
    *_Zvs = Zs - dPiPZs;
    _Zvs->replace_nan(+INF);
    _Zvs->set_name("Zvs");
  }
  return _Zvs;
}//



/*-----------------------------------------*/
/* Compute the theorical number of monomes */
/*-----------------------------------------*/
int SGTELIB::Surrogate_PRS::get_nb_PRS_monomes(const int nvar, const int degree){
  // Return the number of lines in the matrix M computed in get_PRS_monomes()
  int S = 1;
  int v = nvar;
  for (int d = 1 ; d<=degree ; d++){
    S += v;
    v = (v*(nvar+d))/(d+1);
  }
  return S;
}//



/*----------------------------------*/
/*     BUILD THE INDEX MATRICES     */
/*----------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate_PRS::get_PRS_monomes(const int nvar, const int degree){

  double * z = new double [nvar];
  SGTELIB::Matrix M("Monomes",1,nvar);
  bool continuer;

  int i,j,di,ci;

  // Loop on the number of non null terms in z
  // c is the number of non-null terms of the monome.
  // We start with the monomes with only one non-null term.
  for (int c=1 ; c<=std::min(degree,nvar) ; c++){
    for (int d=c ; d<=degree ; d++){
          
      // First monome (c,d)
      z[0] = d-c+1;
      for (i=1 ; i<c ; i++) 
        z[i] = 1;
      for (i=c ; i<nvar ; i++) 
        z[i] = 0;

      continuer = true;
      while (continuer){
        M.add_row(z);
        // Pivot
        i = 0;
        while ( (i<nvar-1) && (z[i]<=z[i+1]) && ( (z[i]<=1) || (z[i+1]>=d-c+1))  )
          i++;
        // Transfert
        if (i < nvar-1){
          z[i+1]++;
          for (j=0; j<=i ; j++){
            z[j] = 0;
          }
          ci = c;
          di = d;
          for (j=i+1 ; j<nvar ; j++){
            ci -= (z[j]!=0);
            di -= static_cast<int>(z[j]);
          }
          if ( (ci==0) && (di>0) ){
            z[i+1] = z[i+1]+di;
          }
          else{
            for (int j=0; j<ci; j++){
              z[j] = 1;
              z[0] -= z[j];
            }
            z[0] += di;
          }
        }
        else{
            continuer = false;
        }
      } // loop while
    }// loop d
  }// loop c
  delete [] z;
  return M;
}//





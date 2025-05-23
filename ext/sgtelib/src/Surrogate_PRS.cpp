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

#include "Surrogate_PRS.hpp"

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS::Surrogate_PRS ( SGTELIB::TrainingSet & trainingset,
                                        const SGTELIB::Surrogate_Parameters& param) :
  SGTELIB::Surrogate ( trainingset , param ),
  _q                 ( 0           ),
  _M                 ( "M",0,0     ),
  _H                 ( "H",0,0     ),
  _Ai                ( "Ai",0,0    ),
  _alpha             ( "alpha",0,0 ),
  _preComputeForJacobianAndHessianDone( false){
  #ifdef SGTELIB_DEBUG
    std::cout << "constructor PRS\n";
  #endif
}//


/*----------------------------*/
/*          destructor        */
/*----------------------------*/
SGTELIB::Surrogate_PRS::~Surrogate_PRS ( void ) {
    
    if (_preComputeForJacobianAndHessianDone)
    {
        for (int i=0 ; i<_n ; i++)
        {
            delete [] _M_dxidxj[i];
            delete [] _alpha_dxidxj[i];
            
        }
        delete [] _M_dxidxj;
        delete [] _alpha_dxidxj;
        
        delete [] _M_dxj;
        delete [] _alpha_dxj;
    }

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
    
    _M = Matrix("M",0,0);
    _H = Matrix( "H",0,0);
    _Ai = Matrix( "Ai",0,0    );
    _alpha = Matrix( "alpha",0,0 );

  // If _q is too big or there is not enough points, then quit
  if (_q>200)
      return false;
    
// Let's try to building model with less points by adding a small ridge
  if ( (_q>pvar) && (_param.get_ridge()==0) )
          _param.set_ridge(0.001);

  // Compute the exponents of the basis functions
  _M = get_PRS_monomes(nvar,_param.get_degree());

  // DESIGN MATRIX H
  _H = compute_design_matrix ( _M , get_matrix_Xs() );

  // Compute alpha
  if ( !  compute_alpha())
      return false;

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
      if (_trainingset.get_X_nbdiff(j)>1)
      {
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

/*-------------------------------------------------*/
/*      Helper to compute PRS grad matrix          */
/*-------------------------------------------------*/
void SGTELIB::Surrogate_PRS::compute_dxi_matrices (
                                            SGTELIB::Matrix& Monomes,
                                            SGTELIB::Matrix& alpha,
                                            const int i ) const
{
    
    int l,k,exponent,d_exponent;
    
    const int nbMonomes = Monomes.get_nb_rows();
    const int m = alpha.get_nb_cols(); // Number of outputs
    
    if (nbMonomes != alpha.get_nb_rows())
    {
        throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                                  "compute_dxi_matrices: monome and alpha are incompatible !" );
    }
    
    // Loop on the monomes
    for (k=0 ; k<nbMonomes ; k++)
    {
        exponent = int(Monomes.get(k,i));
        d_exponent = ((exponent > 0) ? exponent-1:0);
        Monomes.set(k,i,d_exponent);
        for (l = 0 ; l < m ; l++)
        {
            alpha.set(k,l,alpha.get(k,l)*double(exponent));
        }
    }
}//



/*--------------------------------------*/
/*       compute alpha                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_PRS::compute_alpha ( void ){

  const SGTELIB::Matrix   Ht = _H.transpose();
  const SGTELIB::Matrix & Zs = get_matrix_Zs();

  // Ridge
  double r = _param.get_ridge();

    if (_H.has_inf() || _H.has_nan() || Ht.has_inf() || Ht.has_nan())
    {
        return false;
    }
  // COMPUTE COEFS
  if (r>0)
  {
    _Ai = (Ht*_H+r*SGTELIB::Matrix::identity(_q)).SVD_inverse();
    //_Ai = (Ht*_H+r*SGTELIB::Matrix::identity(_q)).cholesky_inverse();
  }
  else
  {
      _Ai = (Ht*_H).SVD_inverse();
      //_Ai = (Ht*_H).cholesky_inverse();
      
      // We may not have enough points to compute all coefficients of the monome
      // Let's try with a small ridge
      if (_Ai.has_nan())
      {
          r = 1E-3;
          _Ai = (Ht*_H+r*SGTELIB::Matrix::identity(_q)).SVD_inverse();
          
      }
  }
  if (_Ai.has_nan())
  {
      return false;
  }
    
  _alpha = _Ai * (Ht * Zs);
    
    SGTELIB::Matrix sing_val = (Ht*_H).get_singular_values();
    double sing_val_min = sing_val.min();
    if (sing_val_min > 0)
    {
        _cond = sing_val.max()/sing_val_min;
    }
    else
    {
        _cond = INF;
    }
    
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
/*       predict obj Z(ZZs only)             */
/*--------------------------------------*/
void SGTELIB::Surrogate_PRS::predict_obj ( const SGTELIB::Matrix & XX,
                                            SGTELIB::Matrix * ZZs,
                                           bool forceOnEmptyTrainingSet )
{

    SGTELIB::Matrix grad("grad",_m,_n);
    SGTELIB::Matrix hess("hess",_n,_n);

    predict_grad(XX, &grad, forceOnEmptyTrainingSet);

    for (int k=0; k < _m; ++k)
    {
      predict_hessian(XX, &hess, k, forceOnEmptyTrainingSet);
      for (int i=0; i < _n; ++i)
      {
        ZZs->set(0, k, grad.get(k, i) * XX.get(0, i));
        for (int j=0; j < _n; ++j)
        {
          ZZs->set(0, k, ZZs->get(0, k) + XX.get(0, i) * hess.get(i, j) * XX.get(0, j));
        }
      }
    }
}//

/*--------------------------------------*/
/*       predict grad Z(ZZs only)             */
/*--------------------------------------*/
void SGTELIB::Surrogate_PRS::predict_grad ( const SGTELIB::Matrix & XX,
                                            SGTELIB::Matrix * ZZs,
                                           bool forceOnEmptyTrainingSet )
{
    if ( ! _ready && forceOnEmptyTrainingSet )
    {
        _M = get_PRS_monomes(_n, _param.get_degree());
    }
    
    // Scale the points
    SGTELIB::Matrix XXs(XX);
    
    if (!forceOnEmptyTrainingSet)
    {
        _trainingset.X_scale(XXs);
    }

    Matrix dFdxi("dFdxi",1,static_cast<int>(_m));
    
    preComputeForJacobianAndHessian();
    
    int ii = 0;  // ii is for index of true variables in the model
    for (int i = 0 ; i < _n ; i++ )
    {
        if ( _trainingset.get_X_nbdiff(i) >1)
        {
            
            SGTELIB::Matrix MX = compute_design_matrix(_M_dxj[i], XXs);
            // UnScaled output (dF/dX)*(sX/sF)
            for(int j = 0 ; j < _m ; j++)
            {
                dFdxi.set_col( MX * _alpha_dxj[i].get_col(j), j); // It should be a single value
                if (! forceOnEmptyTrainingSet)
                {
                    dFdxi.set(0,j,dFdxi.get(0,j)* _trainingset.get_X_scaling_a(ii)/_trainingset.get_Z_scaling_a(j)) ;
                }
            }
            ii++;
        }
        else
        {
            dFdxi.fill(0.0);
        }
        
        // Use UnScaled output
        ZZs->set_col(dFdxi.transpose(),i);
    }
}//

void SGTELIB::Surrogate_PRS::preComputeForJacobianAndHessian()
{
    if (_preComputeForJacobianAndHessianDone)
        return;

    _M_dxidxj = new SGTELIB::Matrix *[_n];
    _alpha_dxidxj = new SGTELIB::Matrix *[_n];
    
    _M_dxj = new SGTELIB::Matrix[_n];
    _alpha_dxj = new SGTELIB::Matrix[_n];
    
    int jj= 0;
    for (int  j = 0 ; j < _n ; ++j )
    {
        _M_dxidxj[j] = new SGTELIB::Matrix[_n];
        _alpha_dxidxj[j] = new SGTELIB::Matrix[_n];
        
        _M_dxj[j] = SGTELIB::Matrix(_M);
        _alpha_dxj[j] = SGTELIB::Matrix(_alpha);
        
        if (_trainingset.get_X_nbdiff(j) > 1)
        {
            compute_dxi_matrices(_M_dxj[j],_alpha_dxj[j],jj);
            
            int ii = 0 ; //index of true variables in the model
            for (int i = 0 ; i < _n ; i++ )
            {
                _M_dxidxj[j][i] =  SGTELIB::Matrix(_M_dxj[j]);
                _alpha_dxidxj[j][i] = SGTELIB::Matrix(_alpha_dxj[j]);

                if (_trainingset.get_X_nbdiff(i) >1)
                {
                    compute_dxi_matrices(_M_dxidxj[j][i],_alpha_dxidxj[j][i],ii);
                    ii++;
                }
            }
            jj++;
        }
    }

    _preComputeForJacobianAndHessianDone = true;
}


/*--------------------------------------*/
/*   predict hessian Z(ZZs only)       */
/*--------------------------------------*/
void SGTELIB::Surrogate_PRS::predict_hessian(const SGTELIB::Matrix & XX,
                                            SGTELIB::Matrix * ZZs,
                                            int k /* output k */,
                                            bool forceOnEmptyTrainingSet )
{
    if ( ! _ready && forceOnEmptyTrainingSet )
    {
        _M = get_PRS_monomes(_n, _param.get_degree());
    }
    
    // Scale the points
    SGTELIB::Matrix XXs(XX);
    
    if (!forceOnEmptyTrainingSet)
    {
        _trainingset.X_scale(XXs);
    }
    
    preComputeForJacobianAndHessian();

    double d2Fkdxidxj=0.0;
    int jj = 0 ;  //index of true variables in the model
    for (int j = 0 ; j < _n ; j++ )
    {
        if (_trainingset.get_X_nbdiff(j) > 1)
        {
            //Matrix M_dxj(_M);
            //Matrix alpha_dxj(_alpha);
            //compute_dxi_matrices(M_dxj,alpha_dxj,jj);
            int ii = 0 ; //index of true variables in the model
            for (int i = 0 ; i < _n ; i++ )
            {
                
                if (_trainingset.get_X_nbdiff(i) >1)
                {
                    
                    //Matrix M_dxidxj(M_dxj);
                    //Matrix alpha_dxidxj(alpha_dxj);

                    //compute_dxi_matrices(M_dxidxj,alpha_dxidxj,ii);
                    
                    // Matrix tmp = compute_design_matrix(M_dxidxj, XXs) * alpha_dxidxj.get_col(k);
                    Matrix tmp = compute_design_matrix(_M_dxidxj[j][i], XXs) * _alpha_dxidxj[j][i].get_col(k);

                    
                    d2Fkdxidxj = tmp.get(0,0);
                    
                    // UnScaled output (d2F/dXi/dXj)*(sXi*sXj/sF^2)
                    if (! forceOnEmptyTrainingSet)
                    {
                        d2Fkdxidxj = d2Fkdxidxj* _trainingset.get_X_scaling_a(ii)*_trainingset.get_X_scaling_a(jj)/_trainingset.get_Z_scaling_a(k) ;
                        
                    }
                    ii++;
                }
                else
                {
                    d2Fkdxidxj = 0.0;
                }
                //  Use UnScaled output
                ZZs->set(i,j,d2Fkdxidxj);
            }
            jj++;
        }
        else
        {
            for (int i = 0 ; i < _n ; i++ )
            {
                //  Use UnScaled output
                ZZs->set(i,j,0.0);
            }
        }
        
    }
}



// Predict only objectives (used in Surrogate Ensemble Stat)
void SGTELIB::Surrogate_PRS::predict_private_objective ( const std::vector<SGTELIB::Matrix *> & XXd,
                                                         SGTELIB::Matrix * ZZsurr_around            ) {
  check_ready(__FILE__,__FUNCTION__,__LINE__);

  const size_t pxx = XXd.size();
  SGTELIB::Matrix _alpha_obj ("alpha_obj", _q, 1);

  // Get only objectives values is alpha
  for (int j=0 ; j<_m ; j++)
  {
    if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ)
    {
      _alpha_obj = _alpha.get_col(j);
      break;
    }
  } // end for j

  // Loop on all pxx points 
  for (int i=0 ; i<static_cast<int>(pxx) ; i++)
  {
    // XXd[i] is of dimension nbd * _n
    ZZsurr_around->set_row( ( compute_design_matrix(_M, *(XXd[i])) * _alpha_obj ).transpose() , i );
  } // end for i

}



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
  for (int d = 1 ; d<=degree ; d++)
  {
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


/*--------------------------------------*/
/*        Get model coefficients        */
/*--------------------------------------*/
SGTELIB::Matrix SGTELIB::Surrogate_PRS::get_alpha() const
{
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  return _alpha;
}

// Get all model outputs on point x
SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelOut(
  const SGTELIB::Matrix & XX,
  bool forceOnEmptyTrainingSet)
{

    // Init the matrices for prediction
    // Creation of matrix for input / output of SGTELIB model
    SGTELIB::Matrix Mpredict (  "M_predict", 1, static_cast<int>(_m));
/*
    if (XX.get_nb_cols() != _n){
    display(std::cout);
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
                 "predict(): dimension error" );
    }
    *ZZ = SGTELIB::Matrix("ZZ",XX.get_nb_rows(),_m);
*/
    if ( ! _ready && forceOnEmptyTrainingSet )
    {
        _M = get_PRS_monomes(_n, _param.get_degree());
    }
    
    // Scale the points
    SGTELIB::Matrix XXs(XX);
    
    if (!forceOnEmptyTrainingSet)
    {
        _trainingset.X_scale(XXs);
    }

    predict_private( XXs , &Mpredict );
  /*
    {
        check_ready(__FILE__,__FUNCTION__,__LINE__);
        
        if ((XX.get_nb_rows() == _n) && (XX.get_nb_cols() == 1) ) {
          this->predict(XX.transpose(), &Mpredict);
        } else {
          this->predict(XX, &Mpredict);
        }
    }
    */
   
    return _trainingset.Z_unscale(Mpredict);
}

double SGTELIB::Surrogate_PRS::getModelObj(
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{

    int j;
    for (j = 0; j <_m ; j++)
    {
        if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ)
        {
            if ((X.get_nb_rows() == _n) && (X.get_nb_cols() == 1))
            {
              return getModelOut(X.transpose(), j, forceOnEmptyTrainingSet);
            }
            else
            {
              return getModelOut(X, j, forceOnEmptyTrainingSet);
            }
        }
    }
    if (j == _m-1)
    {
        throw SGTELIB::Exception ( __FILE__ , __LINE__ , "No obj");
    }
    return 0.0;
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelGrad(
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{
    SGTELIB::Matrix Grad (  "Grad", static_cast<int>(_n), static_cast<int>(1));
    SGTELIB::Matrix Mpredict_grad (  "grad_predict", static_cast<int>(_m), static_cast<int>(_n));
    getModelGrad(&Grad, &Mpredict_grad, X, forceOnEmptyTrainingSet);
    return Grad;
}

void SGTELIB::Surrogate_PRS::getModelGrad(
  SGTELIB::Matrix * g,
  SGTELIB::Matrix * Mpredict_grad, // _m x _n matrix
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{

    if ((X.get_nb_rows() == _n) && (X.get_nb_cols() == 1))
    {
        predict_grad(X.transpose(), Mpredict_grad, forceOnEmptyTrainingSet);
    }
    else
    {
        predict_grad(X, Mpredict_grad, forceOnEmptyTrainingSet);
    }

    int j;
    for ( j = 0 ; j <_m ; j++)
    {
        if (_trainingset.get_bbo(j)==SGTELIB::BBO_OBJ)
        {
            *g = Mpredict_grad->get_row(j).transpose();
        }
    }
    if (j == _m-1)
    {
        throw SGTELIB::Exception ( __FILE__ , __LINE__ , "No obj");
    }
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelHessian(
  const SGTELIB::Matrix & X,
  int j,
  bool forceOnEmptyTrainingSet
)
{
    SGTELIB::Matrix Mpredict_hessian (  "hessian_predict", static_cast<int>(_n), static_cast<int>(_n));
    getModelHessian(&Mpredict_hessian, X, j, forceOnEmptyTrainingSet);
    return Mpredict_hessian;
}

void SGTELIB::Surrogate_PRS::getModelHessian(
  SGTELIB::Matrix * H,
  const SGTELIB::Matrix & X,
  int j,
  bool forceOnEmptyTrainingSet
)
{
    
    if ((X.get_nb_rows() == _n) && (X.get_nb_cols() == 1))
    {
        predict_hessian(X.transpose(), H, j, forceOnEmptyTrainingSet);
    }
    else
    {
        predict_hessian(X, H, j, forceOnEmptyTrainingSet);
    }
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelHessian(
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{
    SGTELIB::Matrix Hessian (  "Hessian", static_cast<int>(_n), static_cast<int>(_n));
    getModelHessian(&Hessian, X, forceOnEmptyTrainingSet);
    return Hessian;
}

void SGTELIB::Surrogate_PRS::getModelHessian(
  SGTELIB::Matrix * H,
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{
    int j;
    for (j=0 ; j < _m ; j++)
    {
        if (_trainingset.get_bbo(j) == SGTELIB::BBO_OBJ)
        {
            getModelHessian(H, X, j, forceOnEmptyTrainingSet) ;
        }
    }
    if (j == _m-1)
    {
        throw Exception(__FILE__, __LINE__, "Assertion error: no objective");
    }
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelCons(
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{

    SGTELIB::Matrix cons("cx", static_cast<int>(_m - 1), 1);
    getModelCons(&cons, X, forceOnEmptyTrainingSet);
    return cons;
}

void SGTELIB::Surrogate_PRS::getModelCons(
  SGTELIB::Matrix * cons,
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{
    int k = 0;
    int n = X.get_nb_cols();
    if (n != _n){
      throw Exception ( __FILE__ , __LINE__ ,
                  "TrainingSet::TrainingSet(): dimension error" );
    }

    for (int j = 0; j <_m ; j++)
    {
      
        if (_trainingset.get_bbo(j)!=SGTELIB::BBO_OBJ)
        {
            cons->set(k, 0, getModelOut(X, j, forceOnEmptyTrainingSet));
            k += 1;
        }
    }
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelJacobian(
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{
    SGTELIB::Matrix Mpredict_grad (  "grad_predict", static_cast<int>(_m), static_cast<int>(_n));
    int ncon = _m - 1; // one objective
    SGTELIB::Matrix J (  "Jx", ncon, static_cast<int>(_n));
    getModelJacobian(&J, &Mpredict_grad, X, forceOnEmptyTrainingSet);
    return J;
}

void SGTELIB::Surrogate_PRS::getModelJacobian(
  SGTELIB::Matrix * J,
  SGTELIB::Matrix * Mpredict_grad,
  const SGTELIB::Matrix & X,
  bool forceOnEmptyTrainingSet
)
{
    if ((X.get_nb_rows() == _n) && (X.get_nb_cols() == 1))
    {
        predict_grad(X.transpose(), Mpredict_grad, forceOnEmptyTrainingSet);
    }
    else
    {
        predict_grad(X, Mpredict_grad, forceOnEmptyTrainingSet);
    }

    int k = 0;
    for (int j=0; j < _m; ++j)
    {
        if (_trainingset.get_bbo(j)!=SGTELIB::BBO_OBJ)
        {
            for (int i=0; i < _n; ++i)
            {
                J->set(k, i, Mpredict_grad->get(j, i));
            }
            k++;
        }
    }
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelLagGrad(
  const SGTELIB::Matrix & X,
  const SGTELIB::Matrix & Y,
  double sigma,
  bool forceOnEmptyTrainingSet
)
{
    SGTELIB::Matrix Mpredict_grad (  "grad_predict", static_cast<int>(_m), static_cast<int>(_n));
    int ncon = _m - 1; // one objective
    SGTELIB::Matrix Jx (  "Jx", ncon, static_cast<int>(_n));
    SGTELIB::Matrix Gx (  "Gx", static_cast<int>(_n), static_cast<int>(1));
    getModelLagGrad(&Gx, &Mpredict_grad, &Jx, X, Y, sigma, forceOnEmptyTrainingSet);
    return Gx;
}

void SGTELIB::Surrogate_PRS::getModelLagGrad(
  SGTELIB::Matrix * Gx,
  SGTELIB::Matrix * Mpredict_grad,
  SGTELIB::Matrix * Jx,
  const SGTELIB::Matrix & X,
  const SGTELIB::Matrix & Y,
  double sigma,
  bool forceOnEmptyTrainingSet
)
{
    getModelGrad(Gx, Mpredict_grad, X);
    Gx->multiply(sigma);
  
    getModelJacobian(Jx, Mpredict_grad, X);
    SGTELIB::Matrix tmp ("JxtY", Jx->get_nb_cols(), 1);
    SGTELIB::Matrix::inplace_product(tmp, Jx->transpose(), Y);
    tmp.multiply(-1);
    Gx->add(tmp);
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::getModelLagHessian(
  const SGTELIB::Matrix & X,
  const SGTELIB::Matrix & Y,
  double sigma,
  bool forceOnEmptyTrainingSet
)
{
    SGTELIB::Matrix lagHessian("lagHessian",static_cast<int>(_n),static_cast<int>(_n));
    getModelLagHessian(&lagHessian, X, Y, sigma, forceOnEmptyTrainingSet);
    return lagHessian;
}

void SGTELIB::Surrogate_PRS::getModelLagHessian(
  SGTELIB::Matrix * lagHessian,
  const SGTELIB::Matrix & X,
  const SGTELIB::Matrix & Y,
  double sigma,
  bool forceOnEmptyTrainingSet
)
{
    SGTELIB::Matrix outHessian("tmp",static_cast<int>(_n),static_cast<int>(_n));

    if (X.get_nb_cols() != static_cast<int>(_n) || X.get_nb_rows() != 1 )
    {
        throw Exception(__FILE__, __LINE__, "X matrix has wrong dimensions!");
    }

    if (Y.get_nb_rows() != static_cast<int>(_m - 1) || Y.get_nb_cols() != 1 )
    {
        throw Exception(__FILE__, __LINE__, "Multipliers matrix has wrong dimensions!");
    }

    lagHessian->fill(0);
    
    int k = 0;
    for (int j = 0; j < _m ; j++)
    {
        getModelHessian(&outHessian, X, j, forceOnEmptyTrainingSet);
        
        if (_trainingset.get_bbo(j) == SGTELIB::BBO_OBJ)
        {
            outHessian.multiply(sigma) ;
        }
        else
        {
            outHessian.multiply(-Y.get(k, 0));
            k += 1;
        }

        lagHessian->add(outHessian);
    }
}

SGTELIB::Matrix SGTELIB::Surrogate_PRS::compute_multiplier(
    const SGTELIB::Matrix & Grad,
    const SGTELIB::Matrix & Jacobian,
    const double rank_tol)
{
  int ncon = Jacobian.get_nb_rows();
  SGTELIB::Matrix multiplier("multiplier",static_cast<int>(ncon),static_cast<int>(1));
  compute_multiplier(&multiplier, Grad, Jacobian, rank_tol);
  return multiplier;
}

void SGTELIB::Surrogate_PRS::compute_multiplier(
    SGTELIB::Matrix * multiplier,
    const SGTELIB::Matrix & Grad,
    const SGTELIB::Matrix & Jacobian,
    const double rank_tol)
{
    const int ncon = Jacobian.get_nb_rows();
    const int nvar = Jacobian.get_nb_cols();

    if (Grad.get_nb_rows() != nvar || Grad.get_nb_cols() != 1)
    {
        throw Exception(__FILE__, __LINE__, "Grad dimensions are not ok!");
    }
    if (Jacobian.get_nb_cols() != nvar)
    {
        throw Exception(__FILE__, __LINE__, "Jacobian dimensions are not ok!");
    }
    if (Jacobian.has_nan())
    {
        throw Exception(__FILE__, __LINE__, "Jacobian contains NaN");
    }

    // The multipliers are the solutions of the least-square equation
    //  Jacobian multiplier = Grad.
    //  We use the SVD to compute the inverse of Jacobian t Jacobian.

    // init matrices for SVD
    double ** U = new double *[nvar];
    double  * W = new double  [ncon];
    double ** V = new double *[ncon];
    for (int i = 0 ; i < nvar ; ++i ) {
        U[i] = new double[ncon];
    }
    for (int i = 0 ; i < ncon ; ++i ) {
        V[i] = new double[ncon];
    }

    // Perform SVD
    std::string error_msg;
    // Compute the SVD of the transpose Jacobian
    Jacobian.transpose().SVD_decomposition(error_msg , U, W, V, 1000000000);
 
    int rank = 0;
    for (int i = 0; i < ncon; i++)
    {
        if (std::fabs(W[i]) > rank_tol)
        {
            rank++;
        }
        else
        {
            W[i] = 0;
        }
    }

    SGTELIB::Matrix Wm = SGTELIB::Matrix("Wm", ncon, ncon);
    for (int i = 0 ; i < ncon ; i++)
    {
        for (int j = 0 ; j < ncon ; j++)
        {
            const double Wmi = (i == j) && W[i] != 0 ? 1.0 / (W[i] * W[i]) : 0;
            Wm.set(i, j, Wmi);
        }
    }
    SGTELIB::Matrix Vm = SGTELIB::Matrix("Vm", ncon, ncon, V);

    *multiplier = SGTELIB::Matrix::product(Wm, Vm.transpose(), Jacobian, Grad);
    *multiplier = SGTELIB::Matrix::product(Vm, *multiplier);

    for (int i=0 ; i < nvar ; i++)
    {
        delete [] U[i];
    }
    delete [] U;
    for (int j=0 ; j < ncon ; j++){
        delete [] V[j];
    }
    delete [] V;
    delete [] W;
}


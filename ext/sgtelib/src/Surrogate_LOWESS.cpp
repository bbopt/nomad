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

#include "Surrogate_LOWESS.hpp"

//#define SGTELIB_LOWESS_DEV
//#define SGTELIB_DEBUG


const int GAMMA_EXP = 2;

/*----------------------------*/
/*         constructor        */
/*----------------------------*/
SGTELIB::Surrogate_LOWESS::Surrogate_LOWESS ( SGTELIB::TrainingSet & trainingset,
                                        SGTELIB::Surrogate_Parameters param) :
  SGTELIB::Surrogate ( trainingset , param ),
  _q                 ( 0        ),
  _q_old             ( 99999999 ),
  _degree            ( 0        ),
  _H                 ( NULL     ),
  _W                 ( NULL     ),
  _A                 ( NULL     ),
  _HWZ               ( NULL     ),
  _u                 ( NULL     ),
  _old_u             ( NULL     ),
  _old_x             ( NULL     ),
  _x_multiple        ( NULL     ),
  _ZZsi              ("ZZsi",0,0){
  #ifdef SGTELIB_DEBUG
    std::cout << "constructor LOWESS\n";
  #endif
}//


/*--------------------------------------*/
/*           delete matrices            */
/*--------------------------------------*/
void SGTELIB::Surrogate_LOWESS::delete_matrices ( void ) {

  if (_W) delete [] _W; 
  _W = NULL;

  if (_u) delete [] _u; 
  _u = NULL;

  if (_x_multiple) delete [] _x_multiple;
  _x_multiple = NULL;

  #ifdef SGTELIB_LOWESS_DEV
    if (_old_u) delete [] _old_u; 
    _old_u = NULL;

    if (_old_x) delete [] _old_x; 
    _old_x = NULL;
  #endif

  const int p = std::min(_p_old,_p);
  if (_H){
    for (int i=0 ; i<p ; i++) delete [] _H[i];
    delete [] _H;
  }
  _H = NULL;

  const int q = std::min(_q_old,_q);
  if (_A){
    for (int i=0 ; i<q ; i++) delete [] _A[i];
    delete [] _A;
  }
  _A = NULL;

  if (_HWZ){
    for (int i=0 ; i<q ; i++) delete [] _HWZ[i];
    delete [] _HWZ;
  }
  _HWZ = NULL;
}//


/*--------------------------------------*/
/*               destructor             */
/*--------------------------------------*/
SGTELIB::Surrogate_LOWESS::~Surrogate_LOWESS ( void ) {
  delete_matrices();
}//




/*----------------------------*/
/*          display           */
/*----------------------------*/
void SGTELIB::Surrogate_LOWESS::display_private ( std::ostream & out ) const {
  out << "q: " << _q << "\n";
}//

/*--------------------------------------*/
/*               init                   */
/*--------------------------------------*/
bool SGTELIB::Surrogate_LOWESS::init_private ( void ) {
  
  const int pvar = _trainingset.get_pvar(); 

  if (pvar<2) return false;

  // Get the number of basis functions.
  const int n1=_trainingset.get_X_nbdiff1();
  const int n2=_trainingset.get_X_nbdiff2();
  const int q10 = 1+n1;
  const int q15 = 1+n1+n2;
  const int q20 = 1+n1+n2*(n2+1)/2;
  const int degree_max = _param.get_degree();
  if ((pvar>q20) & (degree_max>=2)){
    _q = q20;
    _degree = 20;
  }
  else if ((pvar>q15) & (degree_max>=2)){
    _q = q15;
    _degree = 15;
  }
  else if ((pvar>q10) & (degree_max>=1)){
    _q = q10;
    _degree = 10;
  }
  else{
    _q = 1;
    _degree = 0;
  }

  #ifdef SGTELIB_DEBUG
    std::cout << "_q = " << _q << " (degree=" << double(_degree)/10 << ")\n";
  #endif

  // Init matrices for prediction

  delete_matrices();

  if ( ! _W){
    _W = new double [_p];
  } 
  if ( ! _A){
    _A = new double * [_q];
    for (int j=0 ; j<_q ; j++) _A[j] = new double [_q];
  }
  if ( ! _H){
    _H = new double * [_p];
    for (int j=0 ; j<_p ; j++) _H[j] = new double [_q];
  }
  if ( ! _HWZ){
    _HWZ = new double * [_q];
    for (int j=0 ; j<_q ; j++) _HWZ[j] = new double [_m];
  }
  if ( ! _u){
    _u = new double [_q];
    for (int i=0 ; i<_q ; i++) _u[i] = 0.0;
  }
  if ( ! _x_multiple){
    _x_multiple = new bool [_m];
    for (int j=0 ; j<_m ; j++) _x_multiple[j] = (_trainingset.get_X_nbdiff(j)>1);
  }

  #ifdef SGTELIB_LOWESS_DEV
    if ( ! _old_u){
      _old_u = new double [_q];
      for (int i=0 ; i<_q ; i++) _old_u[i] = 0.0;
    }

    if ( ! _old_x){
      _old_x = new double [_n];
      for (int i=0 ; i<_n ; i++) _old_x[i] = 0.0;
    }
  #endif

  _ZZsi = SGTELIB::Matrix("ZZsi",1,_m);
  #ifdef SGTELIB_DEBUG
    std::cout << "Line " << __LINE__ << "(End of private build)\n";
  #endif

  _q_old = _q;
  // C.Tribes jan 17th, 2017 --- update _p_old to prevent memory leak
  _p_old = _p;

  return true;   
}//


/*--------------------------------------*/
/*               build                  */
/*--------------------------------------*/
bool SGTELIB::Surrogate_LOWESS::build_private ( void ) {
  _ready = true;
  return true;  
}


/*--------------------------------------*/
/*       predict (ZZs only)             */
/*--------------------------------------*/
void SGTELIB::Surrogate_LOWESS::predict_private ( const SGTELIB::Matrix & XXs,
                                                     SGTELIB::Matrix * ZZs ) {

  check_ready(__FILE__,__FUNCTION__,__LINE__);
  const int pxx = XXs.get_nb_rows();
  if (pxx>1){
    for (int i=0 ; i<XXs.get_nb_rows() ; i++){
      #ifdef SGTELIB_DEBUG
        std::cout << "============================================\n";
        std::cout << "Prediction of point " << i << "/" << XXs.get_nb_rows() << "\n";
        std::cout << "============================================\n";
      #endif
      predict_private_single ( XXs.get_row(i) );
      ZZs->set_row( _ZZsi , i );
    }
  }
  else{
    predict_private_single ( XXs );
    *ZZs = _ZZsi;
  }
}//

/*--------------------------------------*/
/*       predict (for one point)        */
/*--------------------------------------*/
void SGTELIB::Surrogate_LOWESS::predict_private_single ( const SGTELIB::Matrix XXs , int i_exclude ) {
  if (XXs.get_nb_rows()!=1){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"predict_private_single : XXs must have only one row." );
  }

  int i,j,j1,j2,k;
  double d;

  #ifdef SGTELIB_DEBUG
    std::cout << "i_exclude = " << i_exclude << "\n";
  #endif

  #ifdef SGTELIB_LOWESS_DEV
    int clock_start;
    clock_start = clock();
  #endif

  // Distance Matrix
  // D : distance between points of XXs and other points of the trainingset
  SGTELIB::Matrix D = _trainingset.get_distances(XXs,get_matrix_Xs(),_param.get_distance_type());

  // Preset
  const std::string preset = _param.get_preset();


  // ==================================
  // GAMMA DISTRIBUTION
  // ==================================
  // Number of points taken into account
  // = p if no point is excluded
  // = p-1 if one point is excluded (ie:i_exclude!=-1).
  const double p_divide = double(_p)-double(i_exclude != -1);
  // Empirical mean & variance of the distances
  SGTELIB::Matrix Distances = D;
  if (GAMMA_EXP==2) Distances=SGTELIB::Matrix::hadamard_square(Distances);
  const double mean = Distances.sum()/p_divide;
  double var  = SGTELIB::Matrix::hadamard_square(Distances+(-mean)).sum()/p_divide;
  if (var == 0) {
      var = EPSILON;
  }
  #ifdef SGTELIB_DEBUG
    std::cout << "mean var = " << mean << " " << var << "\n";
  #endif
  if ( (mean<0) || (var<0) ){
    std::cout << "mean: " << mean << "\n";
    std::cout << "var: " << var << "\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Error on computation of mean and var" );
  }
  // Gamma parameters
  const double gamma_shape = mean*mean/var;
  const double gamma_scale = var/mean;  






  #ifdef SGTELIB_LOWESS_DEV
    // Write in a file the values of dq
    if (i_exclude==-1){
      const SGTELIB::Matrix R = D.rank();
      i = 0;
      while (R.get(i)!=_q-1) i++;
      const double dq_emp = D.get(i);
      const double dq_gam = pow(SGTELIB::gammacdfinv(double(_q)/double(_p),gamma_shape,gamma_scale),1./GAMMA_EXP);

      std::ofstream fileout;
      fileout.open ("data_dq.txt" , std::fstream::out | std::fstream::app);
      fileout << dq_emp << " " << dq_gam << " " << _q << " " << _p << " " << gamma_shape << " " << gamma_scale << "\n";
      fileout.close();
    }
  #endif



  if (preset=="D"){
    // ========================================================
    // Distance only
    for (i=0 ; i<_p ; i++) _W[i] = D.get(i);
  }
  else if (preset=="DEN"){
    // ========================================================
    // Distance, normalized with empirical method
    const SGTELIB::Matrix R = D.rank();
    i = 0;
    while (R.get(i)!=_q-1) i++;
    const double dq = 2.0*D.get(i);
    for (i=0 ; i<_p ; i++) _W[i] = D.get(i)/dq;
  }
  else if (preset=="DGN"){
    // ========================================================
    // Distance, normalized with Gamma method
    const double dq = pow(SGTELIB::gammacdfinv(double(_q)/double(p_divide),gamma_shape,gamma_scale),1./GAMMA_EXP);
    for (i=0 ; i<_p ; i++) _W[i] = D.get(i)/dq;
  }
  else if ( (preset=="RE") || (preset=="REN") ){
    // ========================================================
    // Rank, computed with empirical method
    const SGTELIB::Matrix R = D.rank();
    for (i=0 ; i<_p ; i++) _W[i] = R.get(i);
    if (preset=="REN"){
      for (i=0 ; i<_p ; i++) _W[i] /= (double(_p)-1.0);
    }
  }
  else if ( (preset=="RG") || (preset=="RGN") ){
    // ========================================================
    // Rank, computed with gamma method
    for (i=0 ; i<_p ; i++){
      _W[i] = SGTELIB::gammacdf(pow(D.get(i),GAMMA_EXP),gamma_shape,gamma_scale);
    }
    // DE-Normalization
    if (preset=="RG"){
      for (i=0 ; i<_p ; i++) _W[i] *= (double(_p)-1.0);
    }
  }



  double wsum;
  // For Gamma methods, Handle special case where the variance of the distances is null
  if (var==0){
    for (i=0 ; i<_p ; i++){
      _W[i] = 1.0;
    }
    wsum = _p;
  }
  // Normal case
  else{
    // parameters of the gamma distribution
    const double lambda = _param.get_kernel_coef();
    //std::cout << "lambda : " << lambda << "\n";
    const SGTELIB::kernel_t kt = _param.get_kernel_type();
    // Weights
    wsum = 0;
    for (i=0 ; i<_p ; i++){
      _W[i] = kernel(kt,lambda,_W[i]);
      wsum += _W[i];
    }
  }
 
  // If a point must be excluded from the training points, set its weight to 0.
  if (i_exclude != -1){
    wsum -= _W[i_exclude];
    _W[i_exclude] = 0.0;
    #ifdef SGTELIB_DEBUG
      std::cout << "Exclude training point " << i_exclude << "\n";
    #endif
  }

  #ifdef SGTELIB_DEBUG
    std::cout << "XXs = " << XXs[0] << "\n";
    std::cout << "i_exclude = " << i_exclude << "\n";
    std::cout << "wsum : " << wsum << "\n";
    std::cout << "_W = ";
    for (i=0 ; i<_p ; i++){
      std::cout << _W[i] << " ";
    }
    std::cout << "\n";
  #endif


  if (wsum>EPSILON){
    for (i=0 ; i<_p ; i++){
      _W[i] /= wsum;
    }
  }
  else{
    // If all the weights are negligible, put 1 everywhere
    for (i=0 ; i<_p ; i++){
      _W[i] = 1;
    }
  }

  // Ridge
  double ridge = _param.get_ridge();

  // Build matrices
  const SGTELIB::Matrix & Zs = get_matrix_Zs();
  const SGTELIB::Matrix & Xs = get_matrix_Xs();

  double dx1 = 0;

  // Build H
  for (i=0 ; i<_p ; i++){
    k = 0;
    _H[i][k++] = 1;
    if (_W[i]>EPSILON){
      if (_degree>=10){
        // Linear terms
        for (j=0 ; j<_n && k < _q; j++){
          _H[i][k++] = Xs.get(i,j)-XXs.get(0,j);
        }
      }
      if (_degree>=15){
         // Quad and crossed terms
        for (j1=0 ; j1<_n && k < _q; j1++){
          if (_x_multiple[j1]){
            dx1 = (Xs.get(i,j1)-XXs.get(0,j1));
            _H[i][k++] = dx1*dx1;
            if (_degree>=20){
              for (j2=j1+1 ; j2<_n && k < _q; j2++){
                if (_x_multiple[j2]){
                  _H[i][k++] = dx1*(Xs.get(i,j2)-XXs.get(0,j2));
                }
              }
            }
          }
        } 
      }    
    }
    else{
      // Put row to 0 (except first column at 1)
      for (j=1 ; j<_q ; j++){
        _H[i][j] = 0;
      }
    }
  }

  // Reset A and HWZ
  for (i=0 ; i<_q ; i++){
    for (j=i ; j<_q ; j++){
      _A[i][j] = 0;
    }
    for (j=0 ; j<_m ; j++){
      _HWZ[i][j] = 0;
    }
  }
  #ifdef SGTELIB_DEBUG
    int w_count = 0;
  #endif

  // Build A and HWZ
  double w;
  for (k=0 ; k<_p ; k++){
    w = _W[k];
    if (w>EPSILON){
      #ifdef SGTELIB_DEBUG
        w_count++;
      #endif
      for (i=0 ; i<_q ; i++){
        d = _H[k][i]*w;
        for (j=i ; j<_q ; j++){
          _A[i][j] += d*_H[k][j];
        }
        for (j=0 ; j<_m ; j++){
          _HWZ[i][j] += d*Zs.get(k,j);
        }
      }
    }
  }
  #ifdef SGTELIB_DEBUG
    std::cout << "nb of non-zeros in w: " << w_count << " / " << _p << "\n";
  #endif

  // Symmetry of A 
  for (i=0 ; i<_q ; i++){
    for (j=i+1 ; j<_q ; j++){
      _A[j][i] = _A[i][j];
    }
  }

  // Add ridge term (only for i>1)
  for (i=1 ; i<_q ; i++){
    _A[i][i] += ridge;
  }

  #ifdef SGTELIB_DEBUG
    std::cout << "A = \n";
    for (i=0 ; i<_q ; i++){
      for (j=0 ; j<_q ; j++){
        std::cout << _A[j][i] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "u:     ";
     for (i=0 ; i<_q ; i++) std::cout << _u[i] << " ";
    std::cout << "\n";
  #endif

  #ifdef SGTELIB_LOWESS_DEV
    double time_build;
    time_build = double(clock() - clock_start)/double(CLOCKS_PER_SEC);
  #endif

  //=========================//
  //       RESOLUTION        //
  //=========================//

  const double tol = 1e-12;
  int iter_conj = 0;

  #ifdef SGTELIB_LOWESS_DEV
    // Check if e1 is a solution
    {
      bool e1_sol = true;
      if (_A[0][0]!=1) e1_sol = false;
      else{
        for (i=1 ; i<_q ; i++){
          if (_A[0][i]!=0){
            e1_sol = false;
            break;
          }
        }
      }
      if (e1_sol){
        std::cout << "e1 is a solution!\n";
        //_u[0] = 1;
        //for (i=1 ; i<_q ; i++) _u[i] = 0;
      }
    }
  #endif

  // Initial residual error ||Au-b||_2^2
  double res = 0;
  {
    double Au_i = -1;
    for (i=0 ; i<_q ; i++){
      for (j=0 ; j<_q ; j++){
        Au_i += _A[i][j]*_u[j];
      }
      res += Au_i*Au_i;
      Au_i = 0;
    }
    res = sqrt(res);
  }

  // Choice of the method
  bool USE_CHOL = false;
  bool USE_CONJ = false;
  if (res<1e-4) USE_CONJ = true;
  else          USE_CHOL = true;

  #ifdef SGTELIB_DEBUG
    std::cout << "USE CHOL / CONJ : " << USE_CHOL << " " << USE_CONJ << " ( " << res << " )\n";
  #endif
  #ifdef SGTELIB_LOWESS_DEV
    // Test both methods
    USE_CHOL = true;
    USE_CONJ = true;
    // Time
    double time_chol;
    double time_conj;
    // Save old u
    for (i=0 ; i<_q ; i++) _old_u[i] = _u[i];
  #endif


  if (USE_CONJ){
    double * r = new double [_q];
    double * p = new double [_q];
    double * Ap = new double [_q];

    // Use conjugate
    #ifdef SGTELIB_LOWESS_DEV
      clock_start = clock();
    #endif

    // rr = b-Ax // ==================
    double rr = 0;
    d = +1; // Special initialization of the first value of d
    // to take into account the first term of b (which is 1);
    for (i=0 ; i<_q ; i++){
      for (j=0 ; j<_q ; j++){
        d -= _A[i][j]*_u[j];
      }
      r[i] = d;
      p[i] = d;
      rr += d*d;
      d = 0;
    }
    double rr_old,alpha,pAp;
    
    while ( (rr>tol) && (iter_conj < 100)) {
      // Ap // ===================
      for (i=0 ; i<_q ; i++){
        d = 0;
        for (j=0 ; j<_q ; j++){
          d += _A[i][j]*p[j];
        }
        Ap[i] = d;
      }
      // pAp // ===================
      pAp = 0;
      for (i=0 ; i<_q ; i++) pAp += p[i]*Ap[i];
      // Alpha // =================
      alpha = rr/pAp;
      // u // ======================
      for (i=0 ; i<_q ; i++) _u[i] += alpha*p[i];
      // r // ========================
      for (i=0 ; i<_q ; i++) r[i] -= alpha*Ap[i];
      rr_old = rr;
      rr = 0;
      for (i=0 ; i<_q ; i++) rr += r[i]*r[i];
      // Break ?? // =================
      if (rr <= tol) break;
      // p //=========================
      d = rr/rr_old;
      for (i=0 ; i<_q ; i++){
        p[i] *= d;
        p[i] += r[i];
      }  
      iter_conj++;
    }

    #ifdef SGTELIB_DEBUG
      std::cout << "Conj rr = " << rr << "\n";
      std::cout << "Conj iter = " << iter_conj << "\n";
    #endif
    #ifdef SGTELIB_LOWESS_DEV
      time_conj = double(clock() - clock_start)/double(CLOCKS_PER_SEC);
    #endif
    delete [] r;
    delete [] p;
    delete [] Ap;
  }


  if (USE_CHOL){
    // Use cholesky
    #ifdef SGTELIB_LOWESS_DEV
      clock_start = clock();
    #endif
    SGTELIB::Matrix A("A",_q,_q,_A);
    SGTELIB::Matrix b = SGTELIB::Matrix("b",_q,1);
    b.set(0,0,1.0);
    SGTELIB::Matrix u_mat = SGTELIB::Matrix::cholesky_solve(A,b);
    #ifdef SGTELIB_LOWESS_DEV
      time_chol = double(clock() - clock_start)/double(CLOCKS_PER_SEC);
    #endif
    //std::cout << "Clock (CHOL): " << time_chol << "sec\n";
    for (i=0 ; i<_q ; i++){
      _u[i] = u_mat.get(i,0);
    }
  }

  // Check that u has no nan
  bool has_nan = false;
  for (k=0 ; k<_q ; k++){
    if (!isdef(_u[k])){
      has_nan = true;
      break;
    }
  }

  // _ZZsi is the output of this method, but is not returned as it's an attribut of the class.
  if (has_nan){
    // Find the closest point
    int imin = D.get_min_index();
    _ZZsi = Zs.get_row(imin);
  }
  else{
    // Use u to compute the output
    for (j=0 ; j<_m ; j++){
      d = 0;
      for (k=0 ; k<_q ; k++) d += _u[k]*_HWZ[k][j];
      _ZZsi.set(0,j,d);
    }
  }


  #ifdef SGTLEIB_DEBUG
    std::cout << "_u = [ ";
    for (i=0 ; i<_q ; i++){
      std::cout << _u[i] << " ";
    }
    std::cout << "]\n";
    _ZZsi.display(std::cout);
  #endif

  #ifdef SGTELIB_LOWESS_DEV
    // STATISTICS //
    // Compute norm dx
    double dx = 0;
    for (i=0 ; i<_n ; i++){
      d = _old_x[i] - XXs.get(0,i);
      dx += d*d;
    }
    dx = sqrt(dx);

    // Norm du
    double du = 0;
    for (i=0 ; i<_q ; i++){
      d = _old_u[i] - _u[i];
      du += d*d;
    }
    du = sqrt(du);
    // Display stat
    if (i_exclude==-1){
      std::ofstream myfile;
      const std::string file_name = "LOWESS_times_n"+itos(_n)+".txt";
      myfile.open (file_name.c_str(),std::ios::app);
      myfile << _n << " , " << dx << " " 
             << du << " " << res << " , "
             << time_build << " " << time_chol 
             << " " << time_conj << " , " << iter_conj << "\n";
      myfile.close();
    }

    // Save old x
    for (i=0 ; i<_n ; i++){
      _old_x[i] = XXs.get(0,i);
    }
  #endif


}

/*--------------------------------------*/
/*       compute Zvs                    */
/*--------------------------------------*/
const SGTELIB::Matrix * SGTELIB::Surrogate_LOWESS::get_matrix_Zvs (void){
  check_ready(__FILE__,__FUNCTION__,__LINE__);
  #ifdef SGTELIB_DEBUG
    std::cout << "==========================\n";
    std::cout << "Compute Zvs\n";
    std::cout << "==========================\n";
  #endif
  if ( ! _Zvs){
    _Zvs = new SGTELIB::Matrix("Zvs",_p,_m);
    for (int i=0 ; i<_p ; i++){
      predict_private_single( get_matrix_Xs().get_row(i) , i);
      _Zvs->set_row( _ZZsi ,i);
    }
  }
  #ifdef SGTELIB_DEBUG
    std::cout << "==========================\n";
    std::cout << "END Compute Zvs\n";
    std::cout << "==========================\n";
  #endif
  return _Zvs;
}//





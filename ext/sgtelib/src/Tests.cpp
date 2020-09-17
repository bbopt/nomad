/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/
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

#include "Tests.hpp"
#include "Surrogate_Ensemble.hpp"

void SGTELIB::sand_box (void){

  std::cout << "End of sandbox\nEXIT\n";
}//


void SGTELIB::test_LOWESS_times ( void ){

  std::cout << "====================================================================\n";
  std::cout << "START LOWESS TIMES\n";

  const int m = 1;
  const int pp = 6;
  int n = 0;

  for (int in=0 ; in<1 ; in++){
    if (in==0) n = 16;
    if (in==1) n = 8;
    if (in==2) n = 4;
    if (in==3) n = 2;
    if (in==4) n = 1;

    std::cout <<"--------------------\n";

    int NEXP = 20;

    int p = (n+1)*(n+2);
    //if (p<20) p=20;
    double dx;

    for (int id=10 ; id<11 ; id++){
      dx = pow(10.0,-id);
      std::cout << "n=" << n << ", dx=" << dx << "\n";

      SGTELIB::Matrix DX ("DX",1,n);
      SGTELIB::Matrix X ("X",p,n);
      SGTELIB::Matrix XX ("XX",pp,n);
      SGTELIB::Matrix Z ("Z",p,m);
      SGTELIB::Matrix * ZZ  = new SGTELIB::Matrix("ZZ",pp,m);

      for (int i = 0 ; i<NEXP ; i++){
        std::cout <<  "n=" << n << ", dx=" << dx << ", nexp=" << i << "\n";
        // Build data  
        X.set_random(-5,+5,false);
        Z.set_random(-5,+5,false);
        SGTELIB::TrainingSet C0(X,Z);
        SGTELIB::Surrogate * S0;
        // Build model
        S0 = SGTELIB::Surrogate_Factory(C0,"TYPE LOWESS RIDGE 0.001");
        S0->build();
        // Build XX
        DX.set_random(-1,+1,false);
        XX.set_row(DX,0);
        for (int j=1 ; j<pp ; j++){
          DX.set_random(-1,+1,false);
          DX = DX*(dx/DX.norm());
          DX = DX + XX.get_row(j-1);
          XX.set_row(DX,j);
        }
        // Prediction
        S0->predict(XX,ZZ);
      }// end loop i

      delete ZZ;

    }// end loop id
  }// end loop in

  std::cout << "FINISH LOWESS TIMES\n";
  std::cout << "====================================================================\n";

}//




void SGTELIB::test_many_models (        const std::string & output_file ,
                                        const SGTELIB::Matrix & X0 ,
                                        const SGTELIB::Matrix & Z0 ){
  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_many_models\n";

  // Data information
  const int m = Z0.get_nb_cols();

  // TrainingSet creation
  SGTELIB::TrainingSet C0(X0,Z0);
  #ifdef SGTELIB_DEBUG
    C0.display(std::cout);
  #endif
  // Init model
  SGTELIB::Surrogate * S0;

  // LOAD MODEL LIST in model_list_file
  std::vector<std::string> model_list;
  model_list.push_back("TYPE PRS DEGREE 2"); 
  model_list.push_back("TYPE PRS DEGREE OPTIM");
  model_list.push_back("TYPE KS KERNEL_COEF OPTIM");
  model_list.push_back("TYPE LOWESS KERNEL_COEF OPTIM DEGREE 1"); 
  model_list.push_back("TYPE LOWESS KERNEL_COEF OPTIM DEGREE OPTIM KERNEL_TYPE OPTIM");
  model_list.push_back("TYPE RBF PRESET I");


  // BUILD THE METRIC LIST 
  std::vector<SGTELIB::metric_t> metric_list;

  metric_list.push_back( METRIC_RMSE   );
  metric_list.push_back( METRIC_RMSECV  );
  metric_list.push_back( METRIC_OE   );
  metric_list.push_back( METRIC_OECV   );

  // Open an output file
  std::ofstream out;
  out.open (output_file.c_str());

  out << "list_metrics " ;
  for (std::vector<SGTELIB::metric_t>::iterator it_metric = metric_list.begin() ; it_metric != metric_list.end(); ++it_metric){
    out << SGTELIB::metric_type_to_str(*it_metric) << " ";
  }
  out << "\n\n";

  bool ready;
  // Loop on the models
  for (std::vector<std::string>::iterator it_model = model_list.begin() ; it_model != model_list.end(); ++it_model){

    std::cout << "Build " << *it_model << "\n";
    out << *it_model << "\n";

    // Create model
    S0 = SGTELIB::Surrogate_Factory(C0,*it_model);
    ready = S0->build();

    // Loop on the outputs    
    for (int j=0 ; j<m ; j++){
      out << "output " << j << "\n";
      out << "metrics ";
      // Loop on the metrics
      for (std::vector<SGTELIB::metric_t>::iterator it_metric = metric_list.begin() ; it_metric != metric_list.end(); ++it_metric){
        if (ready)
          out << S0->get_metric (*it_metric,j) << " ";
        else
          out << "-999 ";
      }
      out << "\n";
    }
    out << "\n";
    SGTELIB::surrogate_delete(S0);
  }
  out.close();  

}//


/*----------------------------------------------------*/
/*       TEST QUICK                                   */
/*----------------------------------------------------*/
std::string SGTELIB::test_quick (const std::string & s , const SGTELIB::Matrix & X0){


  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_quick\n";
  std::cout << s << "\n";

  // CONSTRUCT DATA
  SGTELIB::Matrix Z0 = test_functions(X0);
  const int m = Z0.get_nb_cols();
  bool ready;

  // CONSTRUCT REFERENCE MODEL
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  ready = S0->build();

  if ( !  ready){
    surrogate_delete(S0);
    std::cout << "test_quick: model ("+s+") is not ready\n";
    return       "test_quick: model ("+s+") is not ready\n";
  }
    
  // Compute the metrics
  double * emax = new double [m];
  double * rmse = new double [m];
  double * rmsecv= new double [m];
  double * oe = new double [m];
  double * oecv = new double [m];
  double * linv = new double [m];
  for (int j=0 ; j<m ; j++){
    emax[j] = S0->get_metric(SGTELIB::METRIC_EMAX,j);
    rmse[j] = S0->get_metric(SGTELIB::METRIC_RMSE,j);
    rmsecv[j]=S0->get_metric(SGTELIB::METRIC_RMSECV,j);
    oe[j]   = S0->get_metric(SGTELIB::METRIC_OE,j);
    oecv[j] = S0->get_metric(SGTELIB::METRIC_OECV,j);
    linv[j] = S0->get_metric(SGTELIB::METRIC_LINV,j);
  }

  // Display
  std::ostringstream oss;
  oss << "test_quick\n";
  oss << "Surrogate string: " << s << "\n";
  oss << "  j|          emax|          rmse|        rmsecv|            oe|          oecv|          linv|\n";
  oss << "---|--------------|--------------|--------------|--------------|--------------|--------------|\n";
  for (int j=0 ; j<m ; j++){
    oss.width(3);
    oss << j << "|";
    oss.width(14);
    oss << emax[j] << "|";
    oss.width(14);
    oss << rmse[j] << "|";
    oss.width(14);
    oss << rmsecv[j] << "|";
    oss.width(14);
    oss << oe[j] << "|";
    oss.width(14);
    oss << oecv[j] << "|";
    oss.width(14);
    oss << linv[j] << "|\n";
  }
  oss << "---|--------------|--------------|--------------|--------------|--------------|--------------|\n";
  std::cout << oss.str();

  for (int j=0 ; j<m ; j++){
    if ( (isnan(emax[j])) || (isnan(rmsecv[j])) || (isnan(oe[j])) || (isnan(oecv[j])) || (isnan(linv[j])) ){
      std::cout << "There is some nan\n";
      std::cout << "EXIT!\n"; 
      exit(0);
    }
    if ( (isinf(emax[j])) || (isinf(rmse[j])) || ( isinf(rmsecv[j])) || (isinf(oe[j])) || (isinf(oecv[j])) || (isinf(linv[j])) ){
      std::cout << "There is some inf\n";
      std::cout << "EXIT!\n"; 
      exit(0);
    }
  }
  delete [] emax;
  delete [] rmse;
  delete [] rmsecv;
  delete [] oe;
  delete [] oecv;
  delete [] linv;

  SGTELIB::surrogate_delete(S0);

  return oss.str();
}//


/*----------------------------------------------------*/
/*       TEST pxx                                    */
/*----------------------------------------------------*/
std::string SGTELIB::test_pxx (const std::string & s , const SGTELIB::Matrix & X0){

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_pxx\n";
  std::cout << s << "\n";

  SGTELIB::Matrix Z0 = test_functions(X0);

  // CONSTRUCT DATA
  const int n = X0.get_nb_cols();
  const int p = X0.get_nb_rows();
  const int m = Z0.get_nb_cols();

  // CONSTRUCT REFERENCE MODEL
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  bool ready;
  ready = S0->build();

  if ( !  ready){
    surrogate_delete(S0);
    std::cout <<  "test_pxx: model ("+s+") is not ready\n";
    return        "test_pxx: model ("+s+") is not ready\n";
  }

  // Init
  SGTELIB::Matrix XX;
  int pxx;

  // Validation predictions
  SGTELIB::Matrix * ZZ  = NULL;
  SGTELIB::Matrix * std = NULL;
  SGTELIB::Matrix * ei  = NULL;
  SGTELIB::Matrix * cdf = NULL;

  // Reference prediction 
  SGTELIB::Matrix * ZZ0  = NULL;
  SGTELIB::Matrix * std0 = NULL;
  SGTELIB::Matrix * ei0  = NULL;
  SGTELIB::Matrix * cdf0 = NULL;

  for (int i=0 ; i<4 ; i++){ 
    switch (i){
      case 0:
        pxx = 1;
        break;
      case 1:
        pxx = 2;
        break;
      case 2:
        pxx = p;
        break;
      case 3:
        pxx = 2*p;
        break;
      default:
        std::cout << "ERROR i = " << i << "\n";
        exit(0);
    }

    // TESTING POINT(S)
    XX = SGTELIB::Matrix("XX",pxx,n);
    XX.set_random(-10,+10,false);

    // Reference output matrices 
    ZZ0  = new SGTELIB::Matrix("ZZ0" ,pxx,m);
    std0 = new SGTELIB::Matrix("std0",pxx,m);
    ei0  = new SGTELIB::Matrix("ei0" ,pxx,m);
    cdf0 = new SGTELIB::Matrix("cdf0",pxx,m);

    S0->predict(XX,ZZ0,std0,ei0,cdf0);

    for (int k=0 ; k<7 ; k++){

      // Output matrices
      ZZ  = new SGTELIB::Matrix("ZZ" ,pxx,m);
      std = new SGTELIB::Matrix("std",pxx,m);
      ei  = new SGTELIB::Matrix("ei" ,pxx,m);
      cdf = new SGTELIB::Matrix("cdf",pxx,m);

      switch (k){
        case 0:
          S0->predict(XX,ZZ);
          check_matrix_diff(ZZ0 ,ZZ);
          break;
        case 1:
          S0->predict(XX,ZZ,std ,NULL,NULL);
          check_matrix_diff(ZZ0 ,ZZ);
          check_matrix_diff(std0,std);
          break;
        case 2:
          S0->predict(XX,ZZ,NULL,ei  ,NULL);
          check_matrix_diff(ZZ0 ,ZZ);
          check_matrix_diff(ei0,ei);
          break;
        case 3:
          S0->predict(XX,ZZ,NULL,NULL,cdf );
          check_matrix_diff(ZZ0 ,ZZ);
          check_matrix_diff(cdf0,cdf);
          break;
        case 4:
          S0->predict(XX,ZZ,NULL,ei  ,cdf );
          check_matrix_diff(ZZ0 ,ZZ);
          check_matrix_diff(ei0 ,ei);
          check_matrix_diff(cdf0,cdf);
          break;
        case 5:
          S0->predict(XX,ZZ,std ,NULL,cdf );
          check_matrix_diff(ZZ0 ,ZZ);
          check_matrix_diff(std0,std);
          check_matrix_diff(cdf0,cdf);
          break;
        case 6:
          S0->predict(XX,ZZ,std ,ei  ,NULL);
          check_matrix_diff(ZZ0 ,ZZ);
          check_matrix_diff(std0,std);
          check_matrix_diff(ei0 ,ei);
          break;
        default:
          std::cout << "ERROR k = " << k << "\n";
          exit(0);

      }// end switch

      delete ZZ;
      delete std;
      delete ei;
      delete cdf;

    }// end loop k

    delete ZZ0;
    delete std0;
    delete ei0;
    delete cdf0;

  }//end loop i

  SGTELIB::surrogate_delete(S0);

  return "test_pxx ok\n";
}//



/*----------------------------------------------------*/
/*       TEST update                                  */
/*----------------------------------------------------*/
std::string SGTELIB::test_update (const std::string & s , const SGTELIB::Matrix & X0){

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_update\n";
  std::cout << s << "\n";

  // CONSTRUCT DATA
  const int p = X0.get_nb_rows();
  const int n = X0.get_nb_cols();
  SGTELIB::Matrix Z0 = test_functions(X0);
  const int m = Z0.get_nb_cols();

  // CONSTRUCT MODEL
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  bool ready;
  ready = S0->build();

  if ( !  ready){
    surrogate_delete(S0);
    std::cout << "test_update: model ("+s+") is not ready\n";
    return       "test_update: model ("+s+") is not ready\n";
  }

  // Testing set
  const int pxx = 3;
  SGTELIB::Matrix XX("XX",pxx,n);  
  XX.set_random(-5,+5,false);

  // Reference prediction 
  SGTELIB::Matrix * ZZ0  = new SGTELIB::Matrix("ZZ0" ,pxx,m);
  SGTELIB::Matrix * std0 = new SGTELIB::Matrix("std0",pxx,m);
  SGTELIB::Matrix * ei0  = new SGTELIB::Matrix("ei0" ,pxx,m);
  SGTELIB::Matrix * cdf0 = new SGTELIB::Matrix("cdf0",pxx,m);
  S0->predict(XX,ZZ0,std0,ei0,cdf0);


  // CONSTRUCT MODEL
  SGTELIB::TrainingSet C1(X0.get_row(0),Z0.get_row(0));
  SGTELIB::Surrogate * S1;
  S1 = SGTELIB::Surrogate_Factory(C1,s);
  S1->build();

  for (int i=1 ; i<p ; i++){
    C1.add_points(X0.get_row(i),Z0.get_row(i));
    S1->build();
  }

  // Validation predictions
  SGTELIB::Matrix * ZZ1  = new SGTELIB::Matrix("ZZ1" ,pxx,m);
  SGTELIB::Matrix * std1 = new SGTELIB::Matrix("std1",pxx,m);
  SGTELIB::Matrix * ei1  = new SGTELIB::Matrix("ei1" ,pxx,m);
  SGTELIB::Matrix * cdf1 = new SGTELIB::Matrix("cdf1",pxx,m);
  S1->predict(XX,ZZ1,std1,ei1,cdf1);

  // Check consistency
  check_matrix_diff(ZZ0 ,ZZ1 );
  check_matrix_diff(std0,std1);
  check_matrix_diff(ei0 ,ei1 );
  check_matrix_diff(cdf0,cdf1);

  // Free space
  SGTELIB::surrogate_delete(S0);
  SGTELIB::surrogate_delete(S1);

  delete ZZ0;
  delete std0;
  delete ei0;
  delete cdf0;

  delete ZZ1;
  delete std1;
  delete ei1;
  delete cdf1;

  return "test_update ok\n";

}//



/*----------------------------------------------------*/
/*       TEST SINGULAR DATA                           */
/*----------------------------------------------------*/
std::string SGTELIB::test_singular_data (const std::string & s ) {

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_singular_data\n";
  std::cout << s << "\n";

  // CONSTRUCT DATA
  const int n = 3;
  const int p = 10;
  SGTELIB::Matrix X0 ("X0",p,n);
  X0.set_random(0,10,false);
  SGTELIB::Matrix Z0 = test_functions(X0);
  const int m = Z0.get_nb_cols();
  Z0.set_name("Z0");

  // Column 0 of X0 is constant (= 0.0);
  X0.set_col(0.0,0);
  // Column 0 of Z0 is constant (= 0.0);
  Z0.set_col(0.0,0);
  // Column 1 of Z0 has some nan
  // Z0.set(2,1,0.0/0.0);
  // Z0.set(5,1,0.0/0.0);
  Z0.set(2,1,SGTELIB::NaN);
  Z0.set(5,1,SGTELIB::NaN);
  // Column 2 of Z0 has some inf
  Z0.set(4,2,SGTELIB::INF);
  Z0.set(7,2,SGTELIB::INF);
  // Column 3 of Z0 has some nan and some inf
  Z0.set(5,3,SGTELIB::INF);
  Z0.set(8,3,SGTELIB::NaN);
  // Z0.set(8,3,0.0/0.0);
  
  #ifdef SGTELIB_DEBUG
    X0.display(std::cout);  
    Z0.display(std::cout);
  #endif

  // CONSTRUCT MODEL
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  bool ready = S0->build();

  if ( !  ready){
    surrogate_delete(S0);
    std::cout << "test_singular_data: model ("+s+") is not ready\n";
    return       "test_singular_data: model ("+s+") is not ready\n";
  }

  // Get the rmse and rmsecv
  double * rmse   = new double [m];
  double * rmsecv = new double [m];
  for (int j=0 ; j<m ; j++){
    rmse[j]   = S0->get_metric(SGTELIB::METRIC_RMSE,j);
    rmsecv[j] = S0->get_metric(SGTELIB::METRIC_RMSECV,j);
  }

  // Display
  std::ostringstream oss;
  oss << "test_singular_data\n";
  oss << "Surrogate string: " << s << "\n";
  oss << "  j|          rmse|        rmsecv|\n";
  oss << "---|--------------|--------------|\n";
  for (int j=0 ; j<m ; j++){
    oss.width(3);
    oss << j << "|";
    oss.width(14);
    oss << rmse[j] << "|";
    oss.width(14);
    oss << rmsecv[j] << "|\n";
  }
  oss << "---|--------------|--------------|\n";

  for (int j=0 ; j<m ; j++){
    if ( ( !  isdef(rmse[j])) || ( !  isdef(rmse[j])) ){
      std::cout << "There are some nan !";
      C0.get_matrix_Xs().display(std::cout);
      exit(0);
    }
  }

  std::cout << oss.str();
  SGTELIB::surrogate_delete(S0);
  delete [] rmse;
  delete [] rmsecv;
  return oss.str();

}//




/*----------------------------------------------------*/
/*       TEST scale                                   */
/*----------------------------------------------------*/
std::string SGTELIB::test_scale (const std::string & s , const SGTELIB::Matrix & X0){

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_scale\n";
  std::cout << s << "\n";

  // CONSTRUCT DATA
  const int p = X0.get_nb_rows();
  const int n = X0.get_nb_cols();
  SGTELIB::Matrix Z0 = test_functions(X0);
  const int m = Z0.get_nb_cols();
  bool ready;

  // CONSTRUCT MODEL
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  ready = S0->build();

  if ( !  ready){
    surrogate_delete(S0);
    std::cout << "test_scale: model ("+s+") is not ready\n";
    return       "test_scale: model ("+s+") is not ready\n";
  }

  // Testing set
  const int pxx = 3;
  SGTELIB::Matrix XX0("XX0",pxx,n);  
  XX0.set_random(-5,+5,false);

  // Reference prediction 
  SGTELIB::Matrix * ZZ0  = new SGTELIB::Matrix("ZZ0" ,pxx,m);
  SGTELIB::Matrix * std0 = new SGTELIB::Matrix("std0",pxx,m);
  SGTELIB::Matrix * ei0  = new SGTELIB::Matrix("ei0" ,pxx,m);
  SGTELIB::Matrix * cdf0 = new SGTELIB::Matrix("cdf0",pxx,m);
  S0->predict(XX0,ZZ0,std0,ei0,cdf0);

  // Build scaling values
  double * ax = new double [n];
  double * bx = new double [n];
  for (int i=0 ; i<n ; i++){
    ax[i] = double(i+2);       // arbitrary values
    bx[i] = 1.5+1/double(i+2); // arbitrary values
  }
  double * az = new double [m];
  double * bz = new double [m];
  for (int i=0 ; i<m ; i++){
    az[i] = 0.0;             // No offset on the constraints (we don't want to change the feasibility)
    bz[i] = 3.0+double(i+2); // arbitrary values
  }

  // CONSTRUCT SCALED DATA
  SGTELIB::Matrix X1 = X0;
  SGTELIB::Matrix Z1 = Z0;  
  for (int i=0 ; i<p ; i++){ 
    for (int j=0 ; j<n ; j++){ 
       X1.set(i,j, ax[j] + bx[j] *  X0.get(i,j) );
    }
    for (int j=0 ; j<m ; j++){ 
       Z1.set(i,j, az[j] + bz[j] *  Z0.get(i,j) );
    }
  }  

  // Build new model
  SGTELIB::TrainingSet C1(X1,Z1);
  SGTELIB::Surrogate * S1;
  S1 = SGTELIB::Surrogate_Factory(C1,s);
  S1->build();

  // Verif prediction 
  SGTELIB::Matrix * ZZ1  = new SGTELIB::Matrix("ZZ1" ,pxx,m);
  SGTELIB::Matrix * std1 = new SGTELIB::Matrix("std1",pxx,m);
  SGTELIB::Matrix * ei1  = new SGTELIB::Matrix("ei1" ,pxx,m);
  SGTELIB::Matrix * cdf1 = new SGTELIB::Matrix("cdf1",pxx,m);

  // CONSTRUCT SCALED PREDICTION POINTS
  SGTELIB::Matrix XX1= XX0;
  for (int i=0 ; i<pxx ; i++){ 
    for (int j=0 ; j<n ; j++){ 
      XX1.set(i,j, ax[j] + bx[j] * XX0.get(i,j) );
    }
  } 
  S1->predict(XX1,ZZ1,std1,ei1,cdf1);

  // Renormalize
  for (int i=0 ; i<pxx ; i++){ 
    for (int j=0 ; j<m ; j++){ 
       ZZ1->set(i,j, (  ZZ1->get(i,j)-az[j] ) / bz[j] );
      std1->set(i,j, ( std1->get(i,j)       ) / bz[j] );
       ei1->set(i,j, (  ei1->get(i,j)       ) / bz[j] );
    }
  }  

  // Check consistency
  std::cout << s << "\n";
  std::cout << "Check ZZ\n";
  check_matrix_diff(ZZ0 ,ZZ1 );
  std::cout << "Check std\n";
  check_matrix_diff(std0,std1);
  std::cout << "Check ei\n";
  check_matrix_diff(ei0 ,ei1 );
  std::cout << "Check cdf\n";
  check_matrix_diff(cdf0,cdf1);

  // Free space
  SGTELIB::surrogate_delete(S0);
  SGTELIB::surrogate_delete(S1);

  delete ZZ0;
  delete std0;
  delete ei0;
  delete cdf0;

  delete ZZ1;
  delete std1;
  delete ei1;
  delete cdf1;

  delete [] ax;
  delete [] bx;
  delete [] az;
  delete [] bz;

  std::cout << "test_scale OK for model " << s << "\n";
  return "test_scale ok\n";
}//




/*----------------------------------------------------*/
/*       TEST dimension                               */
/*----------------------------------------------------*/
std::string SGTELIB::test_dimension (const std::string & s ){

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_dimension\n";
  std::cout << s << "\n";

  // INIT DATA
  int p =0,n=0,m=0,pxx=1;

  bool ready;

  // Training Set and Surrogate
  SGTELIB::TrainingSet * C0=NULL;
  SGTELIB::Surrogate * S0=NULL;

  // Predictions
  SGTELIB::Matrix X0,Z0;
  SGTELIB::Matrix XX, ZZ, STD, EI, CDF;

  // BUILD THE METRIC LIST 
  std::vector<SGTELIB::metric_t> metric_list;
  metric_list.push_back( METRIC_RMSECV );
  metric_list.push_back( METRIC_EMAX   );
  metric_list.push_back( METRIC_EMAXCV );
  metric_list.push_back( METRIC_RMSE   );
  metric_list.push_back( METRIC_RMSECV );
  metric_list.push_back( METRIC_OE     );
  metric_list.push_back( METRIC_OECV   );
  metric_list.push_back( METRIC_LINV   );

  int i_case;
  const int i_case_max = 5;
  for ( i_case = 0 ; i_case < i_case_max ; i_case++ ){
    std::cout << "------------------------------------------------------\n";
    std::cout << "i_case = " << i_case ;
    if (i_case==0){
      std::cout << " (small m,n,p,pxx) ";
      m = 1;
      n = 2;
      p = 3;
      pxx = 4;
    }
    else if (i_case==1){
      std::cout << " (big m) ";
      m = 30;
      n = 1;
      p = 10;
      pxx = 2;
    }
    else if (i_case==2){
      std::cout << " (big n) ";
      m = 1;
      n = 30;
      p = 10;
      pxx = 2;
    }
    else if (i_case==3){
      std::cout << " (big p) ";
      m = 1;
      n = 2;
      p = 30; 
      pxx = 3;
    }
    else if (i_case==4){
      std::cout << " (big pxx) ";
      m = 1;
      n = 2;
      p = 30; 
      pxx = 50;
    }

    std::cout << "m,n,p,pxx = " << m << " " << n << " " << p << " " << pxx << "\n";

    X0 = SGTELIB::Matrix("X0",p,n);
    Z0 = SGTELIB::Matrix("X0",p,m);
    X0.set_random(0,10,false);
    Z0.set_random(0,10,false);

    C0 = new SGTELIB::TrainingSet (X0,Z0);
    S0 = SGTELIB::Surrogate_Factory(*C0,s);
    ready = S0->build();

    // Loop on the metrics
    std::vector<SGTELIB::metric_t>::iterator it;
    int j;
    double v = 0;
    if (ready){
      for (it = metric_list.begin() ; it != metric_list.end(); ++it){
        std::cout << "Metric " << SGTELIB::metric_type_to_str(*it) << "\n"; 
        for (j=0 ; j<m ; j++){
          v = S0->get_metric (*it,j);
          std::cout << "v = " << v << "\n";
        }
      }

      XX  = SGTELIB::Matrix("XX",pxx,n);
      XX.set_random(0,10,false);
      ZZ  = SGTELIB::Matrix("ZZ" ,pxx,m);
      STD = SGTELIB::Matrix("std",pxx,m);
      EI  = SGTELIB::Matrix("ei" ,pxx,m);
      CDF = SGTELIB::Matrix("cdf",pxx,m);

      std::cout << "m,n,p,pxx = " << m << " " << n << " " << p << " " << pxx << "\n";
      std::cout << "predict(XX,&ZZ)...\n";
      S0->predict(XX,&ZZ);  

      std::cout << "m,n,p,pxx = " << m << " " << n << " " << p << " " << pxx << "\n";
      std::cout << "predict(XX,&ZZ,&STD,&EI,&CDF)...\n";
      S0->predict(XX,&ZZ,&STD,&EI,&CDF);
      std::cout << "Finish!\n";
    }
    else{
      std::cout << "Not Ready\n" ;
    }
  }

  // Free space
  SGTELIB::surrogate_delete(S0);
  delete C0;

  std::cout << "test_dimension OK for model " << s << "\n";
  std::cout << "======================================================\n";
  return "test_dimension ok\n";
}//


/*----------------------------------------------------*/
/*       TEST rmse                                     */
/*----------------------------------------------------*/
std::string SGTELIB::test_rmse (const std::string & s , const SGTELIB::Matrix & X0){

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_rmse\n";
  std::cout << s << "\n";

  // CONSTRUCT DATA
  const int p = X0.get_nb_rows();
  SGTELIB::Matrix Z0 = test_functions(X0);
  const int m = Z0.get_nb_cols();
  bool ready;

  // CONSTRUCT MODEL
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  ready = S0->build();

  if ( !  ready){
    surrogate_delete(S0);
    std::cout << "test_rmse: model ("+s+") is not ready\n";
    return       "test_rmse: model ("+s+") is not ready\n";
  }

  // Get the rmse
  double * rmse = new double [m];
  for (int j=0 ; j<m ; j++){
    rmse[j] = S0->get_metric(SGTELIB::METRIC_RMSE,j);
  }

  // GET THE PREDICTION ON THE TP OUTPUT (matrix Zh)
  const SGTELIB::Matrix Zh = S0->get_matrix_Zh();

  // Recompute the prediction on the TP 
  SGTELIB::Matrix Zh_verif("Zh_verif",p,m);
  SGTELIB::Matrix z("z",1,m);
  for (int i=0 ; i<p ; i++){ 
    S0->predict(X0.get_row(i),&z);
    Zh_verif.set_row(z,i);
  }

  // Recompute the rmse
  double * rmse_verif = new double [m];
  double e;
  for (int j=0 ; j<m ; j++){
    e = 0;
    for (int i=0 ; i<p ; i++){
      e += pow(Z0.get(i,j)-Zh_verif.get(i,j),2);
    }
    rmse_verif[j] = sqrt(e/p);
  }

  // Display
  std::ostringstream oss;
  oss << "test_rmse\n";
  oss << "Surrogate string: " << s << "\n";
  oss << "  j|          rmse|    rmse_verif|          diff|\n";
  oss << "---|--------------|--------------|--------------|\n";
  for (int j=0 ; j<m ; j++){
    oss.width(3);
    oss << j << "|";
    oss.width(14);
    oss << rmse[j] << "|";
    oss.width(14);
    oss << rmse_verif[j] << "|";
    oss.width(14);
    oss << rmse[j]-rmse_verif[j] << "|\n";
    if ( fabs(rmse[j]-rmse_verif[j])>1e-6 ){
      oss << "Error! Diff is too big!\n";
    }
  }
  oss << "---|--------------|--------------|--------------|\n";

  std::cout << oss.str();
  SGTELIB::surrogate_delete(S0);
  delete [] rmse;
  delete [] rmse_verif;


  return oss.str();
}//


/*----------------------------------------------------*/
/*       TEST RMSECV                                  */
/*----------------------------------------------------*/
std::string SGTELIB::test_rmsecv (const std::string & s , const SGTELIB::Matrix & X0){

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_rmsecv\n";
  std::cout << s << "\n";

  // CONSTRUCT DATA
  const int p = X0.get_nb_rows();
  const int n = X0.get_nb_cols();
  const SGTELIB::model_t mt = SGTELIB::Surrogate_Parameters::read_model_type(s);
  SGTELIB::Matrix Z0 = test_functions(X0);
  const int m = Z0.get_nb_cols();

  #ifdef SGTELIB_DEBUG
    X0.display(std::cout);
    Z0.display(std::cout);
  #endif

  double dmean0, dmeanv, kc0, kcv , xsa0, xsav;

  // CONSTRUCT REFERENCE MODEL
  bool ready;
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  ready = S0->build();
  dmean0 = C0.get_Ds_mean();
  xsa0 = C0.get_X_scaling_a(0);


  // Get original kernel coefficient
  kc0 = 0.0;
  S0->get_param().get_kernel_coef();

  // Check ready
  if ( !  ready){
    surrogate_delete(S0);
    std::cout << "test_rmsecv: model ("+s+") is not ready\n";
    return       "test_rmsecv: model ("+s+") is not ready\n";
  }

  // Get the RMSECV metric
  double * rmsecv = new double [m];
  double * rmse  = new double [m];
  for (int j=0 ; j<m ; j++){
    rmsecv[j] = S0->get_metric(SGTELIB::METRIC_RMSECV,j);
    rmse[j] = S0->get_metric(SGTELIB::METRIC_RMSE,j);
  }
  // Delete the original model
  SGTELIB::surrogate_delete(S0);

  // INIT THE CROSS VALIDATION MATRICES
  SGTELIB::Matrix X0i ("X0i",p-1,n);   
  SGTELIB::Matrix Z0i ("Z0i",p-1,m); 
  for (int i=1 ; i<p ; i++){
    // Skip the first line
    X0i.set_row(X0.get_row(i),i-1);
    Z0i.set_row(Z0.get_row(i),i-1);
  }

  // Cross-validation models and trainingset
  SGTELIB::Surrogate * Sv;
  SGTELIB::TrainingSet * Cv;

  // Output for one point
  SGTELIB::Matrix Zvi ("Zvi",1,m); 

  // Cross validation output for all the training points:
  // (Supposed to be identical to _Zv of class Surrogate)
  SGTELIB::Matrix Zv_verif  ("Zv_verif",p,m);

  // BUILD THE CV MODELS
  for (int i=0 ; i<p ; i++){
    // Build the trainingset (without point i);
    Cv = new SGTELIB::TrainingSet(X0i,Z0i);
    Cv->build();
    dmeanv = Cv->get_Ds_mean();
    xsav = Cv->get_X_scaling_a(0);
    // Init surrogate
    Sv = SGTELIB::Surrogate_Factory(*Cv,s);
    // Correct ks for RBF and KS
    kcv = kc0*(dmeanv/dmean0)*(xsa0/xsav);
    if (mt==SGTELIB::RBF) static_cast<SGTELIB::Surrogate_RBF*>(Sv)->set_kernel_coef(kcv);
    if (mt==SGTELIB::KS)  static_cast<SGTELIB::Surrogate_KS* >(Sv)->set_kernel_coef(kcv);
    // Build
    Sv->build();
    // Do the prediction on point i
    Sv->predict(X0.get_row(i),&Zvi);
    Zv_verif.set_row(Zvi,i);
    // update the matrices so that they lack the point i+1
    if (i<p-1){
      X0i.set_row(X0.get_row(i),i);
      Z0i.set_row(Z0.get_row(i),i);
    }
    delete Cv;
    surrogate_delete(Sv);
  }


  // Re-Compute the rmsecv
  double * rmsecv_verif = new double [m];  
  double e;
  for (int j=0 ; j<m ; j++){
    e = 0;
    for (int i=0 ; i<p ; i++){
      e += pow(Z0.get(i,j)-Zv_verif.get(i,j),2);
    }
    rmsecv_verif[j] = sqrt(e/p);
  }

  // Display
  double d;
  std::ostringstream oss; 
  oss << "Surrogate string: " << s << "\n";
  oss << "  j|          rmse||        rmsecv|  rmsecv_verif|      rel diff|\n";
  oss << "---|--------------||--------------|--------------|--------------|\n";
  for (int j=0 ; j<m ; j++){
    oss.width(3);
    oss << j << "|";
    oss.width(14);
    oss << rmse[j] << "||";
    oss.width(14);
    oss << rmsecv[j] << "|";
    oss.width(14);
    oss << rmsecv_verif[j] << "|";
    oss.width(14);
    d = 2*fabs(rmsecv[j]-rmsecv_verif[j])/(rmsecv[j]+rmsecv_verif[j]);
    oss << d << "|\n";
    if (d>0.01){
      oss << "Error! Diff is too big!\n";
    }
  }
  oss << "---|--------------||--------------|--------------|--------------|\n";

  delete [] rmse;
  delete [] rmsecv;
  delete [] rmsecv_verif;

  std::cout << oss.str();
  return oss.str();
}//



/*----------------------------------------------------*/
/*       TEST MULTIPLE                                */
/*----------------------------------------------------*/
std::string SGTELIB::test_multiple_occurrences (const std::string & s ){

  std::cout << "======================================================\n";
  std::cout << "SGTELIB::test_multiple_occurences\n";
  std::cout << s << "\n";

  // CONSTRUCT DATA
  const int p = 20;
  const int n = 2;

  // Build X0
  SGTELIB::Matrix X0 ("X0",p,n);
  X0.set_random(-3,+8,false);
  // Build Z0
  SGTELIB::Matrix Z0 = test_functions(X0);
  const int m = Z0.get_nb_cols();
  // Create multiple occurences in X0
  X0.set_row( X0.get_row(0) , 1 );
  X0.set_row( X0.get_row(0) , 2 );
  X0.set_row( X0.get_row(10) , 11 );

  #ifdef SGTELIB_DEBUG
    X0.display(std::cout);
    Z0.display(std::cout);
  #endif

  // CONSTRUCT REFERENCE MODEL
  bool ready;
  SGTELIB::TrainingSet C0(X0,Z0);
  SGTELIB::Surrogate * S0;
  S0 = SGTELIB::Surrogate_Factory(C0,s);
  ready = S0->build();

 // Some data to correct the kernel coef in KS and RBF.
  const SGTELIB::model_t mt = SGTELIB::Surrogate_Parameters::read_model_type(s);
  double dmean0, dmeanv, kc0, kcv , xsa0, xsav;
  dmean0 = C0.get_Ds_mean();
  xsa0 = C0.get_X_scaling_a(0);
  // Get original kernel coefficient
  kc0 = 0.0;
  kc0 = S0->get_param().get_kernel_coef();


  // Check ready
  if ( ! ready){
    surrogate_delete(S0);
    std::cout << "test_rmsecv: model ("+s+") is not ready\n";
    return       "test_rmsecv: model ("+s+") is not ready\n";
  }

  // Get the RMSECV metric
  double * rmsecv = new double [m];
  double * rmse  = new double [m];
  for (int j=0 ; j<m ; j++){
    rmsecv[j] = S0->get_metric(SGTELIB::METRIC_RMSECV,j);
    rmse[j] = S0->get_metric(SGTELIB::METRIC_RMSE,j);
  }
  // Delete the original model
  SGTELIB::surrogate_delete(S0);

  // INIT THE CROSS VALIDATION MATRICES
  SGTELIB::Matrix X0i ("X0i",p-1,n);   
  SGTELIB::Matrix Z0i ("Z0i",p-1,m); 
  for (int i=1 ; i<p ; i++){
    // Skip the first line
    X0i.set_row(X0.get_row(i),i-1);
    Z0i.set_row(Z0.get_row(i),i-1);
  }

  // Cross-validation models and trainingset
  SGTELIB::Surrogate * Sv;
  SGTELIB::TrainingSet * Cv;

  // Output for one point
  SGTELIB::Matrix Zvi ("Zvi",1,m); 

  // Cross validation output for all the training points:
  // (Supposed to be identical to _Zv of class Surrogate)
  SGTELIB::Matrix Zv_verif  ("Zv_verif",p,m);

 
  //std::cout << "SGTELIB::CV\n";
  // BUILD THE CV MODELS
  for (int i=0 ; i<p ; i++){
    std::cout << "BUILD CV MODELS " << i << "\n";
    // Build the trainingset (without point i);
    Cv = new SGTELIB::TrainingSet(X0i,Z0i);
    Cv->build();
    dmeanv = Cv->get_Ds_mean();
    xsav = Cv->get_X_scaling_a(0);
    // Init surrogate
    Sv = SGTELIB::Surrogate_Factory(*Cv,s);

    // Correct ks for RBF and KS
    kcv = kc0*(dmeanv/dmean0)*(xsa0/xsav);
    if (mt==SGTELIB::RBF) static_cast<SGTELIB::Surrogate_RBF*>(Sv)->set_kernel_coef(kcv);
    if (mt==SGTELIB::KS)  static_cast<SGTELIB::Surrogate_KS* >(Sv)->set_kernel_coef(kcv);
    // Build
    Sv->build();
    
    // Do the prediction on point i
    Sv->predict(X0.get_row(i),&Zvi);
    Zv_verif.set_row(Zvi,i);

    // update the matrices so that they lack the point i+1
    if (i<p-1){
      X0i.set_row(X0.get_row(i),i);
      Z0i.set_row(Z0.get_row(i),i);
    }

    delete Cv;
    surrogate_delete(Sv);
  }


  // Re-Compute the rmsecv
  double * rmsecv_verif = new double [m];  
  double e;
  for (int j=0 ; j<m ; j++){
    e = 0;
    for (int i=0 ; i<p ; i++){
      e += pow(Z0.get(i,j)-Zv_verif.get(i,j),2);
    }
    rmsecv_verif[j] = sqrt(e/p);
  }

  // Display
  double d;
  std::ostringstream oss; 
  oss << "Surrogate string: " << s << "\n";
  oss << "  j|          rmse||        rmsecv|  rmsecv_verif|      rel diff|\n";
  oss << "---|--------------||--------------|--------------|--------------|\n";
  for (int j=0 ; j<m ; j++){
     oss.width(3);
    oss << j << "|";
    oss.width(14);
    oss << rmse[j] << "||";
    oss.width(14);
    oss << rmsecv[j] << "|";
    oss.width(14);
    oss << rmsecv_verif[j] << "|";
    oss.width(14);
    d = 2*fabs(rmsecv[j]-rmsecv_verif[j])/(rmsecv[j]+rmsecv_verif[j]);
    oss << d << "|\n";
    if (d>0.01){
      oss << "Error! Diff is too big!\n";
    }
  }
  oss << "---|--------------||--------------|--------------|--------------|\n";

  delete [] rmse;
  delete [] rmsecv;
  delete [] rmsecv_verif;

  std::cout << oss.str();

  return oss.str();
}//






/*----------------------------------------------------*/
/*       test functions                               */
/*----------------------------------------------------*/
SGTELIB::Matrix SGTELIB::test_functions (const SGTELIB::Matrix & X){

  const int n = X.get_nb_cols(); // Input dim
  const int p = X.get_nb_rows(); // Nb of points
  const int m = 6; // There are 6 test functions, ie 6 outputs

  SGTELIB::Matrix T  ("T" ,p,1); // Aggregate input
  SGTELIB::Matrix ZT ("ZT",p,1); // 1 output
  

  const double div = 1.0/double(n);

  // Build matrix SX
  SGTELIB::Matrix SX ("SX",p,1); // Sum of x for each point
  for (int j=0 ; j<n ; j++){
    SX = SX + X.get_col(j);
  }
  SX = SX * div;

  // Build m
  SGTELIB::Matrix Z ("Z",p,m);
  
  for (int j=0 ; j<m ; j++){
    ZT.fill(0.0);
    for (int i=0 ; i<n ; i++){
      if (i==0){
        T = SX;
      }
      else{
        T = SX - X.get_col(i) * 2 * div;
      }
      ZT = ZT + test_functions_1D (T,j);
    }
    ZT = ZT * div;
    Z.set_col(ZT,j);
  }
  return Z;
}//

/*----------------------------------------------------*/
/*       Create some 1D test functions                */
/*----------------------------------------------------*/
double SGTELIB::test_functions_1D (const double t, const int function_index){

  switch (function_index){
    case 0:
      return 6.0*t*t + t - 1.0; // Quad function
    case 1:
      return t/(1.0+fabs(5.0*t)); // Sigmoid
    case 2:
      return 0.5-exp(-10*t*t); // bump
    case 3:
      return 0.5-((t>-0.2) && (t<0.5)); // square
    case 4:
      return 5.0*t-17.0*pow(t,3)+13*pow(t,5); // Oscillations/polynomial
    case 5:
      return sin(6.0*t)+cos(15.0*sqrt(fabs(t))); // Difficult function
    default:
      std::cout << "function_index : " << function_index << "\n";
      throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"test_function_1D : function_index not recognized" );
    }
}//


/*----------------------------------------------------*/
/*       Create some 1D test functions                */
/*----------------------------------------------------*/
SGTELIB::Matrix SGTELIB::test_functions_1D (const SGTELIB::Matrix & T, const int function_index){
  if (T.get_nb_cols()!=1){
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"test_function_1D : only for column vector!" );
  }
  const int p = T.get_nb_rows();
  SGTELIB::Matrix Z ("Z(T)",p,1);
  for (int i=0 ; i<p ; i++){
    Z.set(i,0,test_functions_1D(T.get(i,0),function_index));
  }
  return Z;
}//


/*----------------------------------------------------*/
/*       Check differences between two matrices       */
/*----------------------------------------------------*/
void SGTELIB::check_matrix_diff(const SGTELIB::Matrix * A, const SGTELIB::Matrix * B){
  // Check not NULL
  if ( !  A){
    std::cout << "A is NULL\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"check_matrix_diff : A is NULL" );
  }
  if ( !  B){
    std::cout << "B is NULL\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"check_matrix_diff : B is NULL" );
  }

  // Check dimension
  if (A->get_nb_rows()!=B->get_nb_rows()){
    std::cout << "Different number of rows !! " << A->get_nb_rows() << " " << B->get_nb_rows() << "\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"check_matrix_diff : != nb of rows" );
  }
  const int m = A->get_nb_rows();
  if (A->get_nb_cols()!=B->get_nb_cols()){
    std::cout << "Different number of cols !! " << A->get_nb_cols() << " " << B->get_nb_cols() << "\n";
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"check_matrix_diff : != nb of cols" );
  }
  const int n = A->get_nb_cols();

  double va,vb,dab;
  bool eij = false; // true if there is a problem with value (i,j)
  bool e = false;   // true if there is a problem, anywhere in the matrices

  for (int i=0 ; i<m ; i++){ 
    for (int j=0 ; j<n ; j++){
      va = A->get(i,j);
      vb = B->get(i,j);
      eij = false;
      dab = fabs(va-vb)/std::max( 0.5*(fabs(va)+fabs(vb)) , 1.0);
      if (dab>1e-6){
        eij = true;
        std::cout << "diff is too big !\n";
      }
      if (isnan(va)){
        eij = true;
        std::cout << "va is nan !\n";
      }
      if (isnan(vb)){
        eij = true;
        std::cout << "vb is nan !\n";
      }
      if (isinf(va)){
        eij = true;
        std::cout << "va is inf !\n";
      }
      if (isinf(vb)){
        eij = true;
        std::cout << "vb is inf !\n";
      }
      if (eij){
        e = true;
        std::cout << "A(" << i << "," << j << ") = " << va << "\n";
        std::cout << "B(" << i << "," << j << ") = " << vb << "\n";
        std::cout << "diff = " << fabs(va-vb) << "\n";
        std::cout << "dab  = " << dab << "\n";
      }
    }
  }
  if (e){
    A->display(std::cout);
    B->display(std::cout);
  }


}//









/*----------------------------------------------------*/
/*       build test data                              */
/*----------------------------------------------------*/
void SGTELIB::build_test_data ( const std::string & function_name , 
                                SGTELIB::Matrix & X0 , 
                                SGTELIB::Matrix & Z0 ){

  int p = 0;
  int n = 0;
  int m = 0;

  if ( (function_name=="hartman3") || (function_name=="hartman6") ){
    int q = 0;    
    SGTELIB::Matrix B,D;
    std::string B_str,D_str;
    if (function_name=="hartman3"){
      B_str = "3.0 10.0 30.0 ;"
              "0.1 10.0 35.0 ;"
              "3.0 10.0 30.0 ;"
              "0.1 10.0 35.0 ;";

      D_str = "0.3689 0.1170 0.2673  ;"
              "0.4699 0.4387 0.7470  ;"
              "0.1091 0.8732 0.5547  ;"
              "0.03815 0.5743 0.8828 ;";
    }
    else if (function_name=="hartman6"){
      B_str = "10.0 3.0  17.0 3.5  1.7  8.0  ;"
              "0.05 10.0 17.0 0.1  8.0  14.0 ;"
              "3.0  3.5  1.7  10.0 17.0 8.0  ;"
              "17.0 8.0  0.05 10.0 0.1  14.0 ;";

      D_str = "0.1312 0.1696 0.5569 0.0124 0.8283 0.5886 ;"
              "0.2329 0.4135 0.8307 0.3736 0.1004 0.9991 ;"
              "0.2348 0.1451 0.3522 0.2883 0.3047 0.6650 ;"
              "0.4047 0.8828 0.8732 0.5743 0.1091 0.0381 ;";
    }

    B = SGTELIB::Matrix::string_to_matrix(B_str);
    D = SGTELIB::Matrix::string_to_matrix(D_str);   
    n = B.get_nb_cols();
    q = B.get_nb_rows();
    m = 1;

    p = 100*(n+1);
    X0 = SGTELIB::Matrix("X0",p,n);
    Z0 = SGTELIB::Matrix("Z0",p,m);
    X0.set_random(0.0,1.0,false);

    double zi,eik,ak=0;
    for (int i=0 ; i<p ; i++){
      zi = 0;
      for (int k=0 ; k<q ; k++){
        eik = 0;
        for (int j=0 ; j<n ; j++){
          eik -= B.get(k,j)*pow(X0.get(i,j)-D.get(k,j),2.0);
        }
        switch (k){
          case 0: ak=1.0; break;
          case 1: ak=1.2; break;
          case 2: ak=3.0; break;
          case 3: ak=3.2; break;
        }
        zi -= ak * exp(eik);        
      }
      Z0.set(i,0,zi);
    }
    return;
  }// end hartman


  if ((function_name=="branin-hoo") || (function_name=="braninhoo")) {
    n = 2;
    m = 1;
    p = 100*(n+1);
    X0 = SGTELIB::Matrix("X0",p,n);
    Z0 = SGTELIB::Matrix("Z0",p,m);
    X0.set_random(0.0,1.0,false);
    X0.set_col(-5.0+15.0*X0.get_col(0),0);
    X0.set_col( 0.0+15.0*X0.get_col(1),1);
    #ifdef SGTELIB_DEBUG
      X0.display(std::cout);
    #endif

    double zi,x1,x2;
    for (int i=0 ; i<p ; i++){
      x1 = X0.get(i,0);
      x2 = X0.get(i,1);
      zi = pow(x2-5.1*x1*x1*.25/PI+5.*x1/PI-6.,2.) + 10.*(1.-1./(8.*PI))*cos(x1) + 10.;
      Z0.set(i,0,zi);
    }
    return;
  }// end branin-hoo



  if (function_name=="camelback") {
    n = 2;
    m = 1;
    p = 100*(n+1);
    X0 = SGTELIB::Matrix("X0",p,n);
    Z0 = SGTELIB::Matrix("Z0",p,m);
    X0.set_random(0.0,1.0,false);
    X0.set_col(-3.0+6.0*X0.get_col(0),0);
    X0.set_col(-2.0+4.0*X0.get_col(1),1);
    #ifdef SGTELIB_DEBUG
      X0.display(std::cout);
    #endif

    double zi,x1,x2;
    for (int i=0 ; i<p ; i++){
      x1 = X0.get(i,0);
      x2 = X0.get(i,1);
      zi = (pow(x1,4.)/3.-2.1*x1*x1+4.)*x1*x1 + x1*x2 + (4.*x2*x2-4.)*x2*x2;
      Z0.set(i,0,zi);
    }
    return;
  }// end camelback


  if (function_name=="rosenbrock") {
    n = 2;
    m = 1;
    p = 100*(n+1);
    X0 = SGTELIB::Matrix("X0",p,n);
    Z0 = SGTELIB::Matrix("Z0",p,m);
    X0.set_random(0.0,1.0,false);
    X0 = -5.+15.*X0;
    #ifdef SGTELIB_DEBUG
      X0.display(std::cout);
    #endif

    double zi,xj,xjp;
    for (int i=0 ; i<p ; i++){
      zi = 0;
      for (int j=0 ; j<n-1 ; j++){
        xj = X0.get(i,j);
        xjp = X0.get(i,j+1);
        zi += pow(1-xj,2)+100*pow(xjp-xj*xj,2);
      }
      Z0.set(i,0,zi);
    }
    return;
  }// end camelback


  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"build_test_data : function name not recognized" );

}



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
#include "sgtelib.hpp"
#include "Surrogate_Factory.hpp"
#include "Surrogate_Utils.hpp"
#include <fstream>
#include <string>
using namespace SGTELIB;

/*--------------------------------------*/
/*           main function              */
/*--------------------------------------*/
int main ( int argc , char ** argv ) {

  int i,j,ikw;
  bool error = false;

  //============================================
  // Check action
  //============================================
  std::vector<std::string> keyword;
  // Consider only the action keywords.
  keyword.push_back("-predict");
  keyword.push_back("-help");
  keyword.push_back("-test");
  keyword.push_back("-server");
  keyword.push_back("-best");

  // Check that exactly one action is required.
  std::string action = "null";
  for (i=1 ; i<argc ; i++){
    for (j=0 ; j<static_cast<int>(keyword.size()) ; j++){
      //std::cout << argv[i] << " vs " << keyword.at(j) << "\n";
      if (streq(argv[i],keyword.at(j))){
        // action kw found.
        if (streq(action,"null")){
          action = argv[i];
        }
        else{
          std::cout << "Error! Two actions have been requested: \"" << action << "\" and \"" << argv[i] << "\".\n";
          std::cout << "You can only perform one action at a time.\n";
          error = true;
        }
      }
    }
  }
  if (!strcmp(action.c_str(),"null")){
    std::cout << "Error! No action was found.\n";
    error = true;
  }

  //============================================
  // Add other keywords. 
  //============================================
  keyword.push_back("-model");
  keyword.push_back("-verbose");
  const int NKW = static_cast<int>(keyword.size());
  // Create empty strings to store the information following each keyword.
  std::vector<std::string> info;
  for (i=0 ; i<NKW ; i++){
    info.push_back("");
  }

  //============================================
  // Parse
  //============================================
  i = 1;
  j = -1;
  while (i<argc){
    // Check if the word is one of the keywords.
    bool found_kw = false;
    for  (ikw=0 ; ikw<NKW ; ikw++){
      if (!strcmp( keyword.at(ikw).c_str(),argv[i] )){
        j = ikw;
        found_kw = true;
        //std::cout << "keyword: " <<keyword.at(ikw) << "\n";
        break;
      }
    }
    if (j==-1){
      // If no kw was specified, break.
      error = true;
      break;
    }
    // Then continue ready the words for the current keyword.
    info.at(j) += " ";
    if ( ! found_kw) info.at(j) += std::string(argv[i]);
    i++;
  }

 //============================================ 
  // help
  //============================================ 
  bool verbose=false;
  for (i=0 ; i<NKW ; i++) {
    if (!strcmp(keyword.at(i).c_str(),"-verbose")) break;
  }
  if (info.at(i).size()){
    std::cout << "verbose mode\n";
    verbose = true;
  }


  //============================================
  // Deblank
  //============================================
  for (i=0 ; i<NKW ; i++) {
    if (info.at(i).size()){
      info.at(i) = SGTELIB::deblank(info.at(i));
    }
  }

  //============================================
  // Display error
  //============================================
  if (error){
    // Display results
    for (i=0 ; i<NKW ; i++) {
      if (info.at(i).size()){
        std::cout << keyword.at(i) << " \"" << info.at(i) << "\"\n";
      }
    }
    std::cout << "Could not parse command.\n";
  }

  //============================================ 
  // Search for best model
  //============================================ 
  if (!strcmp(action.c_str(),"-best")){
    for (i=0 ; i<NKW ; i++) {
      if (!strcmp(keyword.at(i).c_str(),"-best")) break;
    }
    SGTELIB::sgtelib_best(info.at(i),verbose);
  }

 
  
  //============================================ 
  // Get model name
  //============================================
  std::string model;
  for (i=0 ; i<NKW ; i++) {
    if (!strcmp(keyword.at(i).c_str(),"-model")) break;
  }
  model = info.at(i);
  if (!strcmp(model.c_str(),"")){
    model = "TYPE ENSEMBLE PRESET SUPER1 METRIC PRESS";
  }


  //============================================ 
  // predict
  //============================================ 
  if (!strcmp(action.c_str(),"-predict")){
    std::cout << "model: " << model << "\n";
    for (i=0 ; i<NKW ; i++) {
      if (!strcmp(keyword.at(i).c_str(),"-predict")) break;
    }
    SGTELIB::sgtelib_predict(info.at(i),model);
  }

  //============================================ 
  // help
  //============================================ 
  if (!strcmp(action.c_str(),"-help") || error){
    std::string model;
    for (i=0 ; i<NKW ; i++) {
      if (!strcmp(keyword.at(i).c_str(),"-help")) break;
    }
    SGTELIB::sgtelib_help(info.at(i));
  }

  //============================================ 
  // server
  //============================================ 
  if (!strcmp(action.c_str(),"-server")){
    std::cout << "model: " << model << "\n";
    SGTELIB::sgtelib_server(model,verbose);
  }

  //============================================ 
  // test
  //============================================ 
  if (!strcmp(action.c_str(),"-test")){
    SGTELIB::sgtelib_test();
  }

  //std::cout << "Quit sgtelib.\n";

  return 0;
}







/*--------------------------------------*/
/*           one prediction             */
/*--------------------------------------*/
void SGTELIB::sgtelib_predict( const std::string & file_list , const std::string & model ){
  bool error = false;
  std::string file;
  SGTELIB::Matrix X,Z,XX,ZZ;
  std::istringstream in_line (file_list);	
  if ( ( ! error) && (in_line >> file) && (SGTELIB::exists(file)) ){
    std::cout << "Read file " << file << "\n";
    X = SGTELIB::Matrix(file);
    //X.display(std::cout);
  }
  else{
    std::cout << "Could not find " << file << "\n";
    error = true;
  }
  if ( ( ! error) && (in_line >> file) && (SGTELIB::exists(file)) ){
    std::cout << "Read file " << file << "\n";
    Z = SGTELIB::Matrix(file);
    //Z.display(std::cout);
  }
  else{
    std::cout << "Could not find " << file << "\n";
    error = true;
  }
  if ( ( ! error) &&  (in_line >> file) && (SGTELIB::exists(file)) ){
    std::cout << "Read file " << file << "\n";
    XX = SGTELIB::Matrix(file);
    //XX.display(std::cout);
  }
  else{
    std::cout << "Could not find " << file << "\n";
    error = true;
  }
  if ( ! (in_line >> file)){
    std::cout << "No zz file (display output in terminal)\n";
    file = "null";
  }
  if ( ! error){
    SGTELIB::TrainingSet TS(X,Z);
    SGTELIB::Surrogate * S = Surrogate_Factory(TS,model);
    S->build();
    ZZ = SGTELIB::Matrix("ZZ",XX.get_nb_rows(),Z.get_nb_cols());
    S->predict(XX,&ZZ);
    ZZ.set_name("ZZ");
    if (strcmp(file.c_str(),"null")){
      std::cout << "Write output matrix in " << file << "\n";
      ZZ.write(file);
    }
    else{
      ZZ.display(std::cout);
    }

  }
  if (error){
    sgtelib_help();
  }
}//






/*--------------------------------------*/
/* find best model for a set of data    */
/*--------------------------------------*/
void SGTELIB::sgtelib_best( const std::string & file_list , const bool verbose){
  bool error = false;
  std::string file;
  SGTELIB::Matrix X,Z;
  std::istringstream in_line (file_list);	
  if ( ( ! error) && (in_line >> file) && (SGTELIB::exists(file)) ){
    std::cout << "Read file " << file << "\n";
    X = SGTELIB::Matrix(file);
  }
  else{
    std::cout << "Could not find " << file << "\n";
    error = true;
  }
  if ( ( ! error) && (in_line >> file) && (SGTELIB::exists(file)) ){
    std::cout << "Read file " << file << "\n";
    Z = SGTELIB::Matrix(file);
  }
  else{
    std::cout << "Could not find " << file << "\n";
    error = true;
  }
  if (error){
    sgtelib_help();
    return;
  }


  std::string m;
  double v;

  SGTELIB::TrainingSet TS(X,Z);
  for (int i=0 ; i<2 ; i++){
    std::cout << "=============================================================\n";
    if (i==0){
      std::cout << "Test models for optimization...\n";
      m = "AOECV";
    }
    else{
      std::cout << "Test models for prediction...\n";
      m = "ARMSECV";
    }
    std::cout << "=============================================================\n";

    std::vector<std::string> model_list;
    model_list.push_back("TYPE PRS DEGREE OPTIM RIDGE OPTIM");
    model_list.push_back("TYPE KS  KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE OPTIM");
    model_list.push_back("TYPE RBF PRESET O KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE OPTIM");
    model_list.push_back("TYPE RBF PRESET R KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE OPTIM RIDGE OPTIM");
    model_list.push_back("TYPE RBF PRESET I KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE OPTIM RIDGE OPTIM");
    model_list.push_back("TYPE LOWESS DEGREE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE OPTIM RIDGE OPTIM");
    model_list.push_back("TYPE KRIGING");
    model_list.push_back("TYPE CN");
    std::vector<double> model_metric;

    double vmin = +INF;
    std::cout << "Testing ";
    if (verbose) std::cout << "...\n";
    for (std::vector<std::string>::iterator it_model = model_list.begin() ; it_model != model_list.end(); ++it_model){
      std::string model_def = *it_model+" METRIC "+m+" BUDGET 500";
      SGTELIB::Surrogate * S = Surrogate_Factory(TS,model_def);
      if (verbose) std::cout << *it_model << "\n";
      else std::cout << model_type_to_str(S->get_type()) << " ";
      S->build();
      if (S->is_ready()){
        v = S->get_metric(SGTELIB::str_to_metric_type(m),0);
        model_metric.push_back(v);
        vmin = std::min(v,vmin);
        if (verbose) std::cout << "===> " << v << "\n";
      }
      else{
        model_metric.push_back(+INF);
        if (verbose) std::cout << "===> no ready\n";
      }
      // Store model definition after parameter optimization
      *it_model = S->get_string();
      // Delete model
      SGTELIB::surrogate_delete(S);
    }
    std::cout << "\n";
    std::cout << "Best metric: " << m << " = " << vmin << "\n";
    std::cout << "Best model(s):\n ";
    for (int k=0 ; k<int(model_list.size()) ; k++){
      if (model_metric.at(k)<=vmin+EPSILON){
        std::cout << model_list.at(k) << "\n";
      }
    }

  }
  std::cout << "=============================================================\n";

    
}//


/*--------------------------------------*/
/*           sgtelib server             */
/*--------------------------------------*/
void SGTELIB::sgtelib_server( const std::string & model , const bool verbose ){

  SGTELIB::TrainingSet * TS = NULL;
  SGTELIB::Surrogate * S = NULL;
  SGTELIB::Matrix X, Z, std, ei, cdf;
  std::string dummy_str;

  std::cout << "========== SERVER ==========================\n";  


  std::cout << "Start server\n";
  std::cout << "Remove all flag files...\n";
  dummy_str = system("rm flag_* 2>/dev/null");
  std::cout << "Ok.\n";

  int dummy_k = 0;
  int pxx = 0;
  int m = 0;
  const bool display = verbose;

  std::ofstream out;
  std::ifstream in;
  
  while (true){

    SGTELIB::wait(0.01);

    if (SGTELIB::exists("flag_new_data_transmit")){
      //------------------------------
      // NEW DATA
      //------------------------------
      std::cout << "============flag: new_data===================\n";
      dummy_k = 0;
      dummy_str = system("mv flag_new_data_transmit flag_new_data_received");
      if (display) std::cout << "Read new data";
      X = SGTELIB::Matrix("new_data_x.txt");
      Z = SGTELIB::Matrix("new_data_z.txt");
      if (display) X.display_short(std::cout);
      if (display) Z.display_short(std::cout);
      std::cout << X.get_nb_rows() << " new data points...\n";

      if ( ! S){
        if (display) std::cout << "First data: Build Trainig Set & Model\n";
        TS = new SGTELIB::TrainingSet(X,Z);
        S = Surrogate_Factory(*TS,model);
        m = TS->get_output_dim();
        dummy_str = system("rm -f flag_not_ready");
      }
      else{
        if (display) std::cout << "Add points to TS\n";
        TS->add_points(X,Z);
      }
      dummy_str = system("rm new_data_x.txt new_data_z.txt");
      if (display) std::cout << "Waiting...\n";

    }
    else if (SGTELIB::exists("flag_predict_transmit")){
      //------------------------------
      // PREDICT
      //------------------------------
      std::cout << "============flag: predict==================" << dummy_k++ << "\n";
      dummy_str = system("mv flag_predict_transmit flag_predict_received");

      bool ready = false;
      if (S){
        if (display) std::cout << "Build Sgte (" << TS->get_nb_points() << " pts)\n";
        ready = S->build();      
        if (display) S->display(std::cout);
      }

      // Read input matrix
      if (display) std::cout << "Read matrix... ";
      X = SGTELIB::Matrix("flag_predict_received");
      pxx = X.get_nb_rows();
      Z = SGTELIB::Matrix("ZZ",pxx,m);
      std = SGTELIB::Matrix("std",pxx,m);
      ei = SGTELIB::Matrix("ei",pxx,m);
      cdf = SGTELIB::Matrix("cdf",pxx,m);
      if (ready){
        // Do prediction
        if (display) std::cout << "Predict... ";
        S->predict(X,&Z,&std,&ei,&cdf);
        //S->predict(X,&Z,NULL,NULL,NULL);
        if (display) X.display_short(std::cout);
        if (display) Z.display_short(std::cout);
      }
      else{
        if (display) std::cout << "Surrogate not ready\n";
        dummy_str = system("touch flag_not_ready");
        Z.fill(+SGTELIB::INF);
      }
      
      // Open stream and write matrices
      if (display) std::cout << "Write output.\n";    
      out.open ("flag_predict_received");
      {
        // Correct matrix names
        Z.set_name("Z");
        std.set_name("std");
        ei.set_name("ei");
        cdf.set_name("cdf");
        // Write matrices 
        Z.display(out);
        std.display(out);
        ei.display(out);
        cdf.display(out);
      }
      out.close(); 

      // Change flag
      dummy_str = system("mv flag_predict_received flag_predict_finished");
      if (display) std::cout << "Waiting...\n";

    }
    else if (SGTELIB::exists("flag_cv_transmit")){
      //------------------------------
      // CV values
      //------------------------------
      std::cout << "============flag: cv values==================" << dummy_k++ << "\n";
      dummy_str = system("mv flag_cv_transmit flag_cv_received");

      bool ready = false;
      if (S) ready = S->build();      
      if ( ! ready){
        if (display) std::cout << "Surrogate not ready\n";
        dummy_str = system("touch flag_not_ready");
      }
      else{
        S->display(std::cout);
      }

      SGTELIB::Matrix Zh, Sh, Zv, Sv;
      if (ready){
        Zh = S->get_matrix_Zh();
        Sh = S->get_matrix_Sh();
        Zv = S->get_matrix_Zv();
        Sv = S->get_matrix_Sv();
      }
      else{
        if (display) std::cout << "Surrogate not ready\n";
        dummy_str = system("touch flag_not_ready");
      }
      
      // Open stream and write matrices
      if (display) std::cout << "Write output.\n";    
      out.open ("flag_cv_received");
      {
        // Correct matrix names
        Zh.set_name("Zh");
        Sh.set_name("Sh");
        Zv.set_name("Zv");
        Sv.set_name("Sv");
        // Write matrices 
        Zh.display(out);
        Sh.display(out);
        Zv.display(out);
        Sv.display(out);
      }
      out.close(); 

      // Change flag
      dummy_str = system("mv flag_cv_received flag_cv_finished");
      if (display) std::cout << "Waiting...\n";

    }
    else if (SGTELIB::exists("flag_metric_transmit")){
      //------------------------------
      // METRIC
      //------------------------------

      std::cout << "============flag: metric===================" << dummy_k++ << "\n";
      dummy_str = system("mv flag_metric_transmit flag_metric_received");

      bool ready = false;
      if (S) ready = S->build();      
      if ( ! ready){
        if (display) std::cout << "Surrogate not ready\n";
        dummy_str = system("touch flag_not_ready");
      }


      std::string metric_string;
      SGTELIB::metric_t mt;
      SGTELIB::Matrix metric_value;

      // Read metric_string
      in.open("flag_metric_received");
      in >> metric_string;
      in.close();
      if (display) std::cout << "metric : " << metric_string << "\n";   
      mt = SGTELIB::str_to_metric_type(metric_string);

      // Open stream and write metric
      if (display) std::cout << "Write output.\n";    
      out.open ("flag_metric_received");
      if (ready){
        metric_value = S->get_metric(mt);
        for (int j=0 ; j<metric_value.get_nb_cols() ; j++){
          out << metric_value[j] << " ";
          if (display){
            std::cout << "output[" << j << "]: " << metric_value[j] << "\n";
          }
        }
      }
      else{
        out << -1;
      }
      out.close(); 

      // Change flag
      dummy_str = system("mv flag_metric_received flag_metric_finished");
      if (display) std::cout << "Waiting...\n";

    }
    else if (SGTELIB::exists("flag_info_transmit")){
      //------------------------------
      // INFO
      //------------------------------

      std::cout << "============flag: info===================" << dummy_k++ << "\n";
      dummy_str = system("mv flag_info_transmit flag_info_received");

      bool ready = false;
      if (S) ready = S->build();      
      if ( ! ready){
        if (display) std::cout << "Surrogate not ready\n";
        dummy_str = system("touch flag_not_ready");
      }

      // Open stream and write metric
      if (ready) S->display(std::cout);
      else std::cout << "Not ready.";
      // Change flag
      dummy_str = system("mv flag_info_received flag_info_finished");
      if (display) std::cout << "Waiting...\n";

    }
    else if (SGTELIB::exists("flag_reset_transmit")){
      //------------------------------
      // RESET
      //------------------------------

      std::cout << "============flag: reset======================" << "\n";
      dummy_str = system("mv flag_reset_transmit flag_reset_received");

      surrogate_delete(S);
      delete TS;
      TS = NULL;
      S = NULL;

      if (display) std::cout << "Pointers: S=" << S << ", TS=" << TS << "\n";

      // Change flag
      dummy_str = system("mv flag_reset_received flag_reset_finished");
      if (display) std::cout << "Waiting...\n";

    }
    else if (SGTELIB::exists("flag_ping")){
      //------------------------------
      // PING 
      //------------------------------
      if (display) std::cout << "ping ";
      // Write state in the ping file

      bool ready = false;
      if (S){
        if (display) std::cout << "Build Sgte (" << TS->get_nb_points() << " pts)\n";
        ready = S->build();
      }
      if (ready){  
        std::cout << "pong: Model is ready.\n";
        dummy_str = system("echo 1 >> flag_ping");
      }
      else{
        std::cout << "pong: Model is not ready.\n";
        dummy_str = system("echo 0 >> flag_ping");
      }
      // Send an answer!
      dummy_str = system("mv flag_ping flag_pong");
    }
    else if (SGTELIB::exists("flag_quit")){
      //------------------------------
      // QUIT
      //------------------------------
      surrogate_delete(S);
      delete TS;
      dummy_str = system("rm flag_quit");
      std::cout << "flag: quit\n";
      break;
    }

  }

  std::cout << "Remove all flag files...";
  dummy_str = system("rm    flag_* 2>/dev/null");
  std::cout << "Ok.\n";

  std::cout << "Quit server\n";

}//





/*--------------------------------------*/
/*           help                       */
/*--------------------------------------*/
void SGTELIB::sgtelib_help( std::string word ){

  int i,j;
  if (!strcmp(word.c_str(),"")) word = "GENERAL";
  //std::cout << "HELP about " << word << "\n";

  std::string ** DATA = SGTELIB::get_help_data();
  const int NL = SGTELIB::dim_help_data();
  bool failedsearch = true;
  bool found = false;

  for (j=0 ; j<3 ; j++){
    // j=0 => search in title
    // j=1 => search in associated keywords
    // j=2 => search in content
    for (i=0 ; i<NL ; i++){
      if ( (SGTELIB::string_find(DATA[i][j],word)) || (streqi(word,"ALL")) ){
        //std::cout << i << " " << j << " " << (SGTELIB::string_find(DATA[i][j],word)) << " " << word << " " << DATA[i][j] << "\n";
        std::cout << "===============================================\n\n";
        std::cout << "  \33[4m" << DATA[i][0] << "\33[0m" << "\n\n";
        //std::cout << "  ***  " << DATA[i][0] << "  ***  " << "\n\n";
        std::cout << DATA[i][2] << "\n\n";
        found = true;
        failedsearch = false;
      }
    }
    // If some data where found for this value of j (i.e. for this depth of search)
    // then we don't search for higher j. 
    if (found) break;
  }

  // Search if the word appears in the keywords of some datum. 
  std::string SeeAlso = "\33[4mSee also\33[0m:";
  found = false;
  for (i=0 ; i<NL ; i++){
    if ( SGTELIB::string_find(DATA[i][1],word) ){
      SeeAlso+=" "+DATA[i][0];
      found = true;
    }
  }
  if (found){
    std::cout << "=======================================\n";
    std::cout << "\n" << SeeAlso << "\n\n";
    std::cout << "=======================================\n";
  }

  if (failedsearch){
      std::cout << "We could not find any information associated to your search.\n";
      SGTELIB::sgtelib_help("MAIN");
  }


}//






/*--------------------------------------*/
/*           test                       */
/*--------------------------------------*/
void SGTELIB::sgtelib_test( void ){

    SGTELIB::sand_box();

    SGTELIB::Matrix X0;
    SGTELIB::Matrix Z0;

    std::cout << "========== TEST MANY MODELS ==========================\n";  

    SGTELIB::build_test_data("hartman6",X0,Z0);
    SGTELIB::test_many_models("output_hartman6.txt",X0,Z0);

    SGTELIB::build_test_data("hartman3",X0,Z0);
    SGTELIB::test_many_models("output_hartman3.txt",X0,Z0);

    SGTELIB::build_test_data("braninhoo",X0,Z0);
    SGTELIB::test_many_models("output_braninhoo.txt",X0,Z0);

    SGTELIB::build_test_data("camelback",X0,Z0);
    SGTELIB::test_many_models("output_camelback.txt",X0,Z0);

    SGTELIB::build_test_data("rosenbrock",X0,Z0);
    SGTELIB::test_many_models("output_rosenbrock.txt",X0,Z0);


    std::cout << "========== END ================================\n"; 
}//

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


#include "sgtelib_help.hpp"


//================================
//  Get dimension of HELP_DATA
//================================
int SGTELIB::dim_help_data (void){
  return 32;
}//

//================================
//  Construct the help data
//================================
std::string ** SGTELIB::get_help_data (void){
  int i;
  const int NL = 32;
  const int NC = 3;
  std::string ** HELP_DATA = new std::string * [NL];
  for (i = 0 ; i<NL ; i++) HELP_DATA[i] = new std::string [NC];
  i = 0;
  //================================
  //      GENERAL
  //================================
  HELP_DATA[i][0] = "GENERAL";
  HELP_DATA[i][1] = "GENERAL MAIN SGTELIB HELP";
  HELP_DATA[i][2] = "sgtelib is a dynamic surrogate modeling library. Given a set of data points [X,z(X)], it allows to estimate the value of z(x) for any x.\n"
" \n"
"sgtelib can be called in 5 modes  \n"
" * -predict: build a model on a set of data points and perform a prediction on a set of prediction points. See PREDICT for more information. This requires the definition of a model with the option -model, see MODEL.\n"
"      sgtelib.exe -model <model description> -predict <input/output files>\n"
"      sgtelib.exe -model TYPE PRS DEGREE 2 -predict x.txt z.txt xx.txt zz.txt\n"
" \n"
" * -server: starts a server that can be interrogated to perform predictions or compute the error metric of a model. The server should be used via the Matlab interface (see SERVER). This requires the definition of a model with the option -model, see MODEL. \n"
"      sgtelib.exe -server -model <model description>\n"
"      sgtelib.exe -server -model TYPE LOWESS SHAPE_COEF OPTIM\n"
" \n"
" * -best: returns the best type of model for a set of data points\n"
"      sgtelib.exe -best <x file name> <z file name>\n"
"      sgtelib.exe -best x.txt z.txt\n"
" \n"
" * -help: allows to ask for some information about some keyword.\n"
"      sgtelib.exe -help keyword\n"
"      sgtelib.exe -help DEGREE\n"
"      sgtelib.exe -help LOWESS\n"
"      sgtelib.exe -help\n"
" \n"
" * -test: runs a test of the sgtelib library.\n"
"      sgtelib.exe -test\n"
" ";
  i++;
  //================================
  //      PREDICT
  //================================
  HELP_DATA[i][0] = "PREDICT";
  HELP_DATA[i][1] = "PREDICT PREDICTION INLINE SGTELIB";
  HELP_DATA[i][2] = "Performs a prediction in command line on a set of data provided through text files. If no ZZ file is provided, the predictions are displayed in the terminal. If no model is provided, the default model is used.\n"
" \n"
"Example\n"
"      sgtelib.exe -predict <x file name> <z file name> <xx file name>\n"
"      sgtelib.exe -predict x.txt z.txt xx.txt -model TYPE PRS DEGREE 2\n"
"      sgtelib.exe -predict x.txt z.txt xx.txt zz.txt";
  i++;
  //================================
  //      BEST
  //================================
  HELP_DATA[i][0] = "BEST";
  HELP_DATA[i][1] = "GENERAL BEST SGTELIB";
  HELP_DATA[i][2] = "Displays the description of the model that best fit the data provided in text files.\n"
" \n"
"Example\n"
"       sgtelib.exe -best x_file.txt z.file.txt ";
  i++;
  //================================
  //      SERVER
  //================================
  HELP_DATA[i][0] = "SERVER";
  HELP_DATA[i][1] = "SERVER MATLAB SGTELIB";
  HELP_DATA[i][2] = "Starts a sgtelib server. See MATLAB_SERVER for more details.\n"
" \n"
"Example\n"
"      sgtelib.exe -server -model TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM";
  i++;
  //================================
  //      MODEL
  //================================
  HELP_DATA[i][0] = "MODEL";
  HELP_DATA[i][1] = "MODEL FIELD DESCRIPTION MODEL_DESCRIPTION DEFINITION MODEL_DEFINITION TYPE";
  HELP_DATA[i][2] = "Models in sgtelib are defined by using a succession of field names (see FIELD for the list of possible fields) and field values. Each field name is made of one single word. Each field value is made of one single word or numerical value. It is good practice to start by the field name TYPE, followed by the model type. \n"
" \n"
"Possible field names  \n"
" * TYPE: mandatory field that specifies the type of model \n"
" * DEGREE: degree of the model for PRS and LOWESS models \n"
" * RIDGE: regularization parameter for PRS, RBF and LOWESS models \n"
" * KERNEL_TYPE: Kernel function for RBF, KS, LOWESS and KRIGING models \n"
" * KERNEL_SHAPE: Shape coefficient for RBF, KS and LOWESS models \n"
" * METRIC_TYPE: Error metric used as criteria for model parameter optimization/selection \n"
" * DISTANCE_TYPE: Metric used to compute the distance between points \n"
" * PRESET: Special information for some types of model \n"
" * WEIGHT_TYPE: Defines how the weights of Ensemble of model are computed \n"
" * BUDGET: Defines the parameter optimization budget \n"
" * OUTPUT: Defines the output text file ";
  i++;
  //================================
  //      FIELD
  //================================
  HELP_DATA[i][0] = "FIELD";
  HELP_DATA[i][1] = "FIELD NAME FIELD_NAME MODEL DEFINITION DESCRIPTION";
  HELP_DATA[i][2] = "A model description is composed of field names and field values.\n"
" \n"
"Example\n"
"      TYPE <model type> FIELD1 <field 1 value> FIELD2 <field 2 value>";
  i++;
  //================================
  //      PRS
  //================================
  HELP_DATA[i][0] = "PRS";
  HELP_DATA[i][1] = "TYPE POLYNOMIAL RESPONSE SURFACE QUADRATIC";
  HELP_DATA[i][2] = "PRS (Polynomial Response Surface) is a type of model. \n"
" \n"
"Authorized fields for this type of model  \n"
" * DEGREE (Can be optimized) \n"
" * RIDGE (Can be optimized) \n"
" \n"
"Example\n"
"      TYPE PRS DEGREE 2\n"
"      TYPE PRS DEGREE OPTIM RIDGE OPTIM";
  i++;
  //================================
  //      PRS_EDGE
  //================================
  HELP_DATA[i][0] = "PRS_EDGE";
  HELP_DATA[i][1] = "TYPE POLYNOMIAL RESPONSE SURFACE QUADRATIC DISCONTINUITY DISCONTINUITIES EDGE";
  HELP_DATA[i][2] = "PRS_EDGE (Polynomial Response Surface EDGE) is a type of model that allows to model discontinuities at 0 by using additional basis functions. \n"
" \n"
"Authorized fields for this type of model  \n"
" * DEGREE (Can be optimized) \n"
" * RIDGE (Can be optimized) \n"
" \n"
"Example\n"
"      TYPE PRS_EDGE DEGREE 2\n"
"      TYPE PRS_EDGE DEGREE OPTIM RIDGE OPTIM";
  i++;
  //================================
  //      PRS_CAT
  //================================
  HELP_DATA[i][0] = "PRS_CAT";
  HELP_DATA[i][1] = "TYPE POLYNOMIAL RESPONSE SURFACE QUADRATIC DISCONTINUITY DISCONTINUITIES";
  HELP_DATA[i][2] = "PRS_CAT (Categorical Polynomial Response Surface) is a type of model that allows to build one PRS model for each different value of the first component of x. \n"
" \n"
"Authorized fields for this type of model  \n"
" * DEGREE (Can be optimized) \n"
" * RIDGE (Can be optimized) \n"
" \n"
"Example\n"
"      TYPE PRS_CAT DEGREE 2\n"
"      TYPE PRS_CAT DEGREE OPTIM RIDGE OPTIM";
  i++;
  //================================
  //      RBF
  //================================
  HELP_DATA[i][0] = "RBF";
  HELP_DATA[i][1] = "TYPE RADIAL BASIS FUNCTION KERNEL";
  HELP_DATA[i][2] = "RBF (Radial Basis Function) is a type of model. \n"
" \n"
"Authorized fields for this type of model  \n"
" * KERNEL_TYPE (Can be optimized) \n"
" * KERNEL_COEF (Can be optimized) \n"
" * DISTANCE_TYPE (Can be optimized) \n"
" * RIDGE (Can be optimized) \n"
" * PRESET: \"O\" for RBF with linear terms and orthogonal constraints, \"R\" for RBF with linear terms and regularization term, \"I\" for RBF with incomplete set of basis functions. This parameter cannot be optimized. \n"
" \n"
"Example\n"
"      TYPE RBF KERNEL_TYPE D1 KERNEL_SHAPE OPTIM DISTANCE_TYPE NORM2";
  i++;
  //================================
  //      KS
  //================================
  HELP_DATA[i][0] = "KS";
  HELP_DATA[i][1] = "TYPE KERNEL SMOOTHING SMOOTHING_KERNEL";
  HELP_DATA[i][2] = "KS (Kernel Smoothing) is a type of model. \n"
" \n"
"Authorized fields for this type of model  \n"
" * KERNEL_TYPE (Can be optimized) \n"
" * KERNEL_COEF (Can be optimized) \n"
" * DISTANCE_TYPE (Can be optimized) \n"
" \n"
"Example\n"
"      TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM";
  i++;
  //================================
  //      KRIGING
  //================================
  HELP_DATA[i][0] = "KRIGING";
  HELP_DATA[i][1] = "TYPE GAUSSIAN PROCESS GP COVARIANCE";
  HELP_DATA[i][2] = "KRIGING is a type of model. \n"
" \n"
"Authorized fields for this type of model  \n"
" * RIDGE (Can be optimized) \n"
" * DISTANCE_TYPE (Can be optimized) \n"
" \n"
"Example\n"
"      TYPE KRIGING";
  i++;
  //================================
  //      LOWESS
  //================================
  HELP_DATA[i][0] = "LOWESS";
  HELP_DATA[i][1] = "TYPE LOCALLY WEIGHTED REGRESSION LOWESS LOWER RIDGE DEGREE KERNEL";
  HELP_DATA[i][2] = "LOWESS (Locally Weighted Regression) is a type of model. \n"
" \n"
"Authorized fields for this type of model  \n"
" * DEGREE: Must be 1 (default) or 2 (Can be optimized). \n"
" * RIDGE (Can be optimized) \n"
" * KERNEL_TYPE (Can be optimized) \n"
" * KERNEL_COEF (Can be optimized) \n"
" * DISTANCE_TYPE (Can be optimized) \n"
" \n"
"Example\n"
"      TYPE LOWESS DEGREE 1      TYPE LOWESS DEGREE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE D1      TYPE LOWESS DEGREE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM DISTANCE_TYPE OPTIM";
  i++;
  //================================
  //      ENSEMBLE
  //================================
  HELP_DATA[i][0] = "ENSEMBLE";
  HELP_DATA[i][1] = "TYPE WEIGHT SELECT SELECTION";
  HELP_DATA[i][2] = "ENSEMBLE is a type of model. \n"
" \n"
"Authorized fields for this type of model  \n"
" * WEIGHT: Defines how the ensemble weights are computed. \n"
" * METRIC: Defines which metric is used to compute the weights. \n"
" \n"
"Example\n"
"      TYPE ENSEMBLE WEIGHT SELECT METRIC OECV      TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV DISTANCE_TYPE NORM2 BUDGET 100";
  i++;
  //================================
  //      TYPE
  //================================
  HELP_DATA[i][0] = "TYPE";
  HELP_DATA[i][1] = "MODEL DESCRIPTION DEFINITION TYPE PRS KS PRS_EDGE PRS_CAT RBF LOWESS ENSEMBLE KRIGING CN";
  HELP_DATA[i][2] = "The field name TYPE defines which type of model is used.\n"
" \n"
"Possible model type  \n"
" * PRS: Polynomial Response Surface \n"
" * KS: Kernel Smoothing \n"
" * PRS_EDGE: PRS EDGE model \n"
" * PRS_CAT: PRS CAT model \n"
" * RBF: Radial Basis Function Model \n"
" * LOWESS: Locally Weighted Regression \n"
" * ENSEMBLE: Ensemble of surrogates \n"
" * KRIGING: Kriging model \n"
" * CN: Closest neighbor \n"
" \n"
"Example\n"
"      TYPE PRS: defines a PRS model.\n"
"      TYPE ENSEMBLE: defines an ensemble of models.";
  i++;
  //================================
  //      DEGREE
  //================================
  HELP_DATA[i][0] = "DEGREE";
  HELP_DATA[i][1] = "PRS LOWESS PRS_CAT PRS_EDGE";
  HELP_DATA[i][2] = "The field name DEGREE defines the degree of a polynomial response surface. The value must be an integer ge 1. \n"
"Allowed for models of type PRS, PRS_EDGE, PRS_CAT, LOWESS. \n"
"Default values  \n"
" * For PRS models, the default degree is 2. \n"
" * For LOWESS models, the degree must be 1 (default) or 2. \n"
" \n"
"Example\n"
"      TYPE PRS DEGREE 3 defines a PRS model of degree 3.\n"
"      TYPE PRS_EDGE DEGREE 2 defines a PRS_EDGE model of degree 2.\n"
"      TYPE LOWESS DEGREE OPTIM defines a LOWESS model where the degree is optimized.";
  i++;
  //================================
  //      RIDGE
  //================================
  HELP_DATA[i][0] = "RIDGE";
  HELP_DATA[i][1] = "PRS LOWESS PRS_CAT PRS_EDGE RBF";
  HELP_DATA[i][2] = "The field name RIDGE defines the regularization parameter of the model. \n"
"Allowed for models of type PRS, PRS_EDGE, PRS_CAT, LOWESS, RBF. \n"
"Possible values Real value ge 0. Recommended values are 0 and 0.001. \n"
"Default values Default value is 0.01. \n"
"Example\n"
"      TYPE PRS DEGREE 3 RIDGE 0 defines a PRS model of degree 3 with no ridge.\n"
"      TYPE PRS DEGREE OPTIM RIDGE OPTIM defines a PRS model where the degree and ridge coefficient are optimized.";
  i++;
  //================================
  //      KERNEL_TYPE
  //================================
  HELP_DATA[i][0] = "KERNEL_TYPE";
  HELP_DATA[i][1] = "KS RBF LOWESS GAUSSIAN BI-QUADRATIC BIQUADRATIC TRICUBIC TRI-CUBIC INVERSE SPLINES POLYHARMONIC";
  HELP_DATA[i][2] = "The field name KERNEL_TYPE defines the type of kernel used in the model. The field name KERNEL is equivalent. \n"
"Allowed for models of type RBF, RBFI, Kriging, LOWESS and KS. \n"
"Possible values  \n"
" * D1: Gaussian kernel (default) \n"
" * D2: Inverse Quadratic Kernel \n"
" * D3: Inverse Multiquadratic Kernel \n"
" * D4: Bi-quadratic Kernel \n"
" * D5: Tri-cubic Kernel \n"
" * D6: Exponential Sqrt Kernel \n"
" * D7: Epanechnikov Kernel \n"
" * I0: Multiquadratic Kernel \n"
" * I1: Polyharmonic splines, degree 1 \n"
" * I2: Polyharmonic splines, degree 2 \n"
" * I3: Polyharmonic splines, degree 3 \n"
" * I4: Polyharmonic splines, degree 4 \n"
" * OPTIM: The type of kernel is optimized \n"
" \n"
"Example\n"
"      TYPE KS KERNEL_TYPE D2 defines a KS model with Inverse Quadratic Kernel\n"
"      TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM defines a KS model with optimized kernel shape and type";
  i++;
  //================================
  //      KERNEL_COEF
  //================================
  HELP_DATA[i][0] = "KERNEL_COEF";
  HELP_DATA[i][1] = "KS RBF LOWESS";
  HELP_DATA[i][2] = "The field name KERNEL_COEF defines the shape coefficient of the kernel function. Note that this field name has no impact for KERNEL_TYPES I1, I2, I3 and I4 because these kernels do not include a shape parameter. \n"
"Allowed for models of type RBF, KS, KRIGING, LOWESS. \n"
"Possible values Real value ge 0. Recommended range is [0.1 , 10]. For KS and LOWESS model, small values lead to smoother models. \n"
"Default values By default, the kernel coefficient is optimized. \n"
"Example\n"
"      TYPE RBF KERNEL_COEF 10 defines a RBF model with a shape coefficient of 10.\n"
"      TYPE KS KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM defines a KS model with optimized kernel shape and type";
  i++;
  //================================
  //      DISTANCE_TYPE
  //================================
  HELP_DATA[i][0] = "DISTANCE_TYPE";
  HELP_DATA[i][1] = "KS RBF CN LOWESS";
  HELP_DATA[i][2] = "The field name DISTANCE_TYPE defines the distance function used in the model. \n"
"Allowed for models of type RBF, RBF, KS, LOWESS. \n"
"Possible values  \n"
" * NORM1: Euclidian distance \n"
" * NORM2: Distance based on norm 1 \n"
" * NORMINF: Distance based on norm infty \n"
" * NORM2_IS0: Tailored distance for discontinuity in 0. \n"
" * NORM2_CAT: Tailored distance for categorical models. \n"
" \n"
"Default values Default value is NORM2. \n"
"Example\n"
"      TYPE KS DISTANCE NORM2_IS0 defines a KS model tailored for VAN optimization.";
  i++;
  //================================
  //      WEIGHT
  //================================
  HELP_DATA[i][0] = "WEIGHT";
  HELP_DATA[i][1] = "ENSEMBLE SELECTION WTA1 WTA2 WTA3 WTA4 WTA";
  HELP_DATA[i][2] = "The field name WEIGHT defines the method used to compute the weights w of the ensemble of models. The keyword WEIGHT_TYPE is equivalent. Allowed for models of type ENSEMBLE.\n"
" \n"
"Possible values  \n"
" * WTA1: w_k propto metric_sum - metric_k (default) \n"
" * WTA3: w_k propto (metric_k + alpha metric_mean)^beta \n"
" * SELECT: w_k propto 1 if metric_k = metric_min \n"
" * OPTIM: w minimizes metric(w) \n"
" \n"
"Example\n"
"      TYPE ENSEMBLE WEIGHT SELECT METRIC RMSECV defines an ensemble of models which selects the model that has the best RMSECV.\n"
"      TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV defines an ensemble of models where the weights w are computed to minimize the RMSECV of the model.";
  i++;
  //================================
  //      OUTPUT
  //================================
  HELP_DATA[i][0] = "OUTPUT";
  HELP_DATA[i][1] = "OUT DISPLAY";
  HELP_DATA[i][2] = "Defines a text file in which informations will be recorded.";
  i++;
  //================================
  //      OPTIM
  //================================
  HELP_DATA[i][0] = "OPTIM";
  HELP_DATA[i][1] = "OPTIM BUDGET PARAMETERS PARAMETER OPTIMIZATION";
  HELP_DATA[i][2] = "The field value OPTIM indicate that the model parameter must be optimized. The default optimization criteria is the AOECV error metric.\n"
" \n"
"Parameters that can be optimized  \n"
" * DEGREE \n"
" * RIDGE \n"
" * KERNEL_TYPE \n"
" * KERNEL_COEF \n"
" * DISTANCE_TYPE \n"
" \n"
"Example\n"
"      TYPE PRS DEGREE OPTIM\n"
"      TYPE LOWESS DEGREE OPTIM KERNEL_TYPE OPTIM KERNEL_SHAPE OPTIM METRIC ARMSECV";
  i++;
  //================================
  //      METRIC
  //================================
  HELP_DATA[i][0] = "METRIC";
  HELP_DATA[i][1] = "PARAMETER OPTIMIZATION CHOICE SELECTION OPTIM BUDGET ENSEMBLE";
  HELP_DATA[i][2] = "The field name METRIC defines the metric used to select the parameters of the model (including the weights of Ensemble models).\n"
" \n"
"Allowed for models of type All types of model.\n"
" \n"
"Possible values  \n"
" * EMAX: Error Max \n"
" * EMAXCV: Error Max with Cross-Validation \n"
" * RMSE: Root Mean Square Error \n"
" * RMSECV: RMSE with Cross-Validation \n"
" * OE: Order Error \n"
" * OECV: Order Error with Cross-Validation \n"
" * LINV: Invert of the Likelihood \n"
" * AOE: Aggregate Order Error \n"
" * AOECV: Aggregate Order Error with Cross-Validation \n"
" \n"
"Default values AOECV.\n"
" \n"
"Example\n"
"      TYPE ENSEMBLE WEIGHT SELECT METRIC RMSECV defines an ensemble of models which selects the model that has the best RMSECV.";
  i++;
  //================================
  //      BUDGET
  //================================
  HELP_DATA[i][0] = "BUDGET";
  HELP_DATA[i][1] = "PARAMETER PARAMETERS OPTIM OPTIMIZATION";
  HELP_DATA[i][2] = "Budget for model parameter optimization. The number of sets of model parameters that are tested is equal to the optimization budget multiplied by the the number of parameters to optimize. \n"
"Default values 20\n"
" \n"
"Example\n"
"      TYPE LOWESS KERNEL_SHAPE OPTIM METRIC AOECV BUDGET 100\n"
"      TYPE ENSEMBLE WEIGHT OPTIM METRIC RMSECV BUDGET 50";
  i++;
  //================================
  //      SGTELIB_SERVER_START
  //================================
  HELP_DATA[i][0] = "SGTELIB_SERVER_START";
  HELP_DATA[i][1] = "Matlab server interface";
  HELP_DATA[i][2] = "Start a sgtelib model in a server from Matlab. \n"
" \n"
"Example\n"
"       sgtelib_server_start('TYPE PRS'); Start a sgtelib server with a PRS model\n"
"       sgtelib_server_start('TYPE LOWESS DEGREE 1'); Start a Lowess model\n"
"       sgtelib_server_start(model_name,true); Start a model defined in model_name and keep the window open";
  i++;
  //================================
  //      SGTELIB_SERVER_NEWDATA
  //================================
  HELP_DATA[i][0] = "SGTELIB_SERVER_NEWDATA";
  HELP_DATA[i][1] = "Matlab server interface data newdata";
  HELP_DATA[i][2] = "Add data points to the sgtelib model from Matlab. \n"
" \n"
"Example\n"
"       sgtelib_server_newdata(X,Z); Add data points [X,Z]";
  i++;
  //================================
  //      SGTELIB_SERVER_PREDICT
  //================================
  HELP_DATA[i][0] = "SGTELIB_SERVER_PREDICT";
  HELP_DATA[i][1] = "Matlab server interface prediction predict";
  HELP_DATA[i][2] = "Perform a prediction from Matlab.\n"
" \n"
"Example\n"
"       [ZZ,std,ei,cdf] = sgtelib_server_predict(XX); Prediction at points XX.";
  i++;
  //================================
  //      SGTELIB_SERVER_INFO
  //================================
  HELP_DATA[i][0] = "SGTELIB_SERVER_INFO";
  HELP_DATA[i][1] = "Matlab server interface";
  HELP_DATA[i][2] = "Command from Matlab. Use sgtelib_server_info to display information about the model.";
  i++;
  //================================
  //      SGTELIB_SERVER_METRIC
  //================================
  HELP_DATA[i][0] = "SGTELIB_SERVER_METRIC";
  HELP_DATA[i][1] = "Matlab server interface RMSE OECV RMSECV OE METRIC";
  HELP_DATA[i][2] = "Command from Matlab. Use sgtelib_server_stop(metric_name) to access the error metric of the model.\n"
" \n"
"Example\n"
"       m = sgtelib_server_metric('OECV'); Return the OECV error metric\n"
"       m = sgtelib_server_metric('RMSE'); Return the RMSE error metric";
  i++;
  //================================
  //      SGTELIB_SERVER_RESET
  //================================
  HELP_DATA[i][0] = "SGTELIB_SERVER_RESET";
  HELP_DATA[i][1] = "Matlab server interface reset";
  HELP_DATA[i][2] = "Reset the model of the sgtelib server from Matlab.";
  i++;
  //================================
  //      SGTELIB_SERVER_STOP
  //================================
  HELP_DATA[i][0] = "SGTELIB_SERVER_STOP";
  HELP_DATA[i][1] = "Matlab server interface stop";
  HELP_DATA[i][2] = "Stop the sgtelib server from Matlab.";
  //================================
  return HELP_DATA;
}//

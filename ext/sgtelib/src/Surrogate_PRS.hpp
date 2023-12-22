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

#ifndef __SGTELIB_SURROGATE_PRS__
#define __SGTELIB_SURROGATE_PRS__

#include "Surrogate.hpp"

namespace SGTELIB {

  /*--------------------------------------*/
  /*         Surrogate_PRS class        */
  /*--------------------------------------*/
  class DLL_API Surrogate_PRS : public SGTELIB::Surrogate {

    /*--------------------------------------------------------*/
    /*  these members are defined in the Surrogate superclass */
    /*--------------------------------------------------------*/
    // int _p; // number of data points in X and Z
    // int _n; // dimension -- number of variables
    // int _m; // number of outputs (includes the objective)

      
  private:
      
      void preComputeForJacobianAndHessian( void ); // Function to precompute matrices for Hessian of output functions.
      bool _preComputeForJacobianAndHessianDone;
      
      SGTELIB::Matrix *_M_dxj;
      SGTELIB::Matrix *_alpha_dxj;
      SGTELIB::Matrix **_M_dxidxj;
      SGTELIB::Matrix **_alpha_dxidxj;
      
      
      double _cond;
      
  protected:

    int _q; // Nb of basis function
    SGTELIB::Matrix _M; // Monomes
    SGTELIB::Matrix _H; // Design matrix
    SGTELIB::Matrix _Ai; // Inverse of Ht*H
    SGTELIB::Matrix _alpha; // Coefficients

    virtual const SGTELIB::Matrix compute_design_matrix ( const SGTELIB::Matrix& Monomes, 
                                                          const SGTELIB::Matrix & Xs );

    virtual void compute_dxi_matrices ( SGTELIB::Matrix& Monomes,
                                        SGTELIB::Matrix& alpha,
                                        const int i ) const ;
      
    // build model (private):
    virtual bool build_private (void) override;

    void predict_private ( const SGTELIB::Matrix & XXs,
                                 SGTELIB::Matrix * ZZs) override;
      
    
      
    virtual void predict_private_objective ( const std::vector<SGTELIB::Matrix *> & XXd,
                                             SGTELIB::Matrix * ZZsurr_around            ) override;


    // Compute metrics
    const SGTELIB::Matrix * get_matrix_Zvs (void) override;
      

    bool compute_alpha ( void );
      

  public:

    // Constructor
    Surrogate_PRS ( SGTELIB::TrainingSet & trainingset ,   
                    SGTELIB::Surrogate_Parameters param) ;

    // destructor:
    virtual ~Surrogate_PRS ( void );

    // Build the monome exponents
    static int get_nb_PRS_monomes(const int nvar, const int degree);
    static SGTELIB::Matrix get_PRS_monomes(const int nvar, const int degree);
    virtual void display_private ( std::ostream & out ) const override;
    
    SGTELIB::Matrix get_alpha() const;
      
    // TEMP. Uncomment for use in unittest. Alloe predict_grad to work
      void set_alpha(const SGTELIB::Matrix & alpha) {_alpha = alpha; } // Set alpha but the model is not ready ( maybe empty _training set)
      
      void set_PRS_monones(const SGTELIB::Matrix & monome ) { _M = monome; _ready = true; }
      
    void predict_obj ( const SGTELIB::Matrix & XX,
                        SGTELIB::Matrix * ZZs,
                       bool forcedOnEmptyTrainingSet = false) ;

    void predict_grad ( const SGTELIB::Matrix & XX,
                        SGTELIB::Matrix * ZZs,
                       bool forcedOnEmptyTrainingSet = false) ;
      
      void predict_hessian ( const SGTELIB::Matrix & XX,
                          SGTELIB::Matrix * ZZs,
                            int j,
                         bool forcedOnEmptyTrainingSet = false) ;
      
    SGTELIB::Matrix getModelOut(const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;
    // Get model output j at point x
    double getModelOut(const SGTELIB::Matrix & XX, int j, bool forceOnEmptyTrainingSet = false ) { return getModelOut(XX, forceOnEmptyTrainingSet).get(0,j); }
    
    double getModelObj(const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;
    SGTELIB::Matrix getModelGrad(const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;
    void getModelGrad(SGTELIB::Matrix * g, SGTELIB::Matrix * Mpredict_grad, const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ; // _m x _n matrix

    // Get hessian of model for output j on point x
    SGTELIB::Matrix getModelHessian(const SGTELIB::Matrix & XX, int j, bool forceOnEmptyTrainingSet = false) ;
    void getModelHessian(SGTELIB::Matrix * H, const SGTELIB::Matrix & XX, int j, bool forceOnEmptyTrainingSet = false) ;
    SGTELIB::Matrix getModelHessian(const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;
    void getModelHessian(SGTELIB::Matrix * H, const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;

    SGTELIB::Matrix getModelCons(const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;
    void getModelCons(SGTELIB::Matrix * cons, const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;
    SGTELIB::Matrix getModelJacobian(const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;
    void getModelJacobian(SGTELIB::Matrix * J, SGTELIB::Matrix * Mpredict_grad, const SGTELIB::Matrix & XX, bool forceOnEmptyTrainingSet = false) ;

    SGTELIB::Matrix getModelLagGrad(const SGTELIB::Matrix & XX, const SGTELIB::Matrix & Y, double sigma = 1.0, bool forceOnEmptyTrainingSet = false) ;
    void getModelLagGrad(SGTELIB::Matrix * Gx, SGTELIB::Matrix * Mpredict_grad, SGTELIB::Matrix * Jx, const SGTELIB::Matrix & XX, const SGTELIB::Matrix & Y, double sigma = 1.0, bool forceOnEmptyTrainingSet = false) ;
    SGTELIB::Matrix getModelLagHessian(const SGTELIB::Matrix & XX, const SGTELIB::Matrix & Y, double sigma = 1.0, bool forceOnEmptyTrainingSet = false) ;    
    void getModelLagHessian(SGTELIB::Matrix * Hx, const SGTELIB::Matrix & XX, const SGTELIB::Matrix & Y, double sigma = 1.0, bool forceOnEmptyTrainingSet = false) ;
      
    double getModelCond() const { return _cond;}

    // KKT conditions for inequality-problem with bounds
    static SGTELIB::Matrix compute_multiplier(const SGTELIB::Matrix & Grad, const SGTELIB::Matrix & Jacobian, const double rank_tol = 1E-15) ;
    static void compute_multiplier(SGTELIB::Matrix * multiplier, const SGTELIB::Matrix & Grad, const SGTELIB::Matrix & Jacobian, const double rank_tol = 1E-15) ;

  };
}

#endif

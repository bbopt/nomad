/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
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
/**
 \file   QuadModelSld.hpp
 \brief  Quadratic regression or MFN interpolation model (headers)
 \author Sebastien Le Digabel
 \date   2010-08-31
 \see    QuadModelSld.cpp
 */
#ifndef __NOMAD_4_2_QUAD_MODEL_SLD__
#define __NOMAD_4_2_QUAD_MODEL_SLD__

#include "../../Algos/Step.hpp"
#include "../../Eval/EvalPoint.hpp"
#include "../../Type/BBOutputType.hpp"

#include "../../nomad_nsbegin.hpp"

/// Quadratic interpolation type
enum interpolation_type
{
    MFN                          , ///< Minimum Frobenius Norm interpolation.
    REGRESSION                   , ///< Regression.
    WP_REGRESSION                , ///< Well-poised regression.
    UNDEFINED_INTERPOLATION_TYPE   ///< Undefined.
};

enum hnorm_type
{
    L1   ,  ///< norm L1
    L2   ,  ///< norm L2
    LINF    ///< norm Linf
};

/// Class for quadratic regression or MFN interpolation model.
class QuadModelSld  
{
    
private:
    
    std::vector<EvalPoint>                    _Y;    ///< Interpolation points.
    
    const BBOutputTypeList _bbot; ///< Blackbox output types.
    
    interpolation_type _interpolation_type;    ///< Interpolation type.
    
    int                       _n;                     ///< Dimension.
    int                       _nfree;                 ///< Number of free variables.
    int                       _n_alpha;               ///< Number of model coefficients.
    bool                    * _fixed_vars;            ///< Fixed variables.
    int                     * _index;                 ///< Var. indexes with fixed var.
    Point           ** _alpha;                 ///< Model coefficients.
    Point              _center;                ///< Model center.
    Point              _ref;                   ///< Reference for scaling.
    Point              _scaling;               ///< Scaling.
    bool                      _error_flag;            ///< Error flag.
    // std::list<Direction>  _dirP;               ///< Directions used for scaling (may be empty)
    Point _delta_m;                           ///< Mesh size used for scaling
    Double _epsilon;                           ///< Offset for direction scaling
    
    
    Double             _cond;                  ///< Condition number.
    
    /// Initialize alpha (model parameters).
    void init_alpha ( void );
    
    /// Check if an unscaled point is in \c B(center,radius) for a given radius.
    /**
     \param x      The unscaled point -- \b IN.
     \param radius The radius         -- \b IN.
     \return \c true is \c x is in \c B(center,radius).
     */
    //    bool is_within_radius ( const NOMAD::Point & x      ,
    //                            const NOMAD::Point & radius   ) const;
    
    /// Check the interpolation set \c Y.
    /**
     \return \c true if the interpolation set is valid.
     */
    bool check_Y ( void ) const;
    
    /// Check outputs before the integration into \c Y.
    /**
     \param bbo The outputs       -- \b IN.
     \param m   Number of outputs -- \b IN.
     return \c true if the \c m outputs are valid.
     */
    bool check_outputs ( const Point & bbo , int m ) const;
    
    /// Reduce the number of interpolation points.
    /**
     The points are sorted accorded to their distance to the model center.
     \param center      Center of the model                -- \b IN.
     \param max_Y_size  Max number of interpolation points -- \b IN.
     */
    //    void reduce_Y  ( const NOMAD::Point & center , int max_Y_size );
    
    /// Compute condition number.
    /**
     \param W   Matrix W given as a vector -- \b IN.
     \param n   Size of \c W               -- \b IN
     \param eps Epsilon                    -- \b IN.
     */
    void compute_cond ( const double * W , int n , double eps );
    
    /// Compute the cumulated error of a model for one output.
    /**
     The errors are computed on the interpolation set \c Y.
     \param bbo_index   Blackbox output index  -- \b IN.
     \param error       Cumulated error        -- \b OUT.
     \param min_rel_err Min relative error     -- \b OUT.
     \param max_rel_err Max relative error     -- \b OUT.
     \param avg_rel_err Average relative error -- \b OUT.
     */
    void compute_model_error ( int             bbo_index   ,
                              Double & error       ,
                              Double & min_rel_err ,
                              Double & max_rel_err ,
                              Double & avg_rel_err   ) const;
    
    /// Compute the maximal relative error of a model.
    /**
     The error is computed for the interpolation set \c Y.
     \return The maximal relative error.
     */
    Double compute_max_rel_err ( void ) const;
    
    /// Compute the element \c (i,j) of the interpolation matrix \c M(phi,Y).
    /**
     \param i Row index    -- \b IN.
     \param j Column index -- \b IN.
     */
    double compute_M ( int i , int j ) const;
    
    /// Construct Minimum Frobenius Norm (MFN) model.
    /**
     - This occurs when \c p+1 \c < \c (n+1)(n+2)/2.
     \param eps        Epsilon                               -- \b IN.
     \param max_mpn    Maximum \c m+n value for SVD matrices -- \b IN.
     \param max_Y_size Maximum number of elements in \c Y    -- \b IN.
     \return \c true if the construction succeeded
     */
    bool construct_MFN_model ( double eps , int max_mpn , int max_Y_size );
    
    /// Construct regression model.
    /**
     - This occurs when \c p+1 \c >= \c (n+1)(n+2)/2.
     \param eps        Epsilon                               -- \b IN.
     \param max_mpn    Maximum \c m+n value for SVD matrices -- \b IN.
     \param max_Y_size Maximum number of elements in \c Y    -- \b IN.
     \return \c true if the construction succeeded
     */
    bool construct_regression_model ( double eps        ,
                                     int    max_mpn    ,
                                     int    max_Y_size   );
    
    /// Construct well-poised (WP) model.
    /**
     \param max_Y_size Maximum number of elements in \c Y -- \b IN.
     \return \c true if the construction succeeded
     */
    bool construct_WP_model ( int max_Y_size );
    
    /// Find interpolation point with max Lagrange polynomial value.
    /**
     \param  li      Lagrange polynomial             -- \b IN.
     \param  Y       Interpolation points candidates -- \b IN.
     \param  i1      Initial index in \c Y           -- \b IN.
     \param  i2      Final index in \c Y             -- \b IN.
     \param  max_lix Absolute value of the max value -- \b OUT.
     \return Index of interpolation point.
     */
    int find_max_lix ( const Point                     & li      ,
                      const EvalPointList & Y       ,
                      int                                      i1      ,
                      int                                      i2      ,
                      Double                          & max_lix   ) const;
    
    /// Resolution of system \c F.[mu alpha_L]'=[f(Y) 0]' for MFN interpolation.
    /**
     \param U         Matrix \c F=U from the SVD decomposition \c U.W.V' -- \b IN.
     \param W         Matrix \c W from the SVD decomposition \c U.W.V'   -- \b IN.
     \param V         Matrix \c V from the SVD decomposition \c U.W.V'   -- \b IN.
     \param bbo_index Blackbox output index                              -- \b IN.
     \param alpha     Model parameters                                   -- \b IN.
     \param eps       Epsilon                                            -- \b IN.
     */
    void solve_MFN_system ( double      ** U         ,
                           double       * W         ,
                           double      ** V         ,
                           int            bbo_index ,
                           Point & alpha     ,
                           double         eps      ) const;
    
    /// Resolution of system \c F.alpha=M'.f(Y) for the regression.
    /**
     \param M         Matrix \c M                                        -- \b IN.
     \param U         Matrix \c F=U from the SVD decomposition \c U.W.V' -- \b IN.
     \param W         Matrix \c W from the SVD decomposition \c U.W.V'   -- \b IN.
     \param V         Matrix \c V from the SVD decomposition \c U.W.V'   -- \b IN.
     \param bbo_index Blackbox output index                              -- \b IN.
     \param alpha     Model parameters                                   -- \b IN.
     \param eps       Epsilon                                            -- \b IN.
     */
    void solve_regression_system ( double      ** M         ,
                                  double      ** U         ,
                                  double       * W         ,
                                  double      ** V         ,
                                  int            bbo_index ,
                                  Point & alpha     ,
                                  double         eps      ) const;
    
    /// Display Lagrange polynomials.
    /**
     \param l Lagrange polynomials -- \b IN.
     \param Y Interpolation set    -- \b IN.
     */
    void display_lagrange_polynomials
    ( const std::vector<Point      *> & l ,
     const EvalPointList & Y   ) const;
    
    
    /// SVD decomposition.
    /**
     - The \c mxn \c M matrix is decomposed into \c M=U.W.V'.
     \param error_msg Error message when the function returns \c false    -- \b OUT.
     \param M         The input \c mxn matrix; Will be replaced by \c U   -- \b IN/OUT.
     \param W         The output \c nxn diagonal matrix                   -- \b OUT.
     \param V         The output \c nxn matrix                            -- \b OUT.
     \param m         Number of rows in M                                 -- \b IN.
     \param n         Number of columns in M                              -- \b IN.
     \param max_mpn   Maximum allowed value for \c m+n; ignored if \c <=0 -- \b IN
     -- \b optional (default = \c 1500).
     \return A boolean equal to \c true if the decomposition worked.
     */
    bool SVD_decomposition ( std::string & error_msg      ,
                            double     ** M              ,
                            double      * W              ,
                            double     ** V              ,
                            int           m              ,
                            int           n              ,
                            int           max_mpn = 1500 );
    
#ifdef MODEL_STATS
    mutable Double _Yw; ///< Width of the interpolation set \c Y.
    
public:
    
    /// Access to the width of the interpolation set \c X (or \c Y).
    /**
     \return The width of the interpolation set.
     */
    const Double & get_Yw ( void ) const { return _Yw; }
#endif
    
    /*-------------------------------------------------------------------------*/
public:
    
    /// Constructor.
    /**
     \param bbot          Output types                -- \b IN.
     */
    QuadModelSld ( const std::vector<BBOutputType> & bbot , size_t n);
    
    /// Destructor.
    virtual ~QuadModelSld ( void );
    
    void setY(const EvalPointList & evalPoints)
    {
        _Y.resize(evalPoints.size());
        std::copy(evalPoints.begin(), evalPoints.end(), _Y.begin());
        _error_flag = false;
    }
    
    
    /// Evaluate a model at a given point.
    /**
     \param x     The point        -- \b IN.
     \param alpha Model parameters -- \b IN.
     \return Model value.
     */
    Double eval ( const Point & x     ,
                 const Point & alpha   ) const;
    
    /// Compute model \c h and \c f values at a point.
    /**
     \param x      The point                 -- \b IN.
     \param h_min  Value of \c h_min         -- \b IN..
     \param h_norm Norm used to compute \c h -- \b IN..
     \param h      Value of \c h             -- \b OUT.
     \param f      Value of \c f             -- \b OUT.
     */
    void eval_hf ( const Point  & x      ,
                  const Double & h_min  ,
                  hnorm_type     h_norm ,
                  Double       & h      ,
                  Double       & f        ) const;
    
    
    /// Access to the interpolation type.
    /**
     \return The interpolation type.
     */
    const interpolation_type & get_interpolation_type ( void ) const
    {
        return _interpolation_type;
    }
    
    /// Access to the center of the model.
    /**
     \return The center.
     */
    const Point & get_center ( void ) const { return _center; }
    
    /// Access to the dimension.
    /**
     \return The dimension \c n.
     */
    int get_n ( void ) const { return _n; }
    
    /// Access to the number of free variables.
    /**
     \return The number of free variables \c n.
     */
    int get_nfree ( void ) const { return _nfree; }
    
    /// Access to the model parameters.
    /**
     \return The model parameters \c alpha.
     */
    Point ** get_alpha ( void ) const { return _alpha; }
    
    /// Check if the model is ready for evaluations.
    /**
     \return A boolean equal to \c true if the model is ready.
     */
    bool check ( void ) const;
    
    /// Access to the fixed variables.
    /**
     \param i Variable index -- \b IN.
     \return \c true if variable \c i is fixed.
     */
    bool variable_is_fixed ( int i ) const { return _fixed_vars[i]; }
    
    /// Access to the number of interpolation points.
    /**
     \return The number of interpolation points \c nY=p+1.
     */
    int get_nY ( void ) const { return static_cast<int> ( _Y.size() ); }
    
    /// Access to the condition number.
    /**
     \return The condition number.
     */
    const Double & get_cond ( void ) const { return _cond; }
    
    /// Access to the error flag.
    /**
     \return The error flag.
     */
    bool get_error_flag ( void ) const { return _error_flag; }
    
    /// Construct the interpolation set \c Y.
    /**
     \param center               Model center                       -- \b IN.
     \param interpolation_radius Interpolation radius               -- \b IN.
     \param max_Y_size           Maximum number of elements in \c Y -- \b IN.
     */
    //    void construct_Y ( const NOMAD::Point & center               ,
    //                       const NOMAD::Point & interpolation_radius ,
    //                       int                  max_Y_size             );
    
    /// Construct \c m models (one by output).
    /**
     \param use_WP     Use or not well-poisedness            -- \b IN.
     \param eps        Epsilon                               -- \b IN.
     \param max_mpn    Maximum \c m+n value for SVD matrices -- \b IN.
     \param max_Y_size Maximum number of elements in \c Y    -- \b IN.
     */
    void construct ( bool   use_WP     ,
                    double eps        ,
                    int    max_mpn    ,
                    int    max_Y_size   );
    
    /// Define scaling to put all coordinates centered in \c [-r;r].
    /**
     - Looks also for fixed variables.
     \param r The \c r parameter corresponds to \c MODEL_RADIUS_FACTOR -- \b IN.
     */
    void define_scaling ( const Double & r );
    
    /// Define scaling based on directions. See paper: Reducing the number of function evaluations in Mesh Adaptive Direct Search algorithms, Audet, Ianni, LeDigabel, Tribes, 2014
    /**
     - Looks also for fixed variables.
     \param dirP    The \c dirP parameter corresponds to set of directions formin a hyper-cube centered on poll center -- \b IN.
     \param delta_m The \c delta_m parameter is the dimension of the mesh -- \b IN.
     \param epsilon The \c epsilon parameter is the hyper-cube offset from the poll center -- \b IN.
     */
    // void define_scaling_by_directions ( const std::list<Direction> & dirP, const Point & delta_m, const Double &epsilon  );
    
    
    /// Scale a point.
    /**
     \param x The point to scale -- \b IN/OUT.
     \return \c true if the scaling worked.
     */
    bool scale ( Point & x ) const;
    
    /// Unscale a point.
    /**
     \param x The point to unscale -- \b IN/OUT.
     \return \c true if the unscaling worked.
     */
    bool unscale ( Point & x ) const;
    
    /// Unscale the gradient at a point.
    /**
     \param x The grad to unscale -- \b IN/OUT.
     \return \c true if the unscaling worked.
     */
    bool unscale_grad ( Point & x ) const;
    
    
    /// Check if a caled point is inside the trust radius.
    /**
     \param x The scaled point -- \b IN.
     \return  \c true is \c x is in the trust radius.
     */
    bool is_within_trust_radius ( const Point & x ) const;
    
    /// Display the model coefficients.
    void display_model_coeffs ( ) const;
    
    /// Display the interpolation set \c Y.
    /**
     \param title Title of the display block -- \b IN
     --\b optional (default="interpolation set Y").
     */
    void display_Y ( const std::string    & title = "interpolation set Y" ) const;
    
    /// Display cumulated error on the interpolation points.
    void display_Y_error ( ) const;
};


#include "../../nomad_nsend.hpp"

#endif

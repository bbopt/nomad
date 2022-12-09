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
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits.h>
using namespace std;

/*-------------*/
/*  constants  */
/*-------------*/
const double PI       = 3.141592653589793;
const double SQRT5    = sqrt(5.0);
const double SQRT10   = sqrt(10.0);
const double C13      = 13.0;
const double C14      = 14.0;
const double C29      = 29.0;
const double C45      = 45.0;
const double DINT_MAX = INT_MAX;

const double Y1[] = { 0.14 , 0.18 , 0.22 , 0.25 , 0.29 ,
		      0.32 , 0.35 , 0.39 , 0.37 , 0.58 ,
		      0.73 , 0.96 , 1.34 , 2.10 , 4.39   };

const double Y2[] = { 0.1957 , 0.1947 , 0.1735 , 0.1600 ,
		      0.0844 , 0.0627 , 0.0456 , 0.0342 ,
		      0.0323 , 0.0235 , 0.0246 };

const double Y3[] = { 34780 , 28610 , 23650 , 19630 ,
		      16370 , 13720 , 11540 ,  9744 ,
		       8261 ,  7030 ,  6005 ,  5147 ,
		       4427 ,  3820 ,  3307 ,  2872   };

const double Y4[] = { 0.844 , 0.908 , 0.932 , 0.936 , 0.925 ,
		      0.908 , 0.881 , 0.850 , 0.818 , 0.784 ,
		      0.751 , 0.718 , 0.685 , 0.658 , 0.628 ,
		      0.603 , 0.580 , 0.558 , 0.538 , 0.522 ,
		      0.506 , 0.490 , 0.478 , 0.467 , 0.457 ,
		      0.448 , 0.438 , 0.431 , 0.424 , 0.420 ,
		      0.414 , 0.411 , 0.406 };

const double Y5[] = { 1.366 , 1.191 , 1.112 , 1.013 , 0.991 ,
		      0.885 , 0.831 , 0.847 , 0.786 , 0.725 ,
		      0.746 , 0.679 , 0.608 , 0.655 , 0.616 ,
		      0.606 , 0.602 , 0.626 , 0.651 , 0.724 ,
		      0.649 , 0.649 , 0.694 , 0.644 , 0.624 ,
		      0.661 , 0.612 , 0.558 , 0.533 , 0.495 ,
		      0.500 , 0.423 , 0.395 , 0.375 , 0.372 ,
		      0.391 , 0.396 , 0.405 , 0.428 , 0.429 ,
		      0.523 , 0.562 , 0.607 , 0.653 , 0.672 ,
		      0.708 , 0.633 , 0.668 , 0.645 , 0.632 ,
		      0.591 , 0.559 , 0.597 , 0.625 , 0.739 ,
		      0.710 , 0.729 , 0.720 , 0.636 , 0.581 ,
		      0.428 , 0.292 , 0.162 , 0.098 , 0.054   };

const double V[] = { 4 , 2 , 1 , 0.5 , 0.25 , 0.167 , 0.125 ,
		     0.1 , 0.0833 , 0.0714 , 0.0625 };

/*-----------------------------------------*/
/*  enumeration type for the problem type  */
/*-----------------------------------------*/
enum problem_type {
  SMOOTH ,  // smooth
  NONDIFF,  // nonsmooth
  WILD3  ,  // deterministic noise
  NOISY3    // random noise
};

problem_type get_pb_type ( int pb_type ) {
  switch ( pb_type ) {
  case 1: return SMOOTH;
  case 2: return NONDIFF;
  case 3: return WILD3;
  case 4: return NOISY3;
  }
  return SMOOTH;
}

void display_pb_type ( problem_type pb_type , ostream & out ) {
  switch ( pb_type ) {
  case SMOOTH : out << "SMOOTH" ; break;
  case NONDIFF: out << "NONDIFF"; break;
  case WILD3  : out << "WILD3"  ; break;
  case NOISY3 : out << "NOISY3" ; break;
  }
}

/*----------------------*/
/*   get problem name   */
/*----------------------*/
string get_problem_name ( int nprob ) {
  switch ( nprob ) {
  case  1: return "LINEAR_FULL";
  case  2: return "LINEAR_R1";
  case  3: return "LINEAR_R1Z";
  case  4: return "ROSENBROCK";
  case  5: return "HELICAL";
  case  6: return "POWELLSG";
  case  7: return "FREUDENSTEIN_ROTH";
  case  8: return "BARD";
  case  9: return "KOWALIK_OSBORNE";
  case 10: return "MEYER";
  case 11: return "WATSON";
  case 12: return "BOX3";
  case 13: return "JENNRICH_SAMPSON";
  case 14: return "BROWN_DENNIS";
  case 15: return "CHEBYQUAD";
  case 16: return "BROWN";
  case 17: return "OSBORNE1";
  case 18: return "OSBORNE2";
  case 19: return "BDQRTIC";
  case 20: return "CUBE";
  case 21: return "MANCINO";
  case 22: return "HEART8";
  }
  return "";
}

/*-----------------------*/
/*   get instance data   */
/*-----------------------*/
void get_instance_data ( int          ninstance ,       // IN: in {1,2,...,53}
			 problem_type pbtype    ,       // IN
			 int        & nprob     ,       // OUT
			 int        & n         ,       // OUT
			 int        & m         ,       // OUT
			 bool       & ns          ) {   // OUT
  switch ( ninstance ) {
  case  1: nprob =  1; n =  9; m = 45; ns = 0; break; // LINEAR_FULL
  case  2: nprob =  1; n =  9; m = 45; ns = 1; break; // LINEAR_FULL
  case  3: nprob =  2; n =  7; m = 35; ns = 0; break; // LINEAR_R1
  case  4: nprob =  2; n =  7; m = 35; ns = 1; break; // LINEAR_R1
  case  5: nprob =  3; n =  7; m = 35; ns = 0; break; // LINEAR_R1Z
  case  6: nprob =  3; n =  7; m = 35; ns = 1; break; // LINEAR_R1Z
  case  7: nprob =  4; n =  2; m =  2; ns = 0; break; // ROSENBROCK
  case  8: nprob =  4; n =  2; m =  2; ns = 1; break; // ROSENBROCK
  case  9: nprob =  5; n =  3; m =  3; ns = 0; break; // HELICAL
  case 10: nprob =  5; n =  3; m =  3; ns = 1; break; // HELICAL
  case 11: nprob =  6; n =  4; m =  4; ns = 0; break; // POWELLSG
  case 12: nprob =  6; n =  4; m =  4; ns = 1; break; // POWELLSG
  case 13: nprob =  7; n =  2; m =  2; ns = 0; break; // FREUDENSTEIN_ROTH
  case 14: nprob =  7; n =  2; m =  2; ns = 1; break; // FREUDENSTEIN_ROTH
  case 15: nprob =  8; n =  3; m = 15; ns = 0; break; // BARD
  case 16: nprob =  8; n =  3; m = 15; ns = 1; break; // BARD
  case 17: nprob =  9; n =  4; m = 11; ns = 0; break; // KOWALIK_OSBORNE
  case 18: nprob = 10; n =  3; m = 16; ns = 0; break; // MEYER
  case 19: nprob = 11; n =  6; m = 31; ns = 0; break; // WATSON
  case 20: nprob = 11; n =  6; m = 31; ns = 1; break; // WATSON
  case 21: nprob = 11; n =  9; m = 31; ns = 0; break; // WATSON
  case 22: nprob = 11; n =  9; m = 31; ns = 1; break; // WATSON
  case 23: nprob = 11; n = 12; m = 31; ns = 0; break; // WATSON
  case 24: nprob = 11; n = 12; m = 31; ns = 1; break; // WATSON
  case 25: nprob = 12; n =  3; m = 10; ns = 0; break; // BOX3
  case 26: nprob = 13; n =  2; m = 10; ns = 0; break; // JENNRICH_SAMPSON
  case 27: nprob = 14; n =  4; m = 20; ns = 0; break; // BROWN_DENNIS
  case 28: nprob = 14; n =  4; m = 20; ns = 1; break; // BROWN_DENNIS
  case 29: nprob = 15; n =  6; m =  6; ns = 0; break; // CHEBYQUAD
  case 30: nprob = 15; n =  7; m =  7; ns = 0; break; // CHEBYQUAD
  case 31: nprob = 15; n =  8; m =  8; ns = 0; break; // CHEBYQUAD
  case 32: nprob = 15; n =  9; m =  9; ns = 0; break; // CHEBYQUAD
  case 33: nprob = 15; n = 10; m = 10; ns = 0; break; // CHEBYQUAD
  case 34: nprob = 15; n = 11; m = 11; ns = 0; break; // CHEBYQUAD
  case 35: nprob = 16; n = 10; m = 10; ns = 0; break; // BROWN
  case 36: nprob = 17; n =  5; m = 33; ns = 0; break; // OSBORNE1
  case 37: nprob = 18; n = 11; m = 65; ns = 0; break; // OSBORNE2
  case 38: nprob = 18; n = 11; m = 65; ns = 1; break; // OSBORNE2
  case 39: nprob = 19; n =  8; m =  8; ns = 0; break; // BDQRTIC
  case 40: nprob = 19; n = 10; m = 12; ns = 0; break; // BDQRTIC
  case 41: nprob = 19; n = 11; m = 14; ns = 0; break; // BDQRTIC
  case 42: nprob = 19; n = 12; m = 16; ns = 0; break; // BDQRTIC
  case 43: nprob = 20; n =  5; m =  5; ns = 0; break; // CUBE
  case 44: nprob = 20; n =  6; m =  6; ns = 0; break; // CUBE
  case 45: nprob = 20; n =  8; m =  8; ns = 0; break; // CUBE
  case 46: nprob = 21; n =  5; m =  5; ns = 0; break; // MANCINO
  case 47: nprob = 21; n =  5; m =  5; ns = 1; break; // MANCINO
  case 48: nprob = 21; n =  8; m =  8; ns = 0; break; // MANCINO
  case 49: nprob = 21; n = 10; m = 10; ns = 0; break; // MANCINO
  case 50: nprob = 21; n = 12; m = 12; ns = 0; break; // MANCINO
  case 51: nprob = 21; n = 12; m = 12; ns = 1; break; // MANCINO
  case 52: nprob = 22; n =  8; m =  8; ns = 0; break; // HEART8
  case 53: nprob = 22; n =  8; m =  8; ns = 1; break; // HEART8
  }
}

/*--------------------------------------------------*/
/*                   dfoxs function                 */
/*--------------------------------------------------*/
/* transcripted in C++ from the matlab version      */
/* dfoxs.m of S. Wild                               */
/* http://www.mcs.anl.gov/~more/dfo/matlab/dfoxs.m  */
/*--------------------------------------------------*/
void dfoxs ( int      n      ,      // IN : dimension of problem
	     int      nprob  ,      // IN : problem index in {1,2,...,22}
	     double   factor ,      // IN : std pt is multiplied by factor
	     double * x        ) {  // OUT: array of size n

  int    i , j;
  double sum , tmp;

  switch ( nprob ) {
  case  1: // Linear function - full rank or rank 1
  case  2: // Linear function - full rank or rank 1
  case  3: // Linear function - full rank or rank 1
  case  8: // Bard function
  case 19: // Bdqrtic
    for ( i = 0 ; i < n ; ++i )
      x[i] = 1.0;
    break;
  case 4: // Rosenbrock function
    x[0] = -1.2;
    x[1] =  1.0;
    break;
  case 5: // Helical valley function
    x[0] = -1;
    for ( i = 1 ; i < n ; ++i )
      x[i] = 0.0;
    break;
  case 6: // Powell singular function
    x[0] =  3;
    x[1] = -1;
    x[2] =  0;
    x[3] =  1;
    break;
  case 7: // Freudenstein and Roth function
    x[0] =  0.5;
    x[1] = -2.0;
    break;
  case 9: // Kowalik and Osborne function
    x[0] = 0.250;
    x[1] = 0.390;
    x[2] = 0.415;
    x[3] = 0.390;
    break;
  case 10: // Meyer function
    x[0] = .02;
    x[1] = 4000;
    x[2] = 250;
    break;
  case 11: // Watson function
  case 16: // Brown almost-linear function
  case 20: // Cube
    for ( i = 0 ; i < n ; ++i )
      x[i] = 0.5;
    break;
  case 12: // Box 3-dimensional function
    x[0] = 0;
    x[1] = 10;
    x[2] = 20;
    break;
  case 13: // Jennrich and Sampson function
    x[0] = 0.3;
    x[1] = 0.4;
    break;
  case 14: // Brown and Dennis function
    x[0] = 25;
    x[1] =  5;
    x[2] = -5;
    x[3] = -1;
    break;
  case 15: // Chebyquad function
    for ( i = 1 ; i <= n ; ++i )
      x[i-1] = i/(n+1.0);
    break;
  case 17: // Osborne 1 function
    x[0] = 0.5;
    x[1] = 1.5;
    x[2] = 1.0;
    x[3] = 0.01;
    x[4] = 0.02;
    break;
  case 18: // Osborne 2 function
    x[ 0] = 1.3;
    x[ 1] = 0.65;
    x[ 2] = 0.65;
    x[ 3] = 0.7;
    x[ 4] = 0.6;
    x[ 5] = 3.0;
    x[ 6] = 5.0;
    x[ 7] = 7.0;
    x[ 8] = 2.0;
    x[ 9] = 4.5;
    x[10] = 5.5;
    break;
  case 21: // Mancino
    for ( i = 1 ; i <= n ; ++i ) {
      sum = 0;
      for ( j = 1 ; j <= n ; ++j ) {
	tmp  = i;
	tmp  = sqrt(tmp/j);
	sum += (tmp*( pow(sin(log(tmp)),5.0)+ pow(cos(log(tmp)),5.0)));
      }
      x[i-1] = -0.0008710996 * (pow(i-50,3.0) + sum);
    }
    break;
  case 22: // Heart8
    x[0] = -0.300;
    x[1] = -0.390;
    x[2] =  0.300;
    x[3] = -0.344;
    x[4] = -1.200;
    x[5] =  2.690;
    x[6] =  1.590;
    x[7] = -1.500;
    break;
  }

  for ( i = 0 ; i < n ; ++i )
    x[i] *= factor;
}

/*---------------------------------------------------*/
/*                  dfovec function                  */
/*---------------------------------------------------*/
/* transcripted in C++ from the matlab version       */
/* dfovec.m of S. Wild                               */
/* http://www.mcs.anl.gov/~more/dfo/matlab/dfovec.m  */
/*---------------------------------------------------*/
bool dfovec ( int            m     ,       // IN : number of outputs
	      int            n     ,       // IN : dimension of problem
	      const double * x     ,       // IN : array of size n
	      int            nprob ,       // IN : problem index in {1,2,...,22}
	      double       * fvec    ) {   // OUT: array of size m

  int    i , j;
  double tmp , tmp1 , tmp2 , tmp3 , tmp4 , den , sum = 0;

  // 1. Linear function - full rank:
  if ( nprob == 1 ) {
    for ( j = 0 ; j < n ; ++j )
      sum += x[j];
    tmp = 2*sum/m + 1;   
    for ( i = 0 ; i < m ; ++i ) {
      fvec[i] = -tmp;
      if ( i < n )
	fvec[i] += x[i];
    }
  }

  // 2. Linear function - rank 1:
  else if ( nprob == 2 ) {
    for ( j = 1 ; j <= n ; ++j )
      sum += j*x[j-1];
    for ( i = 1 ; i <= m ; ++i )
      fvec[i-1] = i*sum - 1;
  }

  // 3. Linear function - rank 1 with zero columns and rows:
  else if ( nprob == 3 ) {
    for ( j = 2 ; j < n ; ++j ) 
      sum += j*x[j-1];
    for ( i = 0 ; i < m-1 ; ++i )
      fvec[i] = i*sum - 1;
    fvec[m-1] = -1;
  }

  // 4. Rosenbrock function:
  else if ( nprob == 4 ) {
    fvec[0] = 10*(x[1] - x[0]*x[0]);
    fvec[1] = 1 - x[0];
  }

  // 5. Helical valley function:
  else if ( nprob == 5 ) {
    if ( x[0] > 0 )
      tmp = atan ( x[1]/x[0] ) / (2*PI);
    else if ( x[0] < 0 )
      tmp = atan ( x[1]/x[0] ) / (2*PI) + .5;
    else
      tmp = .25;
    fvec[0] = 10*(x[2] - 10*tmp);
    fvec[1] = 10*(sqrt(x[0]*x[0]+x[1]*x[1])-1);
    fvec[2] = x[2];
  }

  // 6. Powell singular function:
  else if ( nprob == 6 ) {
    fvec[0] = x[0] + 10*x[1];
    fvec[1] = SQRT5*(x[2] - x[3]);
    fvec[2] = pow(x[1] - 2*x[2],2.0);
    fvec[3] = SQRT10*pow(x[0] - x[3],2.0);
  }

  // 7. Freudenstein and Roth function:
  else if ( nprob == 7 ) {
    fvec[0] = -C13 + x[0] + ((5 - x[1])*x[1] - 2)*x[1];
    fvec[1] = -C29 + x[0] + ((1 + x[1])*x[1] - C14)*x[1];
  }

  // 8. Bard function:
  else if ( nprob == 8 ) {
    for ( i = 1 ; i <= 15 ; ++i ) {
      tmp1 = i;
      tmp2 = 16-i;
      tmp3 = tmp1;
      if ( i > 8 )
	tmp3 = tmp2;

      den = x[1]*tmp2 + x[2]*tmp3;
      
      if ( den == 0 )
	return false;

      fvec[i-1] = Y1[i-1] - (x[0] + tmp1/den);
    }
  }

  // 9. Kowalik and Osborne function:
  else if ( nprob == 9 ) {
    for ( i = 0 ; i < 11 ; ++i ) {
      tmp1 = V[i]*(V[i] + x[1]);
      tmp2 = V[i]*(V[i] + x[2]) + x[3];
      if ( tmp2 == 0.0 )
	return false;
      fvec[i] = Y2[i] - x[0]*tmp1/tmp2;
    }
  }

  // 10. Meyer function:
  else if ( nprob == 10 ) {
    for ( i = 0 ; i < 16 ; ++i ) {
      den = 5*(i+1) + C45 + x[2];
       if ( den == 0 )
	return false;
      fvec[i] = x[0]*exp(x[1]/den) - Y3[i];
    }
  }

  // 11. Watson function:
  else if ( nprob == 11 ) {
    for ( i = 1 ; i <= 29 ; ++i ) {
      tmp = i/C29;
      tmp1 = 0;
      tmp3 = 1;
      for ( j = 2 ; j <= n ; ++j ) {
	tmp1 += (j-1)*tmp3*x[j-1];
	tmp3 *= tmp;
      }
      tmp2 = 0;
      tmp3 = 1;
      for ( j = 1 ; j <= n ; ++j ) {
	tmp2 += tmp3*x[j-1];
	tmp3 *= tmp;
      }
      fvec[i-1] = tmp1 - tmp2*tmp2 - 1;
    }
    fvec[29] = x[0];
    fvec[30] = x[1] - x[0]*x[0] - 1;
  }

  // 12. Box 3-dimensional function:
  else if ( nprob == 12 ) {
    for ( i = 1 ; i <= m ; ++i ) {
      tmp = i;
      tmp1 = tmp/10.0;
      fvec[i-1] = exp(-tmp1*x[0]) - exp(-tmp1*x[1]) +
	          (exp(-tmp) - exp(-tmp1))*x[2];
    }
  }

  // 13. Jennrich and Sampson function:
  else if ( nprob == 13 ) {
    for ( i = 1 ; i <= m ; ++i )
      fvec[i-1] = 2 + 2*i - exp(i*x[0]) - exp(i*x[1]);
  }

  // 14. Brown and Dennis function:
  else if ( nprob == 14 ) {
    for ( i = 1 ; i <= m ; ++i ) {
      tmp = i/5.0;
      tmp1 = x[0] + tmp*x[1] - exp(tmp);
      tmp2 = x[2] + sin(tmp)*x[3] - cos(tmp);
      fvec[i-1] = tmp1*tmp1 + tmp2*tmp2;
    }
  }

  // 15. Chebyquad function:
  else if ( nprob == 15 ) {
    for ( i = 0 ; i < m ; ++i )
      fvec[i] = 0.0;
    for ( j = 0 ; j < n ; ++j ) {
      tmp1 = 1;
      tmp2 = 2*x[j] - 1;
      tmp3 = 2*tmp2;
      for ( i = 0 ; i < m ; ++i ) {
	fvec[i] += tmp2;
	tmp = tmp3*tmp2 - tmp1;
	tmp1 = tmp2;
	tmp2 = tmp;
      }
    }
    tmp = -1;
    for ( i = 1 ; i <= m ; ++i ) {
      fvec[i-1] /= n;
      if ( tmp > 0 )
	fvec[i-1] += 1.0 / (i*i - 1);
      tmp = -tmp;
    }
  }

  // 16. Brown almost-linear function:
  else if ( nprob == 16 ) {
    sum = -n-1;
    tmp = 1;
    for ( j = 0 ; j < n ; ++j ) {
      sum += x[j];
      tmp *= x[j];
    }
    for ( i = 0 ; i < n-1 ; ++i )
      fvec[i] = x[i] + sum;
    fvec[n-1] = tmp - 1;
  }

  // 17. Osborne 1 function:
  else if ( nprob == 17 ) {
    for ( i = 1 ; i <= 33 ; ++i ) {
      tmp = 10*(i-1);
      tmp1 = exp(-x[3]*tmp);
      tmp2 = exp(-x[4]*tmp);
      fvec[i-1] = Y4[i-1] - (x[0] + x[1]*tmp1 + x[2]*tmp2);
    }
  }

  // 18. Osborne 2 function:
  else if ( nprob == 18 ) {
    for ( i = 0 ; i < 65 ; ++i ) {
      tmp = i/10.0;
      tmp1 = exp(-x[4]*tmp);
      tmp2 = exp(-x[5]*pow(tmp-x[8],2.0));
      tmp3 = exp(-x[6]*pow(tmp-x[9],2.0));
      tmp4 = exp(-x[7]*pow(tmp-x[10],2.0));
      fvec[i] = Y5[i] - (x[0]*tmp1 + x[1]*tmp2 +
			 x[2]*tmp3 + x[3]*tmp4);
    }
  }

  // 19. Bdqrtic:
  else if ( nprob == 19 ) {
    for ( i = 1 ; i <= n-4 ; ++i ) {
      fvec[i-1  ] = (-4*x[i-1]+3.0);
      fvec[n-5+i] = (x[i-1]*x[i-1]+2*x[i]*x[i]+3*x[i+1]*x[i+1]
		     +4*x[i+2]*x[i+2]+5*x[n-1]*x[n-1]);
    }
  }

  // 20. Cube:
  else if ( nprob == 20 ) {
    fvec[0] = (x[0]-1.0);
    for ( i = 1 ; i < n ; ++i )
      fvec[i] = 10*(x[i]-pow(x[i-1],3.0));
  }

  // 21. Mancino:
  else if ( nprob == 21 ) {
    for ( i = 1 ; i <= n ; ++i ) {
      tmp1 = 0;
      for ( j = 1 ; j <= n ; ++j ) {
	tmp3 = i;
	tmp2 = sqrt ( x[i-1]*x[i-1] + tmp3/j) ;
	tmp1 += tmp2*( pow(sin(log(tmp2)),5.0) + pow(cos(log(tmp2)),5.0) );
      }
      fvec[i-1]=1400*x[i-1] + pow(i-50,3.0) + tmp1;
    }
  }

  // 22. Heart8:
  else if ( nprob == 22 ) {
    fvec[0] = x[0] + x[1] + 0.69;
    fvec[1] = x[2] + x[3] + 0.044;
    fvec[2] = x[4]*x[0] + x[5]*x[1] - x[6]*x[2] - x[7]*x[3] + 1.57;
    fvec[3] = x[6]*x[0] + x[7]*x[1] + x[4]*x[2] + x[5]*x[3] + 1.31;
    fvec[4] = x[0]*(x[4]*x[4]-x[6]*x[6]) - 2.0*x[2]*x[4]*x[6] +
      x[1]*(x[5]*x[5]-x[7]*x[7]) - 2.0*x[3]*x[5]*x[7] + 2.65;
    fvec[5] = x[2]*(x[4]*x[4]-x[6]*x[6]) + 2.0*x[0]*x[4]*x[6] +
      x[3]*(x[5]*x[5]-x[7]*x[7]) + 2.0*x[1]*x[5]*x[7] - 2.0;
    fvec[6] = x[0]*x[4]*(x[4]*x[4]-3.0*x[6]*x[6]) +
      x[2]*x[6]*(x[6]*x[6]-3.0*x[4]*x[4]) +
      x[1]*x[5]*(x[5]*x[5]-3.0*x[7]*x[7]) +
      x[3]*x[7]*(x[7]*x[7]-3.0*x[5]*x[5]) + 12.6;
    fvec[7] = x[2]*x[4]*(x[4]*x[4]-3.0*x[6]*x[6]) -
      x[0]*x[6]*(x[6]*x[6]-3.0*x[4]*x[4]) +
      x[3]*x[5]*(x[5]*x[5]-3.0*x[7]*x[7]) -
      x[1]*x[7]*(x[7]*x[7]-3.0*x[5]*x[5]) - 9.48;
  }

  return true;
}

/*---------------------------------------------------*/
/*                   calfun function                 */
/*---------------------------------------------------*/
/* transcripted in C++ from the matlab version       */
/* calfun.m of S. Wild                               */
/* http://www.mcs.anl.gov/~more/dfo/matlab/calfun.m  */
/*---------------------------------------------------*/
double calfun ( int            m        ,     // IN : number of outputs.
		int            n        ,     // IN : dimension of problem
		int            nprob    ,     // IN : problem index in {1,...,22}
		problem_type   probtype ,     // IN : problem type
		const double * x        ,     // IN : array of size n
		bool         & error      ) { // OUT: error flag

  error = false;

  int    i;
  double tmp , f = 0.0 , n1 = 0.0 , n2 = 0.0 , ninf = 0.0 ,
       * fvec = new double[m] , * xc = new double[n];

  // restrict domain for some nondiff problems:
  if ( probtype == NONDIFF &&
       ( nprob  ==  8 ||
	 nprob  ==  9 ||
	 nprob  == 13 || 
	 nprob  == 16 ||
	 nprob  == 17 ||
	 nprob  == 18    ) )
    for ( i = 0 ; i < n ; ++i )
      xc[i] = ( x[i] < 0.0 ) ? 0.0 : x[i];
  else if ( probtype == WILD3 ) {
    for ( i = 0 ; i < n ; ++i ) {
      xc[i] = x[i];
      tmp = fabs ( x[i] );
      n1 += tmp;
      n2 += x[i]*x[i];
      if ( tmp > ninf )
	ninf = tmp;
    }
    n2 = sqrt(n2);
  }
  else
    for ( i = 0 ; i < n ; ++i )
      xc[i] = x[i];
  
  // generate the vector:
  if ( !dfovec ( m , n , xc , nprob , fvec ) ) {   
    error = true;
    delete [] xc;
    delete [] fvec;
    return 1e100;
  }

  delete [] xc;

  // calculate the function value:
  switch ( probtype ) {

  case SMOOTH: // smooth
    for ( i = 0 ; i < m ; ++i )
      f += fvec[i]*fvec[i];
    break;
    
  case NONDIFF: // nonsmooth
    for ( i = 0 ; i < m ; ++i )
      f += fabs ( fvec[i] );
    break;
    
  case WILD3: // deterministic noise
    tmp  = 0.9*sin(100*n1)*cos(100*ninf) + 0.1*cos(n2);
    tmp *= 4.0*tmp*tmp - 3.0;
    for ( i = 0 ; i < m ; ++i )
      f += fvec[i]*fvec[i];
    f *= 1.0 + 0.001*tmp;
    break;

  case NOISY3: // random noise
    for ( i = 0 ; i < m ; ++i ) {
      tmp  = rand() / DINT_MAX;
      tmp  = 1.0 + 0.001 * ( tmp * 2.0 - 1.0 );
      tmp *= fvec[i];
      f += tmp*tmp;
    }    
    break;
  }

  delete [] fvec;

  return f;
}

/*--------------------------------------*/
/*            main function             */
/*--------------------------------------*/
int main ( int argc , char ** argv ) {

  if ( argc != 2 ) {
    cout << "error (no input file)" << endl;
    return 1;
  }

  // problem 7777777_WILD3:
  const int          ninstance = 7;
  const problem_type pbtype    = WILD3;
    
  int          nprob;
  int          n    ;
  int          m    ;
  bool ns;
  get_instance_data ( ninstance , pbtype , nprob , n , m , ns );
    
  // bool ns;
  // get_instance_data ( ninstance , pbtype , nprob , n , m , ns );
  //   {
  //     cout << endl << "MORE_WILD_" << ninstance
  // 	 << "_" << get_problem_name ( nprob ) << "_";
  //     display_pb_type ( pbtype , cout );
  //     cout << " (nprob=" << nprob << ",n="
  // 	 << n << ",m=" << m << ",ns=" << ns << ")"
  // 	 << endl;
  //   }

  int      i;
  double * x = new double[n];

  // read the point:
  ifstream in ( argv[1] );

  for ( i = 0 ; i < n ; ++i )
    in >> x[i];

  in.close();

  if ( in.fail() ) {
    cerr << endl << "could not read input file " << argv[1]
	 << endl << endl;
    delete [] x;
    return 1;
  }

//   cout << "x = ( ";
//   for ( i = 0 ; i < n ; ++i )
//     cout << x[i] << " ";
//   cout << ")" << endl;

  bool   error;
  double f = calfun ( m , n , nprob , pbtype , x , error );

  if ( error )
    cout << "error";
  else
    cout << f;

  cout << endl;

  delete [] x;

  return 0;
}

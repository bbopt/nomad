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

#include "Surrogate_Factory.hpp"


/*----------------------------------------------------------*/
SGTELIB::Surrogate * SGTELIB::Surrogate_Factory (SGTELIB::Matrix & X0,
                                                 SGTELIB::Matrix & Z0,
                                                 const std::string & s ){
  SGTELIB::TrainingSet * TS;
  TS = new SGTELIB::TrainingSet(X0,Z0);
  TS->info();
  throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
       "Surrogate_factory: constructor from matrices is forbiden." );
  return SGTELIB::Surrogate_Factory(*TS,s);
}//
/*----------------------------------------------------------*/





/*----------------------------------------------------------*/
SGTELIB::Surrogate * SGTELIB::Surrogate_Factory ( SGTELIB::TrainingSet & TS,
                                                  const std::string & s ) {
/*----------------------------------------------------------*/

  #ifdef SGTELIB_DEBUG
    std::cout << "SGTELIB::Surrogate_Factory (TS,p) begin\n";
    std::cout << "s = " << s << "\n";
    TS.info();
  #endif

  SGTELIB::Surrogate * S;
  SGTELIB::Surrogate_Parameters p ( s );
  


  switch ( p.get_type() ) {

  case SGTELIB::SVN: 
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,
      "Surrogate_Factory: not implemented yet! \""+s+"\"" );

  case SGTELIB::PRS: 
    S = new Surrogate_PRS(TS,p);
    break;

  case SGTELIB::PRS_EDGE: 
    S = new Surrogate_PRS_EDGE(TS,p);
    break;

  case SGTELIB::PRS_CAT: 
    S = new Surrogate_PRS_CAT(TS,p);
    break;

  case SGTELIB::KS: 
    S = new Surrogate_KS(TS,p);
    break;

  case SGTELIB::CN: 
    S = new Surrogate_CN(TS,p);
    break;

  case SGTELIB::RBF: 
    S = new Surrogate_RBF(TS,p);
    break;

  case SGTELIB::LOWESS: 
    S = new Surrogate_LOWESS(TS,p);
    break;

  case SGTELIB::ENSEMBLE: 
    S = new Surrogate_Ensemble(TS,p);
    break;

  case SGTELIB::KRIGING: 
    S = new Surrogate_Kriging(TS,p);
    break;

  default: 
    throw SGTELIB::Exception ( __FILE__ , __LINE__ ,"Undefined type" );
  }


  #ifdef SGTELIB_DEBUG
    std::cout << "SGTELIB::Surrogate_Factory (TS,p) AFTER set param\n";
    std::cout << "TS.info()\n";
    TS.info();
    std::cout << "S->info()\n";
    S->info();
    std::cout << "SGTELIB::Surrogate_Factory (TS,p) RETURN\n";
  #endif
  return S;

}//


/*----------------------------------------------------------*/
void SGTELIB::surrogate_delete ( SGTELIB::Surrogate * S ){
/*----------------------------------------------------------*/
  if (S){
    #ifdef SGTELIB_DEBUG
      std::cout << "Delete surrogate\n";
    #endif
    delete S;
    S = NULL;
  }
}//


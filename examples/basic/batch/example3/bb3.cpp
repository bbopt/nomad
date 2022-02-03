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
//
//  BB_3
//
//  Created by Yassine on 2021-06-22.
//
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;
/* BB optimisation de l'exemple 3.4 p44 du livre, sans contrainte. point de départ [2,2] minimum de la fonction est censé être le point [2/5 , -1/3]
 */
int main(int argc, const char ** argv) {
    
    double f = 1e20; // pas compris pourquoi on prend cette valeur initiale
    double x[2];
    if (argc >=2) { // pas compris le rôle de cette ligne
        
        ifstream in (argv[1]); // associer le file du point de depart à in
        for ( int i = 0 ; i < 2 ; i++ ) {
          in >> x[i];
        }
        f = pow (5 * x[0]-2 , 4) + pow (5 * x[0]-2, 2) * pow( x[1] , 2) +pow ( 3 * x[1] + 1 , 2); // calcul de l'ojectif f

        if ( in.fail() ) // Echec de la lecture du fichier
            f = 1e20; // pas compris pourquoi on prend cette valeur en cas d'echec
    }
    
       cout << f << endl ; // retourne que l'objectif car aucune contrainte dans ce problème
       return 0; // evaluation went well = retourne 0
  }
    

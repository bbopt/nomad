/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
#include "ibex.h"
#include "fstream"

using namespace ibex;
using namespace std;

class SetVolume : public SetVisitor {

private :

    double volume_yes;
    double volume_no;
    double volume_maybe;
    IntervalVector frame;

public:

  SetVolume(IntervalVector box) {volume_yes = 0; volume_no = 0; volume_maybe = 0; frame=box;}

  /**
   * Function that will be called automatically on every boxes (leaves) of the set.
   */
  void visit_leaf(const IntervalVector& box, BoolInterval status) {
  
    // Intersect the box with the frame
    IntervalVector framebox=box&frame;

    //  Associate a color to the box.
    //  - YES (means "inside") is in green
    //  - NO (means "outside") is in red
    //  - MAYBE (means "boundary") is in blue.

    switch (status) {
    case YES:  volume_yes += framebox.volume(); break;
    case NO:   volume_no += framebox.volume(); break;
    case MAYBE : volume_maybe += framebox.volume(); break;
    }
    
    }

 double get_volume_no(){return volume_no;}
 double get_volume_yes(){return volume_yes;}
 double get_volume_maybe(){return volume_maybe;}
};

int main(int argc, char ** argv)
{
	if (argc < 2) {cout << "Il manque la précision" << endl;}
	
	else {
	
	// Création du set
	System sys("system.txt");
  	IntervalVector box = sys.box;
  	const int n = sys.nb_ctr;
  	Array<NumConstraint> constraints = sys.ctrs;
  	Array<Sep> separators(n);
  	for (int i = 0; i<n; i++)
  	{
		separators.set_ref(i, *new SepFwdBwd(constraints[i]));
  	}
  	SepInter sep_poly_in(separators);
  	Set set(box);
  	double precision = atof(argv[1]);
  	sep_poly_in.contract(set, precision);
  	set.save("set.txt");
		
	// Calcul du volume
	SetVolume s_v(box);
	set.visit(s_v);
	
	// Résultats
	fstream f;
	f.open("volume.txt", fstream::out | fstream::app);
	f << "Précision : " << to_string(precision) << endl;
	f << "Volume of the boxes YES : " << s_v.get_volume_yes() << endl;
	f << "Volume of the boxes NO : " << s_v.get_volume_no() << endl;
	f << "Volume of the boxes MAYBE : " << s_v.get_volume_maybe() << endl;
	f << "Volume of MAYBE / Volume of MAYBE + Volume of YES in % : " << s_v.get_volume_maybe() / (s_v.get_volume_yes() + s_v.get_volume_maybe()) *100 << "%" << endl;
	f << endl << endl;	
	}
}

/*
** This software is governed by the CeCILL-C V2 license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
package Examples;

import static java.lang.Math.pow;

import jNomad.AllParameters;
import jNomad.Double;
import jNomad.EvalPoint;
import jNomad.Evaluator;
import jNomad.MainStep;

/*
 * A simple program to test jNomad with the well-known Rosenbrock problem
 */
public class Rosenbrock {
	static {
		String nomadHome = System.getenv("NOMAD_HOME");
		if (nomadHome == null) {
			System.out.println("NOMAD_HOME environment variable not set...");
			System.out.println(" Using default location /opt/Nomad/");
			System.out.println(
					" if Nomad is installed elsewhere, " + "please set NOMAD_HOME environment variable accordingly");
			nomadHome = new String("/opt/Nomad/"); // default location
		}
		System.out.println("NOMAD_HOME=" + nomadHome);
		try {
			String libpath = new String(nomadHome + "/build/release/interfaces/jNomad/");
			String os = System.getProperty("os.name");
			if (os.startsWith("Mac"))
				System.load(libpath + "libjNomad.jnilib");
			else if (os.startsWith("Windows"))
				System.load(libpath + "libjNomad.dll");
			else if (os.startsWith("Linux"))
				System.load(libpath + "libjNomad.so");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load.\n" + e);
			System.exit(1);
		}
	}

	public class My_Evaluator extends Evaluator {
		public My_Evaluator(AllParameters p) {
			super(p.getEvalParams(), EvalType::BB);
		}

		public boolean eval_x(EvalPoint x, Double h_max, boolean[] count_eval) {

			// f(x,y) = (1-x)^2 + 100*(y-x^2)^2
			double f = pow(1.0 - x.get(0).todouble(), 2.0)
					+ 100.0 * pow((x.get(1).todouble() - pow(x.get(0).todouble(), 2.0)), 2.0);
			String bbo = String.valueOf(f);

			// c1(x,y) = (x-1)^3-y+1
			double c1 = pow(x.get(0).todouble() - 1.0, 3.0) - x.get(1).todouble() + 1.0;
			bbo += " " + String.valueOf(c1);

			// c2(x,y) = x+y-2
			double c2 = x.get(0).todouble() + x.get(1).todouble() - 2.0;
			bbo += " " + String.valueOf(c2);

			x.setBBO(bbo);
			count_eval[0] = true;

			return true;
		}
	}

	public Rosenbrock(int argc, String argv[]) {
		try {

			/// The main step is in charge to solve the optimization problem
			MainStep mainstep = new MainStep();

			// The init command must be done first
			// !!!!! Important to do init FIRST (before p reading) to set the Locale to US
			// !!!!
			// Nomad does not support decimal separator other than "."
			mainstep.init();

			AllParameters p = new AllParameters();

			p.setAttributeValueSizeT("DIMENSION", 2);
			// p.readParamLine("DIMENSION 2");

			//  Lower bound
			p.readParamLine("LOWER_BOUND ( -1.5 -0.5 )");

			// Upper bound
			p.readParamLine("UPPER_BOUND ( 1.5 2.5 )");

			// X0
			p.readParamLine("X0 ( 0.1 0.1 )");

			// Max bb eval
			p.setAttributeValueSizeT("MAX_BB_EVAL", 200);

			// BB input type
			p.readParamLine("BB_INPUT_TYPE ( R R )");

			p.readParamLine("BB_OUTPUT_TYPE OBJ PB EB");

			// Display parameters
			p.setAttributeValueInt("DISPLAY_DEGREE", 2);
			p.setAttributeValueBool("DISPLAY_ALL_EVAL", true);
			p.readParamLine("DISPLAY_STATS EVAL ( SOL ) OBJ CONS_H H_MAX");

      // output set parameters
			// System.out.println(p.getSetAttributeAsString());

			p.checkAndComply();

			mainstep.setAllParameters(p);

			// custom evaluator creation:
			My_Evaluator ev = new My_Evaluator(p);
			mainstep.addEvaluator(ev);

			mainstep.start();
			mainstep.run();
			mainstep.end();

		} catch (RuntimeException e) {
			System.err.println("\nNOMAD has been interrupted (" + e.toString() + ")\n\n");
		}

	}

	public static void main(String args[]) {
		System.out.println("Starting test...");
		new Rosenbrock(args.length, args);
	}
}

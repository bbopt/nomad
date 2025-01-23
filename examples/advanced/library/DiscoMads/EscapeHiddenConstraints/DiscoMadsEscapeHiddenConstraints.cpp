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
/*--------------------------------------------------------------------------*/
/*  example of a program that uses DiscoMads to reveal and escape           */
/* hidden constraints regions of the space of variables                     */
/*--------------------------------------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Algos/Mads/MadsMegaIteration.hpp"
#include "Algos/Mads/SearchMethodAlgo.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/EvalSortType.hpp"
#include "Algos/AlgoStopReasons.hpp"
#include "Util/AllStopReasons.hpp"

// # The problem is described sec. 5.3 "5.3. Design of a styrene production process."
// of [1] and more detailed in section 2.5.4, p.61 of [2].
//
// IMPORTANT
//
// PB is Styrene run in hybrid library/batch mode.
// The styrene standalone executable is registered with the function call
// setAttributeValue("BB_EXE",xxxx)
// To run this optimization, the program must be executed in a path where styrene truth executable is available.
// Styrene sources are available at https://github.com/bbopt/styrene and must be compiled prior to run this optimization.

// [1] Escaping Unknown Discontinuous Regions in Blackbox Optimization
// Charles Audet, Alain Batailly, and Solène Kojtych, SIAM Journal on Optimization, 2022
// doi/10.1137/21M1420915
// [2] Contributions à l'optimisation de systèmes mécaniques non réguliers : reconception
// d'aubes de compresseur
// Solène Kojtych, Ph.D. thesis 2022
// doi/10.1137/21M1420915

void initAllParams(const std::shared_ptr<NOMAD::AllParameters>& allParams)
{
    const int n = 8;
    
    // parameters creation
    allParams->setAttributeValue("DIMENSION", n);
    allParams->setAttributeValue("BB_EXE", std::string("./truth.exe"));
    
    //------ General algorithm parameters
    // starting point and bounds
    std::vector<double> x0 = { 100, 87.82280202, 95.36797348, 0, 0, 49.04338841, 42.41599794, 41.01732603};
    allParams->setAttributeValue("X0", NOMAD::Point(x0) );
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, 0.0 ));
    allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(n, 100.0 ));

    // termination
    allParams->setAttributeValue("MAX_BB_EVAL", 1000);
    allParams->setAttributeValue("MIN_MESH_SIZE", NOMAD::ArrayOfDouble(n,1E-7));

    // display
    allParams->setAttributeValue("DISPLAY_DEGREE", 2);
    allParams->set_DISPLAY_ALL_EVAL(true);   
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe obj"));
    allParams->setAttributeValue("EVAL_STATS_FILE", std::string("statsEnd.txt"));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes= {NOMAD::BBOutputType::EB, NOMAD::BBOutputType::EB, NOMAD::BBOutputType::EB, NOMAD::BBOutputType::EB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::OBJ};   
    bbOutputTypes[11].setRevealingStatus(true); // only OBJ function is revealing for DiscoMads for this example
    
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    //------ DiscoMads parameters
    allParams->setAttributeValue("DISCO_MADS_OPTIMIZATION", true);    
    allParams->setAttributeValue("DISCO_MADS_HID_CONST", true);   // use DiscoMads to escape hidden constraints regions

    // remoteness to hidden constraints regions wished  
    NOMAD::Double exclusion_radius = 15;                                           
    allParams->setAttributeValue("DISCO_MADS_EXCLUSION_RADIUS", exclusion_radius);   
                
    // revealing poll
    allParams->setAttributeValue("DISCO_MADS_REVEALING_POLL_RADIUS",  NOMAD::Double(20));
    allParams->setAttributeValue("DISCO_MADS_REVEALING_POLL_NB_POINTS", n);

    // ------- Recommended parameters for DiscoMads

    // quad models are desactivated as they may be slow with DiscoMads
    allParams->getRunParams()->setAttributeValue("QUAD_MODEL_SEARCH", false);
    allParams->getRunParams()->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_2N);
    allParams->getEvaluatorControlParams()->setAttributeValue("EVAL_QUEUE_SORT",NOMAD::EvalSortType::DIR_LAST_SUCCESS);

    // ------- Specific parameters used in thesis [2]
    // Uncomment the following parameters to reproduce the thesis parameters
    //p.set_SEED(3698370); 
    //p.getRunParams()->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_NP1_NEG);  //comment previous "DIRECTION_TYPE" command
    //p.setAttributeValue( "ANISOTROPIC_MESH", false);
    //p.getRunParams()->setAttributeValue("NM_SEARCH", false);
    //p.setAttributeValue("SPECULATIVE_SEARCH", true);

    // Parameters validation
    allParams->checkAndComply();
    

}


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main()
{
    NOMAD::MainStep TheMainStep;
        
    // Set parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initAllParams(params);
    TheMainStep.setAllParameters(params);
    
    try
    {
        // Algorithm creation and execution
        TheMainStep.start();
        TheMainStep.run();
        TheMainStep.end();
    }

    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }
        
    return 0;
}

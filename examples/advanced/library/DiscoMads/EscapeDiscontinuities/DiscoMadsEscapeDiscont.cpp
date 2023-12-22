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
/* discontinuities in user-defined revealing blackbox outputs               */
/*--------------------------------------------------------------------------*/

#include "Nomad/nomad.hpp"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
// The problem is described sec. 5.1 "Results of a typical run" of [1] and more detailled 
// following eq. 2.52, p.53 in [2].

// [1] Escaping Unknown Discontinuous Regions in Blackbox Optimization
// Charles Audet, Alain Batailly, and Solène Kojtych, SIAM Journal on Optimization, 2022
// doi/10.1137/21M1420915
// [2] Contributions à l'optimisation de systèmes mécaniques non réguliers : reconception
// d'aubes de compresseur
// Solène Kojtych, Ph.D. thesis 2022
// doi/10.1137/21M1420915


class My_Evaluator : public NOMAD::Evaluator
{
public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams, NOMAD::EvalType evalType)
        : NOMAD::Evaluator(evalParams, evalType)
    {}

    ~My_Evaluator() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double& hMax, bool &countEval) const override
    {
    bool eval_ok = false;
    NOMAD::Double f = 1e+20, c1 = 1e+20;
    size_t n = x.size();

        try
        {
            double x0_f = 0;
            double y0_f = 10;
            double radius_f = 12;

            // Objective function
            if(pow(x[0].todouble()-x0_f,2)+pow(x[1].todouble()-y0_f,2)>pow(radius_f,2)){
                f=-0.025*x[1].todouble()+3;
            }
            else{
                f= 0.04*x[1].todouble();
            }
                
            // Constraint
            if (x[0]>0){
                c1=2;
            }
            else{
                c1=-0.1;
            }

            std::string bbo = f.tostring();
            bbo += " " + c1.tostring();
            x.setBBO(bbo);

            eval_ok=true;
            countEval=true;
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }
    

    return eval_ok;
    }

};




void initParams(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 2;   // Number of variables
    p.getPbParams()->setAttributeValue("DIMENSION", n);
    p.getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("OBJ-R PB-R"));

    //------ General algorithm parameters
    // starting point and bounds
    NOMAD::Point X0(n, -5.0);
    p.getPbParams()->setAttributeValue("X0", X0);
    p.getPbParams()->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n,-10.0));
    p.getPbParams()->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(n, 10.0));

    // the algorithm terminates after MAX_BB_EVAL black-box evaluations or if the min mesh size is reached
    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_BB_EVAL", 2000);
    p.setAttributeValue("MIN_MESH_SIZE", NOMAD::ArrayOfDouble(n,1E-9));

    // display
    p.getDispParams()->setAttributeValue("DISPLAY_DEGREE", 2);
    p.set_DISPLAY_ALL_EVAL(true);   
    p.setAttributeValue("EVAL_STATS_FILE", std::string("statsEnd.txt"));

    //------ DiscoMads parameters
    p.setAttributeValue("DISCO_MADS_OPTIMIZATION", true);    // use discoMads to reveal and escape discontinuities

    // definition of discontinuities (in the weak sense)
    NOMAD::Double detection_radius = 0.25;
    p.setAttributeValue("DISCO_MADS_DETECTION_RADIUS", detection_radius);
    p.setAttributeValue("DISCO_MADS_LIMIT_RATE", NOMAD::Double(0.3));

    // remoteness to discontinuities wished  
    NOMAD::Double exclusion_radius = 0.25;                                           
    p.setAttributeValue("DISCO_MADS_EXCLUSION_RADIUS", exclusion_radius);   
                
    // revealing poll
    p.setAttributeValue("DISCO_MADS_REVEALING_POLL_RADIUS",  NOMAD::Double(1.01*(detection_radius+exclusion_radius)));
    p.setAttributeValue("DISCO_MADS_REVEALING_POLL_NB_POINTS", n);

    // ------- Recommended parameters for DiscoMads
    // no parallelism
    p.setAttributeValue("NB_THREADS_OPENMP",1); // DiscoMads works with OpenMP but has not been extensively tested

    // quad models are desactivated as they may be slow with DiscoMads
    p.getRunParams()->setAttributeValue("QUAD_MODEL_SEARCH", false);
    p.getRunParams()->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_2N);
    p.getEvaluatorControlParams()->setAttributeValue("EVAL_QUEUE_SORT",NOMAD::EvalSortType::DIR_LAST_SUCCESS);


    // ------- Specific thesis parameters
    // Uncomment the following parameters to reproduce the thesis parameters
    //p.set_SEED(3698370); 
    //p.getRunParams()->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_NP1_NEG);  //comment previous "DIRECTION_TYPE" command
    //p.setAttributeValue( "ANISOTROPIC_MESH", false);
    //p.getRunParams()->setAttributeValue("NM_SEARCH", false);
    //p.setAttributeValue("SPECULATIVE_SEARCH", true);


    // parameters validation
    p.checkAndComply();
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main (int argc, char **argv)
{
    auto TheMainStep = std::make_unique<NOMAD::MainStep>();

    // Initialize all parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initParams(*params);
    TheMainStep->setAllParameters(params);

    // Custom BB evaluator creation
    auto evBB = std::make_unique<My_Evaluator>(params->getEvalParams(),NOMAD::EvalType::BB);
    TheMainStep->addEvaluator(std::move(evBB));

    try
    {
        // Algorithm creation and execution
        TheMainStep->start();
        TheMainStep->run();
        TheMainStep->end();
    }

    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    return 0;
}

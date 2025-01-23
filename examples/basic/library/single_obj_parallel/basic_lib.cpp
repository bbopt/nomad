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
/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "Nomad/nomad.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

// Number of threads to be used for evaluation of a block of points in parallel
#define NUM_THREADS_EVAL  6

// Wrapper of eval_x used for parallel evaluation.
bool wrapper_eval_x(std::shared_ptr<NOMAD::EvalPoint> & x, const NOMAD::Double& hmax, bool & countEval)
{

    // std::cout << "Thread #" << omp_get_thread_num() << std::endl;
    
    NOMAD::Double c1 = 0.0, c2 = 0.0;

    for (int i = 0; i < 5; i++)
    {
        c1 += ((*x)[i]-1).pow2();
        c2 += ((*x)[i]+1).pow2();
    }
    NOMAD::Double f = (*x)[4]; // objective value
    c1 = c1-25;             // constraint 1
    c2 = 25-c2;             // constraint 2
    std::string bbo = f.tostring() + " " + c1.tostring() + " " + c2.tostring();
    x->setBBO(bbo);

    countEval = true; // count a black-box evaluation

    return true;    // Eval ok
}


class My_Evaluator : public NOMAD::Evaluator
{
public:
    explicit My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
        : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() override = default;

    // Implementation
    std::vector<bool> eval_block(std::vector<std::shared_ptr<NOMAD::EvalPoint>> &block,
                const NOMAD::Double& hMax,
                std::vector<bool> &listCountEval) const override

    {
        std::vector<bool> evalOk(block.size(), false);
        listCountEval.resize(block.size(), false);  //  Evaluations are not counted until eval_x is called and sets countEval
        
        int i;
#ifdef _OPENMP
        #pragma omp parallel for num_threads(NUM_THREADS_EVAL) shared(block,evalOk,listCountEval, hMax) private(i)
#endif
        for(i = 0; i < (int)block.size(); i++)
        {
            bool countEval = false;
            evalOk[i] = wrapper_eval_x(block[i],hMax, countEval);
            listCountEval[i] = countEval;
        }

        return evalOk;
    }
};


void initParams(const std::shared_ptr<NOMAD::AllParameters>& params)
{
    params->setAttributeValue("DIMENSION", 5);             // number of variables

    NOMAD::BBOutputTypeList bbot;   // Definition of output types
    bbot.emplace_back(NOMAD::BBOutputType::OBJ);
    bbot.emplace_back(NOMAD::BBOutputType::PB);
    bbot.emplace_back(NOMAD::BBOutputType::EB);
    params->setAttributeValue("BB_OUTPUT_TYPE", bbot);

    // params->setAttributeValue("DISPLAY_ALL_EVAL", true);   // displays all evaluations.
    params->setAttributeValue("DISPLAY_DEGREE", 2);
    params->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("BBE BLK_EVA BLK_SIZE ( SOL ) OBJ"));

    // Disable searches because they generate 1 or 2 points at a time. We want to fill the blocks.
    params->setAttributeValue("SPECULATIVE_SEARCH", false);
    params->setAttributeValue("NM_SEARCH", false);
    params->setAttributeValue("QUAD_MODEL_SEARCH", false);
    // Use poll method n+1 uni -> 6 trial points per poll on a single pass
    params->setAttributeValue("DIRECTION_TYPE",NOMAD::DirectionType::NP1_UNI);

    params->setAttributeValue("X0", NOMAD::Point(5,0.0));  // starting point

    params->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(5, -6.0)); // all var. >= -6
    NOMAD::ArrayOfDouble ub(5);                    // x_4 and x_5 have no bounds
    ub[0] = 5.0;                              // x_1 <= 5
    ub[1] = 6.0;                              // x_2 <= 6
    ub[2] = 7.0;                              // x_3 <= 7
    params->setAttributeValue("UPPER_BOUND", ub);

    params->setAttributeValue("MAX_BB_EVAL", 100);     // the algorithm terminates after
    // 100 black-box evaluations

    // Max number of points to be given as a block for evaluation
    // This option is required to perform parallel evaluations
    // All points in a block will be evaluated in parallel
    params->setAttributeValue("BB_MAX_BLOCK_SIZE", NUM_THREADS_EVAL);
    
    // NOTE: A single thread is used for Nomad "parallel" evaluation queue.

    // parameters validation:
    params->checkAndComply();
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main(int argc, char ** argv)
{
    try
    {

#ifdef _OPENMP
        // This is important for nested parallel regions. Nomad parallel region
        // for evaluation queue and a parallel region for block evaluation.
        omp_set_nested(true);
        omp_set_dynamic(false);
#endif
        
        auto TheMainStep = std::make_unique<NOMAD::MainStep>();

        // Parameters creation
        auto params = std::make_shared<NOMAD::AllParameters>();
        initParams(params);
        TheMainStep->setAllParameters(params);
        
        // Custom Evaluator
        auto ev = std::make_unique<My_Evaluator>(params->getEvalParams());
        TheMainStep->addEvaluator(std::move(ev));

        // Algorithm creation and execution
        TheMainStep->start();
        TheMainStep->run();
        TheMainStep->end();
    }

    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    return EXIT_SUCCESS;
}

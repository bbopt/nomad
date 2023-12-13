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

// Number of threads to be used in parallel
#define NUM_THREADS   6 

// A structure to pass arguments to the evaluation wrapper function
class My_Evaluator;
typedef struct Arg_Eval_tag {
    std::shared_ptr<NOMAD::EvalPoint> x;
    NOMAD::Double hMax;
    bool countEval;
} Arg_Eval_t;


// Wrapper of eval_x used for parallel evaluation.
bool wrapper_eval_x(Arg_Eval_t& eval_arg)
{
    NOMAD::Double c1 = 0.0, c2 = 0.0;
    NOMAD::EvalPoint& x = *(eval_arg.x);

    for (int i = 0; i < 5; i++)
    {
        c1 += (x[i]-1).pow2();
        c2 += (x[i]+1).pow2();
    }
    NOMAD::Double f = x[4]; // objective value
    c1 = c1-25;             // constraint 1
    c2 = 25-c2;             // constraint 2
    std::string bbo = f.tostring() + " " + c1.tostring() + " " + c2.tostring();
    eval_arg.x->setBBO(bbo);

    eval_arg.countEval = true; // count a black-box evaluation

    return true;    // Eval ok
}



class My_Evaluator : public NOMAD::Evaluator
{
public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
        : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() {}

    // Implementation
    std::vector<bool> eval_block(std::vector<std::shared_ptr<NOMAD::EvalPoint>> &block,
                const NOMAD::Double& hMax,
                std::vector<bool> &listCountEval) const override

    {
        std::vector<bool> evalOk(block.size(), false);
        listCountEval.resize(block.size(), false);  //  Evaluations are not counted until eval_x is called and sets countEval

        // Arguments passed to the evaluation wrapper
        Arg_Eval_t* eval_arg = new Arg_Eval_t[block.size()];

        std::vector<std::shared_ptr<NOMAD::EvalPoint>>::iterator it_x = block.begin(), end_x = block.end();
        size_t i = 0;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        for (it_x = block.begin(); it_x != end_x; ++it_x, ++i)
        {
            // A thread is created for each x in the block.
            eval_arg[i].x = (*it_x);
            eval_arg[i].hMax = hMax.todouble();
            evalOk[i] = wrapper_eval_x(eval_arg[i]);
            listCountEval[i] = eval_arg[i].countEval;
        }
        // End of parallel block

        delete[] eval_arg;

        return evalOk;
    }
};


void initParams(std::shared_ptr<NOMAD::AllParameters>& params)
{
    params->setAttributeValue("DIMENSION", 5);             // number of variables

    NOMAD::BBOutputTypeList bbot;   // Definition of output types
    bbot.push_back(NOMAD::BBOutputType::OBJ);
    bbot.push_back(NOMAD::BBOutputType::PB);
    bbot.push_back(NOMAD::BBOutputType::EB);
    params->setAttributeValue("BB_OUTPUT_TYPE", bbot);

    // params->setAttributeValue("DISPLAY_ALL_EVAL", true);   // displays all evaluations.
    params->setAttributeValue("DISPLAY_DEGREE", 2);
    params->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("BBE BLK_EVA BLK_SIZE ( SOL ) OBJ"));

    // Disable searches because they generate 1 or 2 points at a time. We want to fill the blocks.
    params->setAttributeValue("SPECULATIVE_SEARCH", false);
    params->setAttributeValue("NM_SEARCH", false);
    params->setAttributeValue("QUAD_MODEL_SEARCH", false);
    // Use single pass poll method with n+1 points -> 6 trial points in a block per poll 
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
    // This option is required to perform parallel evaluations in eval_block
    // function above
    params->setAttributeValue("BB_MAX_BLOCK_SIZE", NUM_THREADS);
    
    // A single thread is used for Nomad "parallel" evaluation queue.
    params->setAttributeValue("NB_THREADS_OPENMP", 1);

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

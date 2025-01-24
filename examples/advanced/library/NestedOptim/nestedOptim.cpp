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
/*--------------------------------------------*/
/*  Rosenbrock with dimension greater than 2  */
/*  Inner problem has dimension 2             */
/*  Inner problem is solved with SimpleMads   */
/*  Outer problem is solved with regular Mads*/
/*--------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Algos/SimpleMads/SimpleMads.hpp"
#include "Type/DirectionType.hpp"
#include "Type/EvalSortType.hpp"

const size_t Nout = 10; // 10 variables outer problem, 2 additional (hidden) variables are used in the problem

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_EvaluatorOut : public NOMAD::Evaluator
{
private:
    const NOMAD::Step * _genStep;

public:
    explicit My_EvaluatorOut(const std::shared_ptr<NOMAD::EvalParameters>& evalParams, const NOMAD::Step * genStep )
    : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB),
      _genStep(genStep)
    {}

    ~My_EvaluatorOut() override = default;

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override;
};


/*----------------------------------------*/
/*           user-defined eval_x          */
/*----------------------------------------*/
bool My_EvaluatorOut::eval_x(NOMAD::EvalPoint &x,
                          const NOMAD::Double &hMax,
                          bool &countEval) const
{
    if (Nout%2 != 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Dimension N should be an even number");
    }
    
    //
    // Suboptimization on 2 variables. x (dim=10) is passed as a fixed variables to inner evaluation function
    //
    
    // The function to evaluate an eval point.
    // Outer variables ->x
    std::function<bool(std::vector<NOMAD::SimpleEvalPoint>&)> eval_x = [&x](std::vector<NOMAD::SimpleEvalPoint>& allXIn) -> bool
    {
        for (auto & xIn: allXIn)
        {
            double f=0;
            // outer variables (dim 10) contribution to f
            for ( size_t i = 1 ; i <= Nout/2 ; ++i ) {
                f += pow ( 10 * (x[2*i-1].todouble() - pow(x[2*i-2].todouble(),2) ) , 2 );
                f += pow ( 1 - x[2*i-2].todouble() , 2 );
            }
            // inner variables (dim 2) contribution to f
            f += pow ( 10 * (xIn[1].todouble() - pow(xIn[0].todouble(),2) ) , 2 );
            f += pow ( 1 - xIn[0].todouble() , 2 );
            
            xIn.setF(f);
            xIn.setH(0);
        }
        return true;
    };

    // Inner Sub-Pb parameters
    auto pbParams = std::make_shared<NOMAD::PbParameters>();
    pbParams->setAttributeValue("DIMENSION",2);
    pbParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(2,-5)); // lb = (-5,-5)
    pbParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(2,5));   // ub = (5,5)
    NOMAD::ArrayOfPoint x0s{NOMAD::ArrayOfDouble(2,0)}; // A single X0 = (0,0)
    pbParams->setAttributeValue("X0", x0s);
    pbParams->checkAndComply();
    
    // Run parameters (use default)
    auto runParams = std::make_shared<NOMAD::RunParameters>();
    runParams->setAttributeValue("ANISOTROPIC_MESH",false);
    auto evcParams = NOMAD::EvcInterface::getEvaluatorControl()->getEvaluatorControlGlobalParams(); 
    runParams->checkAndComply(evcParams, pbParams);

    NOMAD::BBOutputTypeList bbot = {NOMAD::BBOutputType::Type::OBJ};

    auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

    // Create a simple mads to solve the inner problem on 2 dimension.
    // Outer variables -> x are available as fixed parameters in the evaluation.
    // Simple Mads is not opportunistic. All poll points are passed at once into eval_x.
    NOMAD::SimpleMads mads(_genStep, madsStopReasons, runParams, pbParams, bbot, eval_x /* eval_x for a block of points */, 1600 /* maxEval */);
    
    // No need to display something at the end of the sub-optimization
    mads.setEndDisplay(false);

    mads.start();
    bool runOk = mads.run();
    mads.end();

    double f = 0;
    if (!runOk)
    {
        std::cout << "Pb with inner mads. Let's continue. f=0." << std::endl;
    }
    else
    {
        // Get the best feas solution
        f = mads.getBestSimpleSolution(true).getF().todouble();
    }
    
    x.setBBO(std::to_string(f));

    countEval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
}


void initAllParams(const std::shared_ptr<NOMAD::AllParameters>& allParams)
{
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", Nout);
    // 100 black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 400*Nout);

    // Starting point
    allParams->setAttributeValue("X0", NOMAD::Point(Nout, 0.5) );

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(Nout, -10.0 ));
    allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(Nout, 10.0 ));
 
    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes = {NOMAD::BBOutputType::OBJ};
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    allParams->setAttributeValue("DISPLAY_DEGREE", 2);
    
    NOMAD::ArrayOfString ds("BBE ( SOL ) OBJ");
    allParams->setAttributeValue("DISPLAY_STATS", ds);

    // Parameters validation
    allParams->checkAndComply();
}


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main()
{
    NOMAD::MainStep TheMainStep;

    auto params = std::make_shared<NOMAD::AllParameters>();
    initAllParams(params);
    TheMainStep.setAllParameters(params);

    // Custom Evaluator
    auto ev = std::make_unique<My_EvaluatorOut>(params->getEvalParams(), &TheMainStep);
    TheMainStep.addEvaluator(std::move(ev));

    try
    {
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

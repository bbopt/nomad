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

#include "Nomad/nomad.hpp"

/*----------------------------------------*/
/*          A geometric problem           */
/* Original: min x1^2+...+x5^2            */
/*           st x1+...+x5 = 10            */
/*              0 <= xi <=5               */
/* Solution is xi=2                       */
/* Recall that Nomad cannot manage        */
/* equality constraint.                   */
/*                                        */
/* Modified problem:                      */
/* Pb dimension is set to n=4             */
/* Set an inequality constraint:          */
/*         x1+...+x4<=10                  */
/* If constraint is verified              */
/*   - Pick up   d = 10-(x1+...+x4)       */
/*   - Compute f=-(x1^2+...+x4^2+d^2)     */
/*   - Count eval                         */
/* If constraint is not verified          */
/*   - f=Inf                              */
/*   - Constraint is EB constraint        */
/*   - Geometric constraint is not costly */
/*   - Do not count eval                  */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
public:
    explicit My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
        : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() override = default;

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double& hMax, bool &countEval) const override
    {
        size_t n = x.size();

        NOMAD::Double f = 0.0;  // Objective value
        NOMAD::Double s1 = 0; // Sum of x1,...,x4
        NOMAD::Double c1 = 0.0; // Constraint 1 -> geometric constraint

        try
        {
            
            for (size_t i = 0; i < n; i++)
            {
                // Let's suppose this geometric constraint is not costly to compute but f is!
                s1 += x[i];
                if (s1 > 10)
                {
                    c1 = s1;
                    std::string bbo = f.tostring() + " " + c1.tostring(); // f is not really computed. But the point is infeasible, and we handle the constraint with EB, so it does not matter. The point is simply discarded.
                    
                    x.setBBO(bbo);
                    countEval = false; // DO NOT count as a blackbox evaluation when geometric constraint is not verified
                    
                    // the evaluation succeeded
                    return true;
                }
                
            }

            if (s1 <= 10)
            {
                c1 = s1 - 10 ; // if s1 <= 10 we can get x1+...+x4 - 10 <= 0 -> Feasible !
                NOMAD::Double d = 10.0 - s1; // With this value of d, we have x1+...+x4+d = 10
                
                // Compute the objective. Let's pretend it is costly!
                f =  d*d;
                for (size_t i = 0; i < n; i++)
                {
                    f += x[i]*x[i];
                }
                std::string bbo = f.tostring() + " " + c1.tostring() ;
                // Simple way to set the BBO
                x.setBBO(bbo);
                countEval = true; // count a black-box evaluation
            }
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }

        return true;     // the evaluation succeeded
    }
};


void initParams(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 4;   // Number of variables of the modified problem
    p.setAttributeValue("DIMENSION", n);
    p.setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("OBJ EB")); // EB constraint: If a point is infeasible it is simply discarded. F (costly part) is not computed.

    p.setAttributeValue("X0", NOMAD::Point(n,5.0));  // starting point (0.0 0.0 0.0 0.0 0.0)
    p.setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, 0.0)); // all var. >= 0
    p.setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(n, 5.0)); // all var. >= 0);

    // the algorithm terminates after 100 black-box evaluations,
    // or 10000 total evaluations, including cache hits and evaluations for
    // which countEval was false.
    p.setAttributeValue("MAX_BB_EVAL", 200);
    p.setAttributeValue("MAX_EVAL", 10000);

    p.setAttributeValue("DISPLAY_DEGREE", 2);
    p.setAttributeValue("DISPLAY_ALL_EVAL", true);
    p.setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("EVAL ( SOL ) OBJ CONS_H"));

    // parameters validation
    p.checkAndComply();
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main()
{
    auto TheMainStep = std::make_unique<NOMAD::MainStep>();

    // Initialize all parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initParams(*params);
    TheMainStep->setAllParameters(params);

    // Custom evaluator creation
    auto ev = std::make_unique<My_Evaluator>(params->getEvalParams());
    TheMainStep->addEvaluator(std::move(ev));

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

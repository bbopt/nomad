/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
        : NOMAD::Evaluator(evalParams)
    {}

    ~My_Evaluator() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double& hMax, bool &countEval) const override
    {
        bool eval_ok = false;
        size_t n = x.size();

        NOMAD::Double f = 0.0;  // Objective value
        NOMAD::Double c1 = 0.0; // Constraint 1
        NOMAD::Double c2 = 0.0; // Constraint 2

        try
        {
            for (size_t i = 0; i < n; i++)
            {
                NOMAD::Double xi = x[i];
                c1 += (xi-1).pow2();
                c2 += (xi+1).pow2();
            }
            c1 = c1-25;
            c2 = 25-c2;

            f = x[n-1];
            std::string bbo = f.tostring() + " " + c1.tostring() + " " + c2.tostring();
            x.setBBO(bbo, _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"));
            eval_ok = true;
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }


        countEval = true;  // count a black-box evaluation

        return eval_ok;     // the evaluation succeeded
    }
};




void initParams(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 6;   // Number of variables
    p.getPbParams()->setAttributeValue("DIMENSION", n);
    p.getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("OBJ PB PB"));

    NOMAD::Point X0(n, 0.0);
    X0[n-1] = -6.0; // starting point (0.0 0.0 0.0 0.0 0.0 -6.0)
    p.getPbParams()->setAttributeValue("X0", X0);
    p.getPbParams()->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -6.0)); // all var. >= -6
    NOMAD::ArrayOfDouble ub(n);     // x_4 and x_5 have no bounds
    ub[0] = 5.0;                    // x_1 <= 5
    ub[1] = 6.0;                    // x_2 <= 6
    ub[2] = 7.0;                    // x_3 <= 7
    ub[n-1] = 6.0;                  // x_6 <= 6
    p.getPbParams()->setAttributeValue("UPPER_BOUND", ub);

    // the algorithm terminates after MAX_BB_EVAL black-box evaluations, or MAX_EVAL total evaluations (including cache hits).
    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_BB_EVAL", 1000);
    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_EVAL", 1000);

    p.getDispParams()->setAttributeValue("DISPLAY_DEGREE", 2);
    p.getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("EVAL ( SOL ) OBJ ( BBO ) CONS_H"));

    p.getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    p.getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);

    // Do not read the cache when starting.
    // We do this because it is an example. In larger problems, an input cache file
    // is usually useful.
    p.getCacheParams()->setAttributeValue("CACHE_FILE", std::string(""));

    // parameters validation
    p.checkAndComply();
}


/// An example of custom H component computation.
/// Here the bound for constraint bbo[1] is set to 0.3 instead of the default 0.
/// WARNING: manipulating H component computation for EB constraint must be done carefully. When PhaseOneSearch is performed (EB constraint not feasible at X0), it will try to reduce EB to zero before using the custom H computation for EB.
NOMAD::Double customComputeHComponent( const NOMAD::BBOutputType & bbOutputType , size_t bboIndex, const NOMAD::Double &bbo )
{
    NOMAD::Double hComp = 0.0 ;

    /// Default value for the bound of a constraint (default not relaxed = 0)
    /// Must greater than zero
    NOMAD::Double relaxedConstBound = 0.0;

    /// Set a relaxed bound for constraint given in bbo[1]
    if ( bboIndex == 1 )
    {
        relaxedConstBound = 1;
    }

    if ( NOMAD::BBOutputType::EB == bbOutputType)
    {
        if ( bbo > 0)
        {
            hComp = NOMAD::INF;
        }
    }
    else if (NOMAD::BBOutputType::PB == bbOutputType)
    {
        if ( bbo > relaxedConstBound )
        {
            hComp = (bbo-relaxedConstBound) * (bbo-relaxedConstBound ) ;
        }
    }

    return hComp;
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

    // Custom evaluator creation
    std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
    TheMainStep->setEvaluator(std::move(ev));

    /// Set the use of custom compute of H (infeasibility measure) components (constraint by contraint)
    NOMAD::Eval::setComputeHComponentFunction(customComputeHComponent);

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

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

#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Cache/CacheBase.hpp"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
    : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override
    {
        bool eval_ok = false;
        // Based on G2.
        NOMAD::Double f, c1 = 0, c2 = 0;
        size_t n = x.size();

        try
        {
            for (size_t i = 0; i < 5 ; i++)
            {
                c1 += pow ( x[i].todouble() -1 , 2 );
                    c2 += pow ( x[i].todouble()+1 , 2 );
            }
            f = x[4];
            c1 = c1 - 25.0;
            c2 = 25.0 - c2;

            std::string bbo = f.tostring();
            bbo += " " + c1.tostring();
            bbo += " " + c2.tostring();
            
            x.setBBO(bbo);

            eval_ok = true;
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }

        countEval = true;
        return eval_ok;
    }
};


void initParams1(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 5;   // Number of variables
    p.getPbParams()->setAttributeValue("DIMENSION", n);
    
    // When using batch mode decomment this and remove My_Evaluator
    // p.getEvalParams()->setAttributeValue("BB_EXE", std::string("./u.exe"));

    NOMAD::Point x0(n, 0.0);
    p.getPbParams()->setAttributeValue("X0", x0);  // starting point * 0
    NOMAD::ArrayOfDouble lowerBound(n,-6.0);
    p.getPbParams()->setAttributeValue("LOWER_BOUND", lowerBound); // * -6.0
    NOMAD::ArrayOfDouble upperBound(n);
    upperBound[0] = 5.0;
    upperBound[1] = 6.0;
    upperBound[2] = 7.0;
    p.getPbParams()->setAttributeValue("UPPER_BOUND", upperBound);  // ( 5 6 7 - - )

    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_EVAL", size_t(1000));
    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_BB_EVAL", size_t(1000));

    p.getDispParams()->setAttributeValue("DISPLAY_DEGREE", 2);
    p.getDispParams()->setAttributeValue("DISPLAY_ALL_EVAL", true);
    p.getDispParams()->setAttributeValue("DISPLAY_INFEASIBLE", true);
    p.getDispParams()->setAttributeValue("DISPLAY_UNSUCCESSFUL", false);
    p.getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("BBE ( SOL ) OBJ"));


    p.getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("OBJ PB PB"));

    p.getEvaluatorControlParams()->setAttributeValue("EVAL_OPPORTUNISTIC", false);
    p.getRunParams()->setAttributeValue("NM_OPTIMIZATION",true);


    // parameters validation
    p.checkAndComply();
}


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main (int argc, char **argv)
{
    auto TheMainStep = std::make_unique<NOMAD::MainStep>();

    // Initialize parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initParams1(*params);
    
    TheMainStep->setAllParameters(params);

    std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
    
    TheMainStep->setEvaluator(std::move(ev));
    
    try
    {
        // Algorithm creation and execution
        TheMainStep->start();
        // Here, clear the cache, ensure we start the test fresh.
        NOMAD::CacheBase::getInstance()->clear();
        TheMainStep->run();
        TheMainStep->end();
    }

    catch(std::exception &e)
    {
        std::cerr << "\nRun has been interrupted (" << e.what() << ")\n\n";
    }

    return 0;
}

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
/**
 \file   example1_lib.cpp
 \brief  Library example for nomad
 \author Viviane Rochon Montplaisir
 \date   2017
 */

#include "Nomad/nomad.hpp"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
public:
    explicit My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
    : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB) // Evaluator for true blackbox evaluations only
    {}

    ~My_Evaluator() override = default;

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override
    {
        bool eval_ok = false;
        // Based on G2.
        NOMAD::Double f = 1e+20, g1 = 1e+20, g2 = 1e+20;
        NOMAD::Double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, prod1 = 1.0, prod2 = 1.0;
        size_t n = x.size();

        try
        {
            for (size_t i = 0; i < n ; i++)
            {
                sum1  += pow(cos(x[i].todouble()), 4);
                sum2  += x[i];
                sum3  += (double)(i+1)*x[i]*x[i];
                prod1 *= pow(cos(x[i].todouble()), 2);
                if (prod2 != 0.0)
                {
                    if (x[i] == 0.0)
                    {
                        prod2 = 0.0;
                    }
                    else
                    {
                        prod2 *= x[i];
                    }
                }
            }

            g1 = -prod2 + 0.75;
            g2 = sum2 -7.5 * (double)n;

            f = 10*g1 + 10*g2;
            if (0.0 != sum3)
            {
                f -= ((sum1 -2*prod1) / sum3.sqrt()).abs();
            }
            // Scale
            if (f.isDefined())
            {
                f *= 1e-5;
            }

            NOMAD::Double c2000 = -f-2000;
            // Double::tostring function uses FULL precision
            std::string bbo = g1.tostring();
            bbo += " " + g2.tostring();
            bbo += " " + f.tostring();
            bbo += " " + c2000.tostring();

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

// Uncomment when testing with failed eval simulation (see display parameters below)
//        // Simulate failed evaluations
//        auto bbe = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
//        if(bbe>=9)
//        {
//            //Failed eval
//            eval_ok=false;
//        }

        return eval_ok;
    }
};


void initAllParams(const std::shared_ptr<NOMAD::AllParameters>& allParams)
{
    // Parameters creation
    // Number of variables
    size_t n = 10;
    allParams->setAttributeValue( "DIMENSION", n);
    // The algorithm terminates after
    // this number of black-box evaluations
    allParams->setAttributeValue( "MAX_BB_EVAL", 1000);
    // Starting point
    allParams->setAttributeValue( "X0", NOMAD::Point(n, 7.0) );

    allParams->getPbParams()->setAttributeValue("GRANULARITY", NOMAD::ArrayOfDouble(n, 0.0000001));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes;
    bbOutputTypes.emplace_back(NOMAD::BBOutputType::PB);     // g1
    bbOutputTypes.emplace_back(NOMAD::BBOutputType::PB);     // g2
    bbOutputTypes.emplace_back(NOMAD::BBOutputType::OBJ);    // f
    bbOutputTypes.emplace_back(NOMAD::BBOutputType::EB);     // c2000
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );
    allParams->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_2N);
    allParams->setAttributeValue("DISPLAY_DEGREE", 2);

// Uncomment when testing with failed eval simulation (see eval_x function)
//    allParams->set_DISPLAY_ALL_EVAL(true);
//    allParams->setAttributeValue("DISPLAY_FAILED", true);
//    allParams->setAttributeValue("DISPLAY_UNSUCCESSFUL", true);
//    allParams->setAttributeValue("STATS_FILE", NOMAD::ArrayOfString("stats.txt bbe obj")); //"stats.txt obj mesh_size success_type"
//    allParams->setAttributeValue("EVAL_STATS_FILE", std::string("statsEnd.txt"));

    // Parameters validation
    allParams->checkAndComply();

}


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main()
{
    NOMAD::MainStep TheMainStep;

    try
    {
        // Set parameters
        auto params = std::make_shared<NOMAD::AllParameters>();
        initAllParams(params);
        TheMainStep.setAllParameters(params);

        // Set evaluator
        auto ev = std::make_unique<My_Evaluator>(params->getEvalParams());
        TheMainStep.addEvaluator(std::move(ev));

        // Run optimization
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


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
/*--------------------------------------------*/
#include "Nomad/nomad.hpp"

const size_t N = 50;

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
private:

public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
    : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override;
};


/*----------------------------------------*/
/*           user-defined eval_x          */
/*----------------------------------------*/
bool My_Evaluator::eval_x(NOMAD::EvalPoint &x,
                          const NOMAD::Double &hMax,
                          bool &countEval) const
{
    
    if (N%2 != 0)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Dimension N should be an even number");
    }
    
    double f=0;
    for ( size_t i = 1 ; i <= N/2 ; ++i ) {
        f += pow ( 10 * (x[2*i-1].todouble() - pow(x[2*i-2].todouble(),2) ) , 2 );
        f += pow ( 1 - x[2*i-2].todouble() , 2 );
    }
    NOMAD::BBOutputTypeList bbOutputTypes ={NOMAD::BBOutputType::OBJ};
    x.setBBO(std::to_string(f), bbOutputTypes);

    countEval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
}


void initAllParams(std::shared_ptr<NOMAD::AllParameters> allParams)
{
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", N);
    // 100 black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 100*N);
    // Starting point
    allParams->setAttributeValue("X0", NOMAD::Point(N, 0.5) );

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(N, -10.0 ));
    allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(N, 10.0 ));

    
    allParams->setAttributeValue("NM_SEARCH",false);
    allParams->setAttributeValue("PSD_MADS_OPTIMIZATION",true);
    allParams->setAttributeValue("PSD_MADS_NB_VAR_IN_SUBPROBLEM",2);
    allParams->setAttributeValue("NB_THREADS_OPENMP",4);
    
    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes = {NOMAD::BBOutputType::OBJ};
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    allParams->setAttributeValue("DISPLAY_DEGREE", 0);
    allParams->setAttributeValue("STATS_FILE", NOMAD::ArrayOfString("stats.txt bbe obj"));

    // Parameters validation
    allParams->checkAndComply();

}


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{

    NOMAD::MainStep TheMainStep;

    auto params = std::make_shared<NOMAD::AllParameters>();
    initAllParams(params);
    TheMainStep.setAllParameters(params);

    std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
    TheMainStep.setEvaluator(std::move(ev));

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

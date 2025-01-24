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
/*  Example of a program that makes NOMAD do a Mads custom user search in   */
/*  addition to the default search.                                         */
/*--------------------------------------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Algos/Mads/Mads.hpp"
#include "Algos/Mads/MadsMegaIteration.hpp"
#include "Algos/Mads/SearchMethodAlgo.hpp"
#include "Algos/SubproblemManager.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/EvalSortType.hpp"
#include "Algos/AlgoStopReasons.hpp"
#include "Util/AllStopReasons.hpp"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
const int N=6;

class My_Evaluator : public NOMAD::Evaluator
{
private:

public:
    explicit My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
    : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() override = default;

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
    NOMAD::Double F(f);
    x.setBBO(F.tostring());
    countEval = true;

    return true;       // the evaluation succeeded
}


void initAllParams(const std::shared_ptr<NOMAD::AllParameters>& allParams)
{
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", N);
    // 100 black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 50*N);
    // Starting point
    allParams->setAttributeValue("X0", NOMAD::Point(N, 0.0) );

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(N, -10.0));
    allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(N, 10.0));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes;
    bbOutputTypes.emplace_back(NOMAD::BBOutputType::OBJ);
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    // Algo
    allParams->setAttributeValue("DIRECTION_TYPE",NOMAD::DirectionType::ORTHO_2N);
    allParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    allParams->setAttributeValue("NM_SEARCH", false);

    // Enable the user search method. See below for registering the callback function
    // to generate the search directions. Both are required.
    allParams->setAttributeValue("USER_SEARCH", true);

    // Display
    allParams->setAttributeValue("DISPLAY_DEGREE", 3);
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj"));
    allParams->setAttributeValue("DISPLAY_ALL_EVAL", true);

    // Parameters validation
    allParams->checkAndComply();
}

// The function to generate user search directions. This is registered as a callback below.
bool userSearchMethodCallback(const NOMAD::Step& step, NOMAD::EvalPointSet & trialPoints)
{
    // Important: by default USER_CALLS are disabled when doing quad model optimization
    // -> NO call to this function when doing quad model search.

    auto mads = dynamic_cast<const NOMAD::Mads*>(step.getRootAlgorithm());
    if (nullptr == mads)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"No Mads available.");
    }
    auto barrier = mads->getMegaIterationBarrier();
    if (nullptr == barrier)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"No barrier available.");
    }
    auto frameCenter = barrier->getFirstPoint();

    // Let's work on the Mads pb.
    auto pbParams = mads->getPbParams();

    // Pb parameters
    auto n = pbParams->getAttributeValue<size_t>("DIMENSION");

    // Mesh delta frame size is  used to scale the proposed search direction
    auto mesh = step.getIterationMesh();
    if (nullptr == mesh)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"No mesh available.");
    }
    //  Frame size of the mesh
     NOMAD::ArrayOfDouble frameSize = mesh->getDeltaFrameSize();

    // Let's try a random direction larger than the mesh frame size.
    NOMAD::Direction dir(n,0.0);

    // Compute unit sphere direction
    NOMAD::Direction::computeDirOnUnitSphere(dir);

    for (size_t i = 0 ; i < n ; i++)
    {
        // Let's explore beyond the frame size
            dir[i] *= frameSize[i]*2.0; // Note: resulting pt is projected on the bounds and on the mesh if required.
    }

    trialPoints.clear();

    // insert the trial point in the trial point set
    NOMAD::EvalPoint ep(*frameCenter->getX() + dir);
    ep.setPointFrom(frameCenter, NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(&step));
    trialPoints.insert(ep);

    return true;
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

    // Custom Evaluator
    std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
    TheMainStep.setEvaluator(std::move(ev));

    // Main step start initializes Mads (default algorithm)
    TheMainStep.start();

    // Registering the callback function to generate Mads user search trial points
    // Setting USER_SEARCH yes is also required (see above)
    auto mads = std::dynamic_pointer_cast<NOMAD::Mads>(TheMainStep.getAlgo(NOMAD::StepType::ALGORITHM_MADS));
    if (nullptr == mads)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot access to Mads algorithm");
    }
    mads->addCallback(NOMAD::CallbackType::USER_METHOD_SEARCH, userSearchMethodCallback);

    TheMainStep.run();
    TheMainStep.end();

    return 0;
}

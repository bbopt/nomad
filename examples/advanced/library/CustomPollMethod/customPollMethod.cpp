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
/*  Example of a program that makes NOMAD do a Mads custom/user poll in     */
/*  addition to the Ortho 2n poll method.                                   */
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
#include "Math/MatrixUtils.hpp"

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

    //  black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 400*N);

    // Starting point
    allParams->setAttributeValue("X0", NOMAD::Point(N, 0.0) );

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(N, -10.0));
    allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(N, 10.0));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes;
    bbOutputTypes.emplace_back(NOMAD::BBOutputType::OBJ);
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes);

    // NOTE: USER_POLL directions can be combined with other direction type or not.
    NOMAD::DirectionTypeList dtList = {NOMAD::DirectionType::USER_POLL, NOMAD::DirectionType::ORTHO_2N};

    allParams->setAttributeValue("DIRECTION_TYPE",dtList);

    allParams->setAttributeValue("QUAD_MODEL_SEARCH", false); // Disable QMS and NM for clarity of the display
    allParams->setAttributeValue("NM_SEARCH", false);

    // Display
    allParams->setAttributeValue("DISPLAY_DEGREE", 3);
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj"));

    // Parameters validation
    allParams->checkAndComply();

}

// The function to generate QRMads user poll directions. This is registered as a callback below.
bool userQRPollMethodCallback(const NOMAD::Step& step, std::list<NOMAD::Direction> & dirs, const size_t n)
{
    // Important: by default USER_CALLS are disabled when doing quad model optimization
    // -> NO call to this function when doing quad model search.

    auto mads = dynamic_cast<const NOMAD::Mads*>(step.getRootAlgorithm());
    if (nullptr == mads)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"No Mads available.");
    }
    // Access to the current poll method frame center
    auto callingPoll = dynamic_cast<const NOMAD::PollMethodBase*>(&step);
    if (nullptr == callingPoll)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"No poll method available.");
    }
    auto frameCenter = callingPoll->getFrameCenter();

    // Let's work on the Mads pb. Remark: Mads does not see the fixed variables.
    auto pbParams = mads->getPbParams();

    // Pb parameters
    auto nPb = pbParams->getAttributeValue<size_t>("DIMENSION");
    if (nPb != n)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Dimension pb.");
    }

    // Mesh delta frame size is  used to scale the proposed search direction
    auto mesh = step.getIterationMesh();
    if (nullptr == mesh)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"No mesh available.");
    }
    //  Box size has the dimension of the Mads problem.
     NOMAD::ArrayOfDouble boxSize = mesh->getDeltaFrameSize();

    dirs.clear();
    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);

    while (dirUnit[0] == 0)
    {
        NOMAD::Direction::computeDirOnUnitSphere(dirUnit);
    }

    // Matrix M
    auto ** M = new double*[n];
    for (size_t i = 0; i < n; ++i)
    {
        M[i] = new double [n];
        M[i][0] = dirUnit[i].todouble();
        for (size_t j = 1; j < n; ++j)
        {
            M[i][j] = (i == j)? 1.0:0.0;
        }
    }

    // std::cout << "M matrix for QR:" <<std::endl;
    for (size_t i = 0; i < n; ++i)
    {
        NOMAD::ArrayOfDouble aod(n);
        for (size_t j = 0; j < n; ++j)
        {
            aod[j]=M[i][j];
        }
        // std::cout << aod.display() << std::endl;
    }

    // Matrices Q and R
    auto ** Q = new double*[n];
    auto ** R = new double*[n];
    for (size_t i = 0; i < n; ++i)
    {
        Q[i] = new double [n];
        R[i] = new double [n];
    }

    std::string error_msg;
    bool success = NOMAD::qr_factorization (error_msg,M,Q,R,static_cast<int>(n),static_cast<int>(n));

    if ( !success || !error_msg.empty())
    {
        std::cerr << "QR decomposition for QR 2N poll method has failed" << std::endl;
        return false;
    }

    // std::cout << "Direction after QR decomposition: " << std::endl;
    // Ordering D_k alternates Qk and -Qk instead of [Q_k -Q_k]
    NOMAD::Direction dir(n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            dir[j] = Q[j][i];
        }

        dirs.push_back(dir);
        dirs.push_back(-dir);
    }

    // Delete M, Q and R:
    for ( size_t i = 0 ; i < n ; ++i )
    {
        delete [] M[i];
        delete [] Q[i];
        delete [] R[i];
    }
    delete [] Q;
    delete [] R;
    delete [] M;

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

    // Registering the callback function to generate Mads user poll trial points via poll directions.
    // The user poll also requires to set DIRECTION_TYPE USER_POLL (see above)
    auto mads = std::dynamic_pointer_cast<NOMAD::Mads>(TheMainStep.getAlgo(NOMAD::StepType::ALGORITHM_MADS));
    if (nullptr == mads)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"Cannot access to Mads algorithm");
    }
    mads->addCallback(NOMAD::CallbackType::USER_METHOD_POLL, userQRPollMethodCallback);

    TheMainStep.run();
    TheMainStep.end();

    return 0;
}

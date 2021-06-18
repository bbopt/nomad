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

#include "../../Algos/QuadModel/QuadModelEvaluator.hpp"
#include "../../Output/OutputQueue.hpp"

// Destructor
NOMAD::QuadModelEvaluator::~QuadModelEvaluator()
{
}


void NOMAD::QuadModelEvaluator::init()
{
    _displayLevel = (std::string::npos != _modelDisplay.find("X"))
                        ? NOMAD::OutputLevel::LEVEL_INFO
                        : NOMAD::OutputLevel::LEVEL_DEBUGDEBUG;

    if (nullptr == _model)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Evaluator: a model is required (nullptr)");
    }
}


//*------------------------------------------------------*/
//*       evaluate the quad model at a given point       */
//*------------------------------------------------------*/
std::vector<bool> NOMAD::QuadModelEvaluator::eval_block(NOMAD::Block &block,
                                               const NOMAD::Double &hMax,
                                               std::vector<bool> &countEval) const
{
    std::vector<bool> evalOk;
    countEval.clear();

    // Verify there is at least one point to evaluate
    if (0 == block.size())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Evaluator: eval_block called with an empty block");
    }

    // points were sent to the evaluator in full space.
    // Convert points to subspace, because model is in subspace.
    for (size_t i = 0; i < block.size(); i++)
    {
        block[i] = std::make_shared<NOMAD::EvalPoint>(block[i]->makeSubSpacePointFromFixed(_fixedVariable));
    }

    size_t m = block.size();
    size_t n = block[0]->size();

    const auto bbot = _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");

    size_t nbConstraints = NOMAD::getNbConstraints(bbot);
    size_t nbModels = nbConstraints+1;

    // Init the matrices for prediction
    // Creation of matrix for input / output of SGTELIB model
    SGTELIB::Matrix M_predict (  "M_predict", static_cast<int>(m), static_cast<int>(nbModels));
    SGTELIB::Matrix X_predict("X_predict", static_cast<int>(m), static_cast<int>(n));

    int j = 0;
    for (auto it = block.begin(); it != block.end(); it++, j++)
    {
        if (!(*it)->isComplete())
        {
            throw NOMAD::Exception(__FILE__, __LINE__, "Evaluator: Incomplete point " + (*it)->display());
        }

        std::string s = "X" + itos(j) +" =" + (*it)->display();
        NOMAD::OutputQueue::Add(s, _displayLevel);

        // Set the input matrix
        for (size_t i = 0; i < (*it)->size(); i++)
        {
            X_predict.set(j, static_cast<int>(i), (*(*it))[i].todouble());
        }

        // Reset point outputs
        // By default, set everything to -1
        // Note: Currently NOMAD cannot set a bbo value by index, so we have to
        // work around by constructing a suitable string.
        // Note: Why set some default values on bbo?
        NOMAD::ArrayOfString defbbo(bbot.size(), "-1");
        (*it)->setBBO(defbbo.display(), bbot, NOMAD::EvalType::MODEL);

    }

    // ------------------------- //
    //   Output Prediction    //
    // ------------------------- //
    NOMAD::OutputQueue::Add("Predict with quadratic formulation... ", _displayLevel);


    // Unfortunately, Sgtelib is not thread-safe.
    // For this reason we have to set part of the eval_x code to critical.
#ifdef _OPENMP
#pragma omp critical(SgtelibEvalX)
#endif // _OPENMP
    {
        _model->check_ready(__FILE__,__FUNCTION__,__LINE__);

        _model->predict(X_predict, &M_predict);
        NOMAD::OutputQueue::Add("ok", _displayLevel);
    }

    j = 0;
    // Verify all points are completely defined
    for (auto it = block.begin(); it != block.end(); it++, j++)
    {
        std::string s = "X" + itos(j) +": " + (*it)->display();
        NOMAD::OutputQueue::Add(s, _displayLevel);
        // ====================================== //
        // Output display                    //
        // ====================================== //
        std::string sObj = "F = ";
        std::string sCons = "C = [ ";
        for (size_t i = 0; i < nbModels; i++)
        {
            if (bbot[i] != NOMAD::BBOutputType::OBJ)
                sCons += std::to_string(M_predict.get(j,static_cast<int>(i))) + " ";
            else
                sObj  += std::to_string(M_predict.get(j,static_cast<int>(i))) + " ";
        }
        s = sObj + ((nbConstraints>0 ) ? sCons+" ]":"") ;
        NOMAD::OutputQueue::Add(s, _displayLevel);

        // ====================================== //
        // Application of the formulation         //
        // ====================================== //
        NOMAD::ArrayOfDouble newbbo(bbot.size(), -1);

        // ----------------- //
        //   Set BBO         //
        // ----------------- //
        for (size_t i = 0; i < nbModels; i++)
        {
            newbbo[i] = M_predict.get(j,static_cast<int>(i));
        }
        NOMAD::ArrayOfDouble fullPrecision(bbot.size(), NOMAD::DISPLAY_PRECISION_FULL);
        (*it)->setBBO(newbbo.display(fullPrecision), bbot, NOMAD::EvalType::MODEL);

        // ================== //
        // Exit Status        //
        // ================== //
        countEval.push_back( true );
        evalOk.push_back(true);
        (*it)->setEvalStatus(NOMAD::EvalStatusType::EVAL_OK, NOMAD::EvalType::MODEL);

    }

    // Convert points back to full space.
    for (size_t i = 0; i < block.size(); i++)
    {
        block[i] = std::make_shared<NOMAD::EvalPoint>(block[i]->makeFullSpacePointFromFixed(_fixedVariable));
    }

    return evalOk;
}



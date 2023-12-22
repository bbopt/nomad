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
/*------------------------------------------------------------------------------------*/
/*  example of a program that makes NOMAD restarts after a chosen stopping condition  */
/* enables to restart from current state with different parameter settings (here VNS) */
/*------------------------------------------------------------------------------------*/

#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Algos/Mads/MadsMegaIteration.hpp"
#include "Algos/AlgoStopReasons.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/LHSearchType.hpp"
#include "Eval/SuccessStats.hpp"
#include <iostream>
#include <fstream>


std::shared_ptr<NOMAD::MeshBase> mesh;
auto params = std::make_shared<NOMAD::AllParameters>();
auto madsStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::MadsStopType>>();

// Problem given to the solver
std::string problembName = "./problems/RASTRIGIN";
std::string bbexe = "./problems/RASTRIGIN/bb.exe";

// Counter for the number of successive failure
int compteur = 0;

// Thresholds
int stopConsFailures = 3; // Threshold for the number of successive failure
int stopMeshIndex = -2; // Threshold for the mesh index

// Stopping conditions initialization
bool stopCondition = false;

// Nombre d'entrée dans la méga itération
int nbEnterMegaIter = 0;

// Number of bbeval with VNS activated
int nbBbeBeforeVNS = 0;
int nbBbeWithVNS = 0;

// Run number
int i = 0;

// Number and state of the iteration
int nbIteration = 0;
bool iterationSucess = false;
bool VNS = false;


void initAllParams(std::shared_ptr<NOMAD::AllParameters> allParams, const size_t n)
{
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", n);
    // Blackbox to evaluate
    allParams->setAttributeValue("BB_EXE", bbexe);
    // 100 black-box evaluations
    // allParams->setAttributeValue("MAX_BB_EVAL", 200);
    // Starting point
        // RASTRIGIN
    allParams->setAttributeValue("X0", NOMAD::Point(n, -5.12) );
    /*
        // COCO
    NOMAD::Point x0(n);
    x0[0] = -4.426533;
    x0[1] = 4.238298;
    x0[2] = 1.037814;
    x0[3] = -3.581491;
    x0[4] = -4.835488;
    x0[5] = 1.023338;
    x0[6] = 2.797704;
    x0[7] = -2.753894;
    x0[8] = -2.859614;
    x0[9] = -4.790012;
    allParams->setAttributeValue("X0", x0);
    */
    // Bounds
        // RASTRIGIN
    //allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -5.12));
    //allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(n, 5.12));
        // COCO
    //allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -5.0));
    //allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(n, 5.0));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes;
    bbOutputTypes.push_back(NOMAD::BBOutputType::OBJ);
        // RASTRIGIN
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj"));
    
    
    // COCO
    /*
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj cons_h"));
    */

    allParams->setAttributeValue("DISPLAY_DEGREE", 2);
    allParams->setAttributeValue("DISPLAY_ALL_EVAL", true);

    // Parameters validation
    allParams->checkAndComply();

}


/*--------------------------------------*/
/* Before each MegaIteration, verify if */
/* the algorithm should stop.           */
/*--------------------------------------*/
void userMegaIterationStart(const NOMAD::Step& step,
                            bool &stop)
{
    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>(&step);
    auto bbe = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval(); // Number of blackbox evaluation
    
    //std::cout << "CallBack START no VNS" << std::endl;
    
    stop = false;
    
    
    if (nullptr != megaIter)
    {
        // std::cout << std::endl << "Start mega itération!" << std::endl;
        ++nbIteration;
        // We use the success of the previous megaIteration. The success is reset to UNDEFINED at the defaultStart before this callback function is called.
        auto success = megaIter->getSuccessType(); //NOT_TRIALS, UNSUCCESSFUL, PARTIAL_SUCCESS, FULL_SUCCESS ou UNDEFINED
        
        // The counter is updated regarding the reuslt of the previous iteration
        if (NOMAD::SuccessType::FULL_SUCCESS == success || NOMAD::SuccessType::PARTIAL_SUCCESS == success) // counter reset for success or partial success
        {
            //std::cout << "SUCCESS OR PARTIAL SUCCESS" << std::endl;
            compteur = 0;
            iterationSucess = true;
        }
        else if (NOMAD::SuccessType::UNDEFINED == success) // no event case
        {
            //std::cout << "Default type set at start" << std::endl;
        }
        else // update if unsuccessful or no trial points
        {
            //std::cout << "UNSUCCESSFUL OR NO_TRIALS" << std::endl;
            ++compteur;
        }
        // std::cout << "Nb consecutive fails:" << compteur << std::endl;
        
        // We can also have access to this information through nomad statistics
        // auto nbConsecutiveFail = megaIter->getConstSuccessStats().getStatsNbConsecutiveFail();
        // std::cout << "Compteur NOMAD : " << nbConsecutiveFail << std::endl;
        
        if (i > 0) {
            // Keep information about the value of mesh index from the previous run
            NOMAD::ArrayOfDouble oldMeshIndices = mesh->getMeshIndex();
            //std::cout << "old mesh indices : " << oldMeshIndices << std::endl;
            
            // Let's pass the mesh
            //std::cout << std::endl << "Mesh :" << std::endl;
            mesh = megaIter->getMesh();
            //std::cout << *mesh << std::endl;
            
            // Set the mesh index value
            mesh->setMeshIndex(oldMeshIndices);
        }
        else {
            // Let's pass the mesh
            //std::cout << std::endl << Mesh :" << std::endl;
            mesh = megaIter->getMesh();
            //std::cout << *mesh << std::endl;
        }
        
        // Let's print a parameter on MAX_BB_EVAL
        // int maxBbEval;
        // maxBbEval = params->getAttributeValue<size_t>("MAX_BB_EVAL");
        // std::cout << "MAX_BB_EVAL : " << maxBbEval << std::endl;

        // Collect mesh index information
        NOMAD::ArrayOfDouble meshIndices = mesh->getMeshIndex();
        NOMAD::ArrayOfDouble meshIndexStop(meshIndices.size(), stopMeshIndex);
        //std::cout << "mesh indices: " << meshIndices << std::endl;
       // std::cout << "mesh stop indices: " << meshIndexStop << std::endl;
        
        // Reinitialize the iteration state
        iterationSucess = false;

        // Stopping conditions:
        stopCondition = (compteur >= stopConsFailures);
        //std::cout << std::endl << "Stopping condition: " << stopCondition << std::endl;
        
        if (stopCondition) // Stop motivated by user conditions
        {
            stop = true;
        }
    }
}

/*--------------------------------------------*/
/* Before each MegaIteration after using VNS, */
/* the algorithm should stop.                 */
/*--------------------------------------------*/
void userMegaIterationStartForVNS(const NOMAD::Step& step,
                          bool &stop)
{
    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>(&step);
    auto bbe = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval(); // Récupère le nombre de blackbox evaluations
    
 
    stop = false;
    
    if (nullptr != megaIter)
    {
    
        ++nbEnterMegaIter;
        ++nbIteration;
        // std::cout << std::endl << "Nb mega iter: " << nbEnterMegaIter << std::endl;
        // We use the success of the previous megaIteration. The success is reset to UNDEFINED at the defaultStart before this callback function is called.
        auto success = megaIter->getSuccessType(); //NO_TRIALS, UNSUCCESSFUL, PARTIAL_SUCCESS, FULL_SUCCESS ou UNDEFINED
        
        // The counter is updated regarding the reuslt of the previous iteration
        if (NOMAD::SuccessType::FULL_SUCCESS == success || NOMAD::SuccessType::PARTIAL_SUCCESS == success) // counter reset for success or partial success
        {
            compteur = 0;
            iterationSucess = true;
        }
        else if (NOMAD::SuccessType::UNDEFINED == success) // no event case
        {
            std::cout << "Default type set at start" << std::endl;
        }
        else // update if unsuccessful or no trial points
        {
            ++compteur;
        }
        //std::cout << "Number of successive fails: " << compteur << std::endl;
     
        
        // Keep information about the value of mesh index from the previous run
        NOMAD::ArrayOfDouble oldMeshIndices = mesh->getMeshIndex();
        // std::cout << "old mesh indices: " << oldMeshIndices << std::endl;
        
        // Let's pass the mesh
        //std::cout << std::endl << "Mesh :" << std::endl;
        mesh = megaIter->getMesh();
        //std::cout << *mesh << std::endl;
        
        // Set the mesh index value
        mesh->setMeshIndex(oldMeshIndices);

        // Collect mesh index information
        NOMAD::ArrayOfDouble meshIndices = mesh->getMeshIndex();
        NOMAD::ArrayOfDouble meshIndexStop(meshIndices.size(), stopMeshIndex);
        //std::cout << "mesh indices : " << meshIndices << std::endl;
        //std::cout << "mesh stop indices : " << meshIndexStop << std::endl;
        
        iterationSucess = false;
        
        // Stopping conditions:
        stopCondition = (compteur >= stopConsFailures);
        // std::cout << std::endl << "Stopping condition: " << stopCondition << std::endl;
        
        // Stop motivated by user conditions : after one megaiteration
        if (NOMAD::SuccessType::UNDEFINED != success && nbEnterMegaIter > 1)
        {
            stop = true;
            nbEnterMegaIter = 0;
            nbBbeWithVNS += (bbe - nbBbeBeforeVNS);
        }
    }
}


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    // Dimension (Number of variables)
    size_t n = 12;

    NOMAD::MainStep TheMainStep;

    initAllParams(params, n);
    TheMainStep.setAllParameters(params);
    
    // Custom batch Evaluator
    auto ev = std::make_unique<NOMAD::Evaluator>(params->getEvalParams(),NOMAD::EvalType::BB);
    TheMainStep.addEvaluator(std::move(ev));
    
    std::vector<NOMAD::EvalPoint> bf;
    std::vector<NOMAD::EvalPoint> bi;
    bool stopLoop = false;
    bool stopReasonIsCallback = false;

    // Main run
    try
    {
        i = 0;
        // successive runs:
        do
        {
            std::cout << std::endl << "MADS run #" + NOMAD::itos(i) << std::endl;
            
            // Stats files
            std::string statsFile = problembName;
            statsFile += "/STATS/stats";
            statsFile += std::to_string(i);
            statsFile += ".txt bbe ( sol ) obj";
            params->setAttributeValue("STATS_FILE", NOMAD::ArrayOfString(statsFile));
            
            std::string detailedStatsFile = problembName;
            detailedStatsFile += "/STATS/detailedStats";
            detailedStatsFile += std::to_string(i);
            detailedStatsFile += ".txt";
            params->setAttributeValue("EVAL_STATS_FILE", std::string(detailedStatsFile));
            
            std::string historyFile = problembName;
            historyFile += "/STATS/history";
            historyFile += std::to_string(i);
            historyFile += ".txt";
            params->setAttributeValue("HISTORY_FILE", std::string(historyFile));
            
            // Set callbacks
    
            if (!stopCondition) // MADS default without VNS
            {
                std::cout << "MADS default no VNS" << std::endl;
                TheMainStep.addCallback(NOMAD::CallbackType::MEGA_ITERATION_START, userMegaIterationStart);
                
                if ( i > 0 )
                {
                    // Seuil compteur dynamique
                    stopConsFailures = static_cast<int>(3*std::ceil(i/5.0));
                    
                    // Desactivate VNS MADS search
                    params->getRunParams()->setAttributeValue("VNS_MADS_SEARCH", false);


                    // New starting points
                    NOMAD::ArrayOfPoint x0s;
                    if (bf.size() > 0)
                    {
                        x0s.push_back(bf[0]);
                    }
                    if (bi.size() > 0)
                    {
                        x0s.push_back(bi[0]);
                    }
                    params->getPbParams()->setAttributeValue("X0", x0s);
                    
                    //std::cout << mesh->getDeltaFrameSize() << std::endl;
                    params->getPbParams()->setAttributeValue("INITIAL_FRAME_SIZE", mesh->getDeltaFrameSize());
                    
                    // at least one evaluation is conducted if points bf and bi are null
                    std::string lhSearchStr = (bf.empty() && bi.empty())  ? " 1 0" : "0 0";
                    params->getRunParams()->setAttributeValue("LH_SEARCH", NOMAD::LHSearchType(lhSearchStr));
                    
                    // Parameters validation
                    params->checkAndComply();
                }
                
            }
            else  // MADS default + VNS for 1 mega iteration only
            {
                // Update thresholds for the stopping criterion
                stopConsFailures = static_cast<int>(3*std::ceil(i/5.0));

                TheMainStep.addCallback(NOMAD::CallbackType::MEGA_ITERATION_START, userMegaIterationStartForVNS);
                nbBbeBeforeVNS = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval(); // Keeps the number of blackbox evaluations before using VNS
                
                // Activate VNS MADS serach
                params->getRunParams()->setAttributeValue("VNS_MADS_SEARCH", true);
                params->getRunParams()->setAttributeValue("VNS_MADS_SEARCH_TRIGGER", NOMAD::Double(1));
                
                // New starting points
                NOMAD::ArrayOfPoint x0s;
                if (bf.size() > 0)
                {
                    x0s.push_back(bf[0]);
                }
                if (bi.size() > 0)
                {
                    x0s.push_back(bi[0]);
                }
                params->getPbParams()->setAttributeValue("X0", x0s);
                
                //std::cout << mesh->getDeltaFrameSize() << std::endl;
                params->getPbParams()->setAttributeValue("INITIAL_FRAME_SIZE", mesh->getDeltaFrameSize());
                
                // at least one evaluation is conducted if points bf and bi are null
                std::string lhSearchStr = (bf.empty() && bi.empty())  ? " 1 0" : "0 0";
                params->getRunParams()->setAttributeValue("LH_SEARCH", NOMAD::LHSearchType(lhSearchStr));
                
                // Parameters validation
                params->checkAndComply();
                
            }

            // The run
            TheMainStep.start();
            TheMainStep.run();
            TheMainStep.end();
            
            auto algoStopReason = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get(TheMainStep.getAllStopReasons()); // Récupère le stop reason de mads
            stopLoop = algoStopReason->testIf(NOMAD::MadsStopType::MESH_PREC_REACHED) || algoStopReason->testIf(NOMAD::MadsStopType::MIN_MESH_INDEX_REACHED) ; // Test si le stop reason correspond MESH_PREC_REACHED
            std::cout << algoStopReason->getStopReasonAsString() << std::endl; // Juste affiche le stop reason.
            
            stopReasonIsCallback = (algoStopReason->getStopReasonAsString() == "User-stopped in a callback function (Base)");
       

            bf.clear();
            bi.clear();
            NOMAD::CacheBase::getInstance()->findBestFeas(bf, NOMAD::Point(n), NOMAD::EvalType::BB,NOMAD::ComputeType::STANDARD);
            NOMAD::CacheBase::getInstance()->findBestInf(bi, NOMAD::INF, NOMAD::Point(n), NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD);
            
            ++i;
            
        } while (!stopLoop); // run until one NOMAD termination criterion is met
       
}

    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    return EXIT_SUCCESS;
}

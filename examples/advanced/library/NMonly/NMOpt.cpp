
#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Cache/CacheBase.hpp"


void initParams1(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 5;   // Number of variables
    p.getPbParams()->setAttributeValue("DIMENSION", n);
    p.getEvalParams()->setAttributeValue("BB_EXE", std::string("./u.exe"));

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

    p.getEvalParams()->setAttributeValue("TMP_DIR", std::string("/tmp"));
    p.getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("OBJ PB EB"));

    p.getRunParams()->setAttributeValue("NM_OPTIMIZATION",true);

    p.getRunParams()->setAttributeValue("HOT_RESTART_ON_USER_INTERRUPT", false);
    p.getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    p.getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);


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
    // Part 1: FIXED_VARIABLE 0-2
    auto params = std::make_shared<NOMAD::AllParameters>();
    initParams1(*params);
    TheMainStep->setAllParameters(params);

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

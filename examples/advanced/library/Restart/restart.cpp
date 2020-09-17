/*--------------------------------------------------------------------------*/
/*  example of a program that makes NOMAD restarts after failed iterations  */
/*--------------------------------------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Algos/MegaIteration.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/LHSearchType.hpp"

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
    NOMAD::Double c1 = 0.0 , c2 = 0.0;
    for ( int i = 0 ; i < 5 ; i++ )
    {

        c1 += (x[i]-1).pow2();
        c2 += (x[i]+1).pow2();
    }
    NOMAD::Double constr1 = c1-25;
    NOMAD::Double constr2 = 25-c2;
    auto bbOutputType = _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    std::string bbo = x[4].tostring();
    bbo += " " + constr1.tostring();
    bbo += " " + constr2.tostring();
    x.setBBO(bbo, bbOutputType, getEvalType());

    countEval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
}


void initAllParams(std::shared_ptr<NOMAD::AllParameters> allParams, const size_t n)
{
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", n);
    // 100 black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 100);
    // Starting point
    allParams->setAttributeValue("X0", NOMAD::Point(n, 0.0) );

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -6.0 )); // all var. >= -6
    NOMAD::ArrayOfDouble ub(n);
    ub[0] = 5.0;    // x_1 <= 5
    ub[1] = 6.0;    // x_2 <= 6
    ub[2] = 7.0;    // x_3 <= 7
    allParams->setAttributeValue("UPPER_BOUND", ub);
    allParams->getPbParams()->setAttributeValue("GRANULARITY", NOMAD::ArrayOfDouble(n, 0.0000001));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes;
    bbOutputTypes.push_back(NOMAD::BBOutputType::OBJ);
    bbOutputTypes.push_back(NOMAD::BBOutputType::EB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::EB);
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    allParams->setAttributeValue("DISPLAY_DEGREE", 2);
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj"));
    allParams->setAttributeValue("DISPLAY_ALL_EVAL", true);

    allParams->getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    allParams->getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);


    // Parameters validation
    allParams->checkAndComply();

}


/*-------------------------------------*/
/* After each MegaIteration, verify if */
/* the algorithm should stop.          */
/*-------------------------------------*/
void userMegaIterationEnd(const NOMAD::Step& step,
                          bool &stop)
{
    auto megaIter = dynamic_cast<const NOMAD::MegaIteration*>(&step);
    auto bbe = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
    auto success = megaIter->getSuccessType();
    if (NOMAD::SuccessType::UNSUCCESSFUL == success && bbe > 10)
    {
        // Stop motivated by user conditions
        stop = true;
    }
}


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    // Dimension (Number of variables)
    size_t n = 5;

    NOMAD::MainStep TheMainStep;

    // Set callbacks
    TheMainStep.addCallback(NOMAD::CallbackType::MEGA_ITERATION_END, userMegaIterationEnd);

    // Set parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initAllParams(params, n);
    TheMainStep.setAllParameters(params);

    std::vector<NOMAD::EvalPoint> bf;
    std::vector<NOMAD::EvalPoint> bi;

    // Main run
    try
    {
        // successive runs:
        for ( int i = 0 ; i < 6 ; ++i )
        {
            std::cout << std::endl << "MADS run #" + NOMAD::itos(i) << std::endl;

            // Custom Evaluator
            std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
            TheMainStep.setEvaluator(std::move(ev));

            // not for the first run:
            if ( i > 0 )
            {
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

            bf.clear();
            bi.clear();
            NOMAD::CacheBase::getInstance()->findBestFeas(bf, NOMAD::Point(n), NOMAD::EvalType::BB);
            NOMAD::CacheBase::getInstance()->findBestInf(bi, NOMAD::INF, NOMAD::Point(n), NOMAD::EvalType::BB);
        }
    }

    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    return EXIT_SUCCESS;
}


#include "../Param/EvaluatorControlGlobalParameters.hpp"


/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::EvaluatorControlGlobalParameters::init()
{
    _typeName = "EvaluatorControl";

    try
    {
        #include "../Attribute/evaluatorControlGlobalAttributesDefinition.hpp"
        registerAttributes(_definition);

        // Note: we cannot call checkAndComply() here, the default values
        // are not valid, for instance DIMENSION, X0, etc.

    }
    catch (NOMAD::Exception & e)
    {
        std::string errorMsg = "Attribute registration failed: ";
        errorMsg += e.what();
        throw NOMAD::Exception(__FILE__,__LINE__, errorMsg);
    }

}

/*----------------------------------------*/
/*            check the parameters        */
/*----------------------------------------*/
void NOMAD::EvaluatorControlGlobalParameters::checkAndComply()
{
    checkInfo();

    if (!toBeChecked())
    {
        // Early out
        return;
    }

    if (getAttributeValueProtected<size_t>("MAX_BB_EVAL",false) <= 0)
    {
        setAttributeValue("MAX_BB_EVAL", NOMAD::INF_SIZE_T);
    }

    if (getAttributeValueProtected<size_t>("MAX_EVAL",false) <= 0)
    {
        setAttributeValue("MAX_EVAL", NOMAD::INF_SIZE_T);
    }

    auto sgteBlockSize = getAttributeValueProtected<size_t>("SGTE_MAX_BLOCK_SIZE", false);
    if (0 == sgteBlockSize)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Parameter SGTE_MAX_BLOCK_SIZE must be positive");
    }
    else if (sgteBlockSize > getAttributeValueProtected<size_t>("MAX_SGTE_EVAL", false))
    {
        setAttributeValue("SGTE_MAX_BLOCK_SIZE", getAttributeValueProtected<size_t>("MAX_SGTE_EVAL", false));
    }

    _toBeChecked = false;

}
// End checkAndComply()





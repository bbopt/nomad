
#include "../Param/AllParameters.hpp"

// Do we need to call checkAndComply() ?
bool NOMAD::AllParameters::toBeChecked() const
{
    bool check =    (!_pbParams          || _pbParams->toBeChecked())
                 || (!_evalParams        || _evalParams->toBeChecked())
                 || (!_evaluatorControlGlobalParams || _evaluatorControlGlobalParams->toBeChecked())
                 || (!_evaluatorControlParams || _evaluatorControlParams->toBeChecked())
                 || (!_runParams         || _runParams->toBeChecked())
                 || (!_cacheParams       || _cacheParams->toBeChecked())
                 || (!_dispParams        || _dispParams->toBeChecked());
    return check;
}



void NOMAD::AllParameters::resetToDefaultValues() noexcept
{
    _runParams->resetToDefaultValues();
    _pbParams->resetToDefaultValues();
    _cacheParams->resetToDefaultValues();
    _dispParams->resetToDefaultValues();
    _evalParams->resetToDefaultValues();
    _evaluatorControlGlobalParams->resetToDefaultValues();
    _evaluatorControlParams->resetToDefaultValues();
}

/*----------------------------------------------------------------*/
/*          read a parameters file and interpret attributes       */
/*----------------------------------------------------------------*/
void NOMAD::AllParameters::read(const std::string &paramFile, bool overwrite , bool resetAllEntries )
{

    // Read all entries
    NOMAD::Parameters::readParamFileAndSetEntries(paramFile, overwrite ,resetAllEntries );

    // Read entries and set attribute values for each type of parameters
    _runParams->readEntries();
    _pbParams->readEntries();
    _evalParams->readEntries();
    _evaluatorControlGlobalParams->readEntries();
    _evaluatorControlParams->readEntries();
    _cacheParams->readEntries();
    _dispParams->readEntries();

}


/*----------------------------------------------------------------*/
/*   read a parameter line, interpret it, and update the links    */
/*----------------------------------------------------------------*/
void NOMAD::AllParameters::readParamLine(const std::string &line)
{
    // Note: Here we have the line, but we do not know which parameter type
    // it applies to. For now, the workaround is to create a temporary
    // ParameterEntry and to extract the parameter name from it.
    // This could be cleaner, since we will create another ParameterEntry
    // in readParamLine().
    auto pe = std::make_unique<NOMAD::ParameterEntry>(line);
    std::string name = pe->getName();

    bool overwrite = true;
    if ( _cacheParams->isRegisteredAttribute(name) )
    {
        _cacheParams->readParamLine(line, overwrite);
    }
    else if ( _dispParams->isRegisteredAttribute(name) )
    {
        _dispParams->readParamLine(line, overwrite);
    }
    else if (_evalParams->isRegisteredAttribute(name))
    {
        _evalParams->readParamLine(line, overwrite);
    }
    else if ( _evaluatorControlParams->isRegisteredAttribute(name) )
    {
        _evaluatorControlParams->readParamLine(line, overwrite);
    }
    else if ( _evaluatorControlGlobalParams->isRegisteredAttribute(name) )
    {
        _evaluatorControlGlobalParams->readParamLine(line, overwrite);
    }
    else if ( _pbParams->isRegisteredAttribute(name) )
    {
        _pbParams->readParamLine(line, overwrite);
    }
    else if ( _runParams->isRegisteredAttribute(name) )
    {
        _runParams->readParamLine(line, overwrite);
    }
    else
    {
        std::string err = "Unknown parameter: " + name;
        std::cerr << err << std::endl;
    }

}

bool NOMAD::AllParameters::isAlgoCompatible(const NOMAD::AllParameters& allP_tmp) const
{
    return _pbParams->isAlgoCompatible(allP_tmp.getPbParams().get()) &&
           _dispParams->isAlgoCompatible(allP_tmp.getDispParams().get()) &&
           _runParams->isAlgoCompatible(allP_tmp.getRunParams().get()) &&
           _evaluatorControlGlobalParams->isAlgoCompatible(allP_tmp.getEvaluatorControlGlobalParams().get()) &&
           _evaluatorControlParams->isAlgoCompatible(allP_tmp.getEvaluatorControlParams().get()) &&
           _evalParams->isAlgoCompatible(allP_tmp.getEvalParams().get()) ;
}


/*----------------------------------------------------------------*/
/*          read a parameters file and interpret attributes       */
/*----------------------------------------------------------------*/
void NOMAD::AllParameters::displayHelp(const std::string &helpSubject , bool devHelp , std::ostream & os )
{

    std::ostringstream ossBasic,ossAdvanced;

    _evalParams->displayHelp(helpSubject, devHelp, ossBasic , ossAdvanced);
    _evaluatorControlGlobalParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _evaluatorControlParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _runParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _pbParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);
    _cacheParams->displayHelp(helpSubject,devHelp,ossBasic , ossAdvanced);
    _dispParams->displayHelp(helpSubject,devHelp, ossBasic , ossAdvanced);

    if ( !devHelp )
    {
        if (ossBasic.str().empty() && ossAdvanced.str().empty())
        {
            os << "No help found for " << helpSubject << std::endl << std::endl;
        }

        if (!ossBasic.str().empty())
        {
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << "-------------------------------- BASIC PARAMETERS -----------------------------" << std::endl;
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << std::endl << ossBasic.str() << std::endl << std::endl;
        }

        if (!ossAdvanced.str().empty())
        {
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << "------------------------------ ADVANCED PARAMETERS ----------------------------" << std::endl;
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << std::endl << ossAdvanced.str() << std::endl << std::endl;
        }
    }
    else
    {
        if ( ossBasic.str().empty() )
        {
            os << "No help found for " << helpSubject << std::endl << std::endl;
        }
        else
        {
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << "----------------------------- DEVELOPPER PARAMETERS ---------------------------" << std::endl;
            os << "-------------------------------------------------------------------------------" << std::endl;
            os << ossBasic.str() << std::endl << std::endl;
        }
    }
}

/*----------------------------------------*/
/*            check the parameters        */
/*----------------------------------------*/
void NOMAD::AllParameters::checkAndComply()
{
    if (!toBeChecked())
    {
        // Early out
        return;
    }

    _evaluatorControlGlobalParams->checkAndComply();
    _evaluatorControlParams->checkAndComply();
    _pbParams->checkAndComply();
    _runParams->checkAndComply(_evaluatorControlGlobalParams, _pbParams);
    _evalParams->checkAndComply(_runParams);
    _cacheParams->checkAndComply(_runParams);
    _dispParams->checkAndComply(_runParams, _pbParams);

}
// End checkAndComply()


std::string NOMAD::AllParameters::getSetAttributeAsString() const
{
    std::string str = _runParams->getSetAttributeAsString()
           + _pbParams->getSetAttributeAsString()
           + _evalParams->getSetAttributeAsString()
           + _evaluatorControlGlobalParams->getSetAttributeAsString()
           + _evaluatorControlParams->getSetAttributeAsString()
           + _cacheParams->getSetAttributeAsString()
           + _dispParams->getSetAttributeAsString() ;
    if ( str.empty() )
        str = "All attributes have default value";
    return str;
}


void NOMAD::AllParameters::display(std::ostream & os, bool flagHelp )
{
    if (toBeChecked())
    {
        std::cerr << "Warning: AllParameters::display(): Parameters are not checked." << std::endl;
    }

    os << "----- RUN PARAMETERS -----" << std::endl;
    _runParams->display(os,flagHelp);

    os << "----- PROBLEM PARAMETERS -----" << std::endl;
    _pbParams->display(os,flagHelp);

    os << "----- EVAL PARAMETERS -----" << std::endl;
    _evalParams->display(os,flagHelp);

    os << "----- EVALUATOR CONTROL PARAMETERS (GLOBAL) -----" << std::endl;
    _evaluatorControlGlobalParams->display(os,flagHelp);

    os << "----- EVALUATOR CONTROL PARAMETERS (BY MAIN THREADK)-----" << std::endl;
    _evaluatorControlParams->display(os,flagHelp);

    os << "----- CACHE PARAMETERS -----" << std::endl;
    _cacheParams->display(os,flagHelp);

    os << "----- DISPLAY PARAMETERS -----" << std::endl;
    _dispParams->display(os,flagHelp);
}



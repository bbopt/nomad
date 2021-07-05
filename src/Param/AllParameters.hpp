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
#ifndef __NOMAD_4_0_ALLPARAMETERS__
#define __NOMAD_4_0_ALLPARAMETERS__


#include "../Param/CacheParameters.hpp"
#include "../Param/DeprecatedParameters.hpp"
#include "../Param/DisplayParameters.hpp"
#include "../Param/EvalParameters.hpp"
#include "../Param/EvaluatorControlGlobalParameters.hpp"
#include "../Param/EvaluatorControlParameters.hpp"
#include "../Param/PbParameters.hpp"
#include "../Param/RunParameters.hpp"
#include "../Type/BBInputType.hpp"
#include "../Type/BBOutputType.hpp"

#include "../nomad_nsbegin.hpp"

/// Container class for all NOMAD parameters.
/**
 Currently we have six classes of parameters: for run/main execution, problem definition, evaluation, evaluation control, cache and for display. \n

 The AllParameters::read function reads the entries for all parameters given in a single parameter file. After reading the entries, the parameters can be checked for inter-value compliance by AllParameters::checkAndComply as a sanity check of parameters. The function AllParameters::checkAndComply can also set some parameters according to other parameters (for example, the granularity of integer variables are set to 1). \n

 When changing the value of a parameter, the class where it belongs is tagged "toBeChecked". Before accessing ANY parameter in ANY class, the AllParameters::checkAndComply function must be called.

 Some parameters have a tag to indicate that they control algorithm execution. To verify that two sets of parameters are compatible we can compare all "algo" tagged parameters with the function AllParameters::isAlgoCompatible. The tagged parameters are the one having "ALGO_COMPATIBILITY_CHECK yes" in their attibute definition file.


\todo add the NOMAD parameters in a container and use foreach. This should prevent many modifications in the header when add a new type of NOMAD parameter.
 */
class AllParameters
{
private:


    // Developer: When adding a new type of NOMAD parameters update the code
    std::shared_ptr<DeprecatedParameters>        _deprecatedParams;
    std::shared_ptr<RunParameters>               _runParams;
    std::shared_ptr<PbParameters>                _pbParams;
    std::shared_ptr<CacheParameters>             _cacheParams;
    std::shared_ptr<DisplayParameters>           _dispParams;
    std::shared_ptr<EvalParameters>              _evalParams;
    std::shared_ptr<EvaluatorControlGlobalParameters>  _evaluatorControlGlobalParams;
    std::shared_ptr<EvaluatorControlParameters>  _evaluatorControlParams;


public:
    /// Constructor
    explicit AllParameters()
      : _deprecatedParams(std::make_shared<DeprecatedParameters>()),
        _runParams(std::make_shared<RunParameters>()),
        _pbParams(std::make_shared<PbParameters>()),
        _cacheParams(std::make_shared<CacheParameters>()),
        _dispParams(std::make_shared<DisplayParameters>()),
        _evalParams(std::make_shared<EvalParameters>()),
        _evaluatorControlGlobalParams(std::make_shared<EvaluatorControlGlobalParameters>()),
        _evaluatorControlParams(std::make_shared<EvaluatorControlParameters>())
    {
    }

    /**
     Do not allow copy constructor.
     Copy constructors are not defined for Parameters class, and
     we need a deep copy here: it does not make sense to copy
     the smart pointers.
     */
    AllParameters(const AllParameters &allParams) = delete;

    /**
     Do not allow copy assignement.
     Copy constructors are not defined for Parameters class, and
     we need a deep copy here: it does not make sense to copy
     the smart pointers.
     */
    AllParameters& operator=(const AllParameters& params) = delete;


    virtual ~AllParameters() {}

    /*-------------------*/
    /* setAttributeValue */
    /*-------------------*/
    // This template function implementation must be in the header to be available in the library


    /**
     Search for a registered attribute with the given name and set its value. If there is no corresponding attribute, an exception is triggered. If the templated type does not correspond to the registered type for the attribute, an exception is triggered.
     */
    template<typename T>
    void setAttributeValue(std::string name, T value)
    {

        if (_evalParams->isRegisteredAttribute(name))
        {
            _evalParams->setAttributeValue<T>(name, value);
        }
        else if (_evaluatorControlGlobalParams->isRegisteredAttribute(name))
        {
            _evaluatorControlGlobalParams->setAttributeValue<T>(name, value);
        }
        else if (_evaluatorControlParams->isRegisteredAttribute(name))
        {
            _evaluatorControlParams->setAttributeValue<T>(name, value);
        }
        else if (_runParams->isRegisteredAttribute(name))
        {
            _runParams->setAttributeValue<T>(name, value);
        }
        else if (_pbParams->isRegisteredAttribute(name))
        {
            _pbParams->setAttributeValue<T>(name, value);
        }
        else if (_dispParams->isRegisteredAttribute(name))
        {
            _dispParams->setAttributeValue<T>(name, value);
        }
        else if (_cacheParams->isRegisteredAttribute(name))
        {
            _cacheParams->setAttributeValue<T>(name, value);
        }
        else if (_deprecatedParams->isRegisteredAttribute(name))
        {
            std::string err = "setAttributeValue: attribute " + name + " is  deprecated";
            throw Exception(__FILE__, __LINE__, err);
        }
        else
        {
            // At this point, Verify att is non-null.
            std::string err = "setAttributeValue: attribute " + name + " is not registered";
            throw Exception(__FILE__, __LINE__, err);
        }
    }

    /*--------------------------*/
    /* getAttributeValue        */
    /*--------------------------*/
    // This template function implementation must be in the header to be available in the library

    /**
     Search for a registered attribute with the given name and return its value. If there is no corresponding attribute, an exception is triggered. If the templated type does not correspond to the registered type for the attribute, an exception is triggered.
     */
    template<typename T> const T&
    getAttributeValue(const std::string &name) const
    {
        if (_evalParams->isRegisteredAttribute(name))
        {
            return _evalParams->getAttributeValue<T>(name);
        }
        else if (_evaluatorControlGlobalParams->isRegisteredAttribute(name))
        {
            return _evaluatorControlGlobalParams->getAttributeValue<T>(name);
        }
        else if (_evaluatorControlParams->isRegisteredAttribute(name))
        {
            return _evaluatorControlParams->getAttributeValue<T>(name);
        }
        else if (_runParams->isRegisteredAttribute(name))
        {
            return _runParams->getAttributeValue<T>(name);
        }
        else if (_pbParams->isRegisteredAttribute(name))
        {
            return _pbParams->getAttributeValue<T>(name);
        }
        else if (_dispParams->isRegisteredAttribute(name))
        {
            return _dispParams->getAttributeValue<T>(name);
        }
        else if (_cacheParams->isRegisteredAttribute(name))
        {
            return _cacheParams->getAttributeValue<T>(name);
        }
        else
        {
            // At this point, Verify att is non-null.
            std::string err = "getAttributeValue: attribute " + name + " is not registered";
            throw Exception(__FILE__, __LINE__, err);
        }
        return _evalParams->getAttributeValue<T>(name);
    }

    const std::shared_ptr<CacheParameters>&             getCacheParams() const { return _cacheParams; }
    const std::shared_ptr<DeprecatedParameters>&        getDeprecatedParams() const { return _deprecatedParams; }
    const std::shared_ptr<DisplayParameters>&           getDispParams() const { return _dispParams; }
    const std::shared_ptr<EvalParameters>&              getEvalParams() const { return _evalParams; }
    const std::shared_ptr<EvaluatorControlParameters>&  getEvaluatorControlParams() const { return _evaluatorControlParams; }
    const std::shared_ptr<EvaluatorControlGlobalParameters>& getEvaluatorControlGlobalParams() const { return _evaluatorControlGlobalParams; }
    const std::shared_ptr<PbParameters>&                getPbParams() const { return _pbParams; }
    const std::shared_ptr<RunParameters>&               getRunParams() const { return _runParams; }


    /// Perform checkAndComply() on all parameters.
    void checkAndComply();

    /// Verify if we need to call checkAndComply().
    bool toBeChecked() const;

    /// All registered attributes are reset to their default value
    void resetToDefaultValues() noexcept;

    /// Read a parameters file into entries.
    void read(const std::string &paramFile, bool overwrite = false , bool resetAllEntries = false);

    /**
     Try readParamLine for each class of parameters until it works.
     If the parameter is not found, throw an exception.
     */
    void readParamLine(const std::string &line);

    static void eraseAllEntries();

    /**
     Compare the compatibility of the current set of parameters with a given set of parameters. The compatibility concerns only parameters influencing the execution of the algorithms (that is those with Attribute::_algoCompatibilityCheck == true). This function is used by the Runner.
     */
    bool isAlgoCompatible(const AllParameters& allP_tmp) const;

    std::string getSetAttributeAsString() const;


    /// Display all attributes (if flagHelp == true, display all help info)
    void display(std::ostream &os , bool flagHelp = false);

    /// Display all attributes
    void displayHelp(const std::string &helpSubject , bool devHelp, std::ostream &os);

    // Include set and get methods from NOMAD 3 for backwards compatibility
    #include "../Param/ParametersNomad3.hpp"

};

#include "../nomad_nsend.hpp"


#endif // __NOMAD_4_0_ALLPARAMETERS__

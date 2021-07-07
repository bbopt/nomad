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
#ifndef __NOMAD_4_0_PARAMETERS__
#define __NOMAD_4_0_PARAMETERS__

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <typeindex>
#include <typeinfo>

#include "../Math/Double.hpp"
#include "../Math/Point.hpp"
#include "../Math/ArrayOfPoint.hpp"
#include "../Type/ListOfVariableGroup.hpp"
#include "../Type/DirectionType.hpp"
#include "../Param/AttributeFactory.hpp"
#include "../Param/ParameterEntries.hpp"
#include "../nomad_platform.hpp"


#include "../nomad_nsbegin.hpp"

typedef std::shared_ptr<Attribute> SPtrAtt;

/// Comparator of shared_ptr<Attribute>
struct lessThanAttribute : public std::binary_function<SPtrAtt, SPtrAtt, bool>
{
    bool operator()(SPtrAtt lhs, SPtrAtt rhs) const
    {
        return (lhs->getName() < rhs->getName());
    }
};

/**
 An attribute set contains shared pointers to Attribute and a lessThanAttribute comparison/binary function for ordering.
 */
typedef std::set<SPtrAtt, lessThanAttribute> AttributeSet;

/**
 The attribute definition structure contains only strings to be interpreted, registered (see Parameters::_definition and Parameters::registerAttributes) and used to set Nomad parameters attributes.
 Each element of the structure except _type and _defaultValue encodes a meta data information of an Attribute as a string.
 */
struct AttributeDefinition
{
    std::string _name; ///< The name is used to get the value getAttributeValue(name)
    std::string _type;  ///< The types of attributes are "registered" by Parameters::registerAttributes
    std::string _defaultValue; ///< The default value must be compatible with the type
    std::string _shortInfo;   ///< Brief description
    std::string _helpInfo;   ///< Info for displaying help
    std::string _keywords;  ///< Keyword entries (separated words) used by help
    std::string _algoCompatibilityCheck;  ///< Encode Attribute::_algoCompatibilityCheck as a string.
    std::string _restartAttribute; ///< Encode Attribute::_restartAttribute as a string.
    std::string _uniqueEntry; ///< Encode Attribute::_uniqueEntry as a string.
};


/// Exception class for an invalid parameter.
class InvalidParameter : public Exception
{
public:
    /// Constructor.
    InvalidParameter(const std::string& file,
                     int line,
                     const std::string& msg)
      : Exception(file, line, msg)
    {
        _typeMsg = "Invalid Parameter.";
    }
};

/// Exception class for a parameter that has not been checked.
class ParameterToBeChecked : public Exception
{
public:
    /// Constructor.
    ParameterToBeChecked(const std::string& file,
                     int line,
                     const std::string& msg)
      : Exception(file, line, msg)
    {
        _typeMsg = "Parameter to be checked.";
    }
};


/// Abstract class for the NOMAD parameters.
/**
 Several types of parameters control the execution of Nomad: RunParameters, PbParameters, CacheParameters, DisplayParameters, EvalParameters, EvaluatorControlGlobalParameters and EvaluatorControlParameters. \n

 All the parameters to control a NOMAD run can be obtained from a single parameter file (see Parameters::readParamFileAndSetEntries) or set using the attribute name and value (see the templated function Parameters::setAttributeValue). \n

 After calling one of specific checkAndComply functions, an attribute value can be obtained using the templated Parameters::getAttributeValue with the name as argument. \n
 */
class Parameters
{
private:


    std::ostringstream _streamedAttribute; ///< The attributes in a format ready to be printed with Parameters::getSetAttributeAsString.


protected:
    /*---------*/
    /* Members */
    /*---------*/
    DLL_UTIL_API static ParameterEntries _paramEntries; ///< The set of entries obtained when reading a parameter file.

    std::string _typeName; ///< The type of parameters: ex. Problem, Run


    bool _toBeChecked; ///< Does checkAndComply() need to be called?

    /**
     The set of attributes (no duplicates) and the binary comparison function
     */
    AttributeSet _attributes;

    /**
     Map of attribute names and type name as string.
     Static to the class.
     */
    DLL_UTIL_API static std::map<std::string,std::string> _typeOfAttributes;

    /// Constructors
    /**
     Attributes are assigned by derived object constructors.
     */
    explicit Parameters()
      : _typeName("Unknown"),
        _toBeChecked(true),
        _attributes()
    {
    }

    /**
     Delete: cannot copy all types of parameters (only the selected derived parameters have their copy constructor)
     */
    Parameters(const Parameters& params) = delete ;

    /**
     Delete: cannot copy all types of parameters (only the selected derived parameters have their copy assignement constructor)
     */
    Parameters& operator=(const Parameters& params) = delete ;

    /// Deep copy of parameters
    /**
     This function is used by the copy constructor of selected derived type of parameters
     */
    void copyParameters ( const Parameters& params );

    virtual ~Parameters() {}


    /**
     Called by derived objects to register attribute
     See Parameters::_definition for some details.
     */
    void registerAttributes ( const std::vector<AttributeDefinition> & attributeDef );

    // Get/Set

    /// Get all the attributes that are set on the current parameters instance
    const AttributeSet getAttributes() const;

    /// Get non interpreted entries
    std::shared_ptr<ParameterEntry> getNonInterpretedParamEntry() const { return _paramEntries.findNonInterpreted(); }
    std::vector<std::shared_ptr<ParameterEntry>> getAllNonInterpretedParamEntries() const { return _paramEntries.findAllNonInterpreted(); }

    /**
     \return \c true if this paramName is found in _paramEntries, \c false otherwise.
     */
    bool isSetByUser(const std::string& paramName) const;

    /// The definition of attributes as a series of string.
    /**
     An attribute is defined by some meta data provided in the structure AttributeDefinition. This structure contains a series of string that are translated into Attribute -s- when calling Parameters::registerAttributes. \n

     This definition is provided as an "initializer list" of a vector of several strings: _definition = { {"attribute1_name","attribute1_type","attribute1_defaultValue",....},{"attribute2_name","attribute2_type",....},...} in the header files for specific attribute definition files: cacheAttributesDefinition.hpp, displayAttributesDefinition.hpp, evalAttributesDefinition.hpp, evaluatorControlAttributesDefinition.hpp, evaluatorControlGlobalAttributesDefinition.hpp, pbAttributesDefinition.hpp and runAttributesDefinition.hpp. These files are programmatically created from their equivalent text files (*.txt) by the WriteAttributeDefinition.exe binary. This pre-processing task is automatically performed (by makefile) if an attribute definition is modified in a text file.\n

     */
    std::vector<AttributeDefinition> _definition;


public:

    /**
     Helper for Parameters::readParamLine OR called by AllParameters::read for each specific type of Parameters such as RunParameters, PbParameters, ....
     */
    void readEntries(const bool overwrite = false, std::string problemDir="");

    /**
     Read a single parameter given in a single line.
     Helper for AllParameters::readParamLine
     */
    void readParamLine(const std::string &line,
                       bool overwrite = false);

    /**
     Read a single parameter given in a single line.
     Helper for Parameters::readParamFileAndSetEntries
     */
    static void readParamLine(const std::string &line,
                              const std::string &paramFile,
                              const int line_number,
                              bool overwrite = false);


    /**
     \return \c true if we need to call checkAndComply().
     */
    bool toBeChecked() const;

    /// Get the names of all the Attributes for this Parameter instance
    std::vector<std::string> getAttributeNames() const;

    /**
     Return the name of a parameter type. Ex. Problem, Run
     */
    std::string getTypeName() const { return _typeName; }

    void resetToDefaultValues() noexcept ;
    void resetToDefaultValue(const std::string& paramName);

    /**
     Read the parameter file. Each line of the file is sent to Parameters::readParamLine.
     Called by AllParameters::read.
     */
    static void readParamFileAndSetEntries(const std::string &paramFile, bool overwrite = false , bool resetAllEntries = false );

    /// Erase all entries
    static void eraseAllEntries() { _paramEntries.eraseAll() ; }


    /// Display all attributes
    void display(std::ostream &os ,
                 bool helpInfo ) ;

    // Display help on given subject
    void displayHelp(const std::string & helpSubject,
                     bool devHelp,
                     std::ostringstream & ossBasic,
                     std::ostringstream & ossAdvanced ) ;




private:
    /**
     Must be implemented by derived object. Used by constructor.
     */
    virtual void init() = 0 ; // Pure virtual

    /// Helper for read
    void readValuesAsArray(const ParameterEntry &pe, ArrayOfDouble &array);

    /// Helper for read
    size_t readValuesForArrayOfPoint(const ParameterEntry &pe, Point &point);

    /// Helper for read
    NOMAD::ArrayOfPoint readPointValuesFromFile(const std::string& pointFile);

    /// Helper for read
    size_t readValuesForVariableGroup(const ParameterEntry &pe, VariableGroup &vg);

protected:
    /*-------------------*/
    /* registerAttribute   */
    /*-------------------*/
    // This template function implementation must be in the header to be available in the library
    /**
     Insert a new attribute by its name and type (a duplicate attribute triggers exception). Call for the AttributeFactor::Create function to create an Attribute with a templated type. Register the type (as a string) and the name of an Attribute into Parameters::_typeOfAttributes.
     See Parameters::_definition for details.
     */
    template<typename T, typename ... ARGS>
    void registerAttribute( std::string name,
                         T default_value,
                         bool algoCompatibilityCheck,
                         bool restartAttribute,
                         bool uniqueEntry,
                         ARGS&&... infoArgs)
    {

        if ( std::is_reference<T>::value
            || std::is_pointer<T>::value
            || std::is_const<T>::value )
        {
            std::string err = "Attribute " + name;
            err += " must be of regular type (no pointer or reference).";
            throw Exception(__FILE__,__LINE__, err);
        }

        NOMAD::toupper(name);

        auto attribute = AttributeFactory{}.Create<T>(name, default_value,
                                                      algoCompatibilityCheck, restartAttribute, uniqueEntry,
                                                      std::forward<ARGS>(infoArgs)...);

        auto ret = _attributes.insert(attribute);

        if ( !ret.second )
        {
            std::string err = "Attribute " + name + " is already in set of attributes.";
            throw Exception(__FILE__,__LINE__, err);
        }

        std::string typeName = typeid(T).name();
        auto nameTypePair = std::pair<std::string, std::string>(name, typeName);
        auto ret_t = _typeOfAttributes.insert(nameTypePair);

        if ( !ret_t.second )
        {
            // Attribute is already in the map of attribute types.
            // Verify the type is the same.
            if ( _typeOfAttributes[name] != typeName )
            {
                std::string err = "Trying to add attribute " + name;
                err += " with type " + typeName;
                err += " which is different from registered type " + _typeOfAttributes[name];
                throw Exception(__FILE__,__LINE__, err);
            }
        }

    }

    const std::string& getAttributeType(const std::string& name)
    {
        auto namecaps = name;
        NOMAD::toupper(namecaps);
        return _typeOfAttributes[namecaps];
    }

    SPtrAtt getAttribute(std::string name) const;


    // getSpValue: value is of correct type for parameter name.
    template<typename T> const T&
    getSpValue(const std::string &name, bool flagCheckException, bool flagDefault = false) const
    {
        // Get attribute from which to get value
        SPtrAtt att;
        att = getAttribute(name);

        if (nullptr == att)
        {
            // At this point, Verify att is non-null.
            std::string err = "getAttributeValue: attribute " + name + " does not exist";
            throw Exception(__FILE__, __LINE__, err);
        }

        // Verify attribute type
        // Must use map access with "at" (not []) because the function is const
        const std::string &typeTName = typeid(T).name();
        if (typeTName != _typeOfAttributes.at(name))
        {
            std::string err = "In getAttributeValue<T> the attribute ";
            err += name + " is not of type T = " + typeTName;
            throw Exception(__FILE__,__LINE__, err);
        }

        // Note: we use getAttributeValue in init() and checkAndComply().
        // We cannot verify toBeChecked() here. It has to be verified at
        // another level.

        // Get value from attribute
        // Dynamic cast to the selected TypeAttribute
        std::shared_ptr<TypeAttribute<T>> sp = std::dynamic_pointer_cast<TypeAttribute<T>>(att);

        if (flagDefault)
        {
            // Get initial value
            return sp->getInitValue();
        }
        else
        {
            // All attributes except DIMENSION must be checked before accessing the value
            if ( _toBeChecked && flagCheckException && name != "DIMENSION" )
            {
                std::string err = "In getAttributeValue<T> the attribute ";
                err += name + " has not been checked";
                throw ParameterToBeChecked(__FILE__,__LINE__, err);
            }
            return sp->getValue();
        }
    }


    /*----------------------------*/
    /* getAttributeValueProtected */
    /*----------------------------*/
    template <typename T>
    struct type{};

    /**
     Intermediate step when calling getAttributeValue
     */
    template<typename T> const T&
    getAttributeValueProtected(const std::string &name, type<T>, bool flagCheckException, bool flagDefault = false) const
    {
        // Generic case
        return getSpValue<T>(name, flagCheckException, flagDefault);
    }

    /**
     Template specialization when getting attribute value for Point.
     This handles the special case where the type is an ArrayOfPoint, but user asks to
     return a Point. So we get the "value" for an ArrayOfPoint and return its first element.
     */
    template<typename T> const Point&
    getAttributeValueProtected(const std::string &name, type<Point>, bool flagCheckException, bool flagDefault = false) const
    {
        auto namecaps = name;
        NOMAD::toupper(namecaps);
        if (typeid(ArrayOfPoint).name() == _typeOfAttributes.at(namecaps))
        {
            // Special case: Attribute type is an ArrayOfPoint, but user asks to
            // return a Point.
            // Get the ArrayOfPoint and return its first element.
            const ArrayOfPoint & aop = getSpValue<ArrayOfPoint>(namecaps, flagCheckException, flagDefault);
            if (aop.size() >= 1)
            {
                return aop[0];
            }
            else
            {
                std::string err = "In getAttributeValue<Point> : the attribute " + name;
                err += " contains no point.";
                throw Exception(__FILE__,__LINE__, err);
            }
        }

        // Default behaviour
        return getSpValue<T>(namecaps, flagCheckException, flagDefault);
    }


    // This template function implementation must be in the header to be
    // available in the library
    template<typename T> const T&
    getAttributeValueProtected(const std::string &name, bool flagCheckException, bool flagDefault = false) const
    {
        // Call method that has type as argument.
        // This will differ the cal for T = Point and other types.
        // Ref: https://www.fluentcpp.com/2017/08/11/how-to-do-partial-template-specialization-in-c/
        return getAttributeValueProtected<T>(name, type<T>{}, flagCheckException, flagDefault);
    }


public:


    // This template function implementation must be in the header to be
    // available in the library
    /**
     This template function is used to obtain the attribute value for any type of Nomad specific Parameters (RunParameters, PbParameters, CacheParameters) if the parameters have been checked. Otherwise an exception is triggered.
     */
    template<typename T> const T&
    getAttributeValue(const std::string &name, bool flagDefault = false) const
    {
        auto namecaps = name;
        NOMAD::toupper(namecaps);
        return getAttributeValueProtected<T>(namecaps,true,flagDefault);
    }

    /**
     This function is called by AllParameters::readParamLine to identify the type of an attribute given in a file. If the attribute has not been registered, an error message is displayed but no exception is triggered and the execution continues.
     */
    bool isRegisteredAttribute(const std::string &name) const
    {
        auto att = getAttribute(name);
        if ( att == nullptr )
            return false;
        else
            return true;

    }

    /**
     Test if algo compatible registered attributes for *this and p are equal
     */
    bool isAlgoCompatible (const Parameters* p);


    std::string getSetAttributeAsString ( void ) const
    {
        return _streamedAttribute.str();
    }


    /*
     This template function implementation must be in the header to be available in the library
     */
    /**
     \param name    The attribute's name -- \b IN.
     \return        \c true if the attribute identified by its name has been registered and the current value is the default value. It returns \c false otherwise.
     */
    template<typename T> bool isAttributeDefaultValue(const std::string &name) const
    {

        std::string typeTName = typeid(T).name();
        auto namecaps = name;
        NOMAD::toupper(namecaps);
        auto att = getAttribute(name);

        // Must use map access with "at" (not []) because the function is const
        if (typeTName != _typeOfAttributes.at(namecaps))
        {
            std::string err = "In isAttributeDefaultValue<T> : the attribute " + name;
            err += " is not of type T = " + typeTName;
            throw Exception(__FILE__,__LINE__, err);
        }

        // Dynamic cast to the selected TypeAttribute
        std::shared_ptr<TypeAttribute<T>> sp = std::dynamic_pointer_cast<TypeAttribute<T>>( att );

        return sp->isDefaultValue();
    }



    /**
     Function to manage different set cases.
     Default: value is of the correct type for parameter designed by name.
     In other cases, value is of a different type, and a type conversion
     is needed. With this approach, we can prevent some type conversion (for example, Point to ArrayOfDouble).
     */
    template<typename T>
    void setSpValueDefault(const std::string& name, T value)
    {
        auto att = getAttribute(name);

        if (nullptr == att)
        {
            // At this point, Verify att is non-null.
            std::string err = "setSpValueDefault: attribute " + name + " does not exist";
            throw Exception(__FILE__, __LINE__, err);
        }
        // Dynamic cast to the selected TypeAttribute
        std::shared_ptr<TypeAttribute<T>> sp = std::dynamic_pointer_cast<TypeAttribute<T>>( att );

        const std::string &typeTName = typeid(T).name();
        if (typeTName != _typeOfAttributes[name])
        {
            std::string err = "setSpValueDefault<T> : the attribute " + name;
            err += " is of type " + _typeOfAttributes[name];
            err += " and not of type T = " + typeTName;
            throw Exception(__FILE__,__LINE__, err);
        }
        if (!sp->uniqueEntry())
        {
            if (typeid(ArrayOfString).name() == _typeOfAttributes.at(name))
            {
                // Special case for DISPLAY parameter
                ArrayOfString* valueAsAos = (ArrayOfString*)(&value);
                ArrayOfString* aos = (ArrayOfString*)(&sp->getValue());
                for (size_t i = 0; i < valueAsAos->size(); i++)
                {
                    aos->add((*valueAsAos)[i]);
                }
                // Convert back to T
                value = *((T*)(aos));
            }
        }
        sp->setValue(value);

        // Keep a track of all the non default set (used by runner to display attributes)
        if ( ! sp->isDefaultValue() )
        {
            _streamedAttribute << " [ ";
            sp->display( _streamedAttribute , false /* no short info*/ );
            _streamedAttribute << " ] " ;
        }

    }

    /**
     Called by setAttributeValue. Generic template function for attribute of type T. The generic version of the function is called when no template specialization is available (partial specialization).
     */
    template<typename T>
    void setSpValue(const std::string& name, T value)
    {
        setSpValueDefault(name, value);
    }

    /**
     Overload of setSpValue for Point -> ArrayOfPoint case.
     Value is of type Point, and it might need to be converted to an ArrayOfPoint for parameter name.
     */
    void setSpValue(const std::string& name, Point value)
    {
        if (typeid(ArrayOfPoint).name() == _typeOfAttributes.at(name))
        {
            // Special case: Attribute type is an ArrayOfPoint, but user sets
            // a Point.
            // Create an ArrayOfPoint and set its first element to the
            // given point value.
            ArrayOfPoint aop;
            aop.push_back(value);

            setSpValue(name, aop);
        }
        else
        {
            // Use default behaviour
            setSpValueDefault(name, value);
        }
    }

    /**
     Overload of setSpValue for DirectionType -> DirectionTypeList case.
     Value is of type DirectionType, and it might need to be converted to an DirectionTypeList for parameter name.
     */
    void setSpValue(const std::string& name, DirectionType value)
    {
        if (typeid(DirectionTypeList).name() == _typeOfAttributes.at(name))
        {
            // Special case: Attribute type is an DirectionTypeList, but user sets
            // a DirectionType.
            // Create an DirectionTypeList and set its first element to the
            // given point value.
            DirectionTypeList dirTypeList;
            dirTypeList.push_back(value);

            setSpValue(name, dirTypeList);
        }
        else
        {
            // Use default behaviour
            setSpValueDefault(name, value);
        }
    }

    /**
     Overload of setSpValue for std::string -> ArrayOfString case.
     Value is of type std::string, and it might need to be converted to an ArrayOfString for parameter name.
     */
    void setSpValue(const std::string& name, std::string value)
    {
        if (typeid(ArrayOfString).name() == _typeOfAttributes.at(name))
        {
            // Special case: Attribute type is an ArrayOfString, but user sets
            // a string.
            // Create an ArrayOfString and set its first element to the
            // given string value.
            ArrayOfString aos;
            aos.add(value);

            setSpValue(name, aos);
        }
        else
        {
            // Use default behaviour
            setSpValueDefault(name, value);
        }
    }

    /**
     Overload of setSpValue for int -> size_t case.
     Value is of type int, and it might need to be converted to a size_t for parameter name.
     */
    void setSpValue(const std::string& name, int value)
    {
        if (typeid(size_t).name() == _typeOfAttributes.at(name))
        {
            // Special case: Attribute type is a size_t, but user sets
            // an int.
            // If the int is negative, use INF.
            size_t st = (value < 0) ? INF_SIZE_T : size_t(value);
            setSpValue(name, st);
        }
        else
        {
            // Use default behaviour
            setSpValueDefault(name, value);
        }
    }


    // This template function implementation must be in the header to be
    // available in the library
    /**
     Function called to set the value of a registered attribute. To have access to the updated value the checkAndComply function must be called (for example: RunParameters::checkAndComply).
     */
    template<typename T>
    void setAttributeValue(const std::string& name, T value)
    {
        // This check should not be necessary because setAttributeValue<T>
        // cannot access T as a reference or a pointer

        if (std::is_reference<T>::value || std::is_pointer<T>::value
            || std::is_const<T>::value)
        {
            std::string err = "setAttributeValue: attribute " + name;
            err += " must be of regular type (no pointer or reference).";
            throw Exception(__FILE__, __LINE__, err);
        }

        auto namecaps = name;
        NOMAD::toupper(namecaps);
        setSpValue(namecaps, value);

        _toBeChecked = true;
    }

    auto getAttributeShortInfo(const std::string &name)
    -> decltype(getAttribute(name)->getShortInfo())
    {
        return getAttribute(name)->getShortInfo();
    }

    auto getAttributeHelpInfo(const std::string &name) const
        -> decltype(getAttribute(name)->getHelpInfo())
    {
        return getAttribute(name)->getHelpInfo();
    }

    /// Access to the number of attributes
    size_t getAttributeSetSize() const { return _attributes.size() ;}

protected:
    // Helpers
    void checkFormat1(const std::shared_ptr<ParameterEntry> pe) const;
    void checkFormatNbEntries(const std::shared_ptr<ParameterEntry> pe, const size_t nbEntries) const;
    void checkFormatBool(const std::shared_ptr<ParameterEntry> pe) const;
    void checkFormatSizeT(const std::shared_ptr<ParameterEntry> pe, size_t &sz) const;
    void checkFormatAllSizeT(const std::shared_ptr<ParameterEntry> pe) const;
    void checkFormatInt(const std::shared_ptr<ParameterEntry> pe, int &i) const;
    void checkFormatString(const std::shared_ptr<ParameterEntry> pe) const;
    void checkFormatArrayOfString(const std::shared_ptr<ParameterEntry> pe) const;
    void checkFormatDouble(const std::shared_ptr<ParameterEntry> pe, Double &d) const;

    void checkInfo() const ;

};

#include "../nomad_nsend.hpp"


#endif // __NOMAD_4_0_PARAMETERS__

//
//  TypeAttribute.hpp
//  nomad
//
//  Created by Christophe Tribes on 2017-12-08.
//  Copyright (c) 2017 GERAD. All rights reserved.
//

#ifndef __NOMAD400_TYPEATTRIBUTE__
#define __NOMAD400_TYPEATTRIBUTE__

#include "../Param/Attribute.hpp"
#include "../Util/defines.hpp"

#include "../nomad_nsbegin.hpp"

/// Class for an Attribute of a templated type.
/**
 TypeAttribute is the templated class that derives from an Attribute and completes the values (current and initial) and type for a Nomad parameter. A specific type of attribute is obtained by calling the templated AttributeFactory::Create function.
 */
template <typename T>
class TypeAttribute : public Attribute {
private:
    T _value;
    const T _initValue;

public:

    TypeAttribute ( ) = delete;

    /// Constructor call by the AttributeFactory::Create
    template<typename ... ARGS>TypeAttribute(std::string Name,
                                             T initValue,
                                             bool algoCompatibilityCheck,
                                             bool restartAttribute,
                                             bool uniqueEntry,
                                             ARGS&& ... infoArgs):
    Attribute(Name,
              algoCompatibilityCheck,
              restartAttribute,
              uniqueEntry,
              std::forward<ARGS>(infoArgs)...),
    _value(initValue),
    _initValue(initValue)
    {}


    virtual ~TypeAttribute() {}

    const T& getValue() const { return _value; }
    const T& getInitValue() const { return _initValue; }
    void setValue(T v) { _value=v; }

    /**
     \return \c true if current value equals initial value, \c false otherwise.
     */
    bool isDefaultValue() const { return (_value==_initValue) ; }

    void display( std::ostream & out , bool flagShortInfo = true ) const
    {

        out << _name << " " << _value;
        if ( flagShortInfo && _shortInfo.size() > 0 )
            out << " (" << _shortInfo << ")";
    }

    void resetToDefaultValue() noexcept { _value = _initValue;}

};


#include "../nomad_nsend.hpp"

#endif // __NOMAD400_TYPEATTRIBUTE__

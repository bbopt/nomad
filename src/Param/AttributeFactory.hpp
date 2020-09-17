//
//  AttributeFactory.hpp
//  nomad
//
//  Created by Christophe Tribes on 2017-12-08.
//  Copyright (c) 2017 GERAD. All rights reserved.
//

#ifndef __NOMAD400_ATTRIBUTEFACTORY__
#define __NOMAD400_ATTRIBUTEFACTORY__

#include "../Param/TypeAttribute.hpp"

#include "../nomad_nsbegin.hpp"


/// Factory to make Attribute with a variable number of arguments passed to the Create function.
/**
 The Parameters::registerAttribute is in charge of calling the Create function.
 */
struct AttributeFactory {

public:

    template <typename T, typename ... ARGS> std::shared_ptr<Attribute> Create(std::string Name,
               T initValue,
               bool algoCompatibilityCheck,
               bool restartAttribute,
               bool uniqueEntry,
               ARGS && ... infoArgs )
    {
        return std::make_shared<TypeAttribute<T>>(Name,
                                                  initValue,
                                                  algoCompatibilityCheck,
                                                  restartAttribute,
                                                  uniqueEntry,
                                                  std::forward<ARGS>(infoArgs)...);

    }

};

#include "../nomad_nsend.hpp"
#endif  // __NOMAD400_ATTRIBUTEFACTORY__

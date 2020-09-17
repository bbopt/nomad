/**
 \file   ParameterEntries.cpp
 \brief  Parameter entries (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-05
 \see    ParameterEntries.hpp
 */
#include "../Param/ParameterEntries.hpp"

/*--------------------------------------------*/
/*                 destructor                 */
/*--------------------------------------------*/
NOMAD::ParameterEntries::~ParameterEntries ( void )
{
    // No need to delete anything since smart pointers are used.
}

/*--------------------------------------------*/
/*      finds a specific entry in the set     */
/*--------------------------------------------*/
std::shared_ptr<NOMAD::ParameterEntry> NOMAD::ParameterEntries::find ( const std::string & name ) const
{
    auto it = _entries.find ( std::make_shared<NOMAD::ParameterEntry>(name) );
    if ( it != _entries.end() )
        return (*it);
    return nullptr;
}

/*----------------------------------------*/
/*      inserts an entry into the set     */
/*----------------------------------------*/
void NOMAD::ParameterEntries::insert(std::shared_ptr<NOMAD::ParameterEntry> entry )
{
    std::shared_ptr<NOMAD::ParameterEntry> cur = find ( entry->getName() );
    if ( cur )
    {
        entry->setUnique ( false );
        cur->setUnique   ( false );
        while ( cur->getNext() )
            cur = cur->getNext();
        cur->setNext ( entry );
    }
    _entries.insert ( entry );
}


/*----------------------------------------*/
/*       erase an entry from the set      */
/*----------------------------------------*/
void NOMAD::ParameterEntries::erase(std::shared_ptr<NOMAD::ParameterEntry> entry )
{
    _entries.erase ( entry );
}

/*----------------------------------------*/
/*       erase all entries                */
/*----------------------------------------*/
void NOMAD::ParameterEntries::eraseAll ( )

{
    _entries.clear();
}

/*----------------------------------------*/
/*       find a non-interpreted entry     */
/*----------------------------------------*/
std::shared_ptr<NOMAD::ParameterEntry> NOMAD::ParameterEntries::findNonInterpreted ( void ) const
{
    for (auto it : _entries)
    {
        if ( !it->hasBeenInterpreted() )
        {
            return it;
        }
    }
    return nullptr;
}


std::vector<std::shared_ptr<NOMAD::ParameterEntry>> NOMAD::ParameterEntries::findAllNonInterpreted() const
{
    std::vector<std::shared_ptr<NOMAD::ParameterEntry>> allNonInterp;
    for (auto it : _entries)
    {
        if ( !it->hasBeenInterpreted() )
        {
            allNonInterp.push_back(it);
        }
    }
    return allNonInterp;
}


/*--------------------------------------------*/
/*                   display                  */
/*--------------------------------------------*/
void NOMAD::ParameterEntries::display(std::ostream &out) const
{
    auto end = _entries.end();
    for (auto it = _entries.begin() ; it != end ; ++it )
    {
        out << **it << std::endl;
    }
}

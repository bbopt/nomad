/**
 \file   ParameterEntries.hpp
 \brief  Parameter entries (headers)
 \author Sebastien Le Digabel
 \date   2010-04-05
 \see    ParameterEntries.cpp
 */
#ifndef __NOMAD400_PARAMETERENTRIES__
#define __NOMAD400_PARAMETERENTRIES__

#include <set>
#include <vector>
#include "../Param/ParameterEntry.hpp"
#include "../Util/Uncopyable.hpp"

#include "../nomad_nsbegin.hpp"

/// Store and manage ParameterEntry.
/**
 - All the ParameterEntry objects are stored as a multiset with a comparison operator based on name.
 - ParameterEntries stores all the parameters provided in a file.
 */
class ParameterEntries : private Uncopyable {

private:

    /// List of ParameterEntry objects (the entries).
    std::multiset<std::shared_ptr<ParameterEntry>, ParameterEntryComp> _entries;

public:

    /// Constructor.
    explicit ParameterEntries ( void ) {}

    /// Destructor.
    virtual ~ParameterEntries ( void );

    /// Find a specific entry in a set.
    /**
     \param  name The name of the wanted ParameterEntry object -- \b IN.
     \return      A pointer to the ParameterEntry object if it
     has been found in the list of entries,
     or \c nullptr otherwise.
     */
    std::shared_ptr<ParameterEntry> find ( const std::string & name ) const;

    /// Insert a new entry in the list of entries.
    /**
     \param entry A pointer to the new ParameterEntry object -- \b IN.
     */
    void insert(std::shared_ptr<ParameterEntry> entry);

    /// Erase an entry from the list of entries.
    /**
     \param entry A pointer to the ParameterEntry object to erase -- \b IN.
     */
    void erase(std::shared_ptr<ParameterEntry> entry);

    /// Erase all entries.
    void eraseAll();

    /// Find a non-interpreted entry.
    /**
     \return A pointer to the first ParameterEntry that has not been
     interpreted so far,
     or \c nullptr if all entries have already been interpreted.
     */
    std::shared_ptr<ParameterEntry> findNonInterpreted ( void ) const;

    /// Find all non-interpreted entries
    std::vector<std::shared_ptr<ParameterEntry>> findAllNonInterpreted() const;

    /// Display.
    /**
     \param out The std::ostream object -- \b IN.
     */
    void display(std::ostream &out) const;
};

/// Display a ParameterEntries object.
/**
 \param out The std::ostream object -- \b IN.
 \param e   The ParameterEntries object to be displayed -- \b IN.
 \return    The std::ostream object.
 */
inline std::ostream& operator<< (std::ostream &out,
                                 const ParameterEntries &e)
{
    e.display(out);
    return out;
}

#include "../nomad_nsend.hpp"

#endif // __NOMAD400_PARAMETERENTRIES__

/**
 \file   ParameterEntry.cpp
 \brief  Parameter entry (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-05
 \see    ParameterEntry.hpp
 */

#include <algorithm>    // for for_each
#include "../Param/ParameterEntry.hpp"

/*-----------------------------------*/
/*  . constructor (from a string)    */
/*  . ignores all entries after '#'  */
/*-----------------------------------*/
NOMAD::ParameterEntry::ParameterEntry ( const std::string & entry           ,
                                         bool                removeComments   )
: _ok                 ( false ) ,
_unique               ( true  ) ,
_next                 ( nullptr  ) ,
_paramFile            ( "" ),
_line                 ( 0 ) ,
_hasBeenInterpreted   ( false )
{
    //int                i , idst;
    std::string        s;
    std::istringstream in ( entry );
    in >> _name;

    if (_name.size()==0)
        return;

    if ( removeComments && _name[0] == '#' )
        _name.clear();
    else
    {


        NOMAD::toupper ( _name );

//         bool stats_file_name_read = false;

        while ( true )
        {
            in >> s;

            if ( in.fail() )
                break;

            // comment:
            if ( removeComments && s[0]=='#' )
                break;

            // string with quotes:
            //bool had_quotes = false;
            if ( s[0] == '\"' || s[0] == '\'' )
            {

                //had_quotes = true;
                char quote = s[0];

                s.erase ( s.begin() );

                if ( s[s.size()-1] == quote )
                    s.resize ( s.size() - 1 );

                else
                {

                    std::string ss;
                    getline ( in , ss , quote );

                    if ( in.fail() || !in.good())
                    {
                        _ok = false;
                        return;
                    }

                    s = s + ss;
                }
            }

            // vector:
            if ( s.size() > 1 && ( s[0] == '[' || s[0] == '(' ) )
            {
                _values.push_back ( s[0]=='[' ? "[" : "(" );
                s.erase(s.begin());
            }
            size_t  sm1 = s.size() - 1;
            char c   = s[sm1];
            if ( s.size() > 1 && ( c == ']' || c == ')' ) )
            {
                s.resize(sm1);
                _values.push_back (s);
                _values.push_back ( c==']' ? "]" : ")" );
                continue;
            }

            // other values:
            _values.push_back ( s );
        }

        if ( !_values.empty() )
            _ok = true;
    }
}

/*------------------------------*/
/*             display          */
/*------------------------------*/
void NOMAD::ParameterEntry::display(std::ostream &out) const
{
    if ( _ok )
    {
        out << _name << ": ";
        std::list<std::string>::const_iterator end = _values.end();
        for ( std::list<std::string>::const_iterator it = _values.begin() ; it != end ; ++it )
            out << "[" << *it << "] ";
    }
}

std::string NOMAD::ParameterEntry::getAllValues( void ) const
{
    std::string all;

    std::for_each ( _values.begin() , _values.end() ,[&](const std::string & s){ all+=s+" ";} );

    return all;

}


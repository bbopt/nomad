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

    all.pop_back(); // Remove trailing white space
    return all;

}


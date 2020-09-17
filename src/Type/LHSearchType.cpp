#include "../Type/LHSearchType.hpp"
#include "../Util/ArrayOfString.hpp"
#include "../Util/Exception.hpp"
#include "../Util/utils.hpp"

NOMAD::LHSearchType::LHSearchType(const std::string& entries)
: _enable(false),
  _lhsearch0(0),
  _lhsearch1(0)
{
    if (!entries.empty())
    {
        NOMAD::ArrayOfString aos(entries);
        if (aos.size() != 2)
        {
            std::string err = "LHSearchType must have 2 entries, got ";
            err += std::to_string(aos.size());
            err += "( " + entries + " )";
            throw Exception(__FILE__, __LINE__, err);
        }
        int lhsearch0, lhsearch1;
        std::string s0 = aos[0];
        std::string s1 = aos[1];
        NOMAD::atoi(s0, lhsearch0);
        NOMAD::atoi(s1, lhsearch1);
        _lhsearch0 = lhsearch0;
        _lhsearch1 = lhsearch1;

        _enable = (_lhsearch0 != 0 || _lhsearch1 != 0);
    }
}

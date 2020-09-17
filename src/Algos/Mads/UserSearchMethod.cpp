
#include "../../Algos/Mads/UserSearchMethod.hpp"

void NOMAD::UserSearchMethod::init()
{
    _name = "User Search Method";

    // Query the enabling parameter here
    //auto userSearch = _runParams->getAttributeValue<bool>("USER_SEARCH");
    bool userSearch = false;
    setEnabled(userSearch);

}

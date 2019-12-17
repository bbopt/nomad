//////////// THIS FILE MUST BE CREATED BY EXECUTING WriteAttributeDefinitionFile ////////////
//////////// DO NOT MODIFY THIS FILE MANUALLY ///////////////////////////////////////////////

#ifndef __NOMAD400_CACHEATTRIBUTESDEFINITION__
#define __NOMAD400_CACHEATTRIBUTESDEFINITION__

_definition = {
{ "MAX_CACHE_SIZE",  "size_t",  "INF",  " Termination criterion on the number of evaluation points stored in the cache ",  " \n  \n . The program terminates as soon as the cache reaches this size. \n  \n . Argument: one positive integer (expressed in number of evaluation points). \n  \n . Example: MAX_CACHE_SIZE 10000 \n  \n . Default: INF\n\n",  "  advanced termination cache  "  , "false" , "false" , "true" },
{ "CACHE_FILE",  "std::string",  "",  " Cache file name ",  " \n  \n . Cache file. If the specified file does not exist, it will be created. \n  \n . Argument: one string. \n  \n . If the string is empty, no cache file will be created. \n  \n . Points already in the cache file will not be reevaluated. \n  \n . Example: CACHE_FILE cache.txt \n  \n . Default: Empty string.\n\n",  "  basic cache file  "  , "false" , "false" , "true" } };

#endif

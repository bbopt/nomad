//////////// THIS FILE MUST BE CREATED BY EXECUTING WriteAttributeDefinitionFile ////////////
//////////// DO NOT MODIFY THIS FILE MANUALLY ///////////////////////////////////////////////

#ifndef __NOMAD_4_6_CACHEATTRIBUTESDEFINITION__
#define __NOMAD_4_6_CACHEATTRIBUTESDEFINITION__

_definition = {
{ "CACHE_FILE",  "std::string",  "",  " Cache file name ",  " \n  \n . Cache file. If the specified file does not exist, it will be created. \n  \n . Argument: one string. \n  \n . If the string is empty, no cache file will be created. \n  \n . Points already in the cache file will not be reevaluated. \n  \n . Example: CACHE_FILE cache.txt \n  \n . Default: Empty string.\n\n",  "  basic cache file  "  , "false" , "false" , "true" },
{ "CACHE_SIZE_MAX",  "size_t",  "INF",  " Maximum number of evaluation points to be stored in the cache ",  " \n  \n . The cache will be purged from older points if it reaches this number \n   of evaluation points. \n  \n . Argument: one positive integer (expressed in number of evaluation points). \n  \n . Example: CACHE_SIZE_MAX 10000 \n  \n . Default: INF\n\n",  "  advanced cache  "  , "false" , "false" , "true" } };

#endif

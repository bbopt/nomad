//////////// THIS FILE MUST BE CREATED BY EXECUTING WriteAttributeDefinitionFile ////////////
//////////// DO NOT MODIFY THIS FILE MANUALLY ///////////////////////////////////////////////

#ifndef __NOMAD_4_6_RUNATTRIBUTESDEFINITIONIBEX__
#define __NOMAD_4_6_RUNATTRIBUTESDEFINITIONIBEX__

_definition = {
{ "USE_IBEX",  "bool",  "false",  " Boolean to determine if we want to use the functionnalities of IBEX ",  " \n  \n . Argument : bool \n  \n . Determine if you want to use the fonctionnalities of IBEX \n  \n . Default: false\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" },
{ "SYSTEM_FILE_NAME",  "string",  "-",  " File with the constraints  ",  " \n  \n . Minibex file name, describing the system (i.e constraints, variables...) of the problem. \n  \n . See the documentation here for more detail : http://www.ibex-lib.org/doc/minibex.html on how to create it. \n  \n . No default value.\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" },
{ "SET_FILE",  "bool",  "false",  " Boolean to determine if the file of the set is already created ",  " \n  \n . Argument : bool \n  \n . Determine if the Set of the problem is already created. \n  \n . Default: false\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" },
{ "SET_FILE_NAME",  "string",  "-",  " File to load with the set  ",  " \n  \n . Argument : string \n  \n . Name of the Set file. \n  \n . No need to be provided if SET_FILE = false. \n  \n . No default value.\n\n",  "  advanced project algorithm ibex snap  "  , "true" , "true" , "true" } };

#endif

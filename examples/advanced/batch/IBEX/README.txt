This file explains how to use NOMAD with IBEX thanks to an exemple.

First, when you have the BB file of the problem you want to solve you have to create a "system" file that defines the problems i.e. defines the variables, the constraints, the initial domain etc.... In this case, the system file is "system_anneau.txt". 

Then, you have to use the files "create_set" or "create_set_with_volume" to create the Set (i.e. the caracterization of the domain with boxes) from the system file. Then you can give these information to NOMAD thanks to the param.txt file. In this case, in the param.txt file, the line 15 claims that you want to use IBEX, the line 16 claims that you already created the Set and the line 17 provides the name of the file that contains the Set.

However, if you don't want to create the Set, you can only provide the system file name thanks to the attribute SYSTEM_FILE_NAME and the Set will be created automatically in NOMAD. The argument SET_FILE (bool) defines if you already created the Set or if you only want to pass the system file.

Two different use of NOMAD with IBEX are available. You can either call the Set in the BB file and project directly in the this file (see bb.cpp) or you can project the point in NOMAD, and in this case you give the informations needed thanks to the param.txt file, such as explained above.

IMPORTANT: Modify the path to IBEX in bb.cpp before building the example

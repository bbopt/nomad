#!/bin/sh

# Show txt files except for deprecated parameters;
# Extract informations in CSV format;
# Sort (by parameter name).
ls $NOMAD_HOME/src/Attribute/*.txt | grep -v deprecated | xargs cat | extractParamInfo.sh | sort > $NOMAD_HOME/doc/user_guide/source/allParameters.csv


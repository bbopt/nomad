#!/bin/sh

# Show txt files except for deprecated parameters;
# Extract informations in CSV format;
# Sort (by parameter name).
csvFile="$NOMAD_HOME/doc/user_guide/source/allParameters.csv"
echo "Name , Type , Argument , Short description , Default" > $csvFile
ls $NOMAD_HOME/src/Attribute/*.txt | grep -v deprecated | xargs cat | extractParamInfo.sh | sort >> $csvFile


#!/bin/sh

# Show txt files except for deprecated parameters;
# Extract informations in CSV format;
# Sort (by parameter name).
ls *.txt | grep -v deprecated | xargs cat | extractParamInfo.sh | sort


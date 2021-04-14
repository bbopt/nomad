#!/bin/bash

egrep -v "^#|^ALGO_COMPATIBILITY_CHECK|^RESTART_ATTRIBUTE" $1 | \
    awk 'BEGIN{ biginfopassed = 0;}
    {
        curline = $0;
        lastline = 0;
        if ("" == $0) {
            # do nothing
        }
        else if ("" == name)
        {
            name = $0;
        }
        else if ("" == type)
        {
            type = $0;
        }
        else if ("" == defval)
        {
            defval = $0;
            if ("\"\"" == defval || "-" == defval)
            {
                defval = "No default";
            }
        }
        else if ("" == info)
        {
            info = $0;
            sub(/\\\( */,"",info);
            sub(/ *\\\)/,"",info);
        }
        else if ("\\(" == $0) {
            # Flag used to ignore long info.
            biginfopassed = 0;
        } else if ("\\)" == $0) {
            biginfopassed = 1;
        } else if (1 == biginfopassed && "" == argument) {
            lastline = 1;
            keywords = $0;
            if ($keywords ~ "developer" || $keywords ~ "developper")
            {
                argument = "developer";
            }
            else if ($keywords ~ "advanced")
            {
                argument = "advanced";
            }
            else
            {
                argument = "basic";
            }
        }
        if (1 == lastline) {
            # Here we print the found fields, with commas for CSV format.
            print name","type","argument",\""info"\","defval;

            # Reset values.
            biginfopassed = 0;
            name = "";
            type = "";
            defval = "";
            info = "";
            argument = "";
        }
    }'

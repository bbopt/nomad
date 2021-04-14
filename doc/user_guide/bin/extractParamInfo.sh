#!/bin/bash

egrep -v "^#|^ALGO_COMPATIBILITY_CHECK|^RESTART_ATTRIBUTE|UNIQUE_ENTRY" $1 | \
    awk 'BEGIN{ biginfopassed = 0;}
    {
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
            sub(/\*/,"\\*",defval);
            #sub(/*/,"toto",defval);
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
        } else if (1 == biginfopassed && 0 == lastline) {
            keywords = keywords" "$0;
            lastline = ($keywords ~ "\\\\)" );  # It took me a while to figure this one out!
        }
        if (1 == lastline) {
            # Compute argument (basic / advanced / developer)
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

            # Here we print the found fields, with commas for CSV format.
            print name","type","argument",\""info"\","defval;

            # Reset values.
            biginfopassed = 0;
            name = "";
            type = "";
            defval = "";
            info = "";
            keywords = "";
            lastline = 0;
        }
    }'

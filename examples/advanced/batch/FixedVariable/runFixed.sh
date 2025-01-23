#!/bin/bash
#Batch executable to run a problem using a series of fixed values.

# Argument #1: path to Nomad install dir

rm -f cache.txt
# Fix 0-2
echo "Fix 0-2"
$1/bin/nomad param1.txt
echo "Cache size: "; wc -l cache.txt
# Fix 2-3
echo "Fix 2-3"
$1/bin/nomad param2.txt
echo "Cache size: "; wc -l cache.txt
# Fix 3-4
echo "Fix 3-4"
$1/bin/nomad param3.txt
echo "Cache size: "; wc -l cache.txt
# Fix nothing
echo "Fix nothing"
$1/bin/nomad param10.txt
echo "Cache size: "; wc -l cache.txt

#Batch executable to run a problem using a series of fixed values.

nomad_exe=$1

rm -f cache.txt
# Fix 0-2
echo "Fix 0-2"
$nomad_exe param1.txt
echo "Cache size: "; wc -l cache.txt
# Fix 2-3
echo "Fix 2-3"
$nomad_exe param2.txt
echo "Cache size: "; wc -l cache.txt
# Fix 3-4
echo "Fix 3-4"
$nomad_exe param3.txt
echo "Cache size: "; wc -l cache.txt
# Fix nothing
echo "Fix nothing"
$nomad_exe param10.txt
echo "Cache size: "; wc -l cache.txt

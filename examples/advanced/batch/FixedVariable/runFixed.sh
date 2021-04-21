#Batch executable to run a problem using a series of fixed values.

rm -f cache.txt out1 out2 out3 out10
# Fix 0-2
echo "Fix 0-2"
../../../../bin/nomad param1.txt >& out1
echo "Cache size: "; wc -l cache.txt
# Fix 2-3
echo "Fix 2-3"
../../../../bin/nomad param2.txt >& out2
echo "Cache size: "; wc -l cache.txt
# Fix 3-4
echo "Fix 3-4"
../../../../bin/nomad param3.txt >& out3
echo "Cache size: "; wc -l cache.txt
# Fix nothing
echo "Fix nothing"
../../../../bin/nomad param10.txt >& out10
echo "Cache size: "; wc -l cache.txt

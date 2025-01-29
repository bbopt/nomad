import sys
import os
import re
input_file_name=sys.argv[1]
with open(input_file_name,'r') as openfile:
    line = openfile.read().strip()
    X = re.sub(r'\s+',' ',line).strip().split()
    openfile.close()

# Standard output is grabbed by Nomad evaluator
dim = 5
f = float(X[4])
g1 = sum([(float(X[i])-1)**2 for i in range(dim)])-25
g2 = 25-sum([(float(X[i])+1)**2 for i in range(dim)])
print(f,g1,g2)

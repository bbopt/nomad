import sys
import os
import re
import numpy as np


input_file_name=sys.argv[1]
X=np.genfromtxt(input_file_name)

#print(X)

# Standard output is grabbed by Nomad evaluator
dim = len(X)
assert (dim%2==0), "Dimension must be an even number"

f = 0
for i in range(1,int(dim/2)):
    f += pow( 10*X[2*i-1] - pow(X[2*i-2],2), 2) + pow(1-X[2*i-2],2)
print(f)
exit(0)

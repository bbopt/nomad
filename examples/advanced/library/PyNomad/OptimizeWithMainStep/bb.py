#!/usr/bin/python
import sys

def bb():
    filename = sys.argv[1]
    file = open(filename, 'r')
    x = file.read().split()
    file.close()

    dim = len(x)
    x = [float(x[i]) for i in range(dim)]

    f = sum([x[i]**2 for i in range(dim)])    
    print(f)

bb()

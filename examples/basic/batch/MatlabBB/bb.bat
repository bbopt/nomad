@echo off
@copy %1 inputMatlab.txt >NUL
matlab -nojvm -batch "run('fun'); exit" 

# Loop Suggest and Observe

This project runs a loop of suggestion and observation using the executables `suggest.exe`, `bbr.exe`, and `observe.exe`.

## Description

The `loopSuggestAndObserve.sh` script executes a series of iterations where points are suggested, evaluated, and observed. The process continues until no new points are suggested or the maximum number of iterations is reached.

## Usage

To run the script, ensure that the executables `suggest.exe`, `bbr.exe`, and `observe.exe` are present in the same directory as the script. Then, execute the script with the following command:

./loopSuggestAndObserve.sh

## Generated Files

. cache.txt, cache[1-9]*.txt: Cache files used and updated in each iteration.

. x*.txt: Files containing the suggested points for each iteration.

. f*.txt: Files containing the evaluation results for each iteration.

. param[1-9]*.txt: Parameter files updated in each iteration.

## Notes

Ensure that the executables suggest.exe, bbr.exe, and observe.exe are properly compiled and accessible by using the CMake build procedure. 

The script is configured for a maximum of 10 iterations, but this number can be adjusted by modifying the condition in the until loop.



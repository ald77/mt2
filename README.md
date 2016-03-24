mt2_calc
========
mT2 calculator with ability to return test vectors

Compilation
-----------
Compile all scripts with

    ./compile.sh

Header files should have file extension ".hpp" and should be place in the "inc" directory. Other c++ code goes in the "src" directory. The file extension ".cxx" is used for executables and ".cpp" for other c++ files.

Testing
-------
A test script for validation is provided and can be run with

    ./run/test.exe

The script generates random decays of a pair produced particle and computes mt2 using both the included calculator (src/mt2.cpp) and the standard bisection calculator as a reference (src/mt2_bisect.cpp). It outputs a file "mt2.root" containing the computed mt2 values and input four-momenta for each decay.

# parallel_fac_int64
This is a simple program for calculating the factorial of large numbers, such as 100M and more. This task is not as simple as it might seem at first glance. In order to calculate the factorial of a large number, the capacity of simple built-in C++ types is not enough. Long arithmetic must be used. And if you try to solve this problem head-on, you will understand that the calculations in the resulting program will take years. In order to effectively implement long arithmetic, namely long multiplication, the program uses the Fast Fourier Transform (FFT) algorithm.

This is a pretty interesting example of a program to understand how FFT works. The program also uses parallel computing, it can use all the computing cores that are in your system, not just one.

There is a compiled executable for Windows in the bin/Release folder.
To build the program yourself use:
make -f Makefile

To run the program you need to pass it the number whose factorial you want to calculate as a parameter, for example:
parallel_fac_int64 1000000

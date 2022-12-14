MATLAB files:

The file "algo1057.m" computes the coefficient matrix [A_{i,j}] that determines the system of inequalities given by (5) in [1]. The other MATLAB files "rempwr2.m" and "trc.m" are used by the main code in "algo1057.m" to compute (2^s mod 1057), where s is a nonnegative number, and trace of a finite field element, respectively. The code in "algo1057.m" generates the file "CFs.txt" which contains the coefficient matrix.

C file:

The file "pw1057.cpp" obtains the truth tables of Boolean functions given by Table 1 in [1] and computes their nonlinearities. It generates the file "TTs.txt" which contains all the corresponding truth tables. The other files "ADK.txt", "INEQs_1057.txt", and "EC_73_15.txt" are used by the code in "pw1057.cpp" and their explanations are available within the same code.

[1] Kavut, S. Boolean Functions Generated from (1057, 31)-Interleaved Sequences. Submitted to TBV Journal of Computer Science and Engineering.

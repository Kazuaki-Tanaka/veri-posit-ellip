witten 11/11/2020 K.Tanaka

This repository provides the set of codes related with the article
"A posteriori verification of the positivity of solutions to elliptic problems"
by Kazuaki Tanaka, Taisei Asai
arXiv:2011.04510

=====================================
verify_positivity_elliptic.m
verifies the positivity of a exact solution of exlliptic problems,
while assuming an H10-error estimation given some numerical approximation and an explicit error bound.
This program requires some constants that can be calculated by the approximate solution.
The constants should be calculated before using this program.


fixedpoint_tanx.cc
calculates a rigorous enclosure of the first positive zero of x = tan(x).
The zero coincides with the first positive zero of the Bessel function with order 1.5.
The computation relies on the function "allsol" packaged in the kv library.


besselzero_order0.cc
calculates a rigorous enclosure of the first positive zero of the Bessel function with order 0.


besselzero_order1.cc
calculates a rigorous enclosure of the first positive zero of the Bessel function with order 1.
=====================================







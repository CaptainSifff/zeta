zeta
====

A repository for the implementation of the Hurwitz Zeta, PolyLog and related functions.
In the hope that it will make its way into gcc.

The series expansions in this code are based on the excellent paper by 
David C. Wood,
The computation of Polylogarithms

Other sources are the wikipedia.

ideas:
- Use Kahan summation for improved stability. The runtime will likely be dominated by function evluations.
- Use a Table of values of Riemann Zeta at Integer values.
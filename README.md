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
- Use a table of values of Riemann Zeta at Integer values.
- Generally a function for zeta at the integers.
- norm(z) vs. abs(z) abs on complex numbers involves a sqrt. whereas abs on reals is just an if().
- replace checks against zeroes as checks against sth. like: smallest FPValue * epsilon
- Better stability for large index s
- We're evaluating the Zeta function on positive integers. Try to use an efficient series for the Bernoulli Numbers for that.
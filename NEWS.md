# jacobi 3.0.0.9000

- Elliptic alpha function.

- Rogers-Ramanujan functions.

- Jacobi theta function with characteristics.

- Allows a negative nome.


# jacobi 3.0.0

- Changed the expression of the `kleinj` function in order that its factors 
avoid a possible float overflow. 

- Major changes in the implementation of the Jacobi theta functions, following 
the new Fortran implementation by Mikael Fremling.


# jacobi 2.3.1

- More unit tests.

- The `halfPeriods` function did not work for a pair of real (`numeric`) 
numbers. This has been fixed with the help of `as.complex`. 


# jacobi 2.3.0

- Lemniscate elliptic functions.

- Dixon elliptic functions.


# jacobi 2.2.0

- Some values of the Jacobi theta functions were wrong as of version 2.1.0.

- Added some unit tests.

- New function `halfPeriods`, computing the half-periods from the elliptic 
invariants.

- New function `ellipticInvariants`, computing the elliptic invariants from 
the half-periods.


# jacobi 2.1.0

- The case when the elliptic invariant `g2` is zero is now handled.

- The method computing the half-periods ratio when the elliptic invariants are 
given led a wrong sign sometimes.


# jacobi 2.0.1

- Minor fix in the C++ code.


# jacobi 2.0.0

- Weierstrass sigma function.

- Weierstrass zeta function.

- Costa surface.

- Vectorization.

- Better accuracy.


# jacobi 1.0.0

First release.

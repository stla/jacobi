## Testing environments

* local R 4.3.2, Windows 10
* Ubuntu, via Github action
* win-builder devel


## R CMD check results

OK


## Issues with Valgrind

Many errors are found by the check under the new version of Valgrind. This is 
surprising. The code is rather elementary. The same algorithms work fine with 
Julia, Fortran, and Haskell. I am doubting of the reliability of this new 
version of Valgrind.
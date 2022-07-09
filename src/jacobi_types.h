#include <Rcpp.h>

typedef std::complex<double> cplx;

cplx dologtheta3(cplx, cplx, unsigned);

bool close(cplx, cplx);

const cplx _i_;

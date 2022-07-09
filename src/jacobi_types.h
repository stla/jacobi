#include <Rcpp.h>

typedef std::complex<double> cplx;

const cplx _i_(0.0, 1.0);

cplx dologtheta3(cplx, cplx, unsigned);

bool close(cplx, cplx);

cplx jtheta1_cpp(cplx, cplx);
cplx jtheta2_cpp(cplx, cplx);
cplx jtheta3_cpp(cplx, cplx);
cplx jtheta4_cpp(cplx, cplx);

cplx ljtheta1_cpp(cplx, cplx);
cplx ljtheta2_cpp(cplx, cplx);
cplx ljtheta3_cpp(cplx, cplx);
cplx ljtheta4_cpp(cplx, cplx);

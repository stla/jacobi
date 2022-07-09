#include "jacobi_types.h"

bool isodd(int n) {
  return n % 2 == 1;
}

bool isreal(cplx z) {
  return z.imag() == 0.0;
}

const cplx _ii_(0.0, 1.0);

cplx _calctheta1_alt1(cplx z, cplx q) {
  int n = -1;
  cplx series(0.0, 0.0);
  int maxiter = 3000;
  const cplx qsq = q * q;
  cplx q_2n(1.0, 0.0);
  cplx q_n_np1(1.0, 0.0);
  while(n < maxiter) {
    n += 1;
    if(n > 0) {
      q_2n *= qsq;
      q_n_np1 *= q_2n;
    }
    double k = 2*n + 1;
    cplx term = q_n_np1 * std::sin(k*z);
    if(isodd(n)) {
      term = -term;
    }
    cplx nextseries = series + term;
    if(n >= 2 && close(nextseries, series)) {
      return 2.0 * std::sqrt(std::sqrt(q)) * series;
    } else {
      series = nextseries;
    }
  }
  Rcpp::stop("Reached 3000 iterations.");
}

cplx _calctheta1_alt2(cplx zopi, cplx topi) {
  int nminus = round(0.5 - zopi.real()) + 1;
  int nplus = nminus - 1;
  cplx series(0.0, 0.0);
  int maxterms = 3000;
  while(nplus - nminus < maxterms) {
    nplus += 1;
    nminus -= 1;
    cplx termm = std::exp(-std::pow((double)(nminus) - 0.5 + zopi, 2) / topi);
    cplx termp = std::exp(-std::pow((double)(nplus) - 0.5 + zopi, 2) / topi);
    if(isodd(nplus)) { 
      termp = -termp;
    } else {
      termm = -termm;
    }
    cplx nextseries = series + (termp + termm);
    if(nplus - nminus > 2 && close(nextseries, series)) {
      return std::sqrt(1.0/(M_PI * topi)) * series;
    }else {
      series = nextseries;
    }
  }
  Rcpp::stop("Reached 3000 terms.");
}

// [[Rcpp::export]]
cplx altjtheta1(cplx z, cplx tau) {
  cplx out;
  if(tau.imag() > 1.3) { // Chosen empirically
    // Large imag(tau) case: compute in terms of q
    cplx q = std::exp(_ii_ * M_PI * tau);
    if(isreal(q)) { 
      if(isreal(z)) {
        // Both inputs are real
        cplx outr = _calctheta1_alt1(z.real(), q.real());
        out = outr.real();
      } else {
        // q is real, but z isn't
        out = _calctheta1_alt1(z, q.real());
      }
    } else {
      // q is not real
      out = _calctheta1_alt1(z, q);
    }
  } else {
    // Small imag(tau) case: compute in terms of t/pi where t = -im * tau
    cplx topi = -_ii_ * (tau/M_PI);
    if (isreal(topi)) {
      if (isreal(z)) {
        // both z and t are real
        cplx outr = _calctheta1_alt2(z.real()/M_PI, topi.real());
        out = outr.real();
      } else {
        // t is real but z isn't
        out = _calctheta1_alt2(z/M_PI, real(topi));
      }
    } else {
      // t is not real.  No point in special casing real z here
      out = _calctheta1_alt2(z/M_PI, topi);
    }
  }
  return out;
}
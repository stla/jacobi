#include "jacobi_types.h"

bool isodd(int n) {
  return n % 2 == 1;
}

bool isreal(cplx z) {
  return z.imag() == 0.0;
}

cplx alpha(cplx z, cplx tau) {
  return std::sqrt(-_i_ * tau) * std::exp(1.0 / tau * _i_ * z * z / M_PI);
}


template <typename T1, typename T2, typename T3>
cplx _calctheta1_alt1(T1 z, T2 q) {
  int n = -1;
  T3 series = 0.0;
  int maxiter = 3000;
  const T2 qsq = q * q;
  T2 q_2n = 1;
  T2 q_n_np1 = 1;
  while(n < maxiter) {
    n += 1;
    if(n > 0) {
      q_2n *= qsq;
      q_n_np1 *= q_2n;
    }
    double k = 2*n + 1;
    T3 term = q_n_np1 * std::sin(k*z);
    if(isodd(n)) {
      term = -term;
    }
    T3 nextseries = series + term;
    if(n >= 3 && close(nextseries, series)) {
      cplx out = 2.0 * std::sqrt(std::sqrt(q)) * series;
      return out;
    } else {
      series = nextseries;
    }
  }
  Rcpp::stop("Reached 3000 iterations.");
}

template <typename T1, typename T2, typename T3>
cplx _calctheta1_alt2(T1 zopi, T2 topi) {
  int nminus = round(0.5 - std::real(zopi)) + 1;
  int nplus = nminus - 1;
  T3 series = 0;
  int maxterms = 3000;
  while(nplus - nminus < maxterms) {
    nplus += 1;
    nminus -= 1;
    T3 termm = std::exp(-std::pow((double)(nminus) - 0.5 + zopi, 2) / topi);
    T3 termp = std::exp(-std::pow((double)(nplus) - 0.5 + zopi, 2) / topi);
    if(isodd(nplus)) { 
      termp = -termp;
    } else {
      termm = -termm;
    }
    T3 nextseries = series + (termp + termm);
    if(nplus - nminus > 2 && close(nextseries, series)) {
      cplx out = std::sqrt(1.0/(M_PI * topi)) * series;
      return out;
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
    cplx q = std::exp(_i_ * M_PI * tau);
    if(isreal(q)) { 
      if(isreal(z)) {
        // Both inputs are real
        cplx outr = _calctheta1_alt1<double, double, double>(z.real(), q.real());
        out = outr.real();
      } else {
        // q is real, but z isn't
        out = _calctheta1_alt1<cplx, double, cplx>(z, q.real());
      }
    } else {
      // q is not real
      out = std::exp(ljtheta1_cpp(z/M_PI, tau)); //_calctheta1_alt2<cplx, cplx, cplx>(z/M_PI, -_i_ * (tau/M_PI));
    }
  } else {
    // Small imag(tau) case: compute in terms of t/pi where t = -im * tau
    cplx topi = -_i_ * (tau/M_PI);
    if (isreal(topi)) {
      if (isreal(z)) {
        // both z and t are real
        out = _calctheta1_alt2<double, double, double>(z.real()/M_PI, topi.real());
      } else {
        // t is real but z isn't
        out = _calctheta1_alt2<cplx, double, cplx>(z/M_PI, topi.real());
      }
    } else {
      // t is not real.  No point in special casing real z here - std::exp(-_i_ * M_PI / tau)
      out = _i_ * _calctheta1_alt1<cplx, cplx, cplx>(z / tau, std::exp(-_i_ * M_PI / tau)) / alpha(z, tau);
    }
  }
  return out;
}

// [[Rcpp::export]]
cplx altjtheta2(cplx z, cplx tau) {
  return altjtheta1(z + M_PI_2, tau);
}

inline cplx expM(cplx z, cplx tau) {
  return std::exp(_i_ * (z + tau * M_PI_4));
}

// [[Rcpp::export]]
cplx altjtheta3(cplx z, cplx tau) {
  return altjtheta2(z - M_PI_2 * tau, tau) * expM(-z, tau);
}

// [[Rcpp::export]]
cplx altjtheta4(cplx z, cplx tau) {
  return altjtheta3(z + M_PI_2, tau);
}

// [[Rcpp::export]]
cplx jtheta1_cpp(cplx z, cplx tau) {
  return altjtheta1(z * M_PI, tau);
}

// [[Rcpp::export]]
cplx jtheta2_cpp(cplx z, cplx tau) {
  return altjtheta2(z * M_PI, tau);
}

// [[Rcpp::export]]
cplx jtheta3_cpp(cplx z, cplx tau) {
  return altjtheta3(z * M_PI, tau);
}

// [[Rcpp::export]]
cplx jtheta4_cpp(cplx z, cplx tau) {
  return altjtheta4(z * M_PI, tau);
}

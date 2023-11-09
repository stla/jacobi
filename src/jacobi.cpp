#include "jacobi_types.h"

template <typename T>
T machinePrecision() {
  const T one = (T)(1.0);
  T epsilon = one;
  do {
    epsilon /= (T)(2.0);
  } while((one + epsilon / 2.0) != 1.0);
  return epsilon;
}

const double Epsilon = 2 / pow(2.0, 52);

bool close(cplx z1, cplx z2) {
  const double mod_z2 = std::abs(z2);
  const double norm = (mod_z2 < Epsilon) ? 1.0 : std::max(std::abs(z1), mod_z2);
  return std::abs(z1 - z2) < Epsilon * norm;
}

bool even(int n) {
  return n % 2 == 0;
}

cplx power(cplx z, int p) {
  if(p == 0) {
    return cplx(1.0);
  }
  if(p == 1) {
    return z;
  }
  cplx tmp = power(z, p / 2);
  if(even(p)) {
    return tmp * tmp;
  } else {
    return z * tmp * tmp;
  }
}

// [[Rcpp::export]]
cplx theta1dash(cplx z, cplx tau) {
  cplx q = std::exp(_i_ * M_PI * tau);
  cplx out(0.0, 0.0);
  double alt = -1.0;
  const cplx qsq = q * q;
  cplx q_2n = 1;
  cplx q_n_np1 = 1;
  for(int n = 0; n < 2000; n++) {
    alt = -alt;
    if(n > 0) {
      q_2n *= qsq;
      q_n_np1 *= q_2n;
    }
    double k = (double)(2 * n + 1);
    cplx outnew = out + alt * q_n_np1 * k * std::cos(k * z);
    if(close(out, outnew)) {
      return 2.0 * std::sqrt(std::sqrt(q)) * out;
    }
    out = outnew;
  }
  Rcpp::stop("Reached 2000 iterations (theta1dash).");
}

double modulo(double a, double p) {
  double i = a > 0 ? std::floor(a / p) : std::ceil(a / p);
  return a - i * p;
}

cplx calctheta3(cplx z, cplx tau) {
  cplx out(1.0, 0.0);
  unsigned n = 0;
  while(true) {
    n++;
    double nn = n;
    cplx qweight = std::exp(nn * _i_ * M_PI * (nn * tau + 2.0 * z)) +
                   std::exp(nn * _i_ * M_PI * (nn * tau - 2.0 * z));
    out += qweight;
    double modulus = std::abs(out);
    if(std::isnan(modulus)) {
      Rcpp::stop("NaN has occured in the summation.");
    } else if(std::isinf(modulus)) {
      Rcpp::stop("Infinity reached during the summation.");
    } else if(n >= 3 && close(out + qweight, out)) {
      break;
    }
  }
  return std::log(out);
}

cplx argtheta3(cplx z, cplx tau, unsigned pass_in, unsigned maxiter) {
  unsigned passes = pass_in + 1;
  if(passes > maxiter) {
    Rcpp::stop("Reached maximal iteration (argtheta3).");
  }
  double z_img = z.imag();
  double tau_img = tau.imag();
  double h = tau_img / 2.0;
  cplx zuse(modulo(z.real(), 1.0), z_img);
  cplx out;
  if(z_img < -h) {
    out = argtheta3(-zuse, tau, passes, maxiter);
  } else if(z_img >= h) {
    double quotient = std::floor(z_img / tau_img + 0.5);
    cplx zmin = zuse - quotient * tau;
    out = -2.0 * M_PI * quotient * _i_ * zmin 
          + argtheta3(zmin, tau, passes, maxiter)
          - _i_ * M_PI * tau * quotient * quotient;
  } else {
    out = calctheta3(zuse, tau);
  }
  return out;
}

cplx dologtheta4(cplx z, cplx tau, unsigned pass_in, unsigned maxiter) {
  return dologtheta3(z + 0.5, tau, pass_in + 1, maxiter);
}

cplx dologtheta3(cplx z, cplx tau, unsigned pass_in, unsigned maxiter) {
  unsigned passes = pass_in + 1;
  cplx tau2;
  double rl = tau.real();
  if(rl >= 0) {
    tau2 = modulo(rl + 1.0, 2) - 1.0 + _i_ * tau.imag();
  } else {
    tau2 = modulo(rl - 1.0, 2) + 1.0 + _i_ * tau.imag();
  }
  rl = tau2.real();
  cplx out;
  if(std::abs(tau2) < 0.98 && tau2.imag() < 0.98) { 
    cplx tauprime = -1.0 / tau2;
    out = _i_ * M_PI * tauprime * z * z
      + dologtheta3(z * tauprime, tauprime, passes, maxiter)
      - std::log(std::sqrt(tau2) / std::sqrt(_i_));
  } else if(rl > 0.6) {
    out = dologtheta4(z, tau2 - 1.0, passes, maxiter);
  } else if(rl <= -0.6) {
    out = dologtheta4(z, tau2 + 1.0, passes, maxiter);
  } else {
    out = argtheta3(z, tau2, 0, maxiter);
  }
  return out;
}


cplx principal_log_branch(cplx z) {
  double y = z.imag();
  if(y > -M_PI && y <= M_PI) {
    return z;
  }
  double twopi = 2.0 * M_PI;
  y = std::fmod(y, twopi);
  if(y > M_PI) {
    y -= twopi;
  }
  cplx out(z.real(), y);
  return out;
}

cplx M(cplx z, cplx tau) {
  return _i_ * M_PI * (z + tau / 4.0);
}

// [[Rcpp::export]]
cplx ljtheta2_cpp(cplx z, cplx tau) {
  return principal_log_branch(
    M(z, tau) + dologtheta3(z + 0.5 * tau, tau, 0, 1000)
  );
}

// [[Rcpp::export]]
cplx jtheta2_cpp(cplx z, cplx tau) {
  return std::exp(ljtheta2_cpp(z, tau));
}

// [[Rcpp::export]]
cplx ljtheta1_cpp(cplx z, cplx tau) {
  return principal_log_branch(ljtheta2_cpp(z - 0.5, tau));
}

// [[Rcpp::export]]
cplx jtheta1_cpp(cplx z, cplx tau) {
  return std::exp(ljtheta1_cpp(z, tau));
}

// [[Rcpp::export]]
cplx ljtheta3_cpp(cplx z, cplx tau) {
  return principal_log_branch(dologtheta3(z, tau, 0, 1000));
}

// [[Rcpp::export]]
cplx jtheta3_cpp(cplx z, cplx tau) {
  return std::exp(ljtheta3_cpp(z, tau));
}

// [[Rcpp::export]]
cplx ljtheta4_cpp(cplx z, cplx tau) {
  return principal_log_branch(dologtheta4(z, tau, 0, 1000));
}

// [[Rcpp::export]]
cplx jtheta4_cpp(cplx z, cplx tau) {
  return std::exp(ljtheta4_cpp(z, tau));
}

// std::complex<double> FMA(std::complex<double> a, std::complex<double> b,
// std::complex<double> c) {
//   double re = std::fma(a.real(), b.real(), c.real()) - std::fma(a.imag(),
//   b.imag(), 0.0); double im = std::fma(a.real(), b.imag(), c.imag()) +
//   std::fma(a.imag(), b.real(), 0.0); std::complex<double> z(re, im); return
//   z;
// }
//
// std::complex<double> logsin(std::complex<double> z) {
//   return _i_*z + std::log(1.0 - std::exp(-2.0*_i_*z)) - std::log(2.0*_i_);
// }

// [[Rcpp::export]]
cplx dlogjtheta1(cplx z, cplx q) {
  int nc = 0;
  const cplx qsq = q * q;
  cplx q2n(1.0, 0.0);
  const cplx Q2 = (1.0 - q) * (1.0 + q);
  cplx Q2n = Q2;
  cplx sum1 = 0;
  const cplx zero(0.0, 0.0);
  for(int n = 1; n < 10000; n++) {
    q2n *= q;
    if(q2n == zero) {
      break;
    }
    const double ndbl = (double)n;
    // const cplx f = std::log(q2n) - std::log(Q2n) + logsin(2.0 * ndbl * z);
    // const cplx term = std::log(q2n) + f;
    // const cplx newsum1 = sum1 + std::exp(term);
    const cplx f = (q2n / Q2n) * std::sin(2.0 * ndbl * z);
    const cplx term = q2n * f;
    const cplx newsum1 = sum1 + term;
    // const cplx newsum1 = FMA(q2n, f, sum1);
    if(newsum1 == sum1) {
      nc++;
      if(nc > 1) {
        break;
      }
    } else {
      nc = 0;
    }
    sum1 = newsum1;
    Q2n = qsq * Q2n + Q2;
    // Q2n = FMA(qsq, Q2n, Q2);
  }
  return 4.0 * sum1 + 1.0 / std::tan(z);
}

////////////////////////////////////////////////////////////////////////////////
Rcomplex toCplx(cplx z) {
  Rcomplex zr;
  zr.r = z.real();
  zr.i = z.imag();
  return zr;
}

cplx fromCplx(Rcomplex zr) {
  cplx z(zr.r, zr.i);
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix JTheta1(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(jtheta1_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix JTheta2(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(jtheta2_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix JTheta2_tau(Rcomplex z0, Rcpp::ComplexMatrix dalet) {
  cplx z = fromCplx(z0);
  int m = dalet.nrow();
  int n = dalet.ncol();
  Rcpp::ComplexMatrix out(m, n);
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector tauj = dalet(Rcpp::_, j);
    Rcpp::ComplexVector outj = out(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      outj(i) = toCplx(jtheta2_cpp(z, fromCplx(tauj(i))));
    }
    out(Rcpp::_, j) = outj;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix JTheta3(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(jtheta3_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix JTheta4(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(jtheta4_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix LJTheta1(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(ljtheta1_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix LJTheta2(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(ljtheta2_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix LJTheta3(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(ljtheta3_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix LJTheta4(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      zj(i) = toCplx(ljtheta4_cpp(fromCplx(zj(i)), tau));
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix dLTheta1(Rcpp::ComplexMatrix z0, Rcomplex dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(z0);
  cplx tau = fromCplx(dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    Rcpp::ComplexVector zoutj(m);
    for(int i = 0; i < m; i++) {
      cplx zij = fromCplx(zj(i));
      if(zij == 0.0) {
        cplx jp0 = std::exp(
          ljtheta2_cpp(0.0, tau) + ljtheta3_cpp(0.0, tau) + 
            ljtheta4_cpp(0.0, tau)
        );
        cplx j1 = jtheta1_cpp(0.0, tau);
        zoutj(i) = toCplx(jp0 / j1);
      } else {
        zoutj(i) = 
          toCplx(theta1dash(zij, tau) / jtheta1_cpp(zij/M_PI, tau));
      }
    }
    z(Rcpp::_, j) = zoutj;
  }
  return z;
}

////////////////////////////////////////////////////////////////////////////////
// cplx alpha(cplx z, cplx tau) {
//   return sqrt(-_i_ * tau) * std::exp(M_PI / tau * _i_ * z * z);
// }
// 
// cplx alpha0(cplx tau) {
//   return sqrt(-_i_ * tau);
// }
// 
// cplx jtheta2transfo0(cplx tau) {
//   return jtheta4_cpp(0.0, -1 / tau) / alpha0(tau);
// }
// 
// cplx jtheta3transfo0(cplx tau) {
//   return jtheta3_cpp(0.0, -1 / tau) / alpha0(tau);
// }



// [[Rcpp::export]]
Rcpp::ComplexMatrix lambda_cpp(Rcpp::ComplexMatrix Dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(Dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      if(Rcpp::ComplexVector::is_na(zj(i))) {
        zj(i) = Rcpp::ComplexVector::get_na();
      } else {
        cplx tau = fromCplx(zj(i));
        zj(i) = toCplx(power(jtheta2_cpp(0.0, tau) / jtheta3_cpp(0.0, tau), 4));
      }
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}

// [[Rcpp::export]]
Rcpp::ComplexMatrix lambda_transfo(Rcpp::ComplexMatrix Dalet) {
  Rcpp::ComplexMatrix z = Rcpp::clone(Dalet);
  int m = z.nrow();
  int n = z.ncol();
  for(int j = 0; j < n; j++) {
    Rcpp::ComplexVector zj = z(Rcpp::_, j);
    for(int i = 0; i < m; i++) {
      if(Rcpp::ComplexVector::is_na(zj(i))) {
        zj(i) = Rcpp::ComplexVector::get_na();
      } else {
        cplx tau = fromCplx(zj(i));
        cplx ratio;
        if(tau.imag() < 0.98 && std::abs(tau) < 0.98) {
          ratio = jtheta2_cpp(0.0, -1.0/ tau) / jtheta3_cpp(0.0, -1.0 / tau);
        } else {
          ratio = jtheta2_cpp(0.0, tau) / jtheta3_cpp(0.0, tau);
        }
        zj(i) = toCplx(power(ratio, 4));
      }
    }
    z(Rcpp::_, j) = zj;
  }
  return z;
}


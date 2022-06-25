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

const double Epsilon = 2.0 * machinePrecision<double>();

bool close(cplx z1, cplx z2) {
  const cplx zero(0.0, 0.0);
  const double norm = (z2 == zero) ? 1.0 : std::max(std::abs(z1), std::abs(z2));
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
cplx theta1dash(cplx z, cplx q) {
  cplx out(0.0, 0.0);
  cplx alt(-1.0, 0.0);
  for(int n = 0; n < 1000; n++) {
    alt = -alt;
    double k = (double)(2 * n + 1);
    cplx outnew = out + alt + power(q, n * (n + 1)) * k * std::cos(k * z);
    if(close(out, outnew)) {
      return 2.0 * std::pow(q, 0.25) * out;
    }
    out = outnew;
  }
  Rcpp::stop("Reached 1000 iterations (theta1dash).");
}

double modulo(double a, double p) {
  double i = a > 0 ? std::floor(a / p) : std::ceil(a / p);
  return a - i * p;
}

const cplx _i_(0.0, 1.0);

cplx calctheta3(cplx z, cplx tau) {
  cplx out(1.0, 0.0);
  unsigned n = 0;
//  bool iterate = true;
  while(n < 1000) {
    n++;
    double nn = n;
    cplx qweight = std::exp(nn * _i_ * M_PI * (nn * tau + 2.0 * z)) +
                   std::exp(nn * _i_ * M_PI * (nn * tau - 2.0 * z));
    out += qweight;
    if(std::abs(out) == 0) {
      Rcpp::stop("log(0)");
    } else if(n >= 3 && close(out + qweight, out)) {
      return std::log(out);
    }
  }
  Rcpp::stop("Reached 1000 iterations.");
}

cplx argtheta3(cplx z, cplx tau, unsigned pass_in) {
  unsigned passes = pass_in + 1;
  if(passes > 1000) {
    Rcpp::stop("passes > 1000 (argtheta3)");
  }
  double z_img = z.imag();
  double h = tau.imag() / 2.0;
  cplx zuse(modulo(z.real(), 1.0), z_img);
  cplx out;
  if(z_img < -h) {
    out = argtheta3(-zuse, tau, passes);
  } else if(z_img >= h) {
    cplx zmin = zuse - tau;
    out = -2.0 * M_PI * _i_ * zmin + argtheta3(zmin, tau, passes) -
          _i_ * M_PI * tau;
  } else {
    out = calctheta3(zuse, tau);
  }
  return out;
}

cplx dologtheta4(cplx z, cplx tau, unsigned pass_in) {
  return dologtheta3(z + 0.5, tau, pass_in + 1);
}

cplx dologtheta3(cplx z, cplx tau, unsigned pass_in) {
  unsigned passes = pass_in + 1;
  if(passes > 1000) {
    Rcpp::stop("passes > 1000 (dologtheta3)");
  }
  cplx tau2;
  double rl = tau.real();
  if(rl >= 0) {
    tau2 = modulo(rl + 1.0, 2) - 1.0 + _i_ * tau.imag();
  } else {
    tau2 = modulo(rl - 1.0, 2) + 1.0 + _i_ * tau.imag();
  }
  rl = tau2.real();
  cplx out;
  if(rl > 0.6) {
    out = dologtheta4(z, tau2 - 1.0, passes);
  } else if(rl <= -0.6) {
    out = dologtheta4(z, tau2 + 1.0, passes);
  } else if(std::abs(tau2) < 0.98 && tau2.imag() < 0.98) {
    cplx tauprime = -1.0 / tau2;
    out = _i_ * M_PI * tauprime * z * z +
          dologtheta3(-z * tauprime, tauprime, passes) -
          std::log(std::sqrt(-_i_ * tau2));
  } else {
    out = argtheta3(z, tau2, 0);
  }
  return out;
}

cplx M(cplx z, cplx tau) {
  return _i_ * M_PI * (z + tau / 4.0);
}

// [[Rcpp::export]]
cplx ljtheta2_cpp(cplx z, cplx tau) {
  return M(z, tau) + dologtheta3(z + 0.5 * tau, tau, 0);
}

// [[Rcpp::export]]
cplx jtheta2_cpp(cplx z, cplx tau) {
  return std::exp(ljtheta2_cpp(z, tau));
}

// [[Rcpp::export]]
cplx ljtheta1_cpp(cplx z, cplx tau) {
  return ljtheta2_cpp(z - 0.5, tau);
}

// [[Rcpp::export]]
cplx jtheta1_cpp(cplx z, cplx tau) {
  return std::exp(ljtheta1_cpp(z, tau));
}

// [[Rcpp::export]]
cplx ljtheta3_cpp(cplx z, cplx tau) {
  return dologtheta3(z, tau, 0);
}

// [[Rcpp::export]]
cplx jtheta3_cpp(cplx z, cplx tau) {
  return std::exp(ljtheta3_cpp(z, tau));
}

// [[Rcpp::export]]
cplx ljtheta4_cpp(cplx z, cplx tau) {
  return dologtheta4(z, tau, 0);
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

std::string rgb2hex(int r, int g, int b, bool with_head) {
  std::stringstream ss;
  if(with_head) {
    ss << "#";
  }
  ss << std::setfill('0') << std::setw(6) << std::hex << (r << 16 | g << 8 | b);
  return ss.str();
}

typedef struct Triplet_tag Triplet;
struct Triplet_tag {
  double a;
  double b;
  double c;
};

/* for RGB */
static const Triplet m[3] = {
    {3.24096994190452134377, -1.53738317757009345794, -0.49861076029300328366},
    {-0.96924363628087982613, 1.87596750150772066772, 0.04155505740717561247},
    {0.05563007969699360846, -0.20397695888897656435, 1.05697151424287856072}};

/* for XYZ */
static const Triplet m_inv[3] = {
    {0.41239079926595948129, 0.35758433938387796373, 0.18048078840183428751},
    {0.21263900587151035754, 0.71516867876775592746, 0.07219231536073371500},
    {0.01933081871559185069, 0.11919477979462598791, 0.95053215224966058086}};

static const double ref_u = 0.19783000664283680764;
static const double ref_v = 0.46831999493879100370;

static const double kappa = 903.29629629629629629630;
static const double epsilon = 0.00885645167903563082;

typedef struct Bounds_tag Bounds;
struct Bounds_tag {
  double a;
  double b;
};

static void get_bounds(double l, Bounds bounds[6]) {
  double tl = l + 16.0;
  double sub1 = (tl * tl * tl) / 1560896.0;
  double sub2 = (sub1 > epsilon ? sub1 : (l / kappa));
  int channel;
  int t;
  for(channel = 0; channel < 3; channel++) {
    double m1 = m[channel].a;
    double m2 = m[channel].b;
    double m3 = m[channel].c;
    for(t = 0; t < 2; t++) {
      double top1 = (284517.0 * m1 - 94839.0 * m3) * sub2;
      double top2 = (838422.0 * m3 + 769860.0 * m2 + 731718.0 * m1) * l * sub2 -
                    769860.0 * t * l;
      double bottom = (632260.0 * m3 - 126452.0 * m2) * sub2 + 126452.0 * t;
      bounds[channel * 2 + t].a = top1 / bottom;
      bounds[channel * 2 + t].b = top2 / bottom;
    }
  }
}

static double ray_length_until_intersect(double theta, const Bounds* line) {
  return line->b / (sin(theta) - line->a * cos(theta));
}

static double max_chroma_for_lh(double l, double h) {
  double min_len = DBL_MAX;
  double hrad = h * 0.01745329251994329577; /* (2 * pi / 360) */
  Bounds bounds[6];
  int i;
  get_bounds(l, bounds);
  for(i = 0; i < 6; i++) {
    double len = ray_length_until_intersect(hrad, &bounds[i]);
    if(len >= 0 && len < min_len) {
      min_len = len;
    }
  }
  return min_len;
}

static double dot_product(const Triplet* t1, const Triplet* t2) {
  return (t1->a * t2->a + t1->b * t2->b + t1->c * t2->c);
}

/* Used for rgb conversions */
static double from_linear(double c) {
  if(c <= 0.0031308) {
    return 12.92 * c;
  } else {
    return 1.055 * pow(c, 1.0 / 2.4) - 0.055;
  }
}

static void xyz2rgb(Triplet* in_out) {
  double r = from_linear(dot_product(&m[0], in_out));
  double g = from_linear(dot_product(&m[1], in_out));
  double b = from_linear(dot_product(&m[2], in_out));
  in_out->a = r;
  in_out->b = g;
  in_out->c = b;
}

static double l2y(double l) {
  if(l <= 8.0) {
    return l / kappa;
  } else {
    double x = (l + 16.0) / 116.0;
    return (x * x * x);
  }
}

static void luv2xyz(Triplet* in_out) {
  if(in_out->a <= 0.00000001) {
    /* Black will create a divide-by-zero error. */
    in_out->a = 0.0;
    in_out->b = 0.0;
    in_out->c = 0.0;
    return;
  }
  double var_u = in_out->b / (13.0 * in_out->a) + ref_u;
  double var_v = in_out->c / (13.0 * in_out->a) + ref_v;
  double y = l2y(in_out->a);
  double x = -(9.0 * y * var_u) / ((var_u - 4.0) * var_v - var_u * var_v);
  double z = (9.0 * y - (15.0 * var_v * y) - (var_v * x)) / (3.0 * var_v);
  in_out->a = x;
  in_out->b = y;
  in_out->c = z;
}

static void lch2luv(Triplet* in_out) {
  double hrad = in_out->c * 0.01745329251994329577; /* (pi / 180.0) */
  double u = cos(hrad) * in_out->b;
  double v = sin(hrad) * in_out->b;
  in_out->b = u;
  in_out->c = v;
}

static void hsluv2lch(Triplet* in_out) {
  double h = in_out->a;
  double s = in_out->b;
  double l = in_out->c;
  double c;
  /* White and black: disambiguate chroma */
  if(l > 99.9999999 || l < 0.00000001) {
    c = 0.0;
  } else {
    c = max_chroma_for_lh(l, h) / 100.0 * s;
  }
  /* Grays: disambiguate hue */
  if(s < 0.00000001) {
    h = 0.0;
  }
  in_out->a = l;
  in_out->b = c;
  in_out->c = h;
}

void hsluv2rgb(double h,
               double s,
               double l,
               double* pr,
               double* pg,
               double* pb) {
  Triplet tmp = {h, s, l};
  hsluv2lch(&tmp);
  lch2luv(&tmp);
  luv2xyz(&tmp);
  xyz2rgb(&tmp);
  *pr = tmp.a;
  *pg = tmp.b;
  *pb = tmp.c;
}

Rcpp::NumericVector hsluv_rgb(Rcpp::NumericVector hsl) {
  Rcpp::NumericVector out(3);
  double r, g, b;
  hsluv2rgb(hsl[0], hsl[1], hsl[2], &r, &g, &b);
  out[0] = r;
  out[1] = g;
  out[2] = b;
  return out;
}

std::string hsluv_hex(double h, double s, double l) {
  Rcpp::NumericVector hsl(3);
  Rcpp::NumericVector out(3);
  Rcpp::IntegerVector outint(3);
  hsl[0] = h;
  hsl[1] = s;
  hsl[2] = l;
  out = hsluv_rgb(hsl);
  for(int i = 0; i < 3; i++) {
    outint[i] = floor(out[i] * 255 + 0.5);
  }
  std::string outhex = (outint[0] == 0 && outint[1] == 0 && outint[2] == 0)
                           ? "#000000"
                           : rgb2hex(outint[0], outint[1], outint[2], true);
  return outhex;
}

////////////////////////////////////////////////////////////////////////////////

std::string colormap1(cplx z) {
  double x = z.real();
  double y = z.imag();
  if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    return "#000000";
  }
  double a = std::arg(z);
  double r = modulo(std::abs(z), 1.0);
  double g = abs(modulo(a, 0.5)) * 2.0;
  double b = abs(modulo(x * y, 1));
  if(std::isnan(b)) {
    return "#000000";
  }
  return rgb2hex((int)lround((1.0 - cos(r - 0.5)) * 8.0 * 255.0),
                 (int)lround((1.0 - cos(g - 0.5)) * 8.0 * 255.0),
                 (int)lround((1.0 - cos(b - 0.5)) * 8.0 * 255.0), true);
}

std::string colormap2(cplx z) {
  double x = z.real();
  double y = z.imag();
  if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) {
    return "#000000";
  }
  double arg = std::arg(z);
  if(arg < 0) {
    arg += 2.0 * M_PI;
  }
  double h = arg * 57.29577951308232087680; /* (180 / pi) */
  double w = 2 * M_PI * log1p(std::abs(z));
  double s = 100.0 * sqrt((1.0 + sin(w)) / 2.0);
  double v = 100.0 * (1.0 + cos(w)) / 2.0;
  return hsluv_hex(h, s, v);
}

void get_mobius_params(cplx gamma,
                       double t,
                       cplx* pa,
                       cplx* pb,
                       cplx* pc,
                       cplx* pd) {
  double mgamma = std::abs(gamma);
  double h = std::sqrt(1.0 - mgamma * mgamma);
  cplx z1(cos(t * M_PI / 2.0), sin(t * M_PI / 2.0));
  cplx d2 = pow(h, t) * z1;
  cplx d1 = std::conj(d2);
  cplx a(d1.real(), -d1.imag() / h);
  cplx b = d2.imag() * gamma / h;
  *pa = a;
  *pb = b;
  *pc = std::conj(b);
  *pd = std::conj(a);
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix Image_eta(Rcpp::NumericVector x, cplx gamma, double t) {
  cplx a, b, c, d;
  get_mobius_params(gamma, t, &a, &b, &c, &d);
  //
  const size_t n = x.size();
  Rcpp::CharacterMatrix Z(n, n);
  for(size_t j = 0; j < n; j++) {
    Rcpp::CharacterVector Zj(n);
    double xj = x(j);
    for(size_t i = 0; i < n; i++) {
      cplx q0(x(i), xj);
      cplx q = (a * q0 + b) / (c * q0 + d);
      if(std::abs(q) > 0.99999 || (q.imag() == 0.0 && q.real() <= 0.0)) {
        Zj(i) = "#15191e";
      } else {
        cplx tau = -_i_ * std::log(q) / M_PI;
        cplx z = std::exp(_i_ * M_PI * tau / 12.0) *
                 jtheta3_cpp((tau + 1.0) / 2.0, 3.0 * tau);
        Zj(i) = colormap1(z);
      }
    }
    Z(Rcpp::_, j) = Zj;
  }
  return Z;
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix Image_E4(Rcpp::NumericVector x, cplx gamma, double t) {
  cplx a, b, c, d;
  get_mobius_params(gamma, t, &a, &b, &c, &d);
  //
  const size_t n = x.size();
  Rcpp::CharacterMatrix Z(n, n);
  for(size_t j = 0; j < n; j++) {
    Rcpp::CharacterVector Zj(n);
    double xj = x(j);
    for(size_t i = 0; i < n; i++) {
      cplx q0(x(i), xj);
      cplx q = (a * q0 + b) / (c * q0 + d);
      if(std::abs(q) > 0.99999 || (q.imag() == 0.0 && q.real() <= 0.0)) {
        Zj(i) = "#15191e";
      } else {
        cplx tau = -_i_ * std::log(q) / M_PI / 2.0;
        cplx z = (std::pow(jtheta2_cpp(0.0, tau), 8.0) +
                  std::pow(jtheta3_cpp(0.0, tau), 8.0) +
                  std::pow(jtheta4_cpp(0.0, tau), 8.0)) /
                 2.0;
        Zj(i) = colormap1(1.0 / z);
      }
    }
    Z(Rcpp::_, j) = Zj;
  }
  return Z;
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix Image_E6(Rcpp::NumericVector x, cplx gamma, double t) {
  cplx a, b, c, d;
  get_mobius_params(gamma, t, &a, &b, &c, &d);
  //
  const size_t n = x.size();
  Rcpp::CharacterMatrix Z(n, n);
  for(size_t j = 0; j < n; j++) {
    Rcpp::CharacterVector Zj(n);
    double xj = x(j);
    for(size_t i = 0; i < n; i++) {
      cplx q0(x(i), xj);
      cplx q = (a * q0 + b) / (c * q0 + d);
      if(std::abs(q) > 0.99999 || (q.imag() == 0.0 && q.real() <= 0.0)) {
        Zj(i) = "#15191e";
      } else {
        cplx tau = -_i_ * std::log(q) / M_PI / 2.0;
        cplx j2 = jtheta2_cpp(0.0, tau);
        cplx j3 = jtheta3_cpp(0.0, tau);
        cplx j4 = jtheta4_cpp(0.0, tau);
        cplx z = (std::pow(j3, 12.0) + std::pow(j4, 12.0) -
                  3.0 * std::pow(j2, 8.0) *
                      (std::pow(j3, 4.0) + std::pow(j4, 4.0))) /
                 2.0;
        Zj(i) = colormap2(z);
      }
    }
    Z(Rcpp::_, j) = Zj;
  }
  return Z;
}

// [[Rcpp::export]]
Rcpp::CharacterMatrix Image_lambda(Rcpp::NumericVector x,
                                   cplx gamma,
                                   double t) {
  cplx a, b, c, d;
  get_mobius_params(gamma, t, &a, &b, &c, &d);
  //
  const size_t n = x.size();
  Rcpp::CharacterMatrix Z(n, n);
  for(size_t j = 0; j < n; j++) {
    Rcpp::CharacterVector Zj(n);
    double xj = x(j);
    for(size_t i = 0; i < n; i++) {
      cplx q0(x(i), xj);
      cplx q = (a * q0 + b) / (c * q0 + d);
      if(std::abs(q) > 0.99999 || (q.imag() == 0.0 && q.real() <= 0.0)) {
        Zj(i) = "#15191e";
      } else {
        cplx tau = -_i_ * std::log(q) / M_PI;
        cplx z = std::pow(jtheta2_cpp(0.0, tau) / jtheta3_cpp(0.0, tau), 4.0);
        Zj(i) = colormap1(z);
      }
    }
    Z(Rcpp::_, j) = Zj;
  }
  return Z;
}

// // [[Rcpp::export]]
// Rcpp::ComplexMatrix MOB(Rcpp::NumericVector x, cplx gamma, double t) {
//   double mgamma = std::abs(gamma);
//   double h = std::sqrt(1.0 - mgamma * mgamma);
//   cplx pz(cos(t * M_PI / 2.0), sin(t * M_PI / 2.0));
//   cplx d2 = pow(h, t) * pz;
//   cplx d1 = std::conj(d2);
//   cplx a(d1.real(), -d1.imag() / h);
//   cplx b = d2.imag() * gamma / h;
//   cplx c = std::conj(b);
//   cplx d = std::conj(a);
//   //
//   const size_t n = x.size();
//   Rcpp::ComplexMatrix Z(n, n);
//   for(size_t j = 0; j < n; j++) {
//     Rcpp::ComplexVector Zj(n);
//     double xj = x(j);
//     for(size_t i = 0; i < n; i++) {
//       cplx q0(x(i), xj);
//       cplx q = (a * q0 + b) / (c * q0 + d);
//       if(std::abs(q) > 0.99 || (q.imag() == 0.0 && q.real() <= 0.0)) {
//         Zj(i) = Rcpp::ComplexVector::get_na();
//       } else {
//         cplx tau = -_i_ * std::log(q) / M_PI;
//         cplx z = (std::pow(jtheta2_cpp(0, tau), 8.0) +
//                   std::pow(jtheta3_cpp(0, tau), 8.0) +
//                   std::pow(jtheta4_cpp(0, tau), 8.0)) /
//                  2.0;
//         Rcomplex zr;
//         zr.r = z.real();
//         zr.i = z.imag();
//         Zj(i) = zr;
//       }
//     }
//     Z(Rcpp::_, j) = Zj;
//   }
//   return Z;
// }
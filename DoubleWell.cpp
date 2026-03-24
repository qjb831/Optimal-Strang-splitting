// doublewell_rcpp.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <complex>
#include <algorithm>
#include <limits>
#include <cmath>
#include <array>
#include <vector>
#include <string>

using namespace Rcpp;

// Build roots of a polynomial using the companion matrix.
// Coefficients must be in increasing order:
// a[0] + a[1] z + ... + a[n] z^n = 0
static arma::cx_vec poly_roots_companion(const arma::vec& a, double tol = 1e-12) {
  int deg = static_cast<int>(a.n_elem) - 1;
  
  // Trim trailing near-zero leading coefficients
  while (deg > 0 && std::abs(a[deg]) < tol) {
    --deg;
  }
  
  if (deg == 0) {
    return arma::cx_vec(); // no roots
  }
  
  if (deg == 1) {
    arma::cx_vec roots(1);
    roots[0] = std::complex<double>(-a[0] / a[1], 0.0);
    return roots;
  }
  
  arma::mat C(deg, deg, arma::fill::zeros);
  
  // Subdiagonal ones
  C.submat(1, 0, deg - 1, deg - 2).eye();
  
  // First row: -a[deg-1:0] / a[deg]
  for (int j = 0; j < deg; ++j) {
    C(0, j) = -a[deg - 1 - j] / a[deg];
  }
  
  return arma::eig_gen(C);
}

static double eval_poly(const arma::vec& coeffs, double z) {
  // coeffs[0] + coeffs[1] z + ... + coeffs[n] z^n
  double val = 0.0;
  for (arma::sword i = static_cast<arma::sword>(coeffs.n_elem) - 1; i >= 0; --i) {
    val = val * z + coeffs[i];
  }
  return val;
}

static double closest_real_root_to_x(const arma::cx_vec& roots,
                                     double x,
                                     double real_tol) {
  bool found = false;
  double best_root = NA_REAL;
  double best_dist = std::numeric_limits<double>::infinity();
  
  for (arma::uword i = 0; i < roots.n_elem; ++i) {
    double re = roots[i].real();
    double im = roots[i].imag();
    
    if (std::abs(im) <= real_tol) {
      double dist = std::abs(re - x);
      if (!found || dist < best_dist) {
        found = true;
        best_dist = dist;
        best_root = re;
      }
    }
  }
  
  return found ? best_root : NA_REAL;
}

static double best_real_root_by_poly_value(const arma::cx_vec& roots,
                                           const arma::vec& coeffs,
                                           double real_tol) {
  bool found = false;
  double best_root = NA_REAL;
  double best_val = std::numeric_limits<double>::infinity();
  
  for (arma::uword i = 0; i < roots.n_elem; ++i) {
    double re = roots[i].real();
    double im = roots[i].imag();
    
    if (std::abs(im) <= real_tol) {
      double val = std::abs(eval_poly(coeffs, re));
      if (!found || val < best_val) {
        found = true;
        best_val = val;
        best_root = re;
      }
    }
  }
  
  return found ? best_root : NA_REAL;
}

// [[Rcpp::export]]
double closest_real_root_rcpp(double x, double epsilon, double y, double sigma,
                              double real_tol = 1e-4, double coeff_tol = 1e-6) {
  if (epsilon == 0.0) {
    Rcpp::stop("epsilon must be nonzero.");
  }
  
  double e2 = epsilon * epsilon;
  double e3 = e2 * epsilon;
  
  double c0 = (3.0*x*x*y)/(8.0*e3) - (x*x*x)/(3.0*e3) - y/(12.0*e3)
    - (sigma*sigma*x*x*x)/e2 - (sigma*sigma*x)/(2.0*e2)
    + (sigma*sigma*y)/(8.0*e2);
    
    double c1 = -(x*y)/(4.0*e3) + (3.0*x*x)/(8.0*e3) + 1.0/(12.0*e3)
      + (x*x*x*x)/(8.0*e3) + (3.0*sigma*sigma)/(8.0*e2);
    
    double c2 = (9.0*sigma*sigma*x)/(4.0*e2) + (3.0*x*x*x)/(2.0*e3)
      + (3.0*y)/(8.0*e3) - x/(4.0*e3) - (9.0*x*x*y)/(8.0*e3);
    
    double c3 = (3.0*x*y)/(4.0*e3) - (9.0*x*x)/(4.0*e3) - 3.0/(8.0*e3)
      - (9.0*sigma*sigma)/(8.0*e2) - (3.0*x*x*x*x)/(8.0*e3);
    
    double c4 = -(3.0*x*x*x)/(2.0*e3) + (3.0*x)/(2.0*e3) - (3.0*y)/(8.0*e3);
    
    double c5 = (27.0*x*x)/(8.0*e3) + 3.0/(8.0*e3);
    
    double c6 = -(9.0*x)/(4.0*e3);
    
    arma::vec coeffs(7);
    coeffs[0] = c0;
    coeffs[1] = c1;
    coeffs[2] = c2;
    coeffs[3] = c3;
    coeffs[4] = c4;
    coeffs[5] = c5;
    coeffs[6] = c6;
    
    arma::cx_vec roots = poly_roots_companion(coeffs, coeff_tol);
    
    if (roots.n_elem == 0) {
      Rcpp::warning("Polynomial degenerated to constant; no roots.");
      return NA_REAL;
    }
    
    // First try: actual real roots of the polynomial, closest to x
    double root = closest_real_root_to_x(roots, x, real_tol);
    if (std::isfinite(root)) {
      return root;
    }
    
    // Fallback: derivative polynomial
    arma::vec dcoeffs(6);
    dcoeffs[0] = 1.0 * coeffs[1];
    dcoeffs[1] = 2.0 * coeffs[2];
    dcoeffs[2] = 3.0 * coeffs[3];
    dcoeffs[3] = 4.0 * coeffs[4];
    dcoeffs[4] = 5.0 * coeffs[5];
    dcoeffs[5] = 6.0 * coeffs[6];
    
    arma::cx_vec droots = poly_roots_companion(dcoeffs, coeff_tol);
    
    // Choose derivative root minimizing |p(b)|
    double best_b = best_real_root_by_poly_value(droots, coeffs, real_tol);
    if (std::isfinite(best_b)) {
      return best_b;
    }
    
    Rcpp::warning("No suitable real root found for polynomial or its derivative.");
    return NA_REAL;
}

// Numerical-Recipes style cubic solver (used to match polyroot behaviour)
static std::array<std::complex<double>,3> cubic_roots_NR(std::array<double,4> a) {
  double a3 = a[0], a2 = a[1], a1 = a[2], a0 = a[3];
  std::array<std::complex<double>,3> res = { std::complex<double>(0.0,0.0),
                                             std::complex<double>(0.0,0.0),
                                             std::complex<double>(0.0,0.0) };
  
  if (std::abs(a3) < 1e-300) {
    if (std::abs(a2) < 1e-300) {
      if (std::abs(a1) < 1e-300) return res;
      res[0] = std::complex<double>(-a0/a1, 0.0);
      return res;
    }
    std::complex<double> disc = std::complex<double>(a1*a1 - 4.0*a2*a0, 0.0);
    res[0] = (-a1 + std::sqrt(disc)) / (2.0*a2);
    res[1] = (-a1 - std::sqrt(disc)) / (2.0*a2);
    return res;
  }
  
  double b = a2 / a3;
  double c = a1 / a3;
  double d = a0 / a3;
  
  double Q = (b*b - 3.0*c) / 9.0;
  double R = (2.0*b*b*b - 9.0*b*c + 27.0*d) / 54.0;
  
  if (R*R < Q*Q*Q) {
    double theta = std::acos(R / std::sqrt(Q*Q*Q));
    double sqrtQ = std::sqrt(Q);
    double x1 = -2.0*sqrtQ*std::cos(theta/3.0) - b/3.0;
    double x2 = -2.0*sqrtQ*std::cos((theta + 2.0*M_PI)/3.0) - b/3.0;
    double x3 = -2.0*sqrtQ*std::cos((theta - 2.0*M_PI)/3.0) - b/3.0;
    res[0] = std::complex<double>(x1,0.0);
    res[1] = std::complex<double>(x2,0.0);
    res[2] = std::complex<double>(x3,0.0);
    return res;
  } else {
    double signR = (R >= 0.0) ? 1.0 : -1.0;
    double inside = std::abs(R) + std::sqrt(R*R - Q*Q*Q);
    double A = -signR * std::pow(inside, 1.0/3.0);
    double B = (A == 0.0) ? 0.0 : Q / A;
    std::complex<double> A_c(A, 0.0), B_c(B, 0.0);
    res[0] = (A_c + B_c) - std::complex<double>(b/3.0, 0.0);
    std::complex<double> tmp = -0.5*(A_c + B_c) - std::complex<double>(b/3.0, 0.0);
    std::complex<double> imagpart = std::complex<double>(0.0, std::sqrt(3.0)/2.0) * (A_c - B_c);
    res[1] = tmp + imagpart;
    res[2] = tmp - imagpart;
    return res;
  }
}

static std::vector<double> cubic_realroots_sorted(double a3, double a2, double a1, double a0) {
  auto cr = cubic_roots_NR({a3,a2,a1,a0});
  std::vector<double> re(3);
  for (int i=0;i<3;++i) re[i] = std::real(cr[i]);
  std::sort(re.begin(), re.end());
  return re;
}

// A(par, center, part) -> scalar
static double A_cpp_scalar(NumericVector par, double b) {
  double eps = par(0);
  return (1.0/eps) * (-3.0 * b * b + 1.0);

}

static double N(double x, const arma::vec &par, double b){
  double eps = par(0);
  double y   = par(1);
  
  double du1 = (1.0/eps) * ( -x*x*x - y + 3.0*b*b*x + b - 3.0*b*b*b );
  return du1;
}

// single-step RK4 for fh (scalar ODE)
static double fh_scalar(double x, double h, const arma::vec &par, double b) {
  double k1 = N(x, par,b);
  double k2 = N(x + 0.5 * h * k1, par,b);
  double k3 = N(x + 0.5 * h * k2, par,b);
  double k4 = N(x + h * k3, par,b);
  double x_new = x + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
  return x_new;
}

// mu scalar
static double mu_scalar(double x, double h, NumericVector par, double b) {
  double A = A_cpp_scalar(par, b);
  double M = std::exp(A * h);
  
  return M * x + (1.0 - M) * b;
}

static double Omega_scalar(double h, NumericVector par, double b) {
  double a = A_cpp_scalar(par, b);
  double sigma12 = par(2) * par(2);
  double val = (-sigma12 / (2.0 * a)) * (1.0 - std::exp(2.0 * a * h));
  return val;
}

// should be equivalent to numDeriv in R
// [[Rcpp::export]]
double jacobian_fh_scalar_rcpp(double x, double h,
                               const arma::vec &par, double b) {
  // Base step
  const double eps = std::numeric_limits<double>::epsilon();
  const double step_base = std::pow(eps, 1.0/3.0);
  
  // Scale step with x
  double scale = std::max(1.0, std::abs(x));
  double h1 = step_base * scale;
  double h2 = h1 / 2.0;
  
  // ---- First central difference (h1)
  double f1p = fh_scalar(x + h1, h, par, b);
  double f1m = fh_scalar(x - h1, h, par, b);
  double D1 = (f1p - f1m) / (2.0 * h1);
  
  // ---- Second central difference (h2)
  double f2p = fh_scalar(x + h2, h, par, b);
  double f2m = fh_scalar(x - h2, h, par, b);
  double D2 = (f2p - f2m) / (2.0 * h2);
  
  // ---- Richardson extrapolation
  // Cancels O(h^2) error → O(h^4)
  double D = (4.0 * D2 - D1) / 3.0;
  return D;
}

// [[Rcpp::export]]
double log_lik_path_rcpp(const Rcpp::NumericMatrix &df_r,
                         const Rcpp::NumericVector &times_r,
                         const Rcpp::NumericVector &theta_drift,
                         const Rcpp::NumericVector &theta_diffusion,
                         Rcpp::RObject method) {
  // Build par vector (combine drift + diffusion)
  int nd = theta_drift.size();
  int ns = theta_diffusion.size();
  NumericVector par(nd + ns);
  for (int i=0;i<nd;i++) par(i) = theta_drift[i];
  for (int i=0;i<ns;i++) par(nd+i) = theta_diffusion[i];
  
  double eps = par(0);
  double y   = par(1);
  double sig = par(2);
  
  // df: expect a single-column matrix or vector of observations (u1)
  arma::mat data = as<arma::mat>(df_r);
  int N = data.n_rows;
  if (N < 2) Rcpp::stop("df must have at least 2 rows.");
  
  double h = times_r[1] - times_r[0];
  
  // data vectors
  arma::vec data_new = data.rows(1, N-1);
  arma::vec data_old = data.rows(0, N-2);
  int L = N - 1;
  
  std::string method_str = Rcpp::as<std::string>(method);
  arma::vec quad_terms(L), log_det_D_terms(L), log_det_Omega_terms(L);
  
  if (method_str == "fix") {
    std::vector<bool> old_cond(L);
    std::vector<bool> new_cond(L);
    
    std::vector<double> roots = cubic_realroots_sorted(-1.0, 0.0, 1.0, -y);
  
    
    double left_Omega = Omega_scalar(h, par, roots.front());
    double right_Omega = Omega_scalar(h, par, roots.back());
    double left_Inv = 1.0 / std::abs(left_Omega) ;
    double right_Inv = 1.0 / std::abs(right_Omega) ;

    for (int i=0;i<L;i++) {
      double x_new = data_new(i);
      double x_old = data_old(i);
      bool cond = (x_old > roots[1]);
      
      double use_b = cond ? roots.back() : roots.front();
      
      double tmp = fh_scalar(x_old, h/2.0, par, use_b);

      double z = fh_scalar(x_new, -h/2.0, par, use_b)-mu_scalar(tmp, h, par, use_b);

      
      double D = jacobian_fh_scalar_rcpp(x_new, -h/2.0, par, use_b);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      log_det_D_terms(i) = std::log(std::abs(D));
      
      double quad, logdetO;
      if (cond) {
        quad = z * right_Inv * z;
        logdetO = std::log(std::abs(right_Omega));
      } else {
        quad = z * left_Inv * z;
        logdetO = std::log(std::abs(left_Omega));
      }
      quad_terms(i) = quad;
      log_det_Omega_terms(i) = logdetO;
    }
    
    double loglik = arma::accu(log_det_Omega_terms) + arma::accu(quad_terms) - 2.0 * arma::accu(log_det_D_terms);
    return loglik;
  }else if (method_str == "unbiased") {
    
    for (int i=0;i<L;i++) {
      double x_new = data_new(i);
      double x_old = data_old(i);
      double use_b = closest_real_root_rcpp(x_old,eps,y,sig);
      double Omega = Omega_scalar(h, par, use_b);
      double Inv_Omega = 1.0 / std::abs(Omega) ;

      double tmp = fh_scalar(x_old, h/2.0, par, use_b);
      double z = fh_scalar(x_new, -h/2.0, par, use_b) - mu_scalar(tmp, h, par, use_b);
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h/2.0, par, use_b);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      log_det_D_terms(i) = std::log(std::abs(D));
      
      double quad, logdetO;
 
      quad = z * Inv_Omega * z;
      logdetO = std::log(std::abs(Omega));

      quad_terms(i) = quad;
      log_det_Omega_terms(i) = logdetO;
    }
    
    double loglik = arma::accu(log_det_Omega_terms) + arma::accu(quad_terms) - 2.0 * arma::accu(log_det_D_terms);
    return loglik;
  }
  
  Rcpp::stop("Unhandled center in log_lik_path_rcpp");
  return NA_REAL;
}

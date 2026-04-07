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

inline bool is_near_singularity(double z, double tol = 0) {
  const double target = 1.0 / std::sqrt(3.0);
  return (std::abs(z - target) < tol) || (std::abs(z + target) < tol);
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
      if (is_near_singularity(re)) continue;
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
double closest_real_root_rcpp(double x, NumericVector par, bool const_term_in_N = false,
                              double real_tol = 1e-12, double coeff_tol = 1e-12) {
  double epsilon = par(0);
  double y       = par(1);
  double sigma   = par(2);
  
  if (epsilon == 0.0) {
    Rcpp::stop("epsilon must be nonzero.");
  }
  
  double e2 = epsilon * epsilon;
  double e3 = e2 * epsilon;
  
  arma::vec coeffs;
  
  if (const_term_in_N) {
    coeffs.set_size(7);
    
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
      
      coeffs[0] = c0;
      coeffs[1] = c1;
      coeffs[2] = c2;
      coeffs[3] = c3;
      coeffs[4] = c4;
      coeffs[5] = c5;
      coeffs[6] = c6;
  } else {
    coeffs.set_size(8);
    
    double c0 = -sigma*sigma*x*x*x/e2 - sigma*sigma*x/(2.0*e2) - x*x*x/(3.0*e3);
    double c1 = (3.0*x*x)/(4.0*e3) + (x*x*x*x)/(8.0*e3) + sigma*sigma/(2.0*e2);
    double c2 = (3.0*x*x*x)/(2.0*e3) + (9.0*sigma*sigma*x)/(4.0*e2) - x/(2.0*e3);
    double c3 = -(15.0*x*x)/(4.0*e3) + 1.0/(12.0*e3) - (5.0*sigma*sigma)/(4.0*e2)
      - (3.0*x*x*x*x)/(8.0*e3);
    double c4 = -(3.0*x*x*x)/(2.0*e3) + (5.0*x)/(2.0*e3);
    double c5 = (9.0*x*x)/(2.0*e3) - 3.0/(8.0*e3);
    double c6 = -3.0*x/e3;
    double c7 = 3.0/(8.0*e3);
    
    coeffs[0] = c0;
    coeffs[1] = c1;
    coeffs[2] = c2;
    coeffs[3] = c3;
    coeffs[4] = c4;
    coeffs[5] = c5;
    coeffs[6] = c6;
    coeffs[7] = c7;
  }
  
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
  if (coeffs.n_elem < 2) {
    Rcpp::warning("Polynomial degree too small for derivative fallback.");
    return NA_REAL;
  }
  
  arma::vec dcoeffs(coeffs.n_elem - 1);
  for (arma::uword i = 1; i < coeffs.n_elem; ++i) {
    dcoeffs[i - 1] = static_cast<double>(i) * coeffs[i];
  }
  
  arma::cx_vec droots = poly_roots_companion(dcoeffs, coeff_tol);
  
  // Choose derivative root minimizing |p(b)|
  double best_b = best_real_root_by_poly_value(droots, coeffs, real_tol);
  if (std::isfinite(best_b)) {
    return best_b;
  }
  
  Rcpp::warning("No suitable real root found for polynomial or its derivative.");
  return NA_REAL;
}

// [[Rcpp::export]]
static std::vector<double> cubic_realroots_sorted(double a3, double a2, double a1, double a0) {
  // Build coefficient vector in R order: constant → highest degree
  Rcpp::NumericVector coeffs = Rcpp::NumericVector::create(a0, a1, a2, a3);
  
  // Call R's polyroot
  Rcpp::Function polyroot("polyroot");
  Rcpp::ComplexVector roots = polyroot(coeffs);
  
  // Extract real parts
  std::vector<double> re;
  re.reserve(roots.size());
  
  for (int i = 0; i < roots.size(); ++i) {
    Rcomplex z = roots[i];
    
    // Keep only (near-)real roots
    if (std::abs(z.i) < 1e-12) {
      re.push_back(z.r);
    }
  }
  
  // If fewer than 3 real roots, still return sorted available ones
  std::sort(re.begin(), re.end());
  
  return re;
}

// A(par, center, part) -> scalar
static double A_cpp_scalar(NumericVector par, double b) {
  double eps = par(0);
  return (1.0/eps) * (-3.0 * b * b + 1.0);

}

static double N(double x, const arma::vec &par, double b, bool const_term_in_N) {
  double eps = par(0);
  double y   = par(1);
  double du1;
  
  if (const_term_in_N) {
    du1 = (1.0 / eps) * (-x*x*x - y + 3.0*b*b*x + b - 3.0*b*b*b);
  } else {
    du1 = (1.0 / eps) * (-x*x*x + 3.0*b*b*x - 2.0*b*b*b);
  }
  
  return du1;
}

// single-step RK4 for fh (scalar ODE)
static double fh_scalar(double x, double h, const arma::vec &par, double b,bool const_term_in_N) {
  double k1 = N(x, par,b, const_term_in_N);
  double k2 = N(x + 0.5 * h * k1, par,b, const_term_in_N);
  double k3 = N(x + 0.5 * h * k2, par,b,const_term_in_N);
  double k4 = N(x + h * k3, par,b,const_term_in_N);
  double x_new = x + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
  return x_new;
}

// mu scalar
static double mu_scalar(double x, double h, NumericVector par, double b ,bool const_term_in_N) {
  double A = A_cpp_scalar(par, b);
  double M = std::exp(A * h);
  
  if (const_term_in_N){
    return M * x + (1.0 - M) * b;  
  }else{
    double Fb = -b*b*b+b-par(1);
    double b_tilde = b - Fb/A;
    return M * x + (1.0 - M) * b_tilde;
  }
}

static double Omega_scalar(double h, NumericVector par, double b) {
  double a = A_cpp_scalar(par, b);
  double sigma12 = par(2);
  double val = (-sigma12 / (2.0 * a)) * (1.0 - std::exp(2.0 * a * h));
  return val;
}

// should be equivalent to numDeriv in R
// [[Rcpp::export]]
double jacobian_fh_scalar_rcpp(double x, double h,
                               const arma::vec &par, double b, bool const_term_in_N) {
  // Base step
  const double eps = std::numeric_limits<double>::epsilon();
  const double step_base = std::pow(eps, 1.0/3.0);
  
  // Scale step with x
  double scale = std::max(1.0, std::abs(x));
  double h1 = step_base * scale;
  double h2 = h1 / 2.0;
  
  // ---- First central difference (h1)
  double f1p = fh_scalar(x + h1, h, par, b, const_term_in_N);
  double f1m = fh_scalar(x - h1, h, par, b, const_term_in_N);
  double D1 = (f1p - f1m) / (2.0 * h1);
  
  // ---- Second central difference (h2)
  double f2p = fh_scalar(x + h2, h, par, b, const_term_in_N);
  double f2m = fh_scalar(x - h2, h, par, b, const_term_in_N);
  double D2 = (f2p - f2m) / (2.0 * h2);
  
  // ---- Richardson extrapolation
  // Cancels O(h^2) error → O(h^4)
  double D = (4.0 * D2 - D1) / 3.0;
  return D;
}

// [[Rcpp::export]]
Rcpp::List log_lik_path_rcpp(const NumericMatrix &df_r,
                             const NumericVector &times_r,
                             const NumericVector &theta_drift,
                             Rcpp::RObject method) {
  
  int nd = theta_drift.size();
  if (nd < 2) stop("theta_drift must have at least 2 elements.");
  
  NumericVector par(nd + 1);
  for (int i = 0; i < nd; i++) par(i) = theta_drift[i];
  par(2) = 0.0; // sigma placeholder
  
  double y = par(1);
  
  arma::mat data = as<arma::mat>(df_r);
  int N = data.n_rows;
  if (N < 2) stop("df must have at least 2 rows.");
  if (times_r.size() < 2) stop("times_r must have at least 2 elements.");
  
  double h = times_r[1] - times_r[0];
  int L = N - 1;
  
  arma::vec data_vec = data.col(0);
  arma::vec data_old = data_vec.rows(0, N - 2);
  arma::vec data_new = data_vec.rows(1, N - 1);
  
  std::string method_str = as<std::string>(method);
  
  std::vector<double> Z(L);
  std::vector<double> B(L);
  
  double sum_sigma = 0.0;
  
  // ===================== FIX =====================
  if (method_str == "fix") {
    
    std::vector<double> roots = cubic_realroots_sorted(-1.0, 0.0, 1.0, -y);
    if (roots.size() == 1) {
      roots = std::vector<double>(3, roots[0]);
    }
    
    double root_left  = roots.front();
    double root_mid   = roots[1];
    double root_right = roots.back();
    
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      bool cond = (x_old > root_mid);
      double use_b = cond ? root_right : root_left;
      B[i] = cond ? 1.0 : 0.0;
      
      double tmp = fh_scalar(x_old, h / 2.0, par, use_b, false);
      double z = fh_scalar(x_new, -h / 2.0, par, use_b, false) - 
        mu_scalar(tmp, h, par, use_b, false);
      
      Z[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      double denom = std::exp(2 * A * h) - 1.0;
      if (std::abs(denom) < 1e-12) denom = 1e-12;
      
      sum_sigma += (2.0 * z * z * A) / denom;
    }
    
    par(2) = sum_sigma / L;
    
    double Omega_left  = Omega_scalar(h, par, root_left);
    double Omega_right = Omega_scalar(h, par, root_right);
    
    double inv_left  = 1.0 / std::abs(Omega_left);
    double inv_right = 1.0 / std::abs(Omega_right);
    
    double log_left  = std::log(std::abs(Omega_left));
    double log_right = std::log(std::abs(Omega_right));
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    
    for (int i = 0; i < L; i++) {
      double x_new = data_new(i);
      double z = Z[i];
      bool cond = (B[i] > 0.5);
      
      double use_b = cond ? root_right : root_left;
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, false);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      
      if (cond) {
        sum_quad += z * inv_right * z;
        sum_logOmega += log_right;
      } else {
        sum_quad += z * inv_left * z;
        sum_logOmega += log_left;
      }
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    
    return List::create(
      Named("ll") = ll,
      Named("sigma") = std::sqrt(par(2))
    );
  }
  
  // ===================== OLD =====================
  else if (method_str == "optimal bias old") {
    
    bool const_term_in_N = true;
    
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      double use_b = closest_real_root_rcpp(x_old, par, const_term_in_N);
      B[i] = use_b;
      
      double tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
      double z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
        mu_scalar(tmp, h, par, use_b, const_term_in_N);
      
      Z[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      double denom = std::exp(2 * A * h) - 1.0;
      if (std::abs(denom) < 1e-12) denom = 1e-12;
      
      sum_sigma += (2.0 * z * z * A) / denom;
    }
    
    par(2) = sum_sigma / L;
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    
    for (int i = 0; i < L; i++) {
      double x_new = data_new(i);
      double z = Z[i];
      double use_b = B[i];
      
      double Omega = Omega_scalar(h, par, use_b);
      double absOmega = std::abs(Omega);
      if (absOmega < 1e-12) absOmega = 1e-12;
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      sum_quad += z * (1.0 / absOmega) * z;
      sum_logOmega += std::log(absOmega);
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    
    return List::create(
      Named("ll") = ll,
      Named("sigma") = std::sqrt(par(2)),
      Named("used_bs") = B
    );
  }
  
  // ===================== NEW =====================
  else if (method_str == "optimal bias new") {
    
    bool const_term_in_N = false;
    
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      double use_b;
      if (i == 0){
        use_b = x_old;
      }else{
        use_b = data_old(i-1);
      }
      B[i] = use_b;
      
      double tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
      double z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
        mu_scalar(tmp, h, par, use_b, const_term_in_N);
      
      Z[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      double denom = std::exp(2 * A * h) - 1.0;
      if (std::abs(denom) < 1e-12) denom = 1e-12;
      
      sum_sigma += (2.0 * z * z * A) / denom;
    }
    
    par(2) = sum_sigma / L;
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    
    for (int i = 0; i < L; i++) {
      double x_new = data_new(i);
      double z = Z[i];
      double use_b = B[i];
      
      double Omega = Omega_scalar(h, par, use_b);
      double absOmega = std::abs(Omega);
      if (absOmega < 1e-12) absOmega = 1e-12;
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      sum_quad += z * (1.0 / absOmega) * z;
      sum_logOmega += std::log(absOmega);
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    
    return List::create(
      Named("ll") = ll,
      Named("sigma") = std::sqrt(par(2)),
      Named("used_bs") = B
    );
  }
  
  // ===================== FALLBACK =====================
  else {
    stop("Unknown method in log_lik_path_rcpp");
  }
}
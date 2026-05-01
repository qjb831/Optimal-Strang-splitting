// doublewell_rcpp.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;



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
// [[Rcpp::export]]
static double fh_scalar(double x, double h, const arma::vec &par, double b,bool const_term_in_N) {
  double k1 = N(x, par,b, const_term_in_N);
  double k2 = N(x + 0.5 * h * k1, par,b, const_term_in_N);
  double k3 = N(x + 0.5 * h * k2, par,b,const_term_in_N);
  double k4 = N(x + h * k3, par,b,const_term_in_N);
  double x_new = x + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
  return x_new;
}

// mu scalar
// [[Rcpp::export]]
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
// [[Rcpp::export]]
static double Omega_scalar(double h, NumericVector par, double b) {
  double a = A_cpp_scalar(par, b);
  double sigma12 = par(2)*par(2);
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
  
  double D = (4.0 * D2 - D1) / 3.0;
  return D;
}
// [[Rcpp::export]]
NumericVector sigma_MLE(NumericVector Z, NumericVector A, const NumericVector &times_r) {
  int n = Z.size();
  NumericVector result(n);
  double h = times_r[1] - times_r[0];
  
  for (int i = 0; i < n; i++) {
    result[i] = (2.0 * Z[i] * Z[i] * A[i]) /
      (std::exp(2.0 * A[i] * h) - 1.0);
  }
  
  return result;
}

// [[Rcpp::export]]
Rcpp::List log_lik_path_rcpp(const NumericMatrix &df_r,
                             const NumericVector &times_r,
                             const NumericVector &theta_drift,
                             const NumericVector &theta_diffusion,
                             Rcpp::RObject method,
                             const NumericVector &b_vec, double uns_fix) {
  
  int nd = theta_drift.size();
  if (nd < 2) stop("theta_drift must have at least 2 elements.");
  
  NumericVector par(nd + 1);
  for (int i = 0; i < nd; i++) par(i) = theta_drift[i];
  par(2) = theta_diffusion(0); 
  
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
  
  std::vector<double> Z_list(L);
  std::vector<double> b_list(L);
  std::vector<double> A_list(L);

  
  // ===================== one b =====================
  if (method_str == "avg bias" || method_str == "godambe"|| method_str == "average N") {
    
    bool const_term_in_N = false;
    double use_b =b_vec[0];
    double Omega  = Omega_scalar(h, par, use_b);
    
    double inv_Omega  = 1.0 / std::abs(Omega);
    
    double log_Omega  = std::log(std::abs(Omega));
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      double tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
      double z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
        mu_scalar(tmp, h, par, use_b, const_term_in_N);
      Z_list[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      A_list[i] = A;
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      

      sum_quad += z * inv_Omega * z;
      sum_logOmega += log_Omega;
  
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    return List::create(
      Named("ll") = ll,
      Named("A_list") = A_list,
      Named("Z_list") = Z_list
    );
  }
  // ===================== two b's =====================
  else if (method_str == "fix" || method_str == "fix penalized" || method_str == "negative fix" || method_str == "MLE splitting" || method_str == "theoretical bias (taylor)"|| method_str == "empirical bias (first-order)" || method_str == "theoretical bias (first-order)") {
    
    bool const_term_in_N;
    if (method_str == "fix" || method_str == "fix penalized" || method_str == "MLE splitting"  || method_str == "theoretical bias (taylor)"){
      const_term_in_N = false;
    }else if (method_str == "empirical bias (first-order)" || method_str == "theoretical bias (first-order)"){
      const_term_in_N = true;
    }
    double left_b  = b_vec[0];
    double right_b = b_vec[1];
    
    
    double Omega_left  = Omega_scalar(h, par, left_b);
    double Omega_right = Omega_scalar(h, par, right_b);
    
    double inv_left  = 1.0 / std::abs(Omega_left);
    double inv_right = 1.0 / std::abs(Omega_right);
    
    double log_left  = std::log(std::abs(Omega_left));
    double log_right = std::log(std::abs(Omega_right));
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;

    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      bool cond = (x_old > uns_fix);
      double use_b = cond ? right_b : left_b;

      double tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
      double z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
        mu_scalar(tmp, h, par, use_b, const_term_in_N);
      
      Z_list[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      A_list[i] = A;
      // double denom = std::exp(2 * A * h) - 1.0;
      // if (std::abs(denom) < 1e-12) denom = 1e-12;
      // 
      // sum_sigma += (2.0 * z * z * A) / denom;
   

      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N);
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
      Named("A_list") = A_list,
      Named("Z_list") = Z_list
    );
  }
  // three b's
  else if (method_str == "zero" || method_str == "uns" ) {
    

    bool const_term_in_N = false;

    double left_b  = b_vec[0];
    double middle_b = b_vec[1];
    double right_b = b_vec[2];
    
    
    double Omega_left  = Omega_scalar(h, par, left_b);
    double Omega_middle  = Omega_scalar(h, par, middle_b);
    double Omega_right = Omega_scalar(h, par, right_b);
    
    double inv_left  = 1.0 / std::abs(Omega_left);
    double inv_middle  = 1.0 / std::abs(Omega_middle);
    double inv_right = 1.0 / std::abs(Omega_right);
    
    double log_left  = std::log(std::abs(Omega_left));
    double log_middle  = std::log(std::abs(Omega_middle));
    double log_right = std::log(std::abs(Omega_right));
    
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      bool cond1 = (x_old > 1/std::sqrt(3.0)-0.25);
      bool cond2 = (x_old < -1/std::sqrt(3.0)+0.25);
      double use_b = cond1 ? right_b : (cond2 ? left_b : middle_b);
      
      double tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
      double z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
        mu_scalar(tmp, h, par, use_b, const_term_in_N);
      
      Z_list[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      A_list[i] = A;
      // double denom = std::exp(2 * A * h) - 1.0;
      // if (std::abs(denom) < 1e-12) denom = 1e-12;
      // 
      // sum_sigma += (2.0 * z * z * A) / denom;
      
      
      double D = jacobian_fh_scalar_rcpp(x_new, -h / 2.0, par, use_b, const_term_in_N);
      if (std::abs(D) < 1e-12) D = (D >= 0 ? 1e-12 : -1e-12);
      
      sum_logD += std::log(std::abs(D));
      
      if (cond1) {
        sum_quad += z * inv_right * z;
        sum_logOmega += log_right;
      } else if (cond2) {
        sum_quad += z * inv_left * z;
        sum_logOmega += log_left;
      } else {
        sum_quad += z * inv_middle * z;
        sum_logOmega += log_middle;
      }
    }
    
    double ll = sum_logOmega + sum_quad - 2.0 * sum_logD;
    
    return List::create(
      Named("ll") = ll,
      Named("A_list") = A_list,
      Named("Z_list") = Z_list
    );
  }
  // ===================== multiple b's =====================
  else if (method_str == "OU expectation x" || method_str == "OU expectation b" || method_str == "OU expectation x test" || method_str == "OU expectation b test"|| method_str == "One-step expectation"|| method_str == "One-step expectation test") {
    double sum_logD = 0.0;
    double sum_quad = 0.0;
    double sum_logOmega = 0.0;
    bool const_term_in_N = false;
    for (int i = 0; i < L; i++) {
      double x_old = data_old(i);
      double x_new = data_new(i);
      
      //double use_b = closest_real_root_rcpp(x_old, par, const_term_in_N);
      double use_b = b_vec[i];
      double tmp = fh_scalar(x_old, h / 2.0, par, use_b, const_term_in_N);
      double z = fh_scalar(x_new, -h / 2.0, par, use_b, const_term_in_N) - 
        mu_scalar(tmp, h, par, use_b, const_term_in_N);
      Z_list[i] = z;
      
      double A = A_cpp_scalar(par, use_b);
      A_list[i] = A;
      
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
      Named("A_list") = A_list,
      Named("Z_list") = Z_list
    );
  }
  else {
    stop("Unknown method in log_lik_path_rcpp");
  }
}
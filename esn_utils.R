generate_W_gaussian <- function(Nx, density, sigma = 1, ro = 1){
  #Matrix
  W <- Matrix::rsparsematrix(Nx, Nx, density, rand.x = rnorm)*sigma
  vaps <- eigen(W, only.values = TRUE)$values
  W <- ro*W/max(abs(vaps))
  if( any(is.na(as.matrix(W))) ) {
    W = matrix(0, Nx, Nx)
  }
  else{
    W = as.matrix(W)
  }
  return(W)
}


generate_Win_gaussian <- function(Nx, Nu, mu = 0, sigma = 1){
  Win <- matrix(rnorm(Nx * (1 + Nu), mu, sigma), Nx, 1 + Nu)
  return(Win)
}


Rcpp::cppFunction("arma::mat calculate_xc(arma::mat u, arma::mat W, arma::mat Win, double alpha, bool bias = true) {
            if(bias){ 
              int Nx = Win.n_rows;
              int T0 = u.n_cols;
              arma::mat x = arma::zeros(Nx, T0);
              arma::mat x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(0)));
              x.col(0) = alpha*x_bar;
              for(int i = 1; i < T0; ++i){
                x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(i)) + W*x.col(i - 1));
                x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
              }
              return(x);
            }
            else{
              int Nx = Win.n_rows;
              int T0 = u.n_cols;
              arma::mat x = arma::zeros(Nx, T0);
              arma::mat x_bar = arma::tanh(Win*u.col(0));
              x.col(0) = alpha*x_bar;
              for(int i = 1; i < T0; ++i){
                x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1));
                x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
              }
              return(x);
            }
            }", depends='RcppArmadillo')


Rcpp::cppFunction("arma::mat calculate_xc_fb(arma::mat u, arma::mat W, arma::mat Win, arma::mat Wfb, arma::mat y, double alpha, bool bias = true) {
            if(bias){ 
              int Nx = Win.n_rows;
              int T0 = u.n_cols;
              arma::mat x = arma::zeros(Nx, T0);
              arma::mat x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(0)));
              x.col(0) = alpha*x_bar;
              for(int i = 1; i < T0; ++i){
                x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(i)) + W*x.col(i - 1) + Wfb*y.col(i - 1));
                x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
              }
              return(x);
            }
            else{
              int Nx = Win.n_rows;
              int T0 = u.n_cols;
              arma::mat x = arma::zeros(Nx, T0);
              arma::mat x_bar = arma::tanh(Win*u.col(0));
              x.col(0) = alpha*x_bar;
              for(int i = 1; i < T0; ++i){
                x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*y.col(i - 1));
                x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
              }
              return(x);
            }
            }", depends='RcppArmadillo')


Rcpp::cppFunction("arma::mat calculate_x_fb(arma::mat u, arma::mat W, arma::mat Win, arma::mat Wfb, arma::mat y, double alpha, bool bias = true) {
             if(bias){ 
                  int Nx = Win.n_rows;
                  int T0 = u.n_cols;
                  arma::mat x = arma::zeros(Nx, T0);
                  arma::mat x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(0)));
                  x.col(0) = alpha*x_bar;
                  for(int i = 1; i < T0; ++i){
                  x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(i)) + W*x.col(i - 1) + Wfb*y.col(i - 1));
                  x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                  }
                  return(x);
                  }
              else{
                  int Nx = Win.n_rows;
                  int T0 = u.n_cols;
                  arma::mat x = arma::zeros(Nx, T0);
                  arma::mat x_bar = arma::tanh(Win*u.col(0));
                  x.col(0) = alpha*x_bar;
                  for(int i = 1; i < T0; ++i){
                  x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*y.col(i - 1));
                  x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                  }
                  return(x);
                  }
                  }", depends='RcppArmadillo')

Rcpp::cppFunction("arma::mat calculate_xfb_test(arma::mat u, arma::mat W, arma::mat Win, 
                  arma::mat Wfb, arma::mat Wout, double alpha, bool bias = true) {
                  int Nx = Win.n_rows;
                  int T0 = u.n_cols;
                  int Ny = Wout.n_rows;
                  arma::mat x = arma::zeros(Nx, T0);
                  arma::mat y = arma::zeros(Ny, T0);
                  arma::mat x_bar = arma::tanh(Win*u.col(0));
                  x.col(0) = alpha*x_bar;
                  y.col(0) = Wout * arma::join_cols(join_cols(arma::ones(1, 1), u.col(0)), x.col(0));
                  for(int i = 1; i < T0; ++i){
                    x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*y.col(i - 1));
                    x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                    y.col(i) = Wout * arma::join_cols(join_cols(arma::ones(1, 1), u.col(i)), x.col(i));
                  }
                  return(arma::join_cols(x, y));
                  }", depends='RcppArmadillo')            

calculate_x = function(u, W, Win, alpha){
  Nx = nrow(Win)
  T0 = ncol(u)
  x = matrix(0, Nx, T0)
  x[, 1] = tanh(Win %*% u[, 1]) 
  if(T0 > 1){
    for(i in 2:T0){
      x[, i] = tanh(Win %*% u[, i]  + W %*% x[, i - 1])
    }
  }
  return(x) 
}

train_Wout <- function(y, u, W, Win, x = NULL, alpha, beta){
  if(is.null(x)){
    x <- calculate_xc(u, W, Win, alpha)
  } 
  Ny <- nrow(y)
  Nu <- nrow(u)
  Nx <- nrow(x)
  T0 <- ncol(u)
  input <- rbind(1, u, x)
  Wout <- t(solve(input %*% t(input) + beta*diag(Nx + Nu + 1), input %*% t(y), tol = 1e-60))
  return(Wout)
}


make_esn_predictions <- function(u, W, Win, x = NULL, alpha, Wout){
  if(is.null(x)){
    x <- calculate_x(u, W, Win, alpha)
  } 
  T0 <- ncol(u)
  return(Wout %*% rbind(1, u, x)) 
}



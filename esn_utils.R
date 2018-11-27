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


Rcpp::cppFunction("arma::mat calculate_xc_fb(arma::mat u, arma::mat W, arma::mat Win, arma::mat Wfb, 
                                             arma::mat y, double alpha, double nu_noise, bool bias = true) {
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
                x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*(y.col(i - 1) + R::runif(-nu_noise, nu_noise)));
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
                  arma::mat Wfb, arma::mat y, arma::mat Wout, int tf_until, double alpha, bool bias = true) {
                  int Nx = Win.n_rows;
                  int T0 = u.n_cols;
                  int Ny = Wout.n_rows;
                  arma::mat x = arma::zeros(Nx, T0);
                  arma::mat x_bar = arma::tanh(Win*u.col(0));
                  x.col(0) = alpha*x_bar;
                  for(int i = 1; i < tf_until; ++i){
                    x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*y.col(i - 1));
                    x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                  }

                  arma::mat res_y = arma::zeros(Ny, T0 - tf_until);
                  x_bar = arma::tanh(Win*u.col(tf_until) + W*x.col(tf_until - 1) + Wfb*y.col(tf_until - 1));
                  x.col(tf_until) = (1. - alpha)*x.col(tf_until - 1) + alpha*x_bar;
                  res_y.col(0) = Wout * arma::join_cols(join_cols(arma::ones(1, 1), u.col(tf_until)), x.col(tf_until));
                  for(int i = tf_until + 1; i < T0; ++i){
                    x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*res_y.col(i - 1 - tf_until));
                    x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                    res_y.col(i - tf_until) = Wout * arma::join_cols(join_cols(arma::ones(1, 1), u.col(i)), x.col(i));
                  }
                  return(res_y);
                  }", depends='RcppArmadillo')            


Rcpp::cppFunction("arma::mat calculate_xfb_test_no_input(arma::mat u, arma::mat W, arma::mat Win, 
                  arma::mat Wfb, arma::mat y, arma::mat Wout, int tf_until, double alpha, bool bias = true) {
                  int Nx = Win.n_rows;
                  int T0 = u.n_cols;
                  int Ny = Wout.n_rows;
                  arma::mat x = arma::zeros(Nx, T0);
                  arma::mat x_bar = arma::tanh(Win*u.col(0));
                  x.col(0) = alpha*x_bar;
                  for(int i = 1; i < tf_until; ++i){
                    x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*y.col(i - 1));
                    x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                  }

                  arma::mat res_y = arma::zeros(Ny, T0 - tf_until);
                  x_bar = arma::tanh(Win*u.col(tf_until) + W*x.col(tf_until - 1) + Wfb*y.col(tf_until - 1));
                  x.col(tf_until) = (1. - alpha)*x.col(tf_until - 1) + alpha*x_bar;
                  res_y.col(0) = Wout * arma::join_cols(arma::ones(1, 1), x.col(tf_until));
                  for(int i = tf_until + 1; i < T0; ++i){
                    x_bar = arma::tanh(Win*u.col(i) + W*x.col(i - 1) + Wfb*res_y.col(i - 1 - tf_until));
                    x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                    res_y.col(i - tf_until) = Wout * arma::join_cols(arma::ones(1, 1), x.col(i));
                  }
                  return(res_y);
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

train_Wout_no_input <- function(y, W, Win, x = NULL, alpha, beta){
  if(is.null(x)){
    x <- calculate_xc(u, W, Win, alpha)
  } 
  Ny <- nrow(y)
  Nu <- nrow(u)
  Nx <- nrow(x)
  T0 <- ncol(u)
  input <- rbind(1, x)
  Wout <- t(solve(input %*% t(input) + beta*diag(Nx + 1), input %*% t(y), tol = 1e-60))
  return(Wout)
}


make_esn_predictions <- function(u, W, Win, x = NULL, alpha, Wout){
  if(is.null(x)){
    x <- calculate_x(u, W, Win, alpha)
  } 
  T0 <- ncol(u)
  return(Wout %*% rbind(1, u, x)) 
}

make_esn_predictions_no_input <- function(u, W, Win, x = NULL, alpha, Wout){
  return(Wout %*% rbind(1, x)) 
}

train_multy_enet = function(x, y, lambda, alpha_enet){
  Ny = ncol(y)
  Nx = ncol(x)
  W = matrix(0, Nx + 1, Ny)
  for(i in 1:Ny){
    fit = glmnet(x = x, y = y[, i], lambda = lambda, alpha = alpha_enet)
    W[, i] = as.matrix(coef(fit, s = lambda))
  }
  return(W)
}


train_Wout_enet = function(y, u, W, Win, x = NULL, alpha, lambda, alpha_enet) {
  if(is.null(x)){
    x <- calculate_x(u, W, Win, alpha)
  } 
  Ny <- dim(y)[1]
  Nu <- dim(u)[1]
  Nx <- dim(x)[1]
  T0 <- dim(u)[2]
  input = rbind(u, x)
  Wout = t(train_multy_enet(x = t(input), y = t(y), lambda = lambda, alpha_enet = alpha_enet))
  return(Wout)
}



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


Rcpp::cppFunction("arma::mat calculate_xc(arma::mat u, arma::mat W, arma::mat Win, double alpha) {
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
            }", depends='RcppArmadillo')


train_Wout <- function(y, u, W, Win, x = NULL, alpha, beta){
  if(is.null(x)){
    x <- calculate_xc(u, W, Win, alpha)
  } 
  Ny <- Nrow(y)
  Nu <- Nrow(u)
  Nx <- Nrow(x)
  T0 <- Ncol(u)
  input <- rbind(1, u, x)
  Wout <- t(solve(input %*% t(input) + beta*diag(Nx + Nu + 1), input %*% t(y)))
  return(Wout)
}


make_esn_predictions <- function(u, W, Win, x = NULL, alpha, Wout){
  if(is.null(x)){
    x <- calculate_x(u, W, Win, alpha)
  } 
  T0 <- Ncol(u)
  return(Wout %*% rbind(1, u, x))
}



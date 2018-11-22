#' Title
#'
#' @param N Var dimension
#' @param mu mean of gaussian distribution
#' @param sigma variance of gaussian distribution
#' @param rho spectral radius of the matrix
#'
#' @return
#' @export
#'
#' @examples
generate_PHI_normal <- function(N, mu = 0, sigma = 1, rho){
  PHI <- matrix(rnorm(N*N, mu, sigma), N, N)
  vaps <- eigen(PHI, only.values = TRUE)$values
  PHI <- rho*PHI/max(abs(vaps))
  return(PHI)
}


#' Function wich generates a PHI = (phi0, Phi1, ..., Phin) parameters for a VAR(n) model
#'
#' @param n n of the Var(n) model
#' @param N dimension of the model
#' @param rho_array vector of the radius of the different PHI matrix
#'
#' @return
#' @export
#'
#' @examples
generate_PHI_varn = function(n, N, rho_array){
  PHI = matrix(rnorm(N), N, 1)
  for(i in 1:n){
    aux_phi = generate_PHI_normal(N = N, mu = 0, sigma = 1, rho = rho_array[i])
    PHI = cbind(PHI, aux_phi)
  }
  return(PHI)
}


#' Title
#'
#' @param PHI 
#'
#' @return
#' @export
#'
#' @examples
check_stability = function(PHI){
  big_mat = PHI[, 2:ncol(PHI), drop = FALSE]
  N = nrow(PHI)
  n = ncol(big_mat) / nrow(PHI)
  for(i in 1:(n - 1)) {
    aux_mat = diag(N)
    for(j in 1:n) {
      if(j < i) {
        aux_mat = cbind(matrix(0, N, N), aux_mat)
      } else if(j > i){
        aux_mat = cbind(aux_mat, matrix(0, N, N))
      }
    }
    big_mat = rbind(big_mat, aux_mat)
  }
  vaps = eigen(big_mat, only.values = TRUE)$values
  vaps = abs(vaps)
  return(max(vaps) < 1)
}


#' Generate data using a var(n) with compact PHI
#'
#' @param T0 Total of samples
#' @param PHI Compact phi matrix
#' @param sigma variance of the gaussian distribution
#' @param initial_condition random_gaussian or ones
#'
#' @return
#' @export
#'
#' @examples
generate_varn = function(T0, PHI, sigma, initial_condition = 'random_gaussian'){
  n = (ncol(PHI) - 1)/ nrow(PHI)
  N = nrow(PHI)
  r_mat = matrix(0, N, T0) 
  if(initial_condition == 'random_gaussian') {
    aux_vec = matrix(rnorm(n*N), n*N, 1)
  } else if(initial_condition == 'ones') {
    aux_vec = matrix(1, n*N, 1)
  }
  
  for(i in 1:T0){
    r_mat[, i] = PHI %*% rbind(1, aux_vec) + MASS::mvrnorm(1, rep(0, N), diag(sigma, N, N))
    aux_vec = rbind(r_mat[, i, drop = FALSE], aux_vec[1:((n - 1)*N), 1, drop = FALSE])
  }
  return(r_mat)
}


#' Title
#'
#' @param n 
#' @param data 
#' @param regularization 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
estimate_varn_parameters = function(n, data, reg = 0.01, tol = 1e-10){
  T0 = ncol(data)
  input = rbind(1, data[, n:(T0 - 1), drop = FALSE])
  for(i in 1:(n - 1)){
    aux_input = data[, (n - i):(T0 - i - 1)]
    input = rbind(input, aux_input)
  }
  res_data = data[ , (n + 1):T0, drop = FALSE]
  inp_tinp = input %*% t(input)
  if( det(inp_tinp) < tol) {
    inp_tinp = inp_tinp + diag(reg, nrow(input), nrow(input))
  }
  return(t(solve(inp_tinp, input %*% t(res_data))))
}


#' Title
#'
#' @param data 
#' @param PHI 
#'
#' @return
#' @export
#'
#' @examples
make_varn_predictions = function(data, PHI){
  n = (ncol(PHI) - 1)/ nrow(PHI)
  T0 = ncol(data)
  input = rbind(1, data[, n:(T0 - 1), drop = FALSE])
  for(i in 1:(n - 1)){
    aux_input = data[, (n - i):(T0 - i - 1)]
    input = rbind(input, aux_input)
  }
  res <-  PHI%*%input
  return(res)
}


#' Title
#'
#' @param real_values 
#' @param pred_values 
#'
#' @return
#' @export
#'
#' @examples
MSE_error <- function(real_values, pred_values){
  err = (real_values - pred_values)^2
  return(mean(apply(err, 2, sum)))
}
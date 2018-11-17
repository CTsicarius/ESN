seed = 40
N <- 5
Ny <- 5
Nu <- 5
n = 5
PHI_rho = c(0.2, 0.2, 0.2, 0.2, 0.2)
N_test <- 10000
sigma <- 0.1
#ESN PARAMETERS
Nx <- 1000
w_ro <- 0.9
w_dens <- 0.05
alpha <- 1
c_beta = 50000
T_min = 100

set.seed(seed)
number_of_train <- c(25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000)
rho_array <- c(0, 0.25, 0.5, 0.75, 0.9)
Nx_array = c(10, 50, 100, 250, 500, 750, 1000)
it_array = Nx_array

MSE_array_var <- matrix(0, length(number_of_train), 1)
MSE_array_esn <- matrix(0, length(number_of_train), length(it_array))
PHI = generate_PHI_varn(n = n, N =N, rho_array = PHI_rho)
check_s = check_stability(PHI = PHI)
if(check_s) {
  print('Stability OK')
} else {
  print('WARNING stability')
}
global_train_data <- generate_varn(T0 = number_of_train[length(number_of_train)], PHI = PHI, sigma =  diag(sigma, N, N)) #4.5 secs
test_data <- generate_varn(T0 = N_test, PHI = PHI, sigma = diag(sigma, N, N))
time0 <- proc.time()
print('OK1')
for (it_idx in 1:length(it_array)){
  print(it_idx)
  #w_ro = it_array[it_idx]
  Nx = it_array[it_idx]
  W <- generate_W_gaussian(Nx = Nx, density = w_dens, sigma = 1, ro = w_ro)
  Win = generate_Win(Nx = Nx, Nu = nrow(global_train_data))
  train_x <- calculate_x_c(u = global_train_data[, 1:number_of_train[length(number_of_train)]], 
                           W = W, Win = Win, alpha = alpha)
  test_x <- calculate_x_c(u = test_data, W = W, Win = Win, alpha = alpha)
  for (i in 1:length(number_of_train)){
    #CODE FOR VAR ERROR
    beta = c_beta/log(number_of_train[i])
    train_data <- global_train_data[, 1:number_of_train[i]]
    PHI_est <- estimate_varn_parameters(n, train_data)
    var_pred_val <- make_varn_predictions(test_data, PHI_est)
    MSE_array_var[i,  1] <- MSE_error(test_data[, (n + 1):N_test, drop = FALSE], var_pred_val)
    
    #CODE FOR ESN ERROR
    Wout <- train_Wout(y = train_data[, 2:number_of_train[i], drop = FALSE],
                       u = train_data[, 1:(number_of_train[i] - 1), drop = FALSE],
                       W = W, Win = Win, x = train_x[, 1:(number_of_train[i] - 1)],
                       alpha = alpha, beta = beta, print_det = print_det)
    
    esn_pred_val <- make_esn_predictions(u = test_data, W = W, Win = Win, x = test_x, alpha = alpha, Wout = Wout)
    MSE_array_esn[i, it_idx] <- MSE_error(test_data[, 3:N_test, drop = FALSE], esn_pred_val[, 2:(N_test - 1), drop = FALSE])
  }
}
print(proc.time() - time0)
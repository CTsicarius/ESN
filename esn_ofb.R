seed = 2
set.seed(seed)
#LOAD DATA
u_train = readRDS(file = 'Data/u_train_house_rising_sun')
y_train = readRDS(file = 'Data/y_train_house_rising_sun')
u_test = readRDS(file = 'Data/u_test_house_rising_sun')
y_test = readRDS(file = 'Data/y_test_house_rising_sun')

#APLICAMOS FUNCION NO LINEAL AL OUTPUT
#y_train = tan(y_train)
#y_test = tan(y_test)

#PARAMETERS
Nu = nrow(u_test)
Nx = 400
Ny = nrow(y_test)
T_min = 500
T_max = 1500
tf_until = 500
N_test = 300
alpha = 1
beta = 1e-12
u_noise = 1e-4
#ESN

W = 0.4*random_101(Nx, Nx, c(0.00625, 0.9875, 0.00625))
print(max(abs(eigen(W)$values)))
Win = random_101(Nx, Nu, c(0.5, 0, 0.5))
Wfb = unif_matrix(Nx, Ny, -2, 2)
u_train = matrix(0, nrow(u_train), ncol(u_train))
x_train = calculate_xc_fb(u = u_train, W = W, Win = Win, Wfb = Wfb, y = y_train, 
                          alpha = alpha, nu_noise = u_noise, bias = FALSE)

#DISCARD
u = u_train[1:Nu, T_min:T_max, drop = FALSE]
x = x_train[1:Nx, T_min:T_max, drop = FALSE]
y = y_train[1:Ny, T_min:T_max, drop = FALSE]

#TRAIN Wout
Wout = train_Wout_no_input(y = y, W = W, x = x, alpha = alpha, beta = beta)
train_pred_val = make_esn_predictions_no_input(u = u, W = W, Win = Win, x = x,  alpha = alpha, Wout = Wout)
train_mse = MSE_error(y, train_pred_val)
print(c('train_MSE', train_mse))

#TEST MSE
u_test = matrix(0, nrow(u_test), ncol(u_test))
test_pred_val = calculate_xfb_test_no_input(u = u_test, W = W, Win = Win, Wfb = Wfb, Wout = Wout, y = y_test[, 1:tf_until, drop = FALSE], 
                                   tf_until = tf_until, alpha=alpha, bias = FALSE)
test_pred_val = test_pred_val[, 1:N_test, drop = FALSE]
dy_test = y_test[, (tf_until + 1):(N_test + tf_until), drop = FALSE]
test_mse = MSE_error(dy_test, test_pred_val)
print(c('test_MSE', test_mse))
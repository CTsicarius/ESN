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
T_min = 1
T_max = 1500
test_T_min = 30
test_T_max = 100
alpha = 1
beta = 1e-12
#ESN
W = 0.2*random_101(Nx, Nx, c(0.00625, 0.9875, 0.00625))
print(max(abs(eigen(W)$values)))
Win = random_101(Nx, Nu, c(0.5, 0, 0.5))
u_train = matrix(0, nrow(u_train), ncol(u_train))
u_train[1, 1] = 1
x_train = calculate_xc(u = u_train, W = W, Win = Win, alpha = alpha, bias = FALSE)

#DISCARD
u = u_train[1:Nu, T_min:T_max, drop = FALSE]
x = x_train[1:Nx, T_min:T_max, drop = FALSE]
y = y_train[1:Ny, T_min:T_max, drop = FALSE]

#TRAIN Wout
Wout = train_Wout(y = y, u = u, W = W, x = x, alpha = alpha, beta = beta)
train_pred_val = make_esn_predictions(u = u, W = W, Win = Win, x = x,  alpha = alpha, Wout = Wout)
train_mse = MSE_error(y, train_pred_val)
print(c('train_MSE', train_mse))

#TEST MSE
x_test = calculate_xc(u = u_test, W = W, Win = Win, alpha = alpha, bias = FALSE)
dx_test = x_test[, test_T_min:test_T_max, drop = FALSE]
du_test = u_test[, test_T_min:test_T_max, drop = FALSE]
dy_test = y_test[, test_T_min:test_T_max, drop = FALSE]
test_pred_val = make_esn_predictions(u = du_test, W = W, Win = Win, x = dx_test, alpha = alpha, Wout = Wout)
test_mse = MSE_error(dy_test, test_pred_val)
print(c('test_MSE', test_mse ))

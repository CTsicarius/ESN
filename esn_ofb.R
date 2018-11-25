seed = 2
set.seed(seed)
#LOAD DATA
u_train = readRDS(file = 'Data/u_train_sinusoidal')
y_train = readRDS(file = 'Data/y_train_sinusoidal')
u_test = readRDS(file = 'Data/u_test_sinusoidal')
y_test = readRDS(file = 'Data/y_test_sinusoidal')

#APLICAMOS FUNCION NO LINEAL AL OUTPUT
#y_train = tan(y_train)
#y_test = tan(y_test)

#PARAMETERS
Nu = nrow(u_test)
Nx = 100
Ny = nrow(y_test)
T_min = 100
T_max = 300
tf_until = 100
N_test = 1000
alpha = 1
beta = 0
#ESN
W = 0.4*random_101(Nx, Nx, c(0.025, 0.95, 0.025))
#print(max(abs(eigen(W)$values)))
Win = random_101(Nx, Nu, c(0.5, 0, 0.5))
Wfb = random_101(Nx, Ny, c(0.5, 0, 0.5))
#Wfb = matrix(0, nrow(Wfb), ncol(Wfb))
#x_train1 = calculate_xc(u = u_train, W = W, Win = Win, alpha = alpha, bias = FALSE)
#u_train = matrix(0, nrow(u_train), ncol(u_train))
x_train = calculate_xc_fb(u = u_train, W = W, Win = Win, Wfb = Wfb, y = y_train, alpha = alpha, bias = FALSE)

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
#u_test = matrix(0, nrow(u_test), ncol(u_test))
#x_test = calculate_xc_fb(u = u_test, W = W, Win = Win, Wfb = Wfb, y = y_test, alpha = alpha, bias = FALSE)
test_pred_val = calculate_xfb_test(u = u_test, W = W, Win = Win, Wfb = Wfb, Wout = Wout, y = y_test[, 1:tf_until, drop = FALSE], 
                                   tf_until = tf_until, alpha=alpha, bias = FALSE)
test_pred_val = test_pred_val[, 1:N_test, drop = FALSE]
dy_test = y_test[, (tf_until + 1):(N_test + tf_until), drop = FALSE]
test_mse = MSE_error(dy_test, test_pred_val)
print(c('test_MSE', test_mse))
seed = 42
set.seed(seed)
#LOAD DATA
u_train = readRDS(file = 'Data/u_train_sinusoidal')
y_train = readRDS(file = 'Data/y_train_sinusoidal')
u_test = readRDS(file = 'Data/u_test_sinusoidal')
y_test = readRDS(file = 'Data/y_test_sinusoidal')

#APLICAMOS FUNCION NO LINEAL AL OUTPUT
y_train = tan(y_train)
y_test = tan(y_test)

#PARAMETERS
Nu = nrow(u_test)
Nx = 100
Ny = nrow(y_test)
T_min = 100
T_max = 300
alpha = 1
beta = 1e-10
#ESN
W = 0.38*random_101(Nx, Nx, c(0.025, 0.95, 0.025))
Win = random_101(Nx, Nu + 1, c(0.5, 0, 0.5))
x = calculate_xc(u = u_train, W = W, Win = Win, alpha = alpha)

#DISCARD
u = u_train[1:Nu, T_min:T_max, drop = FALSE]
x = x[1:Nx, T_min:T_max, drop = FALSE]
y = y_train[1:Ny, T_min:T_max, drop = FALSE]

#TRAIN Wout
Wout = train_Wout(y = y, u = u, W = W, x = x, alpha = alpha, beta = beta)
pred_val = make_esn_predictions(u = u, W = W, Win = Win, x = x,  alpha = alpha, Wout = Wout)
mse = MSE_error(y, pred_val)
print(mse)

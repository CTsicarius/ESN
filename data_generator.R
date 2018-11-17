create_sinusoidal = function(P, n_samples){
  return(matrix(sin((1:n_samples)/P), 1, n_samples))
}

u_train = create_sinusoidal(5, 300)
y_train = (1/2)*(u_train^7)
u_test = create_sinusoidal(5, 100)
y_test = (1/2)*(u_test^7)
saveRDS(u_train, 'Data/u_train_sinusoidal')
saveRDS(y_train, 'Data/y_train_sinusoidal')
saveRDS(u_test, 'Data/u_test_sinusoidal')
saveRDS(y_test, 'Data/y_test_sinusoidal')
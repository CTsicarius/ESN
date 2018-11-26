create_sinusoidal = function(P, n_samples, delay = 0){
  return(matrix(sin(((1 + delay):(n_samples + delay))/P), 1, n_samples))
}

u_train = create_sinusoidal(5, 300, 0)
#y_train = (1/2)*(u_train^7)
y_train = (1/2)*(create_sinusoidal(5, 300, 0))^7
u_test = create_sinusoidal(5, 100000, 0)
#y_test = (1/2)*(u_test^7)
y_test = (1/2)*(create_sinusoidal(5, 100000, 0))^7
#y_test = matrix(10, 1, 100)
saveRDS(u_train, 'Data/u_train_sinusoidal')
saveRDS(y_train, 'Data/y_train_sinusoidal')
saveRDS(u_test, 'Data/u_test_sinusoidal')
saveRDS(y_test, 'Data/y_test_sinusoidal')

house_of_the_rising_sun = function(n_reps){
  vec = c(0, 0, 0, 2, 3, 3, 7, 5, 5, 0, 3, 3, 14, 14, 14, 14, 
   12, 7, 5, 7, 7, 7, 7, 7, 14, 14, 14, 14, 12, 12,
    7, 5, 5, 0, 3, 3, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0)
  return(rep(vec, n_reps)/28)
}

u_train = house_of_the_rising_sun(32)
y_train = (1/2)*(u_train^7)
u_test = house_of_the_rising_sun(4)
y_test = (1/2)*(u_test^7)
saveRDS(u_train, 'Data/u_train_house_rising_sun')
saveRDS(y_train, 'Data/y_train_house_rising_sun')
saveRDS(u_test, 'Data/u_test_house_rising_sun')
saveRDS(y_test, 'Data/y_test_house_rising_sun')
create_sinusoidal = function(P, n_samples, delay = 0){
  return(matrix(sin(((1 + delay):(n_samples + delay))/P), 1, n_samples))
}

u_train = create_sinusoidal(5, 300, 0)
#y_train = (1/2)*(u_train^7)
y_train = (1/2)*(create_sinusoidal(5, 300, -15))^7
u_test = create_sinusoidal(5, 1e6, 0)
#y_test = (1/2)*(u_test^7)
y_test = (1/2)*(create_sinusoidal(5, 1e6, -15))^7
#y_test = matrix(10, 1, 100)
saveRDS(u_train, 'Data/u_train_sinusoidal')
saveRDS(y_train, 'Data/y_train_sinusoidal')
saveRDS(u_test, 'Data/u_test_sinusoidal')
saveRDS(y_test, 'Data/y_test_sinusoidal')

house_of_the_rising_sun = function(n_reps){
  vec = c(0, 0, 0, 2, 3, 3, 7, 5, 5, 0, 3, 3, 14, 14, 14, 14, 
   12, 7, 5, 7, 7, 7, 7, 7, 14, 14, 14, 14, 12, 12,
    7, 5, 5, 0, 3, 3, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0)
  return(matrix(rep(vec, n_reps)/28, 1, n_reps*48))
}

u_train = house_of_the_rising_sun(32)
y_train = house_of_the_rising_sun(32)
u_test = house_of_the_rising_sun(32)
y_test = house_of_the_rising_sun(32)
saveRDS(u_train, 'Data/u_train_house_rising_sun')
saveRDS(y_train, 'Data/y_train_house_rising_sun')
saveRDS(u_test, 'Data/u_test_house_rising_sun')
saveRDS(y_test, 'Data/y_test_house_rising_sun')


multiple_atractor_input_random = function(m, T0, spike_prob){
  res = matrix(0.0, m, T0)
  for(i in 1:T0){
    for(j in 1:m){
      if(runif(1, 0, 1) < spike_prob && all(res[1:j,i] == 0)){
        res[j, i] = 0.5
      }
    }
  }
  return(res)
}

multiple_atractor_input_fixed = function(m, T0, change_step){
  res = matrix(0.0, m, T0)
  res[1, 1] = 0.5
  for(i in 1:floor(T0/change_step)){
    index = sample(1:m, 1)
    res[index, change_step*i] = 0.5
  }
  return(res)
}
  
multiple_atractor_output = function(u){
  m = nrow(u)
  T0 = ncol(u)
  res = matrix(-0.5, m, T0)
  control_vec = rep(FALSE, m)
  for(i in 1:T0){
    if(any(u[, i] == 0.5)){
      control_vec = (u[, i] == 0.5)
    }
    res[control_vec, i] = 0.5
  }
  return(res)
}  
  
#u_train = multiple_atractor_input_random(20, 4000, 0.005)
#u_test = multiple_atractor_input_random(20, 4000, 0.005)
u_train = multiple_atractor_input_fixed(20, 4000, 200)
u_test = multiple_atractor_input_fixed(20, 4000, 200)

y_train = multiple_atractor_output(u_train)
y_test = multiple_atractor_output(u_test)

saveRDS(u_train, 'Data/u_train_multiple_atractor')
saveRDS(y_train, 'Data/y_train_multiple_atractor')
saveRDS(u_test, 'Data/u_test_multiple_atractor')
saveRDS(y_test, 'Data/y_test_multiple_atractor')  
  
mackey_glass = function(T0, delta, tau, alpha = 0.2, beta = 10, gamma = 0.1){
  y = matrix(0, 1, T0)
  for(i in 1:T0){
    y[i] = y[i - 1] + delta*( (alpha*y[i - tau/delta])/() )
  }
}













  
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
  
generate_L8 = function(T0 = 50, reps = 1){
  s_up = 2*seq(1, T0)/T0 - 1
  s_down = 2*seq(-T0, -1)/T0 + 1
  y1 = sin(pi/2*s_down)
  y2 = cos(pi/2*s_up) + 1
  y2 = rev(y2)
  y3 = -sin(pi/2*s_up)
  y3 = rev(y3)
  y4 = -cos(pi/2*s_up) - 1
  y4 = rev(y4)
  x = c(y1, rev(s_up), y3, rev(s_up))
  y = c(s_down, y2, rev(s_up), y4)
  x = rep(x, reps)
  y = rep(y, reps)
  return(list('x' = x, 'y' = 0.5*y))
}

T0 = 50
reps_train = 15
reps_test = 105

train_L8 = generate_L8(T0 = T0, reps = reps_train)
data_train = matrix(c(train_L8$x, train_L8$y), 2, 4*T0*reps_train, byrow = TRUE)
u_train = data_train[, 1:(4*T0*reps_train - 1)]
y_train = data_train[, 2:(4*T0*reps_train)]

test_L8 = generate_L8(T0 = T0, reps = reps_test)
data_test = matrix(c(test_L8$x, test_L8$y), 2, 4*T0*reps_test, byrow = TRUE)
u_test = data_test[, 1:(4*T0*reps_test - 1)]
y_test = data_test[, 2:(4*T0*reps_test)]

saveRDS(u_train, 'Data/u_train_L8')
saveRDS(y_train, 'Data/y_train_L8')
saveRDS(u_test, 'Data/u_test_L8')
saveRDS(y_test, 'Data/y_test_L8')











  
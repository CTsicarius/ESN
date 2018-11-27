random_101 = function(Nrows, Ncols, prob_array){
  P1 = prob_array[1]
  P2 = prob_array[2]
  P3 = prob_array[2] + prob_array[3] 
  aux = runif(Nrows*Ncols, 0, 1)
  aux[aux < P1] = -1
  aux[aux > P1 & aux < P3] = 0
  aux[aux > P3] = 1
  return(matrix(aux, Nrows, Ncols))
}

unif_matrix = function(Nrows, Ncols, min, max){
  return(matrix(runif(Nrows*Ncols, min, max), Nrows, Ncols))
}
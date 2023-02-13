mean_list <- list()
mean_list[[1]] <- function(x) {sin(pi*x)}
mean_list[[2]] <- function(x) {cos(2*pi*x)}
mean_list[[3]] <- function(x) {-cos(2*pi*x)}
mean_list[[4]] <- function(x) {cos(pi*x)}
mean_list[[5]] <- function(x) {cos(pi*x)}

FDAdata <- generate_random_process_values_error(n=100,m=5, p=5, mean_list =mean_list,
                                     covariancef = function(x,y){exp(-abs(x-y))},
                                     distribution = 'normal', depend = F, trans = F, snr=5)




VCMdata <- VCM_process_values_error(n=100, m=5, p = 5, domain=c(0,1), mean_list = rep(list(function(x) {0*x}), 5),
                         coef_list=coef_list , covariancef =function(x,y){exp(-abs(x-y))},
                         distribution = 'normal', snr = 5,  depend = T,
                         trans = T, transfmatrix = NULL, num_basis = 1000)

# Function to generate a design matrix with polynomial basis functions of t
generate_t_design_matrix <- function(X, t, d=1) {

  # Generate a vector k with values 0, ..., d
  k <- 0:d

  # Apply the function to each value in k
  t_design <- sapply(k, function(i){
    # Calculate (X - t)^i
    (X-t)^i
  })

  # Return the design matrix
  return(t_design)
}


#' Generates the XT design matrix
#'
#' @param X a matrix with n observations and p features
#' @param t a time point
#' @param d the degree of the polynomial
#'
#' @return the XT design matrix

generate_XT_design_matrix <- function(X, t, d=1) {

  # Calculate (X - t)^0, (X - t)^1, ..., (X - t)^d
  t_design <- generate_t_design_matrix(X[,1],t,d = d)

  # Loop through the columns of X, excluding the first and last columns
  p <- 1:(ncol(X) - 2)

  # Calculate X[,1] * (X - t)^0, X[,2] * (X - t)^1, ..., X[,p] * (X - t)^d
  XT_design <-  Reduce(cbind,lapply(p, function(i){
    # Calculate X[,i+1] * (X - t)^i
    t_design * X[,i+1]
  }) )

  return(XT_design)
}

# Define the Epanechnikov kernel function
Epa_K <- function(x) {
  # Use sapply to apply the function to each element of x
  y <- sapply(x, function(xi) {
    # Check if xi is within [-c, c]
    if (abs(xi) <= 1) {
      # If xi is within [-c, c], return the Epanechnikov kernel value
      return(0.75 * (1 - xi^2))
    } else {
      # If xi is outside [-c, c], return 0
      return(0)
    }
  })
  return(y)
}


i_localp_VCM_t <- function(data, data_list , t, kernel='Epa', d=1, h=0.1){

  totaldata <- Reduce(rbind,data_list)

  if(kernel == 'Epa') {
    i_kernel <- lapply(data_list, function(data,t,h){
      Epa_K((data[,1]-t)/h)/h/nrow(data)
    },t=t, h=h)
  }

  total_kernel <- diag(Reduce(c, i_kernel))

  total_XT_design <- generate_XT_design_matrix(totaldata, t=t, d=d)

  total_denominator <- solve(t(total_XT_design) %*% total_kernel %*% total_XT_design)


  i_XT_design <- generate_XT_design_matrix(data, t=t, d=d)
  i_kernel_design <- diag(Epa_K((data[,1]-t)/h)/h/nrow(data))
  i_y <- data[,ncol(data)]

  localp_hat_d <- total_denominator %*% t(i_XT_design) %*% i_kernel_design %*% i_y

  return(localp_hat_d)
}

i_localp_VCM_tpoints <- function(data, data_list, time_points, kernel='Epa', d=1, h=0.1) {
  i_estimate <- sapply(time_points, i_localp_VCM_t, data = data,
                       data_list = data_list, kernel='Epa', d=d, h=h)
  return(i_estimate)
}

localp_VCM_i <- function(data_list, time_points, kernel = 'Epa', d=1, h=0.1) {
  localp_VCM_i_estimate <- lapply(data_list, i_localp_VCM_tpoints, data_list= data_list, time_points= time_points,
                                  kernel= 'Epa', d=d, h=h)
  return(localp_VCM_i_estimate)
}


localp_VCM_i(data_list = VCMdata, time_points = c(0.1,0.2,0.3))

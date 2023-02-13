
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

#' Estimate localpolynomial for VCM model.
#'
#' @param data_list A list of spatial data, where each element of the list is a matrix for the raw data.
#' @param time_points A vector of time points to estimate at.
#' @param kernel Character string specifying the type of kernel to use for the estimation. The default is "Epa".
#' @param d Numeric value specifying the degree for local. The default is 1.
#' @param h Numeric value specifying the bandwidth for the kernel. The default is 0.1.
#' @return A list of estimates for individual curves.
#' @export
localp_VCM_i <- function(data_list, time_points, kernel = 'Epa', d=1, h=0.1) {

  # Apply the i_localp_VCM_tpoints function to each element of data_list
  localp_VCM_i_estimate <- lapply(data_list, i_localp_VCM_tpoints, data_list= data_list, time_points= time_points,
                                  kernel= kernel, d=d, h=h)

  # Return the list of localp estimates
  return(localp_VCM_i_estimate)
}



localp_VCM_t <- function(data_list, t, kernel='Epa', d=1, h=0.1){

  totaldata <- Reduce(rbind,data_list)

  if(kernel == 'Epa') {
    i_kernel <- lapply(data_list, function(data,t,h){
    Epa_K((data[,1]-t)/h)/h/nrow(data)
    },t=t, h=h)
  }

  total_kernel <- diag(Reduce(c, i_kernel))

  total_XT_design <- generate_XT_design_matrix(totaldata, t=t, d=d)

  total_denominator <- solve(t(total_XT_design) %*% total_kernel %*% total_XT_design)

  total_y <- totaldata[,ncol(totaldata)]

  localp_hat_d <- total_denominator %*% t(total_XT_design) %*% total_kernel %*% total_y

  return(localp_hat_d)
}


#' Estimate local partial variogram at specified time points
#' @param data_list A list of spatial data, where each element of the list is a matrix for the raw data.
#' @param time_points A vector of time points to estimate at.
#' @param kernel Character string specifying the type of kernel to use for the estimation. The default is "Epa".
#' @param d Numeric value specifying the degree for local. The default is 1.
#' @param h Numeric value specifying the bandwidth for the kernel. The default is 0.1.
#' @return A vector of local partial variogram estimates at the specified time points.
#' @export
localp_VCM_tpoints <- function(data_list, time_points, kernel='Epa', d=1, h=0.1) {

  # Apply the localp_VCM_t function to each time point in time_points
  estimate <- sapply(time_points, localp_VCM_t,
                     data_list = data_list, kernel=kernel, d=d, h=h)

  # Return the vector of local partial variogram estimates
  return(estimate)
}


centerdata_VCM <- function(data_list, kernel = 'Epa', d= 1, h = 0.1 ) {

  mis <- sapply(data_list, function(x){nrow(x)})

  totaldata <- Reduce(rbind,data_list)

  est_values <- t(localp_VCM_tpoints(data_list = data_list, time_points = totaldata[,1], kernel= kernel, d=d, h=d))

  est_0_order <- est_values[, seq(from=1, by=d+1, length.out= ncol(totaldata)-2)]

  est_y <- totaldata[,2:(ncol(totaldata)-1)] * est_0_order

  est_epsilon <- totaldata[,ncol(totaldata)] - rowSums(est_y)

  centerdata <- totaldata

  centerdata[,ncol(centerdata)]  <- est_epsilon

  centerdata_list <- split_matrix_by_rows(centerdata, mis)

  return(centerdata_list)
}


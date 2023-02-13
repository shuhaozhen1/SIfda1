#' Bootstrap Algorithm for Inference in Varying-Coefficient Models
#'
#' The `gaussmbs` function implements a bootstrap algorithm for inference in varying-coefficient models using local polynomial estimation.
#' It returns the estimated critical value and the estimated standard deviation of the coefficients.
#'
#' @param data_list a list of data frames, where each data frame contains the data for one subject.
#' @param time_points a vector of time points at which the coefficients will be estimated.
#' @param bs_times the number of bootstrap replicates to perform (default is 1000).
#' @param kernel the type of kernel to use for local polynomial estimation (default is "Epa").
#' @param d the order of derivative to estimate (default is 1).
#' @param h the bandwidth for local polynomial estimation (default is 0.1).
#'
#' @return a list containing the following objects:
#' \itemize{
#'   \item hat_q: the estimated critical value.
#'   \item est_sd: the estimated standard deviation of the coefficients.
#' }
#'
#' @examples
#' \dontrun{
#'   data_list <- list(data1, data2, data3) # where data1, data2, and data3 are data frames
#'   time_points <- c(0.5, 1.0, 1.5)
#'   result <- gaussmbs(data_list, time_points)
#'   hat_q <- result$hat_q
#'   est_sd <- result$est_sd
#' }
#'
#' @export

gaussmbs <- function(data_list, time_points, bs_times=1000, kernel='Epa', d=1, h = 0.1) {
  centerdata_list <- centerdata_VCM(data_list = data_list, kernel = kernel, d= d, h = h)

  n <- length(data_list)

  i_est <- localp_VCM_i(data_list = centerdata_list, time_points = time_points, kernel = kernel, d=d, h=h)

  p <- ncol(data_list[[1]])-2

  i_est_0_order <- lapply(i_est, function(x) { x <- c(x[ seq(from=1, by=d+1, length.out= p), ] ) } )

  i_est_mat <- n* t(Reduce(cbind,i_est_0_order))

  est_var <- colSums( i_est_mat^2 ) / n

  est_sd <- sqrt(est_var)

  sd_i_est <- t(t(i_est_mat) / est_sd )

  gaussmp <- function(n,x){ max( abs (colSums(rnorm(n) * x) / sqrt(n) ) ) }

  bs_gaussmp <- replicate(bs_times, gaussmp(n=n, x= sd_i_est))

  hat_q <- quantile(bs_gaussmp, 0.95)

  return(list(hat_q =hat_q, est_sd =matrix(est_sd, nrow = p)))
}




confidence_band <- function(data_list, time_points, bs_times=1000, kernel='Epa', d=1, h=NULL ){
  totalt <- Reduce(rbind, data_list)[,1]


  if(is.null(h)){
    h <- 1.2*bw.bcv(totalt)
  } else {
    h <- h
  }

  n <- length(data_list)
  p <- ncol(data_list[[1]])-2

  beta_est <- localp_VCM_tpoints(data_list = data_list, time_points = time_points, kernel = kernel,
                                 d=d,h=h)[seq(from=1, by=d+1, length.out= p) ,]

  bs <- gaussmbs(data_list = data_list, time_points = time_points, bs_times=bs_times, kernel=kernel, d=d, h =0.8*h)

  scb_low <- beta_est - bs$hat_q * bs$est_sd / sqrt(n)

  scb_up <- beta_est + bs$hat_q * bs$est_sd / sqrt(n)

  return(list(scb_low=scb_low, scb_up =scb_up))
}

confidence_band(VCMdata, time_points = c(0.1,0.2,0.3,0.5))


#' Generate Fourier Basis Functions
#'
#' This function generates the Fourier basis functions for a given set of input values.
#'
#' @param x A vector of input values for which the Fourier basis functions are to be evaluated.
#' @param frequency The fundamental frequency of the Fourier basis functions. This parameter determines the periodicity of the functions.
#' @param num_harmonics The number of harmonics to include in the Fourier basis. This determines the complexity of the representation, with more harmonics resulting in a finer representation of the periodic function.
#'
#' @return A matrix with `length(x)` rows and `2 * num_harmonics` columns, where each column represents a cosine or sine term of the Fourier series.
#'
#' @examples
#' x <- seq(0, 2 * pi, by = 0.1)
#' fourier_basis <- generate_fourier_basis(x, 1, 10)
generate_fourier_basis <- function(x, frequency=1, num_harmonics=100) {
  basis <- matrix(ncol = 2 * num_harmonics, nrow = length(x))
  for (i in 1:num_harmonics) {
    basis[, 2 * i - 1] <- cos(2 * pi * frequency * x * i)
    basis[, 2 * i] <- sin(2 * pi * frequency * x * i)
  }
  basis
}

#' Function to generate values of a random process at specified time points
#'
#' @param time_points: a vector of time points at which to evaluate the random process
#' @param p: number of columns in the output, representing different realizations of the random process,
#'  to be interpreted as dimensions. Default is 5.
#' @param meanf: a scalar value representing the mean of the random process
#' @param covariancef: a covariance function that maps two time points to their covariance
#' @param num_basis: the number of basis functions to use for approximation. Default is 1000.
#' @param distribution: the distribution to use for generating weights for the basis functions. Default is "normal".
#'   Accepted values are "normal", "uniform", and "exponential".
#
#' @return A matrix with `length(time_points)` rows and `p+1` columns, the first column records the time points and
#' each other column represents a realization of the random process at the specified time points.
#' @export
generate_random_process_values <- function(time_points, p=5, meanf=function(x){0},
                                           covariancef, num_basis=1000, distribution = "normal") {

  # Calculate the covariance matrix between each pair of time points
  covariance_matrix <- outer(time_points, time_points, covariancef)

  # Decompose the covariance matrix into eigenvectors and eigenvalues
  eigen <- eigen(covariance_matrix)

  # Determine the number of basis functions to use for approximation
  # If the requested number of basis functions is greater than the number of eigenvectors, use all eigenvectors instead
  num_basis <- min(num_basis, ncol(eigen$vectors))

  # Select the first `num_basis` eigenvectors as the basis functions
  basis_functions <- eigen$vectors[, 1:num_basis]
  eigen_values <- eigen$values[1:num_basis]

  # Generate weights for the basis functions
  # The distribution for generating the weights can be specified using the `distribution` argument
  if (distribution == "normal") {
    weights <- matrix(rnorm(num_basis * p), ncol = num_basis)
  } else if (distribution == "uniform") {
    weights <- matrix(runif(num_basis * p), ncol = num_basis)
  } else if (distribution == "exponential") {
    weights <- matrix(rexp(num_basis * p)-1, ncol = num_basis)
  } else {
    stop("Invalid distribution argument. Choose 'normal', 'uniform', or 'exponential'.")
  }

  # Scale the weights by the square root of the eigenvalues
  weights <- weights %*% sqrt(diag(eigen_values))

  # Calculate the values of the random process at each time point
  # by taking the inner product between the weights and the basis functions
  process_values_matrix <- weights %*% t(basis_functions)
  process_values_matrix <- t(process_values_matrix) + meanf(time_points)

  # Return the matrix of random process values
  return(cbind(time_points, process_values_matrix))
}

#' Function to generate time points for n elements
#'
#' @param n: an integer value representing the number of elements in the list
#' @param m: an integer value that is used to determine the number of time points in each element,
#'   sampled from (m-1, m, m+1) independently.
#' @param domain: the domain for time points. Default is c(0,1).
#'
#' @return A list of length n, where each element contains time points uniformly distributed of length m_i.
#' The time points in each element are sorted in ascending order.

generate_time_points <- function(n, m, domain = c(0,1)) {

  # Generate a list to store the time points for each element
  time_points_list <- list()

  for (i in 1:n) {
    # Sample the number of time points for the i-th element from (m-1, m, m+1)
    m_i <- sample(c(m-1, m, m+1), 1)

    # Generate m_i time points uniformly distributed from 0 to 1
    time_points_i <- sort(runif(m_i, min = domain[1], max = domain[2]))

    # Add the time points for the i-th element to the list
    time_points_list[[i]] <- time_points_i
  }

  # Return the list of time points
  return(time_points_list)
}

#' Function to calculate the L2 moment of a random process
#'
#' @param meanf: a function representing the mean of the random process
#' @param covariancef: a covariance function representing the covariance of the random process
#' @param domain: the interval over which the L2 moment is calculated. Default is c(0,1).
#'
#' @return The L2 moment of the random process over the specified domain.

l2moment_of_random_process <- function(meanf = function(x){0}, covariancef, domain = c(0,1)) {
  integrate(function(x) {covariancef(x, x) + meanf(x)^2}, domain[1], domain[2])$value
}


#' Generate Random Process Values with Error
#'
#' The function generates random process values with error for a given number of
#' time points and functions. The user can specify the mean function, covariance
#' function, and distribution for the random process values. The error can be
#' added as independent or dependent and can be transformed using a transformation
#' matrix if desired.
#'
#' @param n The number of time points to generate.
#' @param m The number of random process values to generate.
#' @param p The number of functions to evaluate at each time point.
#' @param domain The range of time points to generate.
#' @param mean_list A list of length p, each element is a mean function.
#' @param covariancef A covariance function to generate random process values.
#' @param distribution The distribution to use for generating random process values.
#' @param snr Signal-to-noise ratio for generating the error variance.
#' @param sig Variance of the random error.
#' @param depend Logical indicating if the error is dependent or independent.
#' @param trans Logical indicating if the process values should be transformed.
#' @param transfmatrix A transformation matrix to use for transforming process values.
#' @param num_basis The number of basis functions to use for generating random process values.
#'
#' @return A list of random process values with error.
#'
#' @examples
#' \dontrun{
#' mean_list <- rep(list(function(x) {0}), 5)
#' covariancef <- function(t1, t2) {exp(-abs(t1 - t2))}
#' values_error <- generate_random_process_values_error(n = 100, m = 1, mean_list = mean_list, covariancef = covariancef)
#' }
#' @export

generate_random_process_values_error <- function(n, m, p = 5, domain=c(0,1), mean_list, covariancef,
                                                 distribution = 'normal', snr = 5, sig = NULL, depend = FALSE,
                                                 trans = FALSE, transfmatrix = NULL, num_basis = 1000){
  # Determine the variance of random error
  if( is.null(sig) ) {
    var_error <- l2moment_of_random_process(meanf = function(x){0},
                                            covariancef = covariancef, domain = domain)/ snr
  } else {
    var_error <- sig^2
  }

  # Generate time points lists
  time_points_list <- generate_time_points(n, m, domain = domain)

  # Apply the data generation function to each list
  process_values_list <- lapply(time_points_list, generate_random_process_values, p = p, meanf = function(x){0},
                                covariancef = covariancef, num_basis = num_basis, distribution = distribution)

  # Transform the independent p-dimensional process to a dependent one by matrix
  if (trans == T) {
    if (is.null(transfmatrix)) {
      transfmatrix <- toeplitz(seq(from = 1, by = -1/p, length.out = p))
    }
    trans_values <- lapply(process_values_list, function(x){
      x[,2:(p+1)] <- x[,2:(p+1)] %*% transfmatrix
      return(x)
    })
  } else {
    trans_values <- process_values_list
  }

  # Add random error
  if (depend == T) {
    values_with_error <- lapply(trans_values, function(x){
      x[,2:(p+1)] <- x[,2:(p+1)] + x[,2] * rnorm(1, sd =sqrt(var_error))
      return(x)
    })
  } else {
    values_with_error <- lapply(trans_values, function(x){
      x[,2:(p+1)] <- x[,2:(p+1)] + rnorm(nrow(x), sd =sqrt(var_error))
      return(x)
    })
  }

  # Add mean function
  evaluated_mean <- lapply(time_points_list, function(t, mean_list) {
    value <- sapply(mean_list, function(f){f(t)})
    return(value)
  }, mean_list = mean_list)

  # Function for element-wise sum
  element_wise_sum_of_matrices <- function(list1, list2) {
    result_list <- list()
    for (i in 1:length(list1)) {
      result_list[[i]] <- list1[[i]] + cbind(rep(0,nrow(list2[[i]])),list2[[i]])
    }
    return(result_list)
  }

  values_mean_error <- element_wise_sum_of_matrices(values_with_error,evaluated_mean)

  return(values_with_error)
}


#' Generate Random Process Values with Error For VCM
#'
#' The function generates random process values with error for a given number of
#' time points and functions. The user can specify the mean function, covariance
#' function, and distribution for the random process values. The error can be
#' added as independent or dependent and can be transformed using a transformation
#' matrix if desired.
#'
#' @param n The number of time points to generate.
#' @param m The number of random process values to generate.
#' @param p The number of functions to evaluate at each time point.
#' @param domain The range of time points to generate.
#' @param mean_list A list of length p, each element is a mean function.
#' @param coef_list A list of length p, each element is a coefficient function.
#' @param covariancef A covariance function to generate random process values.
#' @param distribution The distribution to use for generating random process values.
#' @param snr Signal-to-noise ratio for generating the error variance.
#' @param sig Variance of the random error.
#' @param depend Logical indicating if the error is dependent or independent.
#' @param trans Logical indicating if the process values should be transformed.
#' @param transfmatrix A transformation matrix to use for transforming process values.
#' @param num_basis The number of basis functions to use for generating random process values.
#'
#' @return A list of random process values with error.
#'
#' @examples
#' \dontrun{
#' mean_list <- rep(list(function(x) {0}), 5)
#' covariancef <- function(t1, t2) {exp(-abs(t1 - t2))}
#' values_error <- generate_random_process_values_error(n = 100, m = 1, mean_list = mean_list, covariancef = covariancef)
#' }
#' @export

VCM_process_values_error <- function(n, m, p = 5, domain=c(0,1), mean_list = rep(list(function(x) {0*x}), 5),
                                     coef_list, covariancef,
                                                 distribution = 'normal', snr = 5, sig = NULL, depend = FALSE,
                                                 trans = FALSE, transfmatrix = NULL, num_basis = 1000){
  # Determine the variance of random error
  if( is.null(sig) ) {
    var_error <- l2moment_of_random_process(meanf = function(x){0},
                                            covariancef = covariancef, domain = domain)/ snr
  } else {
    var_error <- sig^2
  }

  # Generate time points lists
  time_points_list <- generate_time_points(n, m, domain = domain)

  # Apply the data generation function to each list
  process_values_list <- lapply(time_points_list, generate_random_process_values, p = p, meanf = function(x){0},
                                covariancef = covariancef, num_basis = num_basis, distribution = distribution)

  # Transform the independent p-dimensional process to a dependent one by matrix
  if (trans == T) {
    if (is.null(transfmatrix)) {
      transfmatrix <- toeplitz(seq(from = 1, by = -1/p, length.out = p))
    }
    trans_values <- lapply(process_values_list, function(x){
      x[,2:(p+1)] <- x[,2:(p+1)] %*% transfmatrix
      return(x)
    })
  } else {
    trans_values <- process_values_list
  }

  # Add mean function
  evaluated_mean <- lapply(time_points_list, function(t, mean_list) {
    value <- sapply(mean_list, function(f){f(t)})
    return(value)
  }, mean_list = mean_list)

  # Function for element-wise sum
  element_wise_sum_of_matrices <- function(list1, list2) {
    result_list <- list()
    for (i in 1:length(list1)) {
      result_list[[i]] <- list1[[i]] + cbind(rep(0,nrow(list2[[i]])),list2[[i]])
    }
    return(result_list)
  }

  values_mean <- element_wise_sum_of_matrices(trans_values,evaluated_mean)

  # Multiply coefficients
  evaluated_coef <- lapply(time_points_list, function(t, coef_list) {
    value <- sapply(coef_list, function(f){f(t)})
    return(value)
  }, coef_list = coef_list)

  # Function for element-wise product
  element_wise_prod_of_matrices <- function(list1, list2) {
    result_list <- list()
    for (i in 1:length(list1)) {
      result_list[[i]] <- cbind(list1[[i]], rowSums( list1[[i]][,-1] * list2[[i]] ))
    }
    return(result_list)
  }

  values_txy <- element_wise_prod_of_matrices(values_mean, evaluated_coef)

  # Add random error
  if (depend == T) {
    values_with_error <- lapply(values_txy, function(x){
      x[,2:(p+1)] <- x[,2:(p+1)] + x[,2] * rnorm(1, sd =sqrt(var_error))
      return(x)
    })
  } else {
    values_with_error <- lapply(values_txy, function(x){
      x[,2:(p+1)] <- x[,2:(p+1)] + rnorm(nrow(x), sd =sqrt(var_error))
      return(x)
    })
  }


  return(values_with_error)
}


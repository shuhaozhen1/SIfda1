# Tool Functions

split_matrix_by_rows <- function(mat, lengths) {
  cumulative_lengths <- c(0, cumsum(lengths))
  submatrices <- lapply(1:(length(cumulative_lengths) - 1), function(i) mat[(cumulative_lengths[i] + 1):cumulative_lengths[i + 1], ])
  submatrices
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




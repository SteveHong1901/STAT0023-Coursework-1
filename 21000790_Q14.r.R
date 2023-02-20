#QUESTION 14

###########################################################################################################################################################
##
##    Function cmds() to conduct multidimensional scaling on a distance matrix. 
##    
##    
##    The arguments are:
##
##    d: A square matrix that represents the distance between pairs of points.
##    s: An integer that specifies the number of dimensions in the 
##       new coordinate system. Default value is number of rows of d minus 1
##
##    For the cmds function, some edge cases include:
##
##    1) d is not a matrix
##    2) Asymmetric distance matrix: If the distance matrix is not symmetric, the function will give a warning.
##    3) Empty distance matrix: If the distance matrix is empty (i.e., has zero rows or columns), the function will produce an error due to an invalid index or dimension.
##    4) Exceeding number of positive eigenvalues (s > r+) or all s is 0
##    5) s is not a positive integer
##    
##    ==> The function will raise an error with a descriptive message indicating the issue.
##
##    The function returns:
##
##    evals: A vector containing the eigenvalues of B.
##    Y:     A matrix containing the new coordinates in the new coordinate system.
##    dY:    A squared distance matrix obtained from the new coordinates.
##    ssd:   A value representing the sum of squared differences between the original 
##           and the new distance matrix.
##
###########################################################################################################################################################


cmds <- function(d, s = nrow(d) - 1) {
  
  # INPUT CHECKING
  
  # Input checking for d
  
  if (!is.matrix(d)) {
    stop("Error: d must be a matrix")
  }
  
  if (!isSymmetric(d)) {
    stop("Warning: d is not symmetric")
  }
  
  if (nrow(d) == 0 || ncol(d) == 0) {
    stop("Error: d is an empty distance matrix")
  }
  
  
  # Input checking for s
  if (!is.numeric(s) || s < 1) {
    stop("Error: s must be a positive integer")
  }
  
  # COMPUTATION
  
  # Find the centred distance matrix
  n <- nrow(d)
  J <- diag(n) - matrix(1, n, n) / n
  B <- -0.5 * J %*% d^2 %*% J
  
  # Compute eigenvalues and eigenvectors (evecs are sorted in decreasing order by default)
  e <- eigen(B, symmetric = TRUE)
  evals <- e$values
  evecs <- e$vectors

  # Slicing the top s positive eigenvalues and corresponding eigenvectors
  
  # Count the no of positive eigenvalues
  r <- sum(evals > 0)
  
  # Following the instruction in the question, we can double check if s is valid and raise error
  if (s > r) {
    warning("Error: Invalid s, there is negative eigenvalues in the selected top positive eigenvalues")
  }
  if (s == 0) {
    stop("Error: Invalid s, there is no positive eigenvalues, d is not suitable for this function")
  }
  
  # Select the s largest eigenvalues and the corresponsing eigenvectors
  evals_reduced <- evals[1:s]
  evecs_reduced <- evecs[, 1:s]
  
  
  # Calculating new coordinates
  Y <- evecs_reduced %*% sqrt(diag(evals_reduced))
  
  # Calculating implied distance matrix and SSD
  dY <- as.matrix(dist(Y))
  ssd <- sum((d - dY)^2)

  # Returning results in a list
  return(list(Y = Y, d.centre = B, evals = evals, dY = dY, ssd = ssd))
}

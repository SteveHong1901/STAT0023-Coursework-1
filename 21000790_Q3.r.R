#QUESTION 3

#PART 1

#############################################################################################################################################
##
##    Function normprob() calculates the exact value of P(a<X≤b)
##    and its approximation using Simpson's rule
##    
##    
##    The arguments are:
##
##    mu:      a numeric input representing the mean of the normal distribution.
##    sigmasq: a numeric input representing the variance of the normal distribution.
##    a:       a numeric input representing the lower limit of the interval.
##    b:       a numeric input representing the upper limit of the interval.
##    n:       an optional numeric input representing the number of subintervals 
##             to be used in the Simpson's rule approximation. Default value is 100.
##
##    For the normprob function, some edge cases include:
##
##    1) sigmasq is zero or negative
##    2) n is zero, negative or odd
##    3) a is equal to or greater than b
##    4) the approximated probability (p.approx) is not inclusively between 0 and 1
##    ==> The function will raise an error with a descriptive message indicating the issue.
##
##    The function returns:
##
##    p.exact:  a numeric output representing the exact probability of the normal random variable falling between a and b.
##    p.approx: a numeric output representing the Simpson's rule approximation to the probability of the normal random 
##              variable falling between a and b.
##    n:        the number of subintervals used in the Simpson's rule approximation. 
##    p.error:  a numeric output representing the absolute error between the exact probability and the Simpson's rule approximation.
##
##
#############################################################################################################################################


normprob <- function(mu, sigmasq, a, b, n = 100) {
  
  # INPUT CHECKING
  
  # Raise warning for non-positive sigmasq, odd n and a >= b
  
  if (sigmasq <= 0) {
    stop("Error: sigmasq should be positive")
  }
  if (n %% 2 == 1 | n <= 0) {
    stop("Error: n should be even and positive")
  }
  if (a>=b) {
    stop("Error: incorrect limits, a must be smaller than b")
  }
  
  
  # COMPUTATION
  
  # Calculate the exact probability using pnorm()
  p.exact <- pnorm(b, mean = mu, sd = sqrt(sigmasq), lower.tail = TRUE) - 
    pnorm(a, mean = mu, sd = sqrt(sigmasq), lower.tail = TRUE)
  
  # Define the function to be integrated in Simpson's rule
  f <- function(x) {
    dnorm(x, mean = mu, sd = sqrt(sigmasq))
  }
  
  # Calculate the Simpson approximation to the probability
  h <- (b - a) / n
  x <- seq(a, b, by = h)
  y <- f(x)
  
  # The Indices used in this Simpson's Formula only works for n > 2, so we partition into the cases where n = 2 and n > 2
  
  if (n == 2) {
  p.approx <- (h / 3) * (y[1] + 4*y[2] + y[3  ])
  }
  if (n > 2) {
  p.approx <- (h / 3) * (y[1] + 4*sum(y[seq(2, n, by = 2)]) + 2*sum(y[seq(3, n-1, by = 2)]) + y[n+1])
  }
  
  # Check if the approximated probability is inclusively between 0 and 1
  if (p.approx > 1 | p.approx < 0) {
    warning("Warning: Invalid approximation, P(a<X≤b) is not inclusively between 0 and 1")
  }
  
  # Calculate the absolute error
  p.error <- p.approx - p.exact
  
  # Return the results as a list 
  return(list(p.exact = p.exact, p.approx = p.approx, n = n, p.error = p.error))
}



#PART 2

#############################################################################################################################################
##
##    Function SimpsonTest() evaluate the error of the Simpson approximation for P(a < X ≤ b) 
##    for different values of n, thus calculate alpha
##    
##    
##    The arguments are:
##
##    mu:      a numeric input representing the mean of the normal distribution.
##    sigmasq: a numeric input representing the variance of the normal distribution.
##    a:       a numeric input representing the lower limit of the interval.
##    b:       a numeric input representing the upper limit of the interval.
##    n.grid:  a numeric vector representing the different values of n to be used 
##             in the Simpson's rule approximations.
##
##    For the normprob function, some edge cases include:
##
##    1) sigmasq is zero or negative 
##    2) Any value in n.grid is zero, negative or odd
##    3) a is equal to or greater than b
##    4) the approximated probability (p.approx) is not inclusively between 0 and 1
##    ==> The function will raise an error with a descriptive message indicating the issue.
##
##    The function returns:
##
##    n.grid:       a numeric vector representing the different values of n to be used 
##                  in the Simpson's rule approximations.
##    abs.error:    a numeric vector representing the errors of Simpson's rule approximations for each value of n in n.grid.
##    alpha:        a numeric output representing the slope of the regression line
##                  
##
#############################################################################################################################################


SimpsonTest <- function(mu, sigmasq, a, b, n.grid) {
  
  # INPUT CHECKING
  
  # Raise warning for non-positive sigmasq, odd or negative n in n.grid and a >= b
  
  if (sigmasq <= 0 ) {
    stop("Error: sigmasq should be positive ")
  }
  if (any(n.grid %% 2 == 1 & n.grid <= 0)) {
    stop("Error: all n values in n.grid should be even and positive")
  }
  if (a>=b) {
    stop("Error: incorrect limits, a must be smaller than b")
  }
  
  # COMPUTATION
  
  # Evaluate Simpson approximation errors for each value of n in n.grid
  
  # Create an empty vector(p.error.list) of 0 to append the errors in for later
  p.error.list <- rep(0, length(n.grid))
  
  # Append errors one by one, corresponding to each value of n in n.grid
  for (i in seq_along(n.grid)) {
    result <- normprob(mu, sigmasq, a, b, n = n.grid[i])
    p.error.list[i] <- abs(result$p.approx - result$p.exact)
  }
  
  # Regress log of absolute approximation errors against log of n.grid
  regression <- lm(log(p.error.list) ~ log(n.grid))
  alpha <- -1 * regression$coefficients[2]
  
  # Return results as a list
  return(list(n.grid = n.grid, abs.error = p.error.list, alpha = alpha))
}



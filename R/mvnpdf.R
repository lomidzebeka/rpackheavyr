#' Title
#'
#' @param x is a matrix containing the observations in columns
#' @param mean mean vector of the multivariate gaussian distribution
#' @param M variance covariance matrix
#' @param Log logical parameter (default TRUE)
#'
#' @return a list of the results
#' @export
#'
#' @examples
#' X <- matrix(c(-0.5, 1.5, 0, 1, -1, 1), nrow = 2)


mvnpdf <-function (x,mean,M,Log=TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  output <- c()
  for (i in 1:n) {
    W <- x[,i, drop = FALSE]
    output <- c(output, (2*pi)^(-p/2)*det(M)^(-1/2)*exp(-(1/2)%*%t(W-mean)%*%solve(M)%*%(W-mean)))
  }
  return(output)
}


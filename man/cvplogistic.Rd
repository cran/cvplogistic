\name{cvplogistic}
\Rdversion{2.0}
\alias{cvplogistic}
\title{Majorization minimization by coordinate descent for concave
  penalized logistic regression }
\description{Compute solution surface for a high-dimensional logistic
  regression model with concave penalty using MMCD, adaptive rescaling
  or LLA-CD algorithms}
\usage{
cvplogistic(y, x, penalty = "mcp", approach = "mmcd", kappa = 1/2.7,
nlambda = 100, lambda.min = 0.01, epsilon = 1e-3, maxit = 1e+3)
}

\arguments{
  \item{y}{response vector with elements 0 or 1. }
  \item{x}{the design matrix of penalized variables. By default, an
  intercept vector will be added when fitting the model.}
  \item{penalty}{a character specifying the penalty. One of "mcp" or
    "scad" should be specified, with "mcp" being the default. }
  \item{approach}{a character specifying the numerical algorithm. One of
    "mmcd", "adaptive" or "llacda" can be specified, with "mmcd" being the
    default. See following details for more information. }
  \item{kappa}{a value specifying the regulation parameter kappa. The
    proper range for kappa is [0, 1).}
  \item{nlambda}{a integer value specifying the number of grids along the
    penalty parameter lambda.}
  \item{lambda.min}{a value specifying how to determine the minimal value of
    penalty parameter lambda. We define lambda_min=lambda_max*lambda.min.
    We suggest lambda.min=0.0001 if n>p; 0.01 otherwise.}
  \item{epsilon}{a value specifying the converge criterion of algorithm.}
  \item{maxit}{an integer value specifying the maximum number of iterations for
    each coordinate.}
}

\details{
  The computation for logistic model with concave penalties is not easy.
  The MMCD package implements the majorization minimization by coordinate
  descent (MMCD) algorithm for computing the solution path for logistic
  model with SCAD or MCP penalties. The algorithm is very efficient and
  stable for high-dimensional data with p>>n.
  For the MCP penalty, the package also implements the adaptive
  rescaling and the local linear approximation by coordinate descent
  algorithms (LLA-CDA) algorithms. For SCAD, only the MMCD algorithm is 
  implemented.

  The regularization parameter controls the concavity of the penalty,
  with larger value of kappa being more concave. When kappa=0, both the
  MCP and SCAD penalty become Lasso penalty. Hence if zero is specified
  for kappa, the algorithm returns Lasso solutions.
}

\value{
  A list of two elements is returned.
  \item{coef}{a matrix of dimension (p+1)*nlambda, with
    p the number of variables (columns) in x. The 1st row is the
    intercept, which is added by default. }
  \item{lambdas}{a vector of length nlambda for the penalty
    parameter lambda, ranging from the largest to the smallest.}
}
\author{Dingfeng Jiang}
\references{Dingfeng Jiang, Jian Huang. Majorization Minimization by
  Coordinate Descent for Concave Penalized Generalized Linear Models.

  Zou, H., Li, R. (2008). One-step Sparse Estimates in Nonconcave Penalized
  Likelihood Models. \emph{Ann Stat}, 364: 1509-1533.

  Breheny, P., Huang, J. (2011). Coordinate Descent Algorithms for Nonconvex
  Penalized Regression, with Application to Biological Feature
  Selection. \emph{Ann Appl Stat}, 5(1), 232-253.

  Jiang, D., Huang, J., Zhang, Y. (2011). The Cross-validated AUC for
  MCP-Logistic Regression with High-dimensional Data. \emph{Stat Methods
  Med Res}, online first, Nov 28, 2011.
}

\seealso{\code{hybrid.logistic}, \code{cv.cvplogistic},
  \code{cv.hybrid}, \code{path.plot}}
\examples{
set.seed(10000)
n=100
y=rbinom(n,1,0.4)
p=10
x=matrix(rnorm(n*p),n,p)

## MCP penalty by MMCD algorithm
out=cvplogistic(y, x, "mcp", "mmcd")
## MCP by adaptive rescaling algorithm
## out=cvplogistic(y, x, "mcp", "adaptive")
## MCP by LLA-CD algorith,
## out=cvplogistic(y, x, "mcp", "llacd")
## SCAD penalty
## out=cvplogistic(y, x, "scad")

## Lasso penalty
out=cvplogistic(y, x, kappa =0)

}
\keyword{models}
\keyword{regression}

\name{cv.cvplogistic}
\Rdversion{2.0}
\alias{cv.cvplogistic}
\title{Tuning parameter selection by k-fold cross validation for concave
  penalized logistic model}
\description{Using k-fold cross-validated area under ROC curve (CV-AUC) to select
  tuning parameter for high-dimensional logistic model with concave penalty}
\usage{
cv.cvplogistic(y, x, penalty = "mcp", approach = "mmcd", nfold = 5,
kappa = 1/2.7, nlambda = 100, lambda.min = 0.01,
epsilon = 1e-3, maxit = 1e+3, seed = 1000)
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
  \item{nfold}{an integer value for k-fold cross validation.}
  \item{kappa}{a value specifying the regulation parameter kappa. The
    proper range for kappa is [0, 1).}
  \item{nlambda}{an integer value specifying the number of grids along the
    penalty parameter lambda.}
  \item{lambda.min}{a value specifying how to determine the minimal value of
    penalty parameter lambda. We define lambda_min=lambda_max*lambda.min.
    We suggest lambda.min=0.0001 if n>p; 0.01 otherwise.}
  \item{epsilon}{a value specifying the converge criterion of algorithm.}
  \item{maxit}{an integer value specifying the maximum number of iterations for
    each coordinate.}
  \item{seed}{randomization seed for cross validation.}
}

\details{
  The computation for logistic model with concave penalties is not easy.
  The MMCD package implements the majorization minimization by coordinate
  descent (MMCD) algorithm for computing the solution path for logistic
  model with SCAD or MCP penalties. The algorithm is very efficient and
  stable for high-dimensional data with p>>n.
  For MCP penalty, the package also implements the adaptive rescaling
  and the local linear approximation by coordinate descent algorithms
  (LLA-CDA) algorithms. For SCAD, only the MMCD algorithm is
  implemented.

  The regularization parameter controls the concavity of the penalty,
  with larger value of kappa being more concave. When kappa=0, both the
  MCP and SCAD penalty become Lasso penalty. Hence if zero is specified
  for kappa, the algorithm returns Lasso solutions.

  To select an appropriate tuning parameter for prediction, we use k-fold
  cross-validated area under ROC curve (CV-AUC) approach. The CV-AUC
  approach calculated the predictive AUC for each validation set by
  using the coefficients estimated from the corresponding training
  set. As the cross validation proceeds, the average predictive AUC is
  calculated. Then the CV-AUC approach chooses the lambda
  corresponding to the maximum average predictive AUC as the tuning
  parameter.
}

\value{
  A list of three elements is returned.
  \item{scvauc}{the CV-AUC corresponding to the selected lambda.}
  \item{slambda}{the selected lambda.}
  \item{scoef}{the regression coefficients corresponding to the selected
    lambda, with the first element being the intercept.}

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

\seealso{\code{cvplogistic}, \code{hybrid.logistic},  \code{cv.hybrid},
  \code{path.plot}} 
\examples{
set.seed(10000)
n=100
y=rbinom(n,1,0.4)
p=10
x=matrix(rnorm(n*p),n,p)

## MCP penalty by MMCD algorithm
out=cv.cvplogistic(y, x, "mcp", "mmcd")
## MCP by adaptive rescaling algorithm
## out=cv.cvplogistic(y, x, "mcp", "adaptive")
## MCP by LLA-CD algorith,
## out=cv.cvplogistic(y, x, "mcp", "llacd")
## SCAD penalty
## out=cv.cvplogistic(y, x, "scad")

## Lasso penalty
## out=cv.cvplogistic(y, x, kappa =0)
}
\keyword{models}
\keyword{regression}

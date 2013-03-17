\name{cv.hybrid}
\Rdversion{2.0}
\alias{cv.hybrid}
\title{Tuning parameter selection by k-fold cross validation for logistic
  models with Lasso-concave hybrid penalty}
\description{Using k-fold cross-validated area under ROC curve to select
  tuning parameter for high-dimensional logistic model with
  Lasso-concave hybrid penalty}
\usage{
cv.hybrid(y, x, penalty = "mcp", nfold = 5,
kappa = 1/2.7, nlambda = 100, lambda.min = 0.01,
epsilon = 1e-3, maxit = 1e+3, seed = 1000)
}

\arguments{
  \item{y}{response vector with elements 0 or 1. }
  \item{x}{the design matrix of penalized variables. By default, an
  intercept vector will be added when fitting the model.}
  \item{penalty}{a character specifying the penalty. One of "mcp" or
    "scad" should be specified, with "mcp" being the default. }
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
  A Lasso-concave hybrid penalty applies SCAD or MCP penalty only to the
  variables selected by Lasso. The idea is to use Lasso as a screen tool
  to filter variables, then apply the SCAD or MCP penalty to the
  variables selected by Lasso for further selection. The computation for
  the hybrid penalty is faster than the standard concave penalty. The
  risk of using the hybrid penalty is that the variable missed by Lasso
  penalty will also not selected by the SCAD/MCP penalty.

  We also use the CV-AUC approach to select tuning parameter for models
  using the Lasso-concave hybrid penalty.
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

\seealso{\code{cvplogistic}, \code{hybrid.logistic},
  \code{cv.cvplogistic}, \code{path.plot}}
\examples{
set.seed(10000)
n=100
y=rbinom(n,1,0.4)
p=10
x=matrix(rnorm(n*p),n,p)

## Lasso-concave hybrid using MCP penalty
out=cv.hybrid(y, x, "mcp")
## Lasso-concave hybrid using SCAD penalty
## out=cv.hybrid(y, x, "scad")
}
\keyword{models}
\keyword{regression}

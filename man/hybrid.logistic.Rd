\name{hybrid.logistic}
\Rdversion{2.0}
\alias{hybrid.logistic}
\title{A Lasso-concave hybrid penalty for logistic regression}
\description{Compute solution surface for a high-dimensional logistic
  regression model with Lasso-concave hybrid penalty for fast variable selection}
\usage{
hybrid.logistic(y, x, penalty = "mcp", kappa = 1/2.7,
nlambda = 100, lambda.min = 0.01, epsilon = 1e-3, maxit = 1e+3)
}

\arguments{
  \item{y}{response vector with elements 0 or 1. }
  \item{x}{the design matrix of penalized variables. By default, an intercept
    vector will be added when fitting the model.}
  \item{penalty}{a character specifying the penalty. One of "mcp" or
    "scad" should be specified, with "mcp" being the default. }
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
  A Lasso-concave hybrid penalty applies SCAD or MCP penalty only to the
  variables selected by Lasso. The idea is to use Lasso as a screen tool
  to filter variables, then apply the SCAD or MCP penalty to the
  variables selected by Lasso for further selection. The computation for
  the hybrid penalty is faster than the standard concave penalty. The
  risk of using the hybrid penalty is that the variable missed by Lasso
  will also not selected by the SCAD/MCP penalty.
}

\value{
  A list of two elements is returned.
  \item{coef}{A matrix of dimension (p+1)*nlambda, with
    p the number of variables (columns) in x. The 1st row is the
    intercept, which is added by default. }
  \item{lambdas}{A vector of length nlambda for the penalty
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

\seealso{\code{cvplogistic}, \code{cv.cvplogistic}, \code{cv.hybrid},
  \code{path.plot}} 
\examples{
set.seed(10000)
n=100
y=rbinom(n,1,0.4)
p=10
x=matrix(rnorm(n*p),n,p)

## Lasso-concave hybrid using MCP penalty
out=hybrid.logistic(y, x, "mcp")
## Lasso-concave hybrid using SCAD penalty
## out=hybrid.logistic(y, x, "scad")
}
\keyword{models}
\keyword{regression}

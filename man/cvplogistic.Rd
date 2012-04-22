\name{cvplogistic}
\Rdversion{2.0}
\alias{cvplogistic}
\title{Majorization minimization by coordinate descent for concave penalized logistic regression}
\description{Compute solution surface for a high-dimensional logistic
  regression model with concave penalty using MMCD, adaptive rescaling
  or LLA algorithms}
\usage{
cvplogistic(y, x, penalty = "mcp", approach = "mmcd", path = "kappa",
nkappa = 10, maxkappa = 0.249, nlambda = 100, minlambda = 0.01,
epsilon = 1e-3, maxit = 1e+3)
}

\arguments{
  \item{y}{response vector with elements 0 or 1. }
  \item{x}{the design matrix of variables. By default, an intercept
    vector will be added when fitting the model.}
  \item{penalty}{a character specifying the penalty. One of "mcp" or
    "scad" should be specified, with "mcp" being the default. }
  \item{approach}{a character specifying the numerical algorithm. One of
    "mmcd", "adaptive" or "llacda" can be specified, with "mmcd" being the
    default. See following details for more information. }
  \item{path}{a character indicating how the solution surface should be
    computed. Three values are allowed, "kappa","lambda" or "hybrid". See
    following details for more information.}
  \item{nkappa}{a integer value specifying the number of grids along the
    regulation parameter kappa. }
  \item{maxkappa}{a value specifying the maximum of the regulation
    parameter kappa. The proper range for kappa is [0,0.25) for MCP and
    [0,0.2) for SCAD if using the MMCD algorithm.}
  \item{nlambda}{a integer value specifying the number of grids along the
    penalty parameter lambda.}
  \item{minlambda}{a value specifying how to determine the minimal value of
    penalty parameter lambda. We define lambda_min=lambda_max*minlambda.
    By default, minlambda is 0.0001 if n>p; 0.01 otherwise.}
  \item{epsilon}{a value specifying the converge criterion of algorithm.}
  \item{maxit}{an integer value specifying the maximum number of iterations for
    each coordinate.}
}

\details{
  The package implements the majorization minimization by coordinate
  descent (MMCD) algorithm for computing the solution surface of a concave
  penalized logistic regression model in high-dimensional data. The MMCD
  algorithm seeks a closed form solution for each coordinate and
  majorizes the loss function to avoid the computation of scaling
  factors. The algorithm is efficient and stable for high-dimensional
  data with p>>n.
  
  The package provides three ways to compute solution surfaces for a
  concave penalized logistic model. The first one is to compute along the
  regularization parameter kappa. That is for a given penalty parameter
  lambda, using the Lasso solution (kappa=0) to initiate the computation
  for MCP or SCAD solutions along grid values of kappa. The second type
  is to compute along the penalty parameter lambda. That is for a given
  regularization parameter kappa, the MCP or SCAD solutions are computed
  along lambda. The solution surface computed along kappa tends to have
  a better performance in terms of model size and false discovery
  rate. Thus, the solution surface along kappa is recommended.
  
  The third type of solution is called hybrid algorithm. The
  hybrid algorithm is specifically designed for the applications aimed
  to identify leading causal predictors. In most cases, the hybrid
  algorithm achieves the same predictive performance as the other two
  methods. This hybrid algorithm can be viewed as a variant of the
  solution surface along kappa. It uses Lasso solution (kappa=0) as the
  initial values and apply the MMCD algorithm to the variables selected
  by Lasso only. In another word, Lasso is used to pre-process the
  variables before applying the MCP penalty. This practice greatly
  reduces the computation burden. However, if Lasso misses some
  variables, they will necessarily be removed from final model.
  
  For MCP penalty, the package also implements the adaptive rescaling
  and the local linear approximation by coordinate descent algorithms
  (LLA-CDA) algorithms. 
  The three types of solution path is also available for the adaptive
  rescaling approach. 
  For the LLA-CDA algorithm, only the solution path along kappa is
  implemented for the LLA-CDA algorithm.
}

\value{
  A list with five elements is returned.
  \item{lambda}{A vector of length nkappa*nlambda for the penalty
    parameter lambda, ranging from the largest to the smallest with
    block size nkappa. }
  \item{kappa}{A vector of length nkappa*nlambda for the regularization
    parameter kappa, ranging from 0 to maxkappa within a block with size
    nkappa.}
  \item{df}{A vector of length nkappa*nlambda indicating the degree of
    freedom(model size, the number of covariates) for the corresponding
    solution. }
  \item{coef.intecept}{A vector of length nkappa*nlambda corresponding
    to the coefficient of the intercept.}
  \item{coef.covariates}{A matrix of dimension p*(nkappa*nlambda), with
    p the number of variables (columns) in x. }
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

\seealso{\code{aic.cvplogistic}, \code{bic.cvplogistic},
  \code{cvauc.cvplogistic}}
\examples{
set.seed(10000)
n=100
y=rbinom(n,1,0.4)
p=50
x=matrix(rnorm(n*p),n,p)
nkappa=5
maxkappa=0.249
nlambda=20
## MCP penalty
penalty="mcp"
approach="mmcd"
path="kappa"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="lambda"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="hybrid"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
approach="adaptive"
path="kappa"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="lambda"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="hybrid"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
## using LLA approach, path option has no effect.
approach="llacda"
maxkappa=0.99
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
## SCAD penalty
maxkappa=0.19
penalty="scad"
path="kappa"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="lambda"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="hybrid"
out=cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
}
\keyword{models}
\keyword{regression}

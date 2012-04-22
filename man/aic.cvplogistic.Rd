\name{aic.cvplogistic}
\Rdversion{2.0}
\alias{aic.cvplogistic}
\title{Tuning parameter selection by AIC criteria for a concave
  penalized logistic regression}
\description{Tuning parameter (kappa, lambda)  selection using AIC
  criteria for a concave penalized logistic regression. Only the models
  with df<n will be considered in the tuning parameter selection.}
\usage{
aic.cvplogistic(y, x, penalty = "mcp", approach = "mmcd",
path = "kappa", nkappa = 10, maxkappa = 0.249,
nlambda = 100, minlambda = 0.01,
epsilon = 1e-3, maxit = 1e+3)
}
\arguments{
  \item{y}{response vector with elements 0 or 1. }
  \item{x}{the design matrix of variables. By default, an intercept
    vector will be added when fitting the model.}
  \item{penalty}{a character specifying the penalty. One of "mcp" or
    "scad" should be specified, with "mcp" being the default. }
  \item{approach}{a character specifying the numerical algorithm. One of
    "mmcd", "adaptive" or "lla" can be specified, with "mmcd" being the
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
  The package implement the majorization minimization by coordinate
  descent (MMCD) algorithm for computing the solution surface of concave
  penalized logistic regression model in high-dimensional data. The MMCD
  algorithm seeks a closed form solution for each coordinate and
  majorizes the loss function to avoid the computation of scaling
  factors. The algorithm is efficient and stable for high-dimensional
  data with p>>n.

  The package provides three ways to compute solution surfaces for a
  concave penalized logistic model. The first one is compute along the
  regularization parameter kappa. That is  the Lasso solution (kappa=0)
  is used to initiate the computation for MCP or SCAD solutions, for a
  given penalty  parameter lambda. The second type is to compute along
  the penalty parameter lambda. That is for a given regularization
  parameter kappa, the MCP or SCAD solutions are computed along lambda.
  The solution surface computed along kappa tends to have a better
  performance in terms of model size and false discovery rate. Thus, the
  solution surface along kappa is recommended.

  The third type of solution is called hybrid algorithm. The
  hybrid algorithm is specifically designed for the applications which
  aims to identify the leading causal predictors. In most cases,   the
  hybrid algorithm achieves the same predictive performance as the
  solution  surface along kappa. This hybrid algorithm can be viewed as
  an variant of the solution  surface along kappa. In the hybrid
  algorithm, Lasso solution (kappa=0) is used as the initial values. The
  hybrid algorithm, however, only apply the MMCD algorithm to the
  variables selected by Lasso. That is   Lasso is used to pre-process
  the variables, this practice greatly   reduces the computation
  burden. However, if Lasso misses one variable,   it will necessarily
  removed from the final model.

  The tuning parameters including the regularization parameter kappa and
  penalty parameter lambda are determined by the AIC criterion. The
  solutions corresponding to the chosen tuning parameters is output as
  the solution for the model. We only consider the models with df<n as the
  candidate models to be chosen from. The AIC values are computed
  based on the solution surfaces computed as specified by the user.
}
\value{
  A list of four elements is returned.
  \item{saic}{the AIC value corresponding to the chosen tuning parameters.}
  \item{slambda}{the chosen penalty parameter.}
  \item{skappa}{the chosen regularization parameter.}
  \item{scoef}{the coefficients of variables in x
    corresponding to the chosen tuning parameter.}
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
\seealso{\code{cvplogistic}, \code{bic.cvplogistic},
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
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa,
nlambda)
path="lambda"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="hybrid"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
approach="adaptive"
path="kappa"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="lambda"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="hybrid"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
## using LLA approach, path option has no effect.
approach="llacda"
maxkappa=0.99
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
## SCAD penalty
maxkappa=0.19
penalty="scad"
path="kappa"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="lambda"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
path="hybrid"
out=aic.cvplogistic(y, x, penalty, approach, path, nkappa, maxkappa, nlambda)
}
\keyword{models}
\keyword{regression}

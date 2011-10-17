\name{cvplogistic}
\Rdversion{1.1}
\alias{cvplogistic}
\title{Majorization minimization by coordinate descent (MMCD) algorithm
  for concave penalized logistic regression for high-dimensional data}
\description{Using MMCD algorithm to compute solution surface for a
  high-dimensional logistic regression model with concave penalty}
\usage{
cvplogistic(y,x,penalty="mcp",path="kappa",nkappa=20,maxkappa=0.249,
nlambda=100,minlambda=ifelse(n>p,0.0001,0.01),epsilon=1e-3,maxit=1e+4)
}

\arguments{
  \item{y}{response vector with elements 0 or 1. }
  \item{x}{the design matrix of variables. By default, an intercept
    vector will be added into the model fitting.}
  \item{penalty}{a character indicating the penalty. One of "mcp" or
    "scad" should be specified, with "mcp" being the default. }
  \item{path}{a character indicating how the solution surface should be
    computed. Three values are allowed, "kappa","lambda", and
    "hybrid". See following details for more information.}
  \item{nkappa}{a integer value specifies the number of grids along the
    regulation parameter kappa. }
  \item{maxkappa}{a value specifies the maximum of the regulation
    parameter kappa. The proper range for kappa is [0,0.25) for MCP and
    [0,0.2) for SCAD.}
  \item{nlambda}{a integer value specifies the number of grids along the
    penalty parameter lambda.}
  \item{minlambda}{a value specifies how to compute the minimal value of
    penalty parameter, which is defined as lambda_min=lambda_max*minlambda.
    By default, minlambda is 0.0001 if n>p; 0.01 otherwise.}
  \item{epsilon}{a value specifies the converge criterion of algorithm.}
  \item{maxit}{a integer value specifies the maximum of iterations for
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
}

\value{
  A list with five elements is returned.
  \item{lambda}{A vector of length nkappa*nlambda for the penalty
    parameter lambda, ranging from the largest to the smallest with
    block size nkappa. }
  \item{kappa}{A vector of length nkappa*nlambda for the regularization
    parameter kappa, ranging from 0 to maxkappa within a block size
    nkappa.}
  \item{df}{A vector of length nkappa*nlambda indicating the degree of
    freedom(model size, the number of covariates) for the corresponding
    solution. }
  \item{coef.intecept}{A vector of length nkappa*nlambda corresponding
    to the coefficient of the intercept.}
  \item{coef.covariates}{A matrix of dimension p*(nkappa*nlambda), with
    p the number of variables in x. }
}
\author{Dingfeng Jiang}
\references{Dingfeng Jiang, Jian Huang. Majorization Minimization by
  Coordinate Descent for Concave Penalized Generalized Linear Models. }
\seealso{\code{aic.cvplogistic}, \code{bic.cvplogistic}, \code{auc.cvplogistic}}
\examples{
seed=10000
n=100
y=rbinom(n,1,0.4)
p=50
x=matrix(rnorm(n*p),n,p)
penalty="mcp"
nkappa=5
maxkappa=0.249
nlambda=20
path="kappa"
out=cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="lambda"
out=cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="hybrid"
out=cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
penalty="scad"
maxkappa=0.19
path="kappa"
out=cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="lambda"
out=cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="hybrid"
out=cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
}
\keyword{models}
\keyword{regression}

\name{bic.cvplogistic}
\Rdversion{1.1}
\alias{bic.cvplogistic}
\title{Tuning parameter selection by BIC criteria for a concave
  penalized logistic regression}
\description{Tuning parameter (kappa, lambda)  selection using BIC
  criteria for a concave penalized logistic regression. Only the models
  with df<n will be considered in the tuning parameter selection.}
\usage{
bic.cvplogistic(y,x,penalty="mcp",path="kappa",nkappa=20,maxkappa=0.249,
nlambda=100,minlambda=ifelse(n>p,0.0001,0.01),epsilon=1e-3,maxit=1e+4)
}
\arguments{
  \item{y}{response vector with elements 0 or 1. }
  \item{x}{the design matrix of variables. By default, an intercept
    vector will be added into the model fitting.}
  \item{penalty}{a character indicating the penalty. One of "mcp" or "scad"
    should be specified, with "mcp" being the default. }
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
  
  The tuning parameters including the regularization parameter kappa and
  penalty parameter are determined by the BIC criterion. The
  solutions corresponding to the chosen tuning parameters is output as
  the solution for the model. We only consider the models with df<n as the
  candidate models to be chosen from. The BIC values are computed
  based on the solution surfaces computed as specified by the user.
}
\value{
  A list of five elements is returned.
  \item{tuning.BIC}{The BIC value corresponding to chosen tuning
    parameter.}
  \item{tuning.lambda}{The chosen penalty parameter.}
  \item{tuning.kappa}{The chosen regularization parameter.}
  \item{tuning.intercept}{The intercept coefficient corresponding to the
    chosen tuning parameter.}
  \item{tuning.covariates}{The coefficients of variables in x
    corresponding to the chosen tuning parameter.}
}
\author{Dingfeng Jiang}
\references{Dingfeng Jiang, Jian Huang. Majorization Minimization by
  Coordinate Descent for Concave Penalized Generalized Linear Models. }
\seealso{\code{cvplogistic}, \code{aic.cvplogistic}, \code{auc.cvplogistic}}
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
out=bic.cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="lambda"
out=bic.cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="hybrid"
out=bic.cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
penalty="scad"
maxkappa=0.19
path="kappa"
out=bic.cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="lambda"
out=bic.cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
path="hybrid"
out=bic.cvplogistic(y,x,penalty,path,nkappa,maxkappa,nlambda)
}
\keyword{models}
\keyword{regression}

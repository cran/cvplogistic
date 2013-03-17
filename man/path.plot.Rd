\name{path.plot}
\Rdversion{2.0}
\alias{path.plot}
\title{Plot the solution path for the concave penalized logistic models}
\description{Plot the path trajectories for the solutions computed by
  the implemented methods.}
\usage{
path.plot(out)
}

\arguments{
  \item{out}{the object return from function \code{cvplogistic} or \code{hybrid.logistic}.}
}

\details{
  The function plots the trajectories of solutions, with x-axis being
  the grids of lambda, and y-axis being the coefficients profile.
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
\code{cv.cvplogistic}}
\examples{
set.seed(10000)
n=100
y=rbinom(n,1,0.4)
p=10
x=matrix(rnorm(n*p),n,p)

## MCP
out=cvplogistic(y, x)
path.plot(out)
## hybrid penalty
## out=hybrid.logistic(y, x, "mcp")
## path.plot(out)
}
\keyword{models}
\keyword{regression}

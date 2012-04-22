## Following code is used to generate the simulated datasets
## Dingfeng Jiang

library('Matrix')
library('cvplogistic')

AR1 <- function(p, rho){
    sigma <- matrix(0, ncol = p, nrow = p)
    for(i in 1:p){
        for (j in 1:i){
            if (i == j) {
                sigma[i, j] = 1
            } else {
                sigma[i, j] = sigma[j, i] = rho^(abs(i - j))
            }
        }
    }
    sigma
}

CS <- function(p, rho){
    sigma <- matrix(rho, ncol = p, nrow = p)
    for(i in 1:p){
        for (j in 1:i){
            if (i == j) sigma[i, j] = 1
        }
    }
    sigma
}


ybi <- function(seed, n, p, coef, eS){
    set.seed(seed)
    rz <- rnorm(n*p)
    nz <- matrix(rz, n, p) %*% diag(sqrt(eS$values)) %*%t (eS$vectors)
    int <- rep(1, n)
    xz <- cbind(int, nz)
    eta <- xz %*% coef
    ppi <- 1/(1 + exp(-eta))
    ## generate binary outcome
    ybi <- rbinom(n, 1, ppi)
    lsnr <- sqrt( t(eta)%*%eta /n)
    return(list(ybi, int, nz, lsnr, lsnr))
}

ybiin <- function(seed, n, p, coef){
    set.seed(seed)
    nz <- matrix(rnorm(n*p), n, p)
    int <- rep(1, n)
    xz <- cbind(int, nz)
    eta <- xz %*% coef
    ppi <- 1/(1 + exp( -eta))
    ## generate binary outcome
    ybi <- rbinom(n, 1, ppi)
    lsnr <- sqrt( t(eta)%*%eta /n)
    return(list(ybi, int, nz, lsnr, lsnr))
}


## A wrapup function to generate the data.
## This mainly provide a nice way to present how the simulated data are generated.
##
## For simulation purpose, it is suggest no to use it since
## the eigen function will be called repeated. For large p,
## it takes a lot of time.


fsimu <- function(n, p, rho, alpha, beta1, type, seed=1000){
    tn <- n
    qq <- 1
    qp <- qq + p
    non0p <- length(beta1)
    beta2 <- rep(0, p - non0p)
    beta <- c(beta1, beta2)
    coef <- c(alpha, beta)
    betaidx <- beta != 0
    coefidx <- coef != 0
    ## some check
    if (length(beta1) != non0p) stop("Not match for beta1!\n")
    if (length(beta) != p) stop("Not match for beta!\n")
    if (length(coef) != qp) stop("Not match for coef!\n")
    ## different type
    if (type == "IN"){
        sigma <- diag(1, p)
        if (!all(dim(sigma) == c(p, p))) stop("Not match for sigma!\n")
        train <- ybiin(seed, tn, p, coef)
    } else if (type == "SP"){
        sigma1 <- CS(non0p,rho)
        p2 <- p-non0p
        sigma2 <- CS(p2, rho)
        sigma <- bdiag(sigma1,sigma2)
        eS <- eigen(sigma, symmetric = TRUE, EISPACK = TRUE)
        train <- ybi(seed,tn,p,coef,eS)
    } else if (type == "PC") {
        sigma1 <- CS(non0p/2, rho)
        sigma2 <- CS(non0p, rho)
        p3 <- p - 1.5*non0p
        sigma3 <- CS(p3, rho)
        sigma <- bdiag(sigma1, sigma2, sigma3)
        eS <- eigen(sigma,  symmetric  =  TRUE,  EISPACK  =  TRUE)
        train <- ybi(seed, tn, p, coef, eS)
    } else if (type == "AR") {
        sigma <- AR1(p, rho)
        eS <- eigen(sigma,  symmetric  =  TRUE,  EISPACK  =  TRUE)
        if (!all(dim(sigma) == c(p, p))) stop("Not match for sigma!\n")
        train <- ybi(seed, tn, p, coef, eS)
    } else if (type =="CS") {
        sigma <- CS(p, rho)
        eS <- eigen(sigma, symmetric = TRUE, EISPACK = TRUE)
        train <- ybi(seed,tn,p,coef,eS)
    }
    list(train[[1]], train[[3]])
}

n <- 500
p <- 100
rho <- 0.5
seed <- 1000
alpha <- 0
beta1 <- c(0.6, -0.6, 1.2, -1.2, 2.4, -0.6, 0.6, -1.2, 1.2, -2.4)
type <- "IN"
data <- fsimu(n, p, rho, alpha, beta1, type, seed)
type <- "SP"
data <- fsimu(n, p, rho, alpha, beta1, type, seed)
type <- "PC"
data <- fsimu(n, p, rho, alpha, beta1, type, seed)
type <- "AR"
data <- fsimu(n, p, rho, alpha, beta1, type, seed)
type <- "CS"
data <- fsimu(n, p, rho, alpha, beta1, type, seed)






## Following is the actual set up for the simulation
## Codes are broken into pieces depending on the specified structure

tn <- n
qq <- 1
qp <- qq + p
non0p <- 10
alpha <- 0.00
beta1 <- c(0.6, -0.6, 1.2, -1.2, 2.4, -0.6, 0.6, -1.2, 1.2, -2.4)
beta2 <- rep(0, p - non0p)
beta <- c(beta1, beta2)
coef <- c(alpha, beta)
betaidx <- beta != 0
coefidx <- coef != 0
## some check
if (length(beta1) != non0p) stop("Not match for beta1!\n")
if (length(beta) != p) stop("Not match for beta!\n")
if (length(coef) != qp) stop("Not match for coef!\n")



## IN strucutre
sigma <- diag(1, p)
if (!all(dim(sigma) == c(p, p))) stop("Not match for sigma!\n")
train <- ybiin(seed, tn, p, coef)
y <- train[[1]]
x <- train[[3]]

## SP strucutre
sigma1 <- CS(non0p,rho)
p2 <- p-non0p
sigma2 <- CS(p2, rho)
sigma <- bdiag(sigma1,sigma2)
eS <- eigen(sigma, symmetric = TRUE, EISPACK = TRUE)
train <- ybi(seed,tn,p,coef,eS)
y <- train[[1]]
x <- train[[3]]

## PC strucutre
sigma1 <- CS(non0p/2, rho)
sigma2 <- CS(non0p, rho)
p3 <- p - 1.5*non0p
sigma3 <- CS(p3, rho)
sigma <- bdiag(sigma1, sigma2, sigma3)
eS <- eigen(sigma,  symmetric  =  TRUE,  EISPACK  =  TRUE)
train <- ybi(seed, tn, p, coef, eS)
y <- train[[1]]
x <- train[[3]]

## AR strucutre
sigma <- AR1(p, rho)
eS <- eigen(sigma,  symmetric  =  TRUE,  EISPACK  =  TRUE)
if (!all(dim(sigma) == c(p, p))) stop("Not match for sigma!\n")
train <- ybi(seed, tn, p, coef, eS)
y <- train[[1]]
x <- train[[3]]


## CS strucutre
sigma <- CS(p, rho)
eS <- eigen(sigma, symmetric = TRUE, EISPACK = TRUE)
train <- ybi(seed,tn,p,coef,eS)
y <- train[[1]]
x <- train[[3]]

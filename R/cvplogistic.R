## Compute solution surface for a concave penalized logistic regression model
## using majorization minimization by coordinate descent (MMCD) algorithm.
## Two concave penalties are considered: SCAD and MCP.
## Three types of solution surfaces are provided for each penalty,
## 1. solution surface computed along kappa
## 2. solution surface computed along lambda
## 3. solution surface using Lasso-MCP hybrid algorithm
## The adaptive resaling approach and local linear approximation are also
## provided as optional choices.
## April 15, 2012
## dyn.load("../src/cvplogistic.dll")


cvplogistic <- function(y, x, penalty = "mcp", approach = "mmcd", path = "kappa",
                        nkappa = 10, maxkappa = 0.249, nlambda = 100,
                        minlambda = 0.01, epsilon = 1e-3, maxit = 1e+3){
    ## error checking
    if (nrow(x) != length(y)) stop("# of rows in X does not match the length of Y! \n")
    ## penalty
    pen <- pmatch(penalty, c("mcp", "scad"))
    if (is.na(pen)) stop("Penalty need to be either 'mcp' or 'scad'!\n")
    ## penalty and computational approach
    app <- pmatch(approach, c("mmcd", "adaptive", "llacda"))
    if ((pen == 1) & (app == 1) & ((maxkappa >= 0.25) | (maxkappa < 0))) {
        stop("Using MMCD algorithm for MCP penalty, the regulation parameter kappa should be in [0, 0.25)!\n")
    }
    if ((pen == 1) & (app == 2) & ((maxkappa >= 1.0) | (maxkappa < 0))) {
        stop("Using adaptive rescaling algorithm for MCP penalty, the regulation parameter kappa should be in [0, 1.0)!\n")
    }
    if ((pen == 1) & (app == 2) & (maxkappa > 0.25)) {
        warning("Using adaptive rescaling algorithm for MCP penalty, the algorithm may not converge for large kappa!\n")
    }
    if ((pen == 1) & (app == 3) & ((maxkappa >= 1.0) | (maxkappa < 0)))  {
        stop("Using LLA-CDA algorithm for MCP penalty, the regulation parameter kappa should be in [0, 1.0)!\n")
    }
    if ((pen == 2) & ((maxkappa >= 0.2) | (maxkappa < 0))) {
        stop("Using MMCD algorithm for SCAD penalty, the regulation parameter kappa should be in [0, 0.2)!\n")
    }
    ## solution surface
    ss <- pmatch(path, c("kappa", "lambda", "hybrid"))
    if (is.na(ss)) stop("Solution path should be one of 'kappa',  'lambda', 'hybrid'!\n")
    ## space assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq+p
    ## create output space
    olmdas <- rep(0, nkappa*nlambda)
    okas <- rep(0, nkappa*nlambda)
    ocoef <- matrix(0, qp, nkappa*nlambda)
    oaic <- rep(0, nkappa*nlambda)
    obic <- rep(0, nkappa*nlambda)
    oobj <- rep(0, nkappa*nlambda)
    odf <- rep(0, nkappa*nlambda)
    ocvx <- rep(0, nkappa*nlambda)
    ## fit through Fortran depend on pen, app and ss
    if (pen == 1){
        if (app == 1){
            if (ss == 1) {
                out <- try(.Fortran("mcpkapa",
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx),
                                    as.double(y), as.double(int), as.double(x),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon), as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 2) {
                out <- try(.Fortran("mcplmda",
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx),
                                    as.double(y), as.double(int), as.double(x),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 3){
                out <- try(.Fortran("mcpkapa2",
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx),
                                    as.double(y), as.double(int), as.double(x),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            }
        } else if (app == 2) {
            if (ss == 1) {
                out <- try(.Fortran("adpmcpkp",
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx),
                                    as.double(y), as.double(int), as.double(x),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 2) {
                out <- try(.Fortran("adpmcplm",
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx),
                                    as.double(y), as.double(int), as.double(x),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 3){
                out <- try(.Fortran("adpmcpkp2",
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx),
                                    as.double(y), as.double(int), as.double(x),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            }
        } else if (app == 3) {
            out <- try(.Fortran("fllabi",
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx),
                                as.double(y), as.double(int), as.double(x),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        }
    } else {
        if (ss == 1) {
            out <- try(.Fortran("scadkapa",
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx),
                                as.double(y), as.double(int), as.double(x),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (ss == 2) {
            out <- try(.Fortran("scadlmda",
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx),
                                as.double(y), as.double(int), as.double(x),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (ss == 3){
            out <- try(.Fortran("scadkapa2",
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx),
                                as.double(y), as.double(int), as.double(x),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        }
    }
    ## organize output
    if (!(inherits(out, 'try-error'))){
        lambda <- out[[1]]
        kappa <- out[[2]]
        coef <- matrix(out[[3]], qp, nkappa*nlambda)
        df <- out[[7]]
        ## ocvx <- out[[8]]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept", paste("x", 1:p, sep = ""))
        } else {
            rownames(coef) <- c("intercept", colnames(x))
        }
        coef.intercept <- coef[1, ]
        coef.covariates <- coef[-1, ]
        list(lambda, kappa, df, coef.intercept, coef.covariates)
    } else {
        stop("Model fitting fails,  double check!\n")
        NULL
    }
}


## Tuning parameter selection by AIC
aic.cvplogistic <- function(y, x, penalty = "mcp", approach = "mmcd", path = "kappa",
                            nkappa = 10, maxkappa = 0.249, nlambda = 100,
                            minlambda = 0.01,
                            epsilon = 1e-3, maxit = 1e+3){
    ## compute the sollution surface
    ssout <- try(cvplogistic(y, x, penalty, approach, path,
                             nkappa, maxkappa, nlambda, minlambda,
                             epsilon, maxit), TRUE)
    xz <- cbind(rep(1, nrow(x)), x)
    if (!(inherits(ssout, 'try-error'))){
        lambdas <- ssout[[1]]
        kas <- ssout[[2]]
        df <- ssout[[3]]
        coef.int <- ssout[[4]]
        coef.cov <- ssout[[5]]
        coef <- rbind(coef.int, coef.cov)
        ## regular solution if df<n
        uidx <- nrow(x) > df
        ulambdas <- lambdas[uidx]
        ukas <- kas[uidx]
        udf <- df[uidx]
        ucoef <- coef[, uidx]
        ## compute the AIC
        logl <- apply(ucoef, 2, function(v){
            eta <- xz %*% coef
            sum(y * eta - log(1 + exp(eta)))
        })
        aic <- -2 * logl + 2 * udf
        ## selection based on AIC
        aicidx <- aic == min(aic)
        saic <- aic[aicidx][1]
        slambda <- ulambdas[aicidx][1]
        ska <- ukas[aicidx][1]
        scoef <- ucoef[, aicidx]
        if (is.matrix(scoef)) scoef <- scoef[, 1]
        list(saic, slambda, ska, scoef)
    } else {
        NULL
    }
}


## Tuning parameter selection by BIC
bic.cvplogistic <- function(y, x, penalty = "mcp", approach = "mmcd", path = "kappa",
                            nkappa = 10, maxkappa = 0.249, nlambda = 100,
                            minlambda = 0.01,
                            epsilon = 1e-3, maxit = 1e+3){
    ## compute the sollution surface
    ssout <- try(cvplogistic(y, x, penalty, approach, path,
                             nkappa, maxkappa, nlambda, minlambda,
                             epsilon, maxit), TRUE)
    xz <- cbind(rep(1, nrow(x)), x)
    if (!(inherits(ssout, 'try-error'))){
        lambdas <- ssout[[1]]
        kas <- ssout[[2]]
        df <- ssout[[3]]
        coef.int <- ssout[[4]]
        coef.cov <- ssout[[5]]
        coef <- rbind(coef.int, coef.cov)
        n <- nrow(x)
        ## regular solution if n > df
        uidx <- n > df
        ulambdas <- lambdas[uidx]
        ukas <- kas[uidx]
        udf <- df[uidx]
        ucoef <- coef[, uidx]
        ## compute the BIC
        logl <- apply(ucoef, 2, function(v){
            eta <- xz %*% coef
            sum(y * eta - log(1 + exp(eta)))
        })
        bic <- -2 * logl + log(n) * udf
        ## selection based on BIC
        bicidx <- bic == min(bic)
        sbic <- bic[bicidx][1]
        slambda <- ulambdas[bicidx][1]
        ska <- ukas[bicidx][1]
        scoef <- ucoef[, bicidx]
        if (is.matrix(scoef)) scoef <- scoef[, 1]
        list(sbic, slambda, ska, scoef)
    } else {
        NULL
    }
}




## Tuning parameter selection using CV-AUC
cvauc.cvplogistic <- function(cv = 5, stratified = TRUE, y, x, penalty = "mcp", approach = "mmcd",
                              path = "kappa", nkappa = 10, maxkappa = 0.249,
                              nlambda = 100, minlambda = 0.01,
                              epsilon = 1e-3, maxit = 1e+3, seed = 1000){
    ## error checking
    if (nrow(x) != length(y)) stop("# of rows in X does not match the length of Y! \n")
    ## penalty
    pen <- pmatch(penalty, c("mcp", "scad"))
    if (is.na(pen)) stop("Penalty need to be either 'mcp' or 'scad'!\n")
    ## penalty and computational approach
    app <- pmatch(approach, c("mmcd", "adaptive", "llacda"))
    if ((pen == 1) & (app == 1) & ((maxkappa >= 0.25) | (maxkappa < 0))) {
        stop("Using MMCD algorithm for MCP penalty, the regulation parameter kappa should be in [0, 0.25)!\n")
    }
    if ((pen == 1) & (app == 2) & ((maxkappa >= 1.0) | (maxkappa < 0))) {
        stop("Using adaptive rescaling algorithm for MCP penalty, the regulation parameter kappa should be in [0, 1.0)!\n")
    }
    if ((pen == 1) & (app == 2) & (maxkappa > 0.25)) {
        warning("Using adaptive rescaling algorithm for MCP penalty, the algorithm may not converge for large kappa!\n")
    }
    if ((pen == 1) & (app == 3) & ((maxkappa >= 1.0) | (maxkappa < 0)))  {
        stop("Using LLA-CDA algorithm for MCP penalty, the regulation parameter kappa should be in [0, 1.0)!\n")
    }
    if ((pen == 2) & ((maxkappa >= 0.2) | (maxkappa < 0))) {
        stop("Using MMCD algorithm for SCAD penalty, the regulation parameter kappa should be in [0, 0.2)!\n")
    }
    ## solution surface
    ss <- pmatch(path, c("kappa", "lambda", "hybrid"))
    if (is.na(ss)) stop("Path need to be: kappa,  lambda,  or hybrid!")
    ## assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1, n)
    qp <- qq+p
    cvk <- cv
    set.seed(seed)
    ## cross validation index
    if (stratified == FALSE) {
        ## cross validation without stratification
        cvpool <- rep(rep(1:cvk), length = n)
        nindex <- sample(cvpool, replace = FALSE)
        oy <- y
        ox <- x
    } else{
        ## strafified cross validation
        n0 <- length(y[y==0])
        n1 <- n-n0
        cvpool0 <- rep(rep(1:cvk), length = n0)
        cvpool1 <- rep(rep(1:cvk), length = n1)
        nidx0 <- sample(cvpool0, replace = FALSE)
        nidx1 <- sample(cvpool1, replace = FALSE)
        nindex <- c(nidx0, nidx1)
        ## reorder the data
        odridx <- order(y)
        oy <- y[odridx]
        ox <- x[odridx,]
    }
    ## create output space
    oout <- rep(0, 3+qp)
    opauc <- rep(0, nkappa*nlambda)
    olmdas <- rep(0, nkappa*nlambda)
    okas <- rep(0, nkappa*nlambda)
    ocoef <- matrix(0, qp, nkappa*nlambda)
    oaic <- rep(0, nkappa*nlambda)
    obic <- rep(0, nkappa*nlambda)
    oobj <- rep(0, nkappa*nlambda)
    odf <- rep(0, nkappa*nlambda)
    ocvx <- rep(0, nkappa*nlambda)
    ofull <- rep(0, nkappa*nlambda)
    cvcvx <- rep(0, nkappa*nlambda)
    cvfull <- rep(0, nkappa*nlambda)
    ## fit through Fortran depend on ss
    if (pen == 1){
        if (app == 1){
            if (ss == 1) {
                out <- try(.Fortran("cvauckapa", as.double(oout), as.double(opauc),
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                    as.integer(cvcvx), as.integer(cvfull),
                                    as.integer(nindex), as.integer(cvk),
                                    as.double(oy), as.double(int), as.double(ox),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 2) {
                out <- try(.Fortran("cvauclmda", as.double(oout), as.double(opauc),
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                    as.integer(cvcvx), as.integer(cvfull),
                                    as.integer(nindex), as.integer(cvk),
                                    as.double(oy), as.double(int), as.double(ox),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 3){
                out <- try(.Fortran("cvauckapa2", as.double(oout), as.double(opauc),
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                    as.integer(cvcvx), as.integer(cvfull),
                                    as.integer(nindex), as.integer(cvk),
                                    as.double(oy), as.double(int), as.double(ox),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            }
        } else if (app == 2){
            if (ss == 1) {
                out <- try(.Fortran("adpcvauckp", as.double(oout), as.double(opauc),
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                    as.integer(cvcvx), as.integer(cvfull),
                                    as.integer(nindex), as.integer(cvk),
                                    as.double(oy), as.double(int), as.double(ox),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 2) {
                out <- try(.Fortran("adpcvauclm", as.double(oout), as.double(opauc),
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                    as.integer(cvcvx), as.integer(cvfull),
                                    as.integer(nindex), as.integer(cvk),
                                    as.double(oy), as.double(int), as.double(ox),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            } else if (ss == 3){
                out <- try(.Fortran("adpcvauckp2", as.double(oout), as.double(opauc),
                                    as.double(olmdas), as.double(okas), as.double(ocoef),
                                    as.double(oaic), as.double(obic), as.double(oobj),
                                    as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                    as.integer(cvcvx), as.integer(cvfull),
                                    as.integer(nindex), as.integer(cvk),
                                    as.double(oy), as.double(int), as.double(ox),
                                    as.integer(n), as.integer(qq), as.integer(p),
                                    as.integer(nkappa), as.double(maxkappa),
                                    as.integer(nlambda), as.double(minlambda),
                                    as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
            }
        } else if (app ==3){
            out <- try(.Fortran("fcvllabi", as.double(oout), as.double(opauc),
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                as.integer(cvcvx), as.integer(cvfull),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy), as.double(int), as.double(ox),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        }
    } else {
        if (ss == 1) {
            out <- try(.Fortran("cvaucsdka", as.double(oout), as.double(opauc),
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                as.integer(cvcvx), as.integer(cvfull),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy), as.double(int), as.double(ox),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (ss == 2) {
            out <- try(.Fortran("cvaucsdlm", as.double(oout), as.double(opauc),
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                as.integer(cvcvx), as.integer(cvfull),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy), as.double(int), as.double(ox),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        } else if (ss == 3){
            out <- try(.Fortran("cvaucsdka2", as.double(oout), as.double(opauc),
                                as.double(olmdas), as.double(okas), as.double(ocoef),
                                as.double(oaic), as.double(obic), as.double(oobj),
                                as.integer(odf), as.integer(ocvx), as.integer(ofull),
                                as.integer(cvcvx), as.integer(cvfull),
                                as.integer(nindex), as.integer(cvk),
                                as.double(oy), as.double(int), as.double(ox),
                                as.integer(n), as.integer(qq), as.integer(p),
                                as.integer(nkappa), as.double(maxkappa),
                                as.integer(nlambda), as.double(minlambda),
                                as.double(epsilon),as.integer(maxit), PACKAGE = "cvplogistic"), TRUE)
        }
    }
    ## selection based on CV-AUC,  df
    if (!(inherits(out, 'try-error'))){
        cvauc <- out[[2]]
        lambda <- out[[3]]
        kappa <- out[[4]]
        coef <- matrix(out[[5]], qp, nkappa*nlambda)
        df <- out[[9]]
        fullmodel <- out[[11]]
        cvfullmodel <- out[[13]]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept", paste("x", 1:p, sep = ""))
        } else {
            rownames(coef) <- c("intercept", colnames(x))
        }
        ## regular solution if n > df
        uidx <- n > df
        ucvauc <- cvauc[uidx]
        ulambda <- lambda[uidx]
        ukappa <- kappa[uidx]
        ucoef <- coef[, uidx]
        ## minmum cvauc
        aucidx <- ucvauc ==  max(ucvauc)
        scvauc <- ucvauc[aucidx][[1]]
        slambda <- ulambda[aucidx][[1]]
        skappa <- ukappa[aucidx][[1]]
        scoef <- ucoef[, aucidx]
        if (is.matrix(scoef)) scoef <- scoef[, 1]
        list(scvauc, slambda, skappa, scoef)
    } else {
         stop("Model fitting fails, double check!\n")
         NULL
    }
}

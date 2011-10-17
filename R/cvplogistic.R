## Compute solution surface for a concave penalized logistic regression model
## using majorization minimization by coordinate descent (MMCD) algorithm.
## Two concave penalties are considered: SCAD and MCP.
## Three types of solution surfaces are provided for each penalty,
## 1. solution surface computed along kappa
## 2. solution surface computed along lambda
## 3. solution surface using Lasso-MCP hybrid algorithm
## Oct 17, 2011


cvplogistic <- function(y,x,penalty="mcp",path="kappa",
                        nkappa=20,maxkappa=0.249,nlambda=100,
                        minlambda=ifelse(n>p,0.0001,0.01),
                        epsilon=1e-3,maxit=1e+4){
    ## error checking
    if (nrow(x)!=length(y)) stop("Rows of X does not match length of Y!")
    if ((penalty=="mcp") & (maxkappa>=0.25)) stop("Regulation paramter kappa for MCP is in [0,0.25)!")
    if ((penalty=="scad") & (maxkappa>=0.2)) stop("Regulation paramter kappa for SCAD is in [0,0.2)!")
    ## penalty
    pen <- match(penalty,c("mcp","scad"))
    if (is.na(pen)) stop("Penalty need to be: mcp or scad!")
    ## solution surface
    ss <- match(path,c("kappa","lambda","hybrid"))
    if (is.na(ss)) stop("Path need to be: kappa, lambda, or hybrid!")
    ## space assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1,n)
    qp <- qq+p
    ## create output space
    olmdas <- rep(0,nkappa*nlambda)
    okas <- rep(0,nkappa*nlambda)
    ocoef <- matrix(0,qp,nkappa*nlambda)
    oaic <- rep(0,nkappa*nlambda)
    obic <- rep(0,nkappa*nlambda)
    oobj <- rep(0,nkappa*nlambda)
    odf <- rep(0,nkappa*nlambda)
    ocvx <- rep(0,nkappa*nlambda)
    ## fit through Fortran depend on pen and ss
    if (pen==1){
        if (ss==1) {
            out <- try(.Fortran("mcpkapa",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("mcplmda",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("mcpkapa2",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else {
        if (ss==1) {
            out <- try(.Fortran("scadkapa",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("scadlmda",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("scadkapa2",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    }
    if (!(inherits(out,'try-error'))){
        lambda <- out[[1]]
        kappa <- out[[2]]
        coef <- matrix(out[[3]],qp,nkappa*nlambda)
        df <- out[[7]]
        ## ocvx <- out[[8]]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept",paste("x",1:p,sep=""))
        } else {
            rownames(coef) <- c("intercept",colnames(x))
        }
        coef.intercept <- coef[1,]
        coef.covariates <- coef[-1,]
        return(list(lambda,kappa,df,coef.intercept,coef.covariates))
    } else {
        stop("Model fitting fails, double check!\n")
        return(NULL)
    }
}




## Tuning parameter selection using AIC
aic.cvplogistic <- function(y,x,penalty="mcp",path="kappa",
                            nkappa=20,maxkappa=0.249,nlambda=100,
                            minlambda=ifelse(n>p,0.0001,0.01),
                            epsilon=1e-3,maxit=1e+4){
    ## error checking
    if (nrow(x)!=length(y)) stop("Rows of X does not match length of Y!")
    if ((penalty=="mcp") & (maxkappa>=0.25)) stop("Regulation paramter kappa for MCP is in [0,0.25)!")
    if ((penalty=="scad") & (maxkappa>=0.2)) stop("Regulation paramter kappa for SCAD is in [0,0.2)!")
    ## penalty
    pen <- match(penalty,c("mcp","scad"))
    if (is.na(pen)) stop("Penalty need to be: mcp or scad!")
    ## solution surface
    ss <- match(path,c("kappa","lambda","hybrid"))
    if (is.na(ss)) stop("Path need to be: kappa, lambda, or hybrid!")
    ## assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1,n)
    qp <- qq+p
    ## create output space
    olmdas <- rep(0,nkappa*nlambda)
    okas <- rep(0,nkappa*nlambda)
    ocoef <- matrix(0,qp,nkappa*nlambda)
    oaic <- rep(0,nkappa*nlambda)
    obic <- rep(0,nkappa*nlambda)
    oobj <- rep(0,nkappa*nlambda)
    odf <- rep(0,nkappa*nlambda)
    ocvx <- rep(0,nkappa*nlambda)
    ## fit through Fortran depend on pen and ss
    if (pen==1){
        if (ss==1) {
            out <- try(.Fortran("mcpkapa",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("mcplmda",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("mcpkapa2",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else {
        if (ss==1) {
            out <- try(.Fortran("scadkapa",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("scadlmda",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("scadkapa2",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    }
    if (!(inherits(out,'try-error'))){
        lambda <- out[[1]]
        kappa <- out[[2]]
        coef <- matrix(out[[3]],qp,nkappa*nlambda)
        df <- out[[7]]
        aic <- out[[4]]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept",paste("x",1:p,sep=""))
        } else {
            rownames(coef) <- c("intercept",colnames(x))
        }
        ## selection based on AIC, df
        dfidx <- df<=n
        aicidx <- aic==min(aic)
        tuning.AIC <- aic[(aicidx&dfidx)]
        tuning.lambda <- lambda[(aicidx&dfidx)]
        tuning.kappa <- kappa[(aicidx&dfidx)]
        tuning.intercept <- coef[1,(aicidx&dfidx)]
        tuning.covariates <- coef[-1,(aicidx&dfidx)]
        return(list(tuning.AIC,tuning.lambda,tuning.kappa,
                    tuning.intercept,tuning.covariates))
    } else {
        stop("Model fitting fails, double check!\n")
        return(NULL)
    }
}



## Tuning parameter selection using BIC
bic.cvplogistic <- function(y,x,penalty="mcp",path="kappa",
                            nkappa=20,maxkappa=0.249,nlambda=100,
                            minlambda=ifelse(n>p,0.0001,0.01),
                            epsilon=1e-3,maxit=1e+4){
    ## error checking
    if (nrow(x)!=length(y)) stop("Rows of X does not match length of Y!")
    if ((penalty=="mcp") & (maxkappa>=0.25)) stop("Regulation paramter kappa for MCP is in [0,0.25)!")
    if ((penalty=="scad") & (maxkappa>=0.2)) stop("Regulation paramter kappa for SCAD is in [0,0.2)!")
    ## penalty
    pen <- match(penalty,c("mcp","scad"))
    if (is.na(pen)) stop("Penalty need to be: mcp or scad!")
    ## solution surface
    ss <- match(path,c("kappa","lambda","hybrid"))
    if (is.na(ss)) stop("Path need to be: kappa, lambda, or hybrid!")
    ## assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1,n)
    qp <- qq+p
    ## create output space
    olmdas <- rep(0,nkappa*nlambda)
    okas <- rep(0,nkappa*nlambda)
    ocoef <- matrix(0,qp,nkappa*nlambda)
    oaic <- rep(0,nkappa*nlambda)
    obic <- rep(0,nkappa*nlambda)
    oobj <- rep(0,nkappa*nlambda)
    odf <- rep(0,nkappa*nlambda)
    ocvx <- rep(0,nkappa*nlambda)
    ## fit through Fortran depend on pen and ss
    if (pen==1){
        if (ss==1) {
            out <- try(.Fortran("mcpkapa",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("mcplmda",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("mcpkapa2",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else {
        if (ss==1) {
            out <- try(.Fortran("scadkapa",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("scadlmda",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("scadkapa2",
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    }
    if (!(inherits(out,'try-error'))){
        lambda <- out[[1]]
        kappa <- out[[2]]
        coef <- matrix(out[[3]],qp,nkappa*nlambda)
        df <- out[[7]]
        aic <- out[[5]]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept",paste("x",1:p,sep=""))
        } else {
            rownames(coef) <- c("intercept",colnames(x))
        }
        ## selection based on AIC, df
        dfidx <- df<=n
        aicidx <- aic==min(aic)
        tuning.AIC <- aic[(aicidx&dfidx)]
        tuning.lambda <- lambda[(aicidx&dfidx)]
        tuning.kappa <- kappa[(aicidx&dfidx)]
        tuning.intercept <- coef[1,(aicidx&dfidx)]
        tuning.covariates <- coef[-1,(aicidx&dfidx)]
        return(list(tuning.AIC,tuning.lambda,tuning.kappa,
                    tuning.intercept,tuning.covariates))
    } else {
        stop("Model fitting fails, double check!\n")
        return(NULL)
    }
}






## Tuning parameter selection using CV-AUC
auc.cvplogistic <- function(cv=5,y,x,penalty="mcp",path="kappa",
                            nkappa=20,maxkappa=0.249,nlambda=100,
                            minlambda=ifelse(n>p,0.0001,0.01),
                            epsilon=1e-3,maxit=1e+4,seed=1000){
    ## error checking
    if (nrow(x)!=length(y)) stop("Rows of X does not match length of Y!")
    if ((penalty=="mcp") & (maxkappa>=0.25)) stop("Regulation paramter kappa for MCP is in [0,0.25)!")
    if ((penalty=="scad") & (maxkappa>=0.2)) stop("Regulation paramter kappa for SCAD is in [0,0.2)!")
    ## penalty
    pen <- match(penalty,c("mcp","scad"))
    if (is.na(pen)) stop("Penalty need to be: mcp or scad!")
    ## solution surface
    ss <- match(path,c("kappa","lambda","hybrid"))
    if (is.na(ss)) stop("Path need to be: kappa, lambda, or hybrid!")
    ## assignment for FORTRAN
    dimx <- dim(x)
    n <- dimx[1]
    p <- dimx[2]
    qq <- 1
    int <- rep(1,n)
    qp <- qq+p
    cvk <- cv
    cvpool <- rep(rep(1:cvk),length=n)
    set.seed(seed)
    nindex <- sample(cvpool,replace=FALSE)
    ## create output space
    oout <- rep(0,3+qp)
    opauc <- rep(0,nkappa*nlambda)
    olmdas <- rep(0,nkappa*nlambda)
    okas <- rep(0,nkappa*nlambda)
    ocoef <- matrix(0,qp,nkappa*nlambda)
    oaic <- rep(0,nkappa*nlambda)
    obic <- rep(0,nkappa*nlambda)
    oobj <- rep(0,nkappa*nlambda)
    odf <- rep(0,nkappa*nlambda)
    ocvx <- rep(0,nkappa*nlambda)
    ofull <- rep(0,nkappa*nlambda)
    cvcvx <- rep(0,nkappa*nlambda)
    cvfull <- rep(0,nkappa*nlambda)
    ## fit through Fortran depend on ss
    if (pen==1){
        if (ss==1) {
            out <- try(.Fortran("cvauckapa",as.double(oout),as.double(opauc),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),as.integer(ofull),
                                as.integer(cvcvx),as.integer(cvfull),
                                as.integer(nindex),as.integer(cvk),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("cvauclmda",as.double(oout),as.double(opauc),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),as.integer(ofull),
                                as.integer(cvcvx),as.integer(cvfull),
                                as.integer(nindex),as.integer(cvk),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("cvauckapa2",as.double(oout),as.double(opauc),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),as.integer(ofull),
                                as.integer(cvcvx),as.integer(cvfull),
                                as.integer(nindex),as.integer(cvk),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }
    } else {
        if (ss==1) {
            out <- try(.Fortran("cvaucsdka",as.double(oout),as.double(opauc),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),as.integer(ofull),
                                as.integer(cvcvx),as.integer(cvfull),
                                as.integer(nindex),as.integer(cvk),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==2) {
            out <- try(.Fortran("cvaucsdlm",as.double(oout),as.double(opauc),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),as.integer(ofull),
                                as.integer(cvcvx),as.integer(cvfull),
                                as.integer(nindex),as.integer(cvk),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        } else if (ss==3){
            out <- try(.Fortran("cvaucsdka2",as.double(oout),as.double(opauc),
                                as.double(olmdas),as.double(okas),as.double(ocoef),
                                as.double(oaic),as.double(obic),as.double(oobj),
                                as.integer(odf),as.integer(ocvx),as.integer(ofull),
                                as.integer(cvcvx),as.integer(cvfull),
                                as.integer(nindex),as.integer(cvk),
                                as.double(y),as.double(int),as.double(x),
                                as.integer(n),as.integer(qq),as.integer(p),
                                as.integer(nkappa),as.double(maxkappa),
                                as.integer(nlambda),as.double(minlambda),
                                as.double(epsilon),as.integer(maxit)),TRUE)
        }


    }
    ## selection based on CV-AUC, df
    if (!(inherits(out,'try-error'))){
        cvauc <- out[[2]]
        lambda <- out[[3]]
        kappa <- out[[4]]
        coef <- matrix(out[[5]],qp,nkappa*nlambda)
        fullmodel <- out[[11]]
        cvfullmodel <- out[[13]]
        if (is.null(colnames(x))){
            rownames(coef) <- c("intercept",paste("x",1:p,sep=""))
        } else {
            rownames(coef) <- c("intercept",colnames(x))
        }
        ## selection based on CV-AUC, df
        dfidx <- fullmodel!=1
        cvdfidx <- cvfullmodel!=1
        aucidx <- cvauc==min(cvauc)
        tuning.CVAUC <- cvauc[(aucidx&dfidx&cvdfidx)]
        tuning.lambda <- lambda[(aucidx&dfidx&cvdfidx)]
        tuning.kappa <- kappa[(aucidx&dfidx&cvdfidx)]
        tuning.intercept <- coef[1,(aucidx&dfidx&cvdfidx)]
        tuning.covariates <- coef[-1,(aucidx&dfidx&cvdfidx)]
        return(list(tuning.CVAUC,tuning.lambda,tuning.kappa,
                    tuning.intercept,tuning.covariates))
    } else {
         stop("Model fitting fails, double check!\n")
         return(NULL)
    }
}

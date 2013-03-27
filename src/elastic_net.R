# this is to implement the elastic net with pathwise coordinate descent optimization
# probably-to-do: 
#   covariance updates
#   sparse updates
#   weighted updates

library(glmnet)
library(MASS)
library(bigmatrix)

S <- function(z, gamma) {
    if(abs(z) <= gamma) return(0)
    else {
        if(z > 0) return(z - gamma)
        else return(z + gamma)
    }
}

obj <- function(residual, beta, nobs, lambda, alpha)
{
    mean.square.loss <- (1 / (2 * nobs)) * t(residual) %*% residual
    regularization <- lambda * (((1 - alpha) / 2) * t(beta) %*% beta + alpha * sum(abs(beta)))
    
    obj.value <- mean.square.loss + regularization
    
    return(obj.value)
}

elanet <- function(x, y, alpha = 0.2, thresh = 1e-7, lambda = NULL, max.iter = 100, standardize = TRUE)
{
    nobs <- nrow(x)
    nvar <- ncol(x)
    
    if(standardize) {
        x.mean <- colMeans(x)
        x <- x - matrix(rep(x.mean, nobs), nrow=nobs, byrow=TRUE)
        x.sd <- apply(x, 2, sd)
        x <- x / x.sd
    } 

    y.mean <- mean(y)
    y <- scale(y, center=TRUE, scale=FALSE)
    
    if(is.null(lambda)) {
        # calculate the path of lambda
        if(alpha != 0) {
            nlambda <- 100
            lambda.max <- max(abs(t(x) %*% y)) / (alpha * nobs)
            lambda.min <- 0.001 * lambda.max
            lambda <- exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda))
        } else {
            nlambda <- 100
            lambda <- seq(0.2, 0.001, length.out = nlambda)
        }
    } else {
        nlambda <- length(lambda)
        lambda <- rev(sort(lambda))
    }
    
    # beta <- rep(0, nvar)
    beta <- matrix(rep(0, nvar * nlambda), nvar, nlambda)
    if(is.null(rownames(beta)))
        rownames(beta) <- paste("V", seq(nvar), sep = "")
    if(is.null(colnames(beta)))
        colnames(beta) <- paste("lambda", seq(nlambda), sep = "")
    dev <- rep(0, nlambda)
    npasses <- rep(0, nlambda)
    
    r <- y
    nulldev <- sum(r^2) # define the null.deviance as the sum square of initial residual.
    thresh <- thresh * nulldev
    
    for(k in 1:nlambda) {
        max.update <- 0
        
        if(k > 1) {
            beta[,k] <- beta[,k-1] # a warm start, which reflects the term "pathwise"
        }
        
        i <- 1
        while(i <= max.iter) {
            for(j in 1:nvar) {
                beta.prev <- beta[,k]
                beta[j,k] <- S((t(x[,j]) %*% r) / nobs + beta[j,k], lambda[k]*alpha) / (1 + lambda[k]*(1-alpha))
            
                delta.r <- x[,j] * (beta[j,k] - beta.prev[j])
                r.new <- r - delta.r
            
                # change of the L2 norm of residual
                obj.update <- sum(delta.r^2)
                max.update <- ifelse(abs(obj.update) >= max.update, abs(obj.update), max.update)
            
                r <- r.new
            }
            if(max.update < thresh) break
        
            i <- i + 1
            max.update <- 0
        }
        
        npasses[k] <- i
        dev[k] <- sum(r^2)
        # r <- y
    }  
    
    if(standardize) {
        scale.mean <- x.mean
        scale.sd <- x.sd
    } else {
        scale.mean <- 0
        scale.sd <- 1
    }

    result <- list(beta = beta, standardize = standardize, scale.mean = scale.mean, scale.sd = scale.sd, intercept = y.mean, nulldev = nulldev, dev = dev, npasses = npasses, lambda = lambda)
    return(result)
}

dbg.info <- function(object)
{
    cat("nulldev: ", object$nulldev, "\n")
    cat("dev: ", object$dev, "\n")
    cat("lambda: ", object$lambda, "\n")
    print(as.matrix(object$beta))
    cat("npasses: ", object$npasses, "\n")
}

elanet.predict <- function(model, x, y)
{
    if(model$standardize) {
        x <- x - matrix(rep(model$scale.mean, nrow(x)), nrow = nrow(x), byrow=TRUE)
        x <- x / model$scale.sd
    }
    
    response <- x %*% model$beta + model$intercept
    beta.size <- ncol(model$beta)
    residual <- matrix(rep(y, beta.size), ncol=beta.size) - response
    dev <- colSums(residual^2)
    rmse <- sqrt(dev / nrow(x))
    
    return(list(response = response, rmse = rmse, dev = dev))
}

main <- function(n, d, alpha, lambda, thresh)
{
    data.path <- file.path(".", paste("n", n, "_d", d, ".RData", sep = ""))
    if(file.exists(data.path))
        load(file = data.path)
    else {
        D = tiger.generator(n=n, d=d, graph="scale-free", prob=0.03)
        x <- as.matrix(D$data[,2:d])
        y <- as.matrix(D$data[,1])
    }
    
    beg.time <- Sys.time()
    cat("my elastic net >>> \n")
    cv.num <- 5
    rmse <- 0
    beta <- 0
    for(iter in 1:cv.num) {
        cat("iter ", iter, "\n")
        ntr <- sample(1:n, 0.5 * n)
        fit1 <- elanet(x[ntr,], y[ntr], alpha = alpha, lambda = lambda)
        beta <- beta + fit1$beta
        dbg.info(fit1)
        rmse <- rmse + elanet.predict(fit1, x[-ntr,], y[-ntr])$rmse
    }
    rmse <- rmse / cv.num
    beta <- beta / cv.num
    print(rmse)
    save(file = data.path, x, y, n, d)
    plot(x = 1:length(rmse), y = rmse, xlab = "step", ylab = "rmes", ylim = c(min(rmse), max(rmse)), type="b")
    #plot.path(beta)
    
    #cat("-----------------------------------\n")
    #cat("glmnet >>> \n")
    #fit2 <- glmnet(x, y, alpha = alpha, lambda = lambda)
    #fit2$dev <- fit2$nulldev * (1 - fit2$dev.ratio)
    #glmnet.time <- Sys.time()
    #dbg.info(fit2)
    #print(difftime(glmnet.time, elanet.time, unit = "sec"))
    
}

plot.path <- function(beta)
{
    y.upper <- max(beta)
    y.lower <- min(beta)
    x.coords <- 1:ncol(beta)
    plot(x.coords, beta[1,], type="o", ylim=c(y.lower-0.1, y.upper+0.1), pch=16, lwd=2)
    for(i in 2:nrow(beta)) {
        lines(x.coords, beta[i,], type="o", lwd=2, pch=16)
    }
    lines(x.coords, rep(0, ncol(beta)), col="red")
}

#main(n = 200, d = 100, alpha = 0, lambda = seq(0.2, 0.01, length.out=10))
main(n = 200, d = 50, alpha = 1, lambda = NULL)
#main(n = 100, d = 2, alpha = 1, lambda = seq(0.1, 0.001, length.out = 100))


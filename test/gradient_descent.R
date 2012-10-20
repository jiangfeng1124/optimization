source("src/gradient_descent.R")

## example: 1-dimensional
f1 <- function(x) {
	return((x[1] - 5.0)^2)
}

g1 <- function(x) {
	return(2 * (x[1] - 5.0))
}

initial.x <- array(c(0.0))

results <- GradientDescent(f1, g1, initial.x, show.trace = FALSE)
#print(paste("minimum.x: ", results$optim.x))
#print(paste("minimum.f: ", results$optim.f))
print(results)

## example: 2-dimensional
eta = 0.9

f2 <- function(x) {
	return((1.0 / 2.0) * (x[1]^2 + eta * x[2]^2))
}

g2 <- function(x) {
	return(c(x[1], eta * x[2]))
}

initial.x <- array(c(1.0, 1.0))
results <- GradientDescent(f2, g2, initial.x, show.trace = FALSE)
print(results)

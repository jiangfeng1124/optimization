source("src/newton.R")

cat(">>> example: 1-dimensional\n")
f <- function(x)
{
	return((x[1] - 5.0)^4)
}

g <- function(x)
{
	return(array(c(4.0 * (x[1] - 5.0)^3)))
}

h <- function(x)
{
	return(array(c(12.0 * (x[1] - 5.0)^2)))
}

initial.x <- array(c(0.0))
results <- Newton(f, g, h, initial.x, show.trace = TRUE)
print(results)

cat(">>> example: 2-dimensional\n")
eta = 0.9

f2 <- function(x)
{
	return((1.0 / 2.0) * (x[1]^2 + eta * x[2]^2))
}

g2 <- function(x)
{
	return(array(c(x[1], eta * x[2])))
}

h2 <- function(x)
{
	return(matrix(c(1.0, 0.0, 0.0, eta), nrow=2, ncol=2))
}

initial.x <- array(c(127.0, 921.0))
results <- Newton(f2, g2, h2, initial.x, show.trace = TRUE)
print(results)

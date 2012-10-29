source("src/backtracking_line_search.R")

Newton <- function(f, g, h, initial.x, tolerance = 10e-8, max.iterations = 1000, show.trace = TRUE)
{
	x <- initial.x

	dx <- -solve(h(x)) %*% g(x)
	l2 <- (t(g(x)) %*% solve(h(x)) %*% g(x))

	i <- 0

	converged <- FALSE

	if(show.trace) {
		print(paste("Iteration: ", i))
		print(paste("x: ", x))
		print(paste("f(x)", f(x)))
		print(paste("g(x)", g(x)))
		print(paste("h(x)", h(x)))
	}

	while(!converged && i < max.iterations) {
		i <- i + 1

		# matrix should be cast into numerical vector
		dx <- -as.numeric(solve(h(x)) %*% g(x))

		step.size <- BacktrackingLineSearch(f, g, x, dx)

		x <- x + step.size * dx

		l2 <- (t(g(x)) %*% solve(h(x)) %*% g(x))
		if(l2/2 <= tolerance) {
			converged = TRUE
		}

		if(show.trace) {
			print(paste("Iteration: ", i))
			print(paste("x: ", x))
			print(paste("f(x)", f(x)))
			print(paste("g(x)", g(x)))
			print(paste("h(x)", h(x)))
		}
	}

	optimization.result <- list(initial.x = initial.x, optim.x = x, optim.f = f(x), num.iter = i, is.converged = converged)
	return(optimization.result)
}


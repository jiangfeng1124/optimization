GradientDescent <- function(f, g, initial.x, tolerance = 10e-8, max.iterations = 1000, show.trace = TRUE)
{
	x <- initial.x

	i = 0
	if(show.trace) {
		print(paste("Iteration: ", i))
		print(paste("x: ", x))
		print(paste("g(x): ", g(x)))
		print(paste("||g(x)||: ", norm(as.matrix(g(x)), type = "F")))
	}

	converged = FALSE

	while(!converged && i < max.iterations) {
		i = i + 1

		step.size <- BacktrackingLineSearch(f, g, x, -g(x))

		x <- x - step.size * g(x)

		if (norm(as.matrix(g(x)), type = "F") < tolerance) {
			converged = TRUE
		}

		if(show.trace) {
			print(paste("Iteration: ", i))
			print(paste("x: ", x))
			print(paste("g(x): ", g(x)))
			print(paste("||g(x)||: ", norm(as.matrix(g(x)), type = "F")))
			print(paste("step-size: ", step.size))
		}
	}

	optimization.result <- list(initial.x = initial.x, optim.x = x, optim.f = f(x), num.iter = i, is.converged = converged)
	return(optimization.result)
}

BacktrackingLineSearch <- function(f, g, x, dx, alpha = 0.1, beta = 0.8, max.iterations = 1000) {
	i <- 0
	
	f.x <- f(x)
	g.x <- g(x)
	angle <- t(g.x) %*% dx
	
	t <- 1
	while(f(x + t * dx) > f.x + alpha * t * angle) {
		t <- beta * t
		i <- i + 1
		
		if(i > max.iterations) {
			print(paste("too many iterations in BacktrackingLineSearch(alpha: ", alpha, " beta: ", beta))
		}
	}

	return(t)
}


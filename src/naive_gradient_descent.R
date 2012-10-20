NaiveGradientDescent <- function(f, g, initial.x, step.size = 0.1, tolrence = 10e-8, max.iterations = 1000, show.trace = TRUE)
{
	x.old <- initial.x
	x.new <- initial.x
	y.old <- Inf
	y.new <- f(x.new)

	i = 0

	converged = FALSE

	if(show.trace) {
		print(paste("Iteration: ", i))
		print(paste("x.new: ", x.new))
		print(paste("f(x.new): ", f(x.new)))
		print(paste("g(x.new): ", g(x.new)))
	}

	while(!converged && i < max.iterations) {
		i = i + 1

		x.old <- x.new
		x.new <- x.new - step.size * g(x.new)

		y.old <- y.new
		y.new <- f(x.new)

		if(abs(y.new - y.old) <= tolrence) {
			converged = TRUE
		}

		if(show.trace) {
			print(paste("Iteration: ", i))
			print(paste("x.new: ", x.new))
			print(paste("f(x.new): ", f(x.new)))
			print(paste("g(x.new): ", g(x.new)))
		}
	}

	optimization.result <- list(initial.x = initial.x, optim.x = x.new, optim.f = y.new, num.iter = i, is.converged = converged)
	return(optimization.result)
}


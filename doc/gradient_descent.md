## GradientDescent(...)
* Arguments:
	* f: A function that is differentiable
	* g: the gradient of f
	* initial.x: A value in the domain of f, from which to start search for a minimum
	* tolerance: How small the norm of the gradient should be for convergence to be declared
	* max.iterations: The maximum number of iterations
	* show.trace: whether or not to show tracing infomation while running
* Returns:
	* initial.x: Same as it is in the Arguments
	* optim.x: The purported minimum of function to be optimized
	* optim.f: The purported minimal value of function to be optimized
	* num.iter: The number of iterations
	* is.converged: Whether or not the algorithm is converged


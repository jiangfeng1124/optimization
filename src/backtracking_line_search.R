BacktrackingLineSearch <- function(f, g, x, dx, c = 0.1, rho = 0.8, max.iterations = 1000) 
{
	i <- 0
	
	f.x <- f(x)
	g.x <- g(x)
	angle <- t(g.x) %*% dx
	
	alpha <- 1
	while(f(x + alpha * dx) > f.x + c * alpha * angle) {
		alpha <- rho * alpha
		i <- i + 1
		
		if(i > max.iterations) {
			print(paste("too many iterations in BacktrackingLineSearch(c: ", c, " rho: ", rho))
		}
	}

	return(alpha)
}


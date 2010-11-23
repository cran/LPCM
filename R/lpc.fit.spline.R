
lpc.fit.spline <- function(lpcsl,  num.knots=100, add.spline=FALSE, lty=1, col=3){
    nbranch <- length(lpcsl)-1   
    # First of all evaluate the spline at each knot
    knots.pi <- list()
    knots.coords <- list()
    for (r in 0:nbranch) {
      #knots.pi[[r+1]] <- lpcsl[[r+1]]$range[1] + 0:(num.knots-1)/num.knots * (lpcsl[[r+1]]$range[2]-lpcsl[[r+1]]$range[1])
      #print(lpcsl[[r+1]]$range)
      knots.pi[[r+1]] <- round(lpcsl[[r+1]]$range[1] + 0:(num.knots-1)/(num.knots-1) * (lpcsl[[r+1]]$range[2]-lpcsl[[r+1]]$range[1]), digits=6)
      coords <- matrix(ncol=num.knots,nrow=length(lpcsl[[r+1]]$splinefun))
      for (j in 1:length(lpcsl[[r+1]]$splinefun))
        coords[j,] <- lpcsl[[r+1]]$splinefun[[j]](knots.pi[[r+1]])
      knots.coords[[r+1]] <- coords
    }
    if (add.spline==TRUE){ 
      for (r in 0:nbranch){
         lines(knots.coords[[r+1]][1,],knots.coords[[r+1]][2,] ,lty=lty, col=col)
       }
    }

    list(knots.coords=knots.coords, knots.pi=knots.pi)
  }

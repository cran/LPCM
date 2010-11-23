lpc <- function(X,h, t0=mean(h),  x0,   mult=1, way = "two", scaled=TRUE,  weights=1, pen=2,
              depth=1, control=lpc.control()){
  iter    <-  control$iter
  boundary <- control$boundary
  convergence.at <- control$convergence.at
  thresh  <- control$pruning.thresh
  rho0    <- control$rho0
  cross   <- control$cross
  Xi       <- as.matrix(X)
  N        <- dim(Xi)[1]
  d        <- dim(Xi)[2]
  if (N %% length(weights) !=0){
       stop("The length of the vector of weights is not a multiple of the sample size.")
  } else {
      weights <-rep(weights, N %/%   length(weights))
  }     
  X1       <- matrix(0,mult,d) # collects starting points for branches of depth 1.
  X2<-X3   <- matrix(0,0,d)    # collect starting points for branches of depth 2 and 3.
  S2<-S3   <- matrix(0,0,d)    # collects candidates for starting points of depth 2 and 3.
  saveall  <- matrix(0,0,d)    # stores the LPC
  countb   <- 0                # counter for branches 17/07/2007
  Lambda   <- matrix(0,0,5)    # stores the parameters (lambda, countb) as well as  depth, order, and  initialization number of the corresponding branch.
  s1       <- apply(Xi, 2, function(dat){ diff(range(dat))})  # range
  if (missing(h)){ if (!scaled){h   <- s1/10} else {h<- 0.1}} # bandwidth by default: 10% of range
  if (length(h)==1){h <- rep(h,d)}
  if (scaled){        # scales the data to lie by its range
        Xi <- sweep(Xi, 2, s1, "/")
   } 
  if (d==1){
      stop("Data set needs to consist of at least two variables!")
  } else if (d > 2 && depth > 1){
      cat("Dimension of data set enforces depth=1! \n")
      depth <- 1
    }
  if (!missing(x0)){
      x0 <- matrix(x0, ncol=d, byrow=TRUE)
      if (scaled){ x0 <- sweep(x0, 2, s1, "/") } # scales the starting point just as the scaled data
      if (dim(x0)[1] < mult) {n <- runif(mult-dim(x0)[1],1,N+1)%/%1; x0 <- rbind(x0,Xi[n,])}  # creates the appropriate number of random starting points
  } else {
       n <- runif(mult,1,N+1)%/%1;     
       x0 <- matrix(Xi[n,], length(n),d)# corrected 16/10          
  }

   if (depth>1){
       kde2x <-   function(X, x, h){
          l<-dim(X)[1]
          1/l*sum(kern(X[,1],x[1],h[1])*kern(X[,2],x[2],h[2]))
       }
   }
   
  for (j in 1:mult){
    
     xo           <- x0[j,]                       # defines appropraite starting point for the jth curve
     X1[j,]       <- xo                           # adds starting point to list
     curve0       <-  followx(Xi, xo, h, t0, iter, way, weights, pen, phi =1, 0,rho0,  boundary, convergence.at,  cross ) # computes LPC of depth 1
     saveall      <- rbind(saveall,curve0[[1]])	  # stores LPC
     l            <- dim(curve0[[5]])[1]	  # number of candidates for junctions
     Lambda       <- rbind(Lambda, cbind(curve0[[6]], countb,1,1,j))
     countb       <- countb + 1
     dimnames(Lambda)[[2]]<- c("lambda", "branch", "depth", "order", "init")

     if (depth >1 && l>=1 ){							  # constructs branches of depth 2.

       for (s in 1:l){

          x          <- curve0[[5]][s,]							  # candidates for junctions on initial curve.
          center.x   <- cov.wt(Xi, wt= kernd(Xi,x,h)*weights)           # computes local covariance and mean at x
          eigen.cov  <- eigen(center.x[[1]])  				# eigenvalues and eigenvec's of local cov matrix
          eigen.vecd <- eigen.cov[[2]][,2]				# second local eigenvector
          new.x      <- center.x[[2]]+ 2*t0 * eigen.vecd 		# double stepsize to escape from initial curve

          if (kde2x(Xi,new.x, h) >= thresh){ 				# Pruning
             X2     <- rbind(X2, t(new.x))				# adds new starting point to list
             x      <- new.x
              curve1 <- followx(Xi, x, h, t0, iter, way ="one", weights, pen,  phi=2, lasteigenvector= eigen.vecd, rho0, boundary,convergence.at,  cross)
             saveall<- rbind(saveall, curve1[[1]])    # stores LPC
             S2     <- rbind(S2,curve1[[5]])          # stores candidates for further junctions
             Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 2, 2,j))
             countb <- countb + 1
          }
          new.x <- center.x[[2]]- 2*t0 * eigen.vecd 		# go in opposite direction
          if (kde2x(Xi,new.x, h) >= thresh){				  # Pruning
             X2     <- rbind(X2, t(new.x))
             x      <- new.x
             curve1 <- followx(Xi, x, h, t0, iter, way="back", weights, pen, phi=2, lasteigenvector=-eigen.vecd, rho0,boundary, convergence.at,  cross)
             saveall<- rbind(saveall, curve1[[1]])
             S2     <- rbind(S2, curve1[[5]])
             Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 2,2,j))
             countb <- countb + 1     
          } # end if kde2x.....
      } # end for (s)
     }# end if (depth)

     k <- dim(S2)[1]							      # no. of candidates for starting points of depth 3.
     if (depth >2 && k>=1 ){					# constructs branches of depth 3
         for (s in 1:k){

            x         <- S2[s,]			  # candidates for junctions on curve of depth 2.
            center.x  <- cov.wt(Xi, wt= kernd(Xi,x,h)*weights)
            eigen.cov <- eigen(center.x[[1]])
            eigen.vecd<- eigen.cov[[2]][,2]
            new.x     <- center.x[[2]]+ 2*t0 * eigen.vecd
            if (kde2x(Xi,new.x, h)  >= thresh){
                  X3      <- rbind(X3, t(new.x))
                  x       <- new.x
                   curve1  <- followx(Xi, x , h, t0, iter, way ="one", weights, pen, phi=3, lasteigenvector= eigen.vecd, rho0, boundary, convergence.at, cross)
                  saveall <- rbind(saveall, curve1[[1]])
                  S3      <- rbind(S3,curve1[[5]])
                  Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 3,2,j))
                  countb <- countb + 1 
                }
            new.x <- center.x[[2]]- 2*t0* eigen.vecd
            if (kde2x(Xi,new.x, h)   >= thresh){
                  X3     <- rbind(X3, t(new.x))
                  x      <- new.x
                 curve1 <- followx(Xi, x, h, t0, iter, way="back", weights, pen, phi=3, lasteigenvector=-eigen.vecd, rho0, boundary, convergence.at, cross)
                  saveall<- rbind(saveall, curve1[[1]])
                  S3     <- rbind(S3,curve1[[5]])
                  Lambda <- rbind(Lambda, cbind(curve1[[6]], countb, 3,2,j))
                  countb <- countb + 1
            } # end if kde2x....
         } # end for (s)
     } # end if (depth)
  } # end for (j)
  if (d==2){
      starting.points <- rbind(X1,X2,X3) 
      dimnames(starting.points)<-list(c(rep(1, dim(X1)[1]), rep(2, dim(X2)[1]), rep(3, dim(X3)[1]) ), NULL)
      #dimnames(starting.points)[[1]]<-c(rep(1, dim(X1)[1]), rep(2, dim(X2)[1]), rep(3, dim(X3)[1]) )
  } else {
      starting.points <- as.matrix(X1)
      dimnames(starting.points)<-list(c(rep(1, dim(X1)[1])), NULL)
  }


  fit<-  c(LPC=list(saveall),
           Parametrization= list(Lambda),
           h=list(h),
           t0=t0,
           starting.points=list(starting.points),
           data=list(Xi),
           scaled=list(scaled),
           weights=list(weights),
           control=list(control),
           Misc=list(list( rho=curve0[[4]], scaled.by = if (scaled) s1 else rep(1,d), adaptive.band=curve0[[7]] ))
           ) 
   class(fit)<-"lpc"
   fit
  
} # end lpc function

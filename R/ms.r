meanshift<-function(X, x, h){
  
  if(is.null(dim(X))){X<-as.vector(X)}
  if(is.vector(X)){
    X<-matrix(X,nrow=length(X))
  }
  x <- as.numeric(x)
  d <- length(x)
  g <- kernd(X, x, h)
  
  ms <- NULL
  for(j in 1:d){
    ms[j]<-sum(X[,j]*g)/sum(g)
  }
  # Alternative: 
  # ms<- apply(X*g, 2, sum)/sum(g)
  # this is faster by itself, but it gets slower as soon as integrated into ms.rep....
  ms
}


# Mean shift iterative function (until convergence ...)
ms.rep <- function (X, x, h, thresh= 0.0001, iter=200) {    
  s  <- 0
  th <- rep(0,iter)
  M  <-matrix(0, iter, length(x))
  x0 <- x
  Xm<-colMeans(X)
  
  for (j in 1: iter){           
    m     <- meanshift(X, x, h)        
    M[j,] <- m
    th[j] <- enorm(m-x)/enorm(Xm-x)          
    if (th[j]<thresh){
      s<-j;                 
      break
    }
    x     <- m
  }
  return(list("Meanshift.points"=M[1:s,], "Threshold.values"=th[1:s], "iterations"=s, "start"=as.numeric(x0),  "final"=m))
}

# Mean shift clustering

ms<-function (X, h, subset, thr = 0.01, scaled = 1, iter=200, plot = TRUE,  ...) 
{
  if(is.null(dim(X))){X<-as.vector(X)}
  if(is.vector(X)){
    X<-matrix(X,nrow=length(X))
  }
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  
  if (missing(subset)) {
    subset <- 1:n
  }
  
  s1<-rep(1,d)
  
  if (scaled==1){
    s1       <- apply(X, 2, function(dat){ diff(range(dat))})  # range 
  } else if (scaled > 0){
    s1<- apply(X, 2, sd)
  } else if (scaled <0){
    stop("Negative values for scaled not allowed.")
  }
  if (missing(h)){
    if (!scaled){
      h<- apply(X, 2, function(dat){ diff(range(dat))})/20
    }  else if (scaled==1) { 
      h<- 0.05
    } else if (scaled >0) {
      h<-0.2
    }
  } # bandwidth by default: 5% of range or 20% of standard deviation 
  
  if (length(h)==1){h <- rep(h,d)}
  
  if (scaled >0){        # scales the data to lie by its range
    X <- sweep(X, 2, s1, "/")
  } 
  Xm<-colMeans(X)
  
  if (plot && d ==2) {
    plot(X, col = "grey70", ...)
  }
  
  finals <- matrix(0, n, d)
  ncluster <- 0
  savecluster <- matrix(0, 0, d)
  cluster.label <- closest.label <- rep(0, n)
  
  all.ms<-list()
  
  for (i in subset) {
    all.ms[[i]] <- ms.rep(X, X[i, ], h,  thresh = thr^2, iter)
    finals[i, ] <- all.ms[[i]]$final
    cluster.dist <- rep(0, ncluster)
    if (ncluster >= 1) {
      for (j in 1:ncluster) {
        cluster.dist[j] <- enorm(savecluster[j, ] - finals[i,])/enorm(savecluster[j, ]- Xm)
      }
    }
    
    if (ncluster == 0 || min(cluster.dist) > thr) {
      ncluster <- ncluster + 1
      savecluster <- rbind(savecluster, finals[i, ])
      cluster.label[i] <- ncluster
    } else {
      cluster.label[i] <- which(cluster.dist == min(cluster.dist))
    }
    if (plot && d == 2){
      lines(rbind(all.ms[[i]]$start, all.ms[[i]]$Meanshift.points), col = cluster.label[i] + 1)
    }    
  }
  
  if (plot && d == 2){
    points(finals, pch = 15)
  } else if (plot && d > 2 && d<=16) {
    pairs(rbind(as.matrix(X), savecluster), col = c(cluster.label + 
                                                      1, rep(1, dim(savecluster)[1])), pch = c(rep(20, 
                                                                                                   dim(X)[1]), rep(24, dim(savecluster)[1])), ...)
  } else if (plot && d>16){
    cat("Data set has inadequate dimension to be plotted in a pairs plot.\n\n")
  }
  
  for (i in subset){
    closest.label[i] <- mindist(savecluster, X[i,])$closest.item
    #closest.coords[i,]<- object$cluster.center[closest.center,]
  }
  
  dimnames(savecluster) <- list(1:ncluster, NULL)
  
  fit <- list(
    cluster.center = savecluster,
    cluster.label = cluster.label,
    closest.label = closest.label,
    h=h,
    data = X,              
    scaled= scaled, 
    scaled.by = s1,
    all.trajectories=all.ms
  )
  class(fit) <- "ms"
  return(fit)
}


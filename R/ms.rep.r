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
  if(is.null(dim(X))){X<-as.vector(X)}
  if(is.vector(X)){
    X<-matrix(X,nrow=length(X))
  }
  d  <- dim(X)[2]
  if (length(h) == 1) {
    h <- rep(h, d)
  }
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
  return(list("Meanshift.points"=M[1:s,], "Threshold.values"=th[1:s], 
              "iterations"=s, "start"=as.numeric(x0),  "final"=m))
}

# inverse meanshift for antimodes (experimental function)
ms.rep.min <- function (X, x, h, thresh=0.000001, iter=200, 
                        adjust.convergence=FALSE, verbose=TRUE ) {
              
    i <- iter+1
    th <- rep(0,iter)
    M  <-  matrix(0, iter, length(x)) 
    S  <-  matrix(0, iter, length(x)) 
    a <-  rep(1,iter)
    x0 <- x
    if(is.null(dim(X))){X<-as.vector(X)}
    if(is.vector(X)){
                  X<-matrix(X,nrow=length(X))
    }
    d  <- dim(X)[2]
    if (length(h) == 1) {
                  h <- rep(h, d)
    }
    
    Xm<-colMeans(X)
    for (j in 1: iter){
                  h<- h*sqrt(a[j])
                  m     <- meanshift(X, x, h)
                  if (any(is.na(m))){i<- j-1; if (verbose){
                    cat("required threshold not reached"); cat("\n")}; M[i,]<-NA; break}
                  M[j,] <- 2*x- m
                  S[j,]  <- m-x
                  th[j] <- enorm(S[j,])/enorm(Xm-x)  
                  if (th[j]<thresh){i<- j; if (verbose){cat("required threshold reached at iteration "); cat(j); cat("\n")}; break}
                  x     <- M[j,]
    }
                
    if (i>iter && verbose) {i<- iter; cat("required threshold not reached");  
                            M[i,]<-NA; cat("\n")
                           }
    return(list("M"=M[1:i,],"Threshold.values"=th[1:i], 
                            "iterations"=i, 
                            "start"=as.numeric(x0), "final"=M[i,]))
    }
          


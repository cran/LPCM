
coverage.raw <-function(X, vec, tau, weights=1, plot.type="p", print=FALSE, label=NULL,...){
     X<- as.matrix(X)
     p <- dim(vec)[1]
     n <- dim(X)[1]
     d <- dim(X)[2]

     if (n %% length(weights) !=0){
         stop("The length of the vector of weights is not a multiple of the sample size.")
     } else {
      weights <-rep(weights, n %/%   length(weights))
     }   

     min.distance  <- rep(0,n)
     ins.distance <- rep(1,n)

     for (i in 1: n){
        #if (i %%10==0){ print(i)}
        if (is.null(label)){
             min.distance[i] <- mindist(vec, X[i,])$mindist
        } else {
            if (!is.matrix(vec)){vec<-matrix(vec, nrow=1)}
            if (as.character(label[i]) %in%  dimnames(vec)[[1]]){
                  min.distance[i] <- mindist(matrix(vec[as.character(label[i]),],nrow=1), X[i,])$mindist   # 11/10/10 experimental, for clustering
            } else { 
                min.distance[i]<- tau+1   #01/11/10 if data point not allocatable to any centre.
            }
        }    
        ins.distance[i] <- (min.distance[i] <= tau) #indicator for being inside/outside the tube
           }
     ci<- weighted.mean(min.distance <= tau, w=weights) 
     if (plot.type %in% c("p","l")) {
        plot(X, col=ins.distance+1,...)
        if (plot.type=="p"){points(vec, col=3,lwd=2 )} else if (plot.type=="l"){lines(vec, col=3,lwd=2 )}
        }
     if (print){print(c(tau,ci))}
     return(list(tau=tau, coverage= ci, min=min.distance, inside= ins.distance))
 }



coverage<-function(X, vec,  taumin=0.02, taumax,  gridsize=25, weights=1, plot.type="o", print=TRUE,...){
  if (missing(taumax)){
  m <-mean(X)
  Xm <- sweep(X, 1, m, "-")
  Xs<- apply(X,1, sd)
  taumax<-max(Xs)
 }
 all.taus      <- seq(taumin, taumax, length=gridsize)
 all.coverages <- rep(0, gridsize)

 for (j in 1:gridsize){
    all.coverages[j]<- coverage.raw(X, vec, all.taus[j],  weights, plot.type=0, print)$coverage
 }

 if (plot.type!=0){
   plot(all.taus, all.coverages, type=plot.type, xlab=expression(tau), ylab=expression(C(tau)),...)
 }
 return(list(tau=all.taus, coverage=all.coverages))
}


lpc.coverage<-function(object, taumin=0.02, taumax, gridsize=25,  quick=TRUE, plot.type="o", print=TRUE, ...){
 
 if (class(object)=="lpc"){
   X <-object$data
   scaled<-object$scaled
   weights <- object$weights
   if (quick){
       lpc.vec<- object$LPC
   } else {
       lpc.vec<-lpc.spline(object, project=TRUE)$closest.coords
   }
 } else if (class(object)=="lpc.spline"){
    X <-object$lpcobject$data
    scaled <- object$lpcobject$scaled
    weights <- object$lpcobject$weights
    if (object$closest.coords[1] != "none") {
       lpc.vec<- object$closest.coords
     } else if (quick){
       lpc.vec<- object$lpcobject$LPC
     } else {
       lpc.vec<- lpc.spline(object$lpcobject,project=TRUE)$closest.coords
    }
  } else {
    stop("Invalid lpc object.")
  }
    
   if (missing(taumax)){   
      if (scaled){taumax<-0.5} else {
         m <-mean(X)
         Xm <- sweep(X, 2, m, "-")
         Xs<- apply(X,2, sd)
        taumax<-max(Xs)
      } 
    }
 
 
  result <- coverage(X, lpc.vec, taumin, taumax, gridsize,  weights, plot.type=0, print)
  if (plot.type!=0){
   plot(result$tau,result$coverage, type=plot.type, xlab=expression(tau), ylab=expression(C(tau)),...)
 }
  return(result)
}  
  
  

lpc.self.coverage <-
function (X, taumin = 0.02, taumax = 0.5, gridsize = 25, x0, 
    mult = 1, way = "two", scaled = TRUE, weights = 1, pen = 2, 
    depth = 1, control = lpc.control(boundary = 0, cross = FALSE), 
    quick = TRUE, plot.type = "o", print = TRUE, ...) 
{
    if (class(X) %in% c("lpc", "lpc.spline")) {
        stop("Invalid data matrix.")
    }
    Xi <- as.matrix(X)
    N <- dim(Xi)[1]
    d <- dim(Xi)[2]
    if (!missing(x0)) {
        x0 <- matrix(x0, ncol = d, byrow = TRUE)
        if (dim(x0)[1] < mult) {
            rn <- runif(mult - dim(x0)[1], 1, N + 1)%/%1
            x0 <- rbind(x0, Xi[rn, ])
        }
    }
    else {
        rn <- runif(mult, 1, N + 1)%/%1
        x0 <- matrix(Xi[rn, ], length(rn), d)
    }
    if ((!scaled) && taumax < 1) {
        warning("Please adjust the range (taumin, taumax) of tube widths by hand, as the data are not scaled.")
    }
    Pm <- NULL
    h0 <- taumin
    h1 <- taumax
    h <- seq(h0, h1, length = gridsize)
    #n <- gridsize
    cover <- matrix(0, gridsize, 2)
    for (i in 1:gridsize) {
        new.h0 <- h[i]
        fit <- lpc(Xi, h = new.h0, t0 = new.h0, x0 = x0, mult = mult, 
            way = way, scaled = scaled, weights = weights, pen = pen, 
            depth = depth, control)
        if (!quick) {
            fit.spline <- lpc.spline(fit, project = TRUE)
        }
        if (quick) {
            Pm[[i]] <- fit$LPC
        }
        else {
            Pm[[i]] <- fit.spline$closest.coords
        }
        cover[i, ] <- as.numeric(coverage.raw(fit$data, Pm[[i]], 
            new.h0, weights, plot.type = 0, print = print)[1:2])
    }
    select <- select.self.coverage(self = cover, sens = 0.02, 
        from = 2/3, plot.type = plot.type, auto = FALSE)
    result <- list(self.coverage.curve = cover, select = select, 
        type = "lpc")
    class(result) <- "self"
    result
}


select.self.coverage <-
function (self, sens = 0.02, from, plot.type = "o", auto = FALSE) 
{
    if (class(self) == "self") {
        cover <- self$self.coverage.curve
    }
    else {
        cover <- self
    }
    if (missing(from)) {
        if (class(self) == "self") {
            from <- switch(self$type, lpc = 2/3, ms = 1/3)
        }
        else stop("Please specify `from' argument.")
    }
    n <- dim(cover)[1]
    diff1 <- diff2 <- rep(0, n)
    diff1[2:n] <- cover[2:n, 2] - cover[1:(n - 1), 2]
    diff2[2:(n - 1)] <- diff1[3:n] - diff1[2:(n - 1)]
    select <- select.coverage <- select.2diff <- NULL
    if (plot.type != 0) {
        plot(cover, type = plot.type, xlab = "h", ylab = "S(h)", 
            ylim = c(0, 1))
    }
    for (i in (3:(n - 1))) {
        if (diff2[i] < -sens && cover[i, 2] > max(from, cover[1:(i - 
            1), 2])) {
            select <- c(select, cover[i, 1])
            select.coverage <- c(select.coverage, cover[i,2])
            select.2diff <- c(select.2diff, diff2[i])
            if (plot.type != 0) {
                segments(cover[i, 1], 0, cover[i, 1], cover[i, 
                  2], col = 3, lty = 1)
            }
        }
    }
    return(list("select"=select, "select.coverage"=select.coverage, "select.2diff"=select.2diff))
}


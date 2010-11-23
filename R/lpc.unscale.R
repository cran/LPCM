


lpc.unscale <-function(object){

      if (class(object)=="lpc"){ lpcobject<-object
                                 splineobject<-NULL
                                }
      if (class(object)=="lpc.spline"){
                                lpcobject <- object$lpcobject
                                splineobject<-object
                              }

     if (!lpcobject$scaled){stop("The lpcobject was not fitted with scaled data, so it cannot be unscaled!")}
      
      #if (missing(lpcobject)){lpcobject <- splineobject$lpcobject}
      
      LPC            <- sweep(lpcobject$LPC,2, lpcobject$Misc$scaled.by, "*")
      start          <- sweep(lpcobject$starting.points,2, lpcobject$Misc$scaled.by, "*")
      data           <- sweep(lpcobject$data,2, lpcobject$Misc$scaled.by , "*")
      knots.coords   <- list(NULL)
      closest.coords <- list (NULL)

      if (!is.null(splineobject)){
        lk  <- length(splineobject$knots.coords)
        for (j in 1:lk){  
              knots.coords[[j]] <-  sweep(splineobject$knots.coords[[j]],1, lpcobject$Misc$scaled.by, "*")
    }
       }

      if (!is.null(splineobject) && splineobject$closest.coords!="none" ){ 
            closest.coords <- sweep(splineobject$closest.coords,2, lpcobject$Misc$scaled.by , "*")
      }

      
    return(list("LPC"=LPC,"data"=data, "starting.points"=start, "knots.coords"=knots.coords,"closest.coords"=closest.coords))    
 }


# up to v 0.39:
 #lpc.unscale <-function(lpcobject, splineobject=NULL){
#
#      if (missing(lpcobject)){lpcobject <- splineobject$lpcobject}
#      
#      LPC            <- sweep(lpcobject$LPC,2, lpcobject$Misc$scaled.by, "*")
#      start          <- sweep(lpcobject$starting.points,2, lpcobject$Misc$scaled.by, "*")
#      data           <- sweep(lpcobject$data,2, lpcobject$Misc$scaled.by , "*")
#      knots.coords   <- list(NULL)
#      closest.coords <- list (NULL)#
#
#      if (!is.null(splineobject)){
#        lk  <- length(splineobject$knots.coords)
#        for (j in 1:lk){  
#              knots.coords[[j]] <-  sweep(splineobject$knots.coords[[j]],1, lpcobject$Misc$scaled.by, "*")
#    }
#       }
#
#      if (!is.null(splineobject) && splineobject$closest.coords!="none" ){ 
#            closest.coords <- sweep(splineobject$closest.coords,2, lpcobject$Misc$scaled.by , "*")
#      }
# 
#    return(list("LPC"=LPC,"data"=data, "starting.points"=start, "knots.coords"=knots.coords,"closest.coords"=closest.coords#))    
# }

lpc.control <- function (iter =100, boundary = 0.005, convergence.at =0.00001, pruning.thresh=0.0, rho0=0.4, cross=TRUE) 
{
    if (boundary!=0 &&  boundary < convergence.at) {
        warning("The boundary correction will not have any effect if its threshold is set to a smaller value than the convergence criterion.") 
      }
    if (iter <= 9) {
        warning("Please choose iter=10 at least. The default value iter=100 has been used instead.")
        iter <- 100
      }
    if (pruning.thresh<0) {
        warning("This should be a non-negative number. The default 0 has been used instead.") 
    pruning.thresh<-0
    }
    if (rho0<0) {
        warning("This should be a non-negative number. The default 0.4 has been used instead.") 
    rho0<-0.4
    }    
    if (!is.logical(cross)){
       warning("cross needs to be a Boolean. The default FALSE has been used instead.")
     }   
    list(iter = iter, boundary = boundary, convergence.at = convergence.at, 
        pruning.thresh=pruning.thresh, rho0=rho0, cross=cross)
}

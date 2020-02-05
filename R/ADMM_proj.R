
l1proj<-function(v, b){
  
  # Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo
  # Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML
  
  stopifnot(b>0)
  
  u <- sort(abs(v),decreasing=TRUE)
  sv <- cumsum(u)
  rho <- max(which(u>(sv-b)/1:length(u)))
  theta <- max(0, (sv[rho]-b)/rho)
  w <-sign(v) * pmax(abs(v)-theta,0)
  
  return(w)
}

#' ADMM algorithm
#'
#' Finds the nearest positive semi-definite matrix with respect to the max norm
#'
#' @param mat Matrix to be projected
#' @param epsilon Approximation of the space of positive semi-definite matrix at epsilon
#' @param mu Penalty parameter of ADMM algorithm
#' @param it.max Number maximum of iterations
#' @param etol Tolerance parameter for the convergence of primal and dual residual
#' @param etol_distance Tolerance parameter for the convergence of the distance
#' 
#' @return list containing \itemize{
#' \item \code{mat} Projected matrix
#' \item \code{df_ADMM} dataframe containing parameters of the ADMM algorithm for each iteration of the algorithm
#' }
#' 
#' @examples 
#' M = matrix(-1,20,20)
#' mat_proj <- BDcocolasso::ADMM_proj(mat=M)$mat
#' 
#' @seealso \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
#' @export
#' 

ADMM_proj<-function(mat,
                    epsilon=1e-4,
                    mu=10,
                    it.max=1e3,
                    etol=1e-4,
                    etol_distance = 1e-4){
  

  
  p<-nrow(mat)
  
  # Initialization
  R<-diag(mat)
  S<-matrix(0,p,p)
  L<-matrix(0,p,p)
  
  itr<-0
  iteration <- eps_R <- eps_S <- eps_primal <- time <- distance <- NULL
  while (itr<it.max) {
    #print(itr)
    Rp<-R
    Sp<-S
    start <- Sys.time()
    # Subproblem I: R step
    W<-mat+S+mu*L
    W.eigdec<-eigen(W, symmetric=TRUE) 
    W.V<-W.eigdec$vectors
    W.D<-W.eigdec$values
    R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)
    
    # Subproblem II: S step
    M<-R-mat-mu*L     
    S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)    
    for (i in 2:p){
      for (j in 1:(i-1)){
        S[j,i]<-S[i,j]
      }
    }
    
    # L step: update the Lagrange parameter
    L<-L-(R-S-mat)/mu
    end <- Sys.time()
    #Stocking the values of different parameters with the number of iterations
    iteration <- c(iteration, itr)
    eps_R <- c(eps_R,max(abs(R-Rp)))
    eps_S <- c(eps_S,max(abs(S-Sp)))
    eps_primal <- c(eps_primal, max(abs(R-S-mat)))
    time <- c(time, end - start)
    distance <- c(distance,max(abs(R-mat)))
    
    # Stopping Rule                        
    #cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if (((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)) || (abs(max(abs(Rp-mat)) - max(abs(R-mat)))<etol_distance)){
      itr<-it.max
    } else {
      itr<-itr+1
    }
    
    if (itr%%20==0) {
      mu<-mu/2
    }
  }
  df_ADMM <- data.frame(iteration = iteration, eps_R = eps_R, eps_S=eps_S, eps_primal=eps_primal, time=time, distance=distance)
  return(list(mat=R,df_ADMM=df_ADMM))
  
}

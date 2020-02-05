getRmat=function(X,y=NA){
  m=dim(X)[1]; n=dim(X)[2]
  R=matrix(0,n,n); rho_paired=c()
  for (j in (1:n)){ #col
    for (i in (1:j)){#row
      colProd=X[,i]*X[,j]; nij=sum(!is.na(colProd))
      R[i,j]=(nij/m)
      R[j,i]=R[i,j]
    }
    rho_j=(sum(X[,j]*y,na.rm=T))/nij
    rho_paired=c(rho_paired,rho_j)
  }
  return(list(rho_paired,R))
}

updateB_maxNorm_hmLasso=function(Akp1,Lk,S_paired,mu,W){
  C=Akp1-S_paired-mu*Lk
  C_vec=as.vector(C); W_vec=as.vector(W); n=dim(C)[1]; m=dim(C)[2]
  WC=abs(C_vec)*W_vec
  WC_sort_inx=order(WC,decreasing = T);
  W_sort=W_vec[WC_sort_inx]; C_sort=C_vec[WC_sort_inx]
  l=1; frac=0
  while (W_sort[l]*abs(C_sort[l])> frac && l<=length(C_vec)){
    frac=(sum(abs(C_sort[1:l]))-mu/2)/(sum(1/(W_sort[1:l])))
    l=l+1
  }
  d=frac
  
  b_vec=mapply(FUN=function(t) ifelse(W_vec[t]*abs(C_vec[t])>d, d*sign(C_vec[t])/W_vec[t],
                                      C_vec[t]),1:length(W_sort))
  Bkp1=matrix(b_vec,n,m)
  return(Bkp1)
}

#' HM-lasso algorithm
#'
#' Finds the nearest positive semi-definite matrix with respect to the Frobenius norm
#'
#' @param sigmaHat covariance matrix
#' @param R weight matrix
#' @param a hmlasso (1) or cocolasso (0)
#' @param iter_max Number maximum of iterations
#' @param epsilon Stopping criterion
#' @param mu Penalty parameter of ADMM algorithm
#' @param tolerance Tolerance parameter for the convergence of primal and dual residual
#' @param norm Frobenius or max norm
#' 
#' @return list containing \itemize{
#' \item \code{Ak} Projected matrix
#' }
#' 
#' @examples  
#' M = matrix(-1,20,20)
#' mat_proj <- BDcocolasso::HM_proj(sigmaHat=M)
#' 
#' @seealso \url{https://arxiv.org/pdf/1811.00255.pdf}
#' @export
#' 

HM_proj=function(sigmaHat,R=NULL,a=1,iter_max=1000,epsilon=1e-4,mu=10,tolerance=1e-4,norm="F"){
  iter=0
  S_paired=sigmaHat; n=nrow(S_paired)
  if (is.null(R)){
    W=matrix(1,n,n)^a
  } else{ W=R^a}

  Ak=S_paired
  Bk=matrix(0,n,n); Lk=matrix(0,n,n)
  while (iter<iter_max){
    A=Bk+S_paired+mu*Lk
    A_eigdec=eigen(A, symmetric=TRUE) 
    A_V=A_eigdec$vectors
    A_D=A_eigdec$values
    Akp1=A_V%*%diag(pmax(A_D,epsilon))%*%t(A_V)
    
    if (norm =='F'){
      Bkp1 = (Akp1-S_paired-mu*Lk)/(mu*W*W+matrix(1,n,n)) 
    } else {
      Bkp1 = updateB_maxNorm_hmLasso(Akp1,Lk,S_paired,mu,W)
    }
    Lkp1=Lk-(Akp1-Bkp1-S_paired)/mu
    
    if (max(max(abs(Akp1-Ak)),max(abs(Bkp1-Bk)),max(abs(Lkp1-Lk)))<tolerance ){
      iter=iter_max
    } else { iter=iter+1}
    Ak=Akp1; Bk=Bkp1; Lk=Lkp1
  }
  return(Ak)
}

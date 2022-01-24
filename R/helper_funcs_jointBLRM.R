

#'@keywords internal
#determines the maximum dosing step in a vector of (positive) dose levels.
#used to determine default escalation steps
max_step_BLRM <- function(dose){
  doselv <- sort(as.numeric(unique(dose)))
  if(length(doselv)>1){
    return(
      max(sapply(1:length(doselv), FUN=function(i){
                                    if(i<length(doselv)){
                                      return(doselv[i+1]/doselv[i])
                                    }else{
                                      return(0)}}
          )
      )
    )
  }else{
    return(1)
  }
}

#'@keywords internal
#transforms sequence of study names to the conventions
harmonize_vecnames_jointBLRM <- function(x){
  res <- x
  res[which(x=="mono1.a")] <- rep(1, times = length(which(x=="mono1.a")))
  res[which(x=="mono1.b")] <- rep(4, times = length(which(x=="mono1.b")))
  res[which(x=="mono2.a")] <- rep(2, times = length(which(x=="mono2.a")))
  res[which(x=="mono2.b")] <- rep(5, times = length(which(x=="mono2.b")))
  res[which(x=="combi.a")] <- rep(3, times = length(which(x=="combi.a")))
  res[which(x=="combi.b")] <- rep(6, times = length(which(x=="combi.b")))

  fx <- levels(factor(res))
  idx_hist <- which(!fx%in%c(1, 2, 3, 4, 5, 6))
  if(!length(idx_hist)==0){
    for(n in 1:length(idx_hist)){
      res[which(x==fx[idx_hist[n]])] <- rep(6+n, times = length(which(x==fx[idx_hist[n]])))
    }
  }
  return(as.numeric(res))
}

#'@keywords internal
#test for integer vectors (adapted from R documentation)
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5){
    if(!is.numeric(x)){
      return(F)
    }
    if(any(is.na(x))){
      return(F)
    }
    return(all(abs(x - round(x)) < tol))
  }


#'@keywords internal
#test for other numbers
is.num <-
  function(x,
           low = -Inf, up= Inf, len = length(x),
           uB=TRUE, lB=TRUE){
    if(!is.numeric(x)){
      return(F)
    }
    if(any(is.na(x))){
      return(F)
    }
    if(!length(x)==len){
      return(F)
    }
    if(uB & lB){
      return(all(low<=x & x<=up))
    }else if(uB & !lB){
      return(all(low<x & x<=up))
    }else if(!uB & lB){
      return(all(low<=x & x<up))
    }else{
      return(all(low<x & x<up))
    }
  }

#'@keywords internal
#check for format of decision rule
is.dec.rule <- function(x){
  na <- names(x)
  if(!all(c("TARGET.PROB", "PAT.AT.MTD", "MIN.PAT",
            "MIN.DLT", "RULE") %in% na)){
    return(F)
  }
  if(!is.num(x$TARGET.PROB, len=1, low=0, up=1, lB=T, uB=F)){
    return(F)
  }
  if(!is.wholenumber(x$PAT.AT.MTD)){
    return(F)
  }
  if(!length(x$PAT.AT.MTD)==1){
    return(F)
  }
  if(!x$PAT.AT.MTD>=0){
    return(F)
  }

  if(!is.wholenumber(x$MIN.PAT)){
    return(F)
  }
  if(!length(x$MIN.PAT)==1){
    return(F)
  }
  if(!x$MIN.PAT>=0){
    return(F)
  }
  if(!is.wholenumber(x$MIN.DLT)){
    return(F)
  }
  if(!length(x$MIN.DLT)==1){
    return(F)
  }
  if(!x$MIN.DLT>=0){
    return(F)
  }

  if(!is.wholenumber(x$RULE)){
    return(F)
  }
  if(!length(x$RULE)==1){
    return(F)
  }
  if(!x$RULE%in%c(1, 2)){
    return(F)
  }
  return(T)
}


#'@keywords internal
#check for format of prior for hyper means
is.prior.mu <- function(x){
  na <- names(x)
  expec <- paste("mu_", c("a1", "a2", "b1", "b2", "eta"), sep="")
  if(!all(expec %in% na)){
    return(F)
  }
  if(!is.num(x$mu_a1, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$mu_a1[2]>0){
    return(F)
  }
  if(!is.num(x$mu_a2, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$mu_a2[2]>0){
    return(F)
  }
  if(!is.num(x$mu_b1, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$mu_b1[2]>0){
    return(F)
  }
  if(!is.num(x$mu_b2, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$mu_b2[2]>0){
    return(F)
  }
  if(!is.num(x$mu_eta, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$mu_eta[2]>0){
    return(F)
  }
  return(T)
}

#'@keywords internal
test.prior.mu <- function(prior.mu){

  names_pmu <- names(prior.mu)
  if(!(all(c("mu_a1", "mu_a2", "mu_b1", "mu_b2", "mu_eta")%in%names_pmu))){
    stop("`prior.mu` must have entries \"mu_a1\",\"mu_b1\",\"mu_a2\",\"mu_b2\", and \"mu_eta\".")
  }
  if(!is.numeric(prior.mu$mu_a1)){

    stop("`prior.mu$mu_a1` must be numeric.")
  }
  if(any(is.na(prior.mu$mu_a1))){

    stop("`prior.mu$mu_a1` must not be NA.")
  }
  if(!length(prior.mu$mu_a1)==2){

    stop("`prior.mu$mu_a1` must have length 2.")
  }
  if(!prior.mu$mu_a1[2]>0){

    stop("`prior.mu$mu_a1[2]` is the SD and cannot be negative.")
  }
  if(!is.numeric(prior.mu$mu_b1)){

    stop("`prior.mu$mu_b1` must be numeric.")
  }
  if(any(is.na(prior.mu$mu_b1))){

    stop("`prior.mu$mu_b1` must not be NA.")
  }
  if(!length(prior.mu$mu_b1)==2){

    stop("`prior.mu$mu_b1` must have length 2.")
  }
  if(!prior.mu$mu_b1[2]>0){

    stop("`prior.mu$mu_b1[2]` is the SD and cannot be negative.")
  }
  if(!is.numeric(prior.mu$mu_a2)){

    stop("`prior.mu$mu_a2` must be numeric.")
  }
  if(any(is.na(prior.mu$mu_a2))){

    stop("`prior.mu$mu_a2` must not be NA.")
  }
  if(!length(prior.mu$mu_a2)==2){

    stop("`prior.mu$mu_a2` must have length 2.")
  }
  if(!prior.mu$mu_a2[2]>0){

    stop("`prior.mu$mu_a2[2]` is the SD and cannot be negative.")
  }

  if(!is.numeric(prior.mu$mu_b2)){

    stop("`prior.mu$mu_b2` must be numeric.")
  }
  if(any(is.na(prior.mu$mu_b2))){

    stop("`prior.mu$mu_b2` must not be NA.")
  }
  if(!length(prior.mu$mu_b2)==2){

    stop("`prior.mu$mu_b2` must have length 2.")
  }
  if(!prior.mu$mu_b2[2]>0){

    stop("`prior.mu$mu_b2[2]` is the SD and cannot be negative.")
  }
  if(!is.numeric(prior.mu$mu_eta)){

    stop("`prior.mu$mu_eta` must be numeric.")
  }
  if(any(is.na(prior.mu$mu_eta))){

    stop("`prior.mu$mu_eta` must not be NA.")
  }
  if(!length(prior.mu$mu_eta)==2){

    stop("`prior.mu$mu_eta` must have length 2.")
  }
  if(!prior.mu$mu_eta[2]>0){

    stop("`prior.mu$mu_eta[2]` is the SD and cannot be negative.")
  }

}

#'@keywords internal
test.prior.mu.covar <- function(prior.mu){

  names_pmu <- names(prior.mu)
  if(!(all(c("mu_g1", "mu_g2")%in%names_pmu))){
    stop("`prior.mu.covar` must have entries \"mu_g1\" and \"mu_g2\".")
  }
  if(!is.numeric(prior.mu$mu_g1)){

    stop("`prior.mu.covar$mu_g1` must be numeric.")
  }
  if(any(is.na(prior.mu$mu_g1))){

    stop("`prior.mu.covar$mu_g1` must not be NA.")
  }
  if(!length(prior.mu$mu_g1)==2){

    stop("`prior.mu.covar$mu_g1` must have length 2.")
  }
  if(!prior.mu$mu_g1[2]>0){

    stop("`prior.mu.covar$mu_g1[2]` is the SD and cannot be negative.")
  }


  if(!is.numeric(prior.mu$mu_g2)){

    stop("`prior.mu.covar$mu_g2` must be numeric.")
  }
  if(any(is.na(prior.mu$mu_g2))){

    stop("`prior.mu.covar$mu_g2` must not be NA.")
  }
  if(!length(prior.mu$mu_g2)==2){

    stop("`prior.mu.covar$mu_g2` must have length 2.")
  }
  if(!prior.mu$mu_g2[2]>0){

    stop("`prior.mu.covar$mu_g2[2]` is the SD and cannot be negative.")
  }


}

#'@keywords internal
test.prior.tau.covar <- function(prior.tau){

  names_ptau <- names(prior.tau)
  if(!(all(c("tau_g1", "tau_g2")%in%names_ptau))){
    stop("`prior.tau.covar` must have entries \"tau_g1\" and \"tau_g2\".")
  }
  if(!is.numeric(prior.tau$tau_g1)){

    stop("`prior.tau.covar$tau_g1` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_g1))){

    stop("`prior.tau.covar$tau_g1` must not be NA.")
  }
  if(!length(prior.tau$tau_g1)==2){

    stop("`prior.tau.covar$tau_g1` must have length 2.")
  }
  if(!prior.tau$tau_g1[2]>0){

    stop("`prior.tau.covar$tau_g1[2]` is the SD and cannot be negative.")
  }


  if(!is.numeric(prior.tau$tau_g2)){

    stop("`prior.tau.covar$tau_g2` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_g2))){

    stop("`prior.tau.covar$tau_g2` must not be NA.")
  }
  if(!length(prior.tau$tau_g2)==2){

    stop("`prior.tau.covar$tau_g2` must have length 2.")
  }
  if(!prior.tau$tau_g2[2]>0){

    stop("`prior.tau.covar$tau_g2[2]` is the SD and cannot be negative.")
  }


}


#'@keywords internal
#check for format of prior for between-trial heterogeneity
is.prior.tau <- function(x){
  na <- names(x)
  expec <- paste("tau_", c("a1", "a2", "b1", "b2", "eta"), sep="")
  if(!all(expec %in% na)){
    return(F)
  }
  if(!is.num(x$tau_a1, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$tau_a1[2]>0){
    return(F)
  }
  if(!is.num(x$tau_a2, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$tau_a2[2]>0){
    return(F)
  }
  if(!is.num(x$tau_b1, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$tau_b1[2]>0){
    return(F)
  }
  if(!is.num(x$tau_b2, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$tau_b2[2]>0){
    return(F)
  }
  if(!is.num(x$tau_eta, len=2, lB=F, uB=F)){
    return(F)
  }
  if(!x$tau_eta[2]>0){
    return(F)
  }
  return(T)
}


#'@keywords internal
test.prior.mu.mono <- function(prior.mu){
  names_pmu <- names(prior.mu)
  if(!(all(c("mu_a", "mu_b")%in%names_pmu))){

    stop("`prior.mu.mono` must have entries \"mu_a\" and \"mu_b\".")
  }
    if(!is.numeric(prior.mu$mu_a)){

      stop("`prior.mu.mono$mu_a` must be numeric.")
    }
    if(any(is.na(prior.mu$mu_a))){

      stop("`prior.mu.mono$mu_a` must not be NA.")
    }
    if(!length(prior.mu$mu_a)==2){

      stop("`prior.mu.mono$mu_a` must have length 2.")
    }
    if(!prior.mu$mu_a[2]>0){

      stop("`prior.mu.mono$mu_a[2]` is the SD and cannot be negative.")
    }
    if(!is.numeric(prior.mu$mu_b)){

      stop("`prior.mu.mono$mu_b` must be numeric.")
    }
    if(any(is.na(prior.mu$mu_b))){

      stop("`prior.mu.mono$mu_b` must not be NA.")
    }
    if(!length(prior.mu$mu_b)==2){

      stop("`prior.mu.mono$mu_b` must have length 2.")
    }
    if(!prior.mu$mu_b[2]>0){

      stop("`prior.mu.mono$mu_b[2]` is the SD and cannot be negative.")
    }

}


#'@keywords internal
test.prior.tau.mono <- function(prior.tau){
  names_ptau <- names(prior.tau)
  if(!(all(c("tau_a", "tau_b")%in%names_ptau))){

    stop("`prior.tau.mono` must have entries \"tau_a\" and \"tau_b\".")
  }
  if(!is.numeric(prior.tau$tau_a)){

    stop("`prior.tau.mono$tau_a` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_a))){

    stop("`prior.tau.mono$tau_a` must not be NA.")
  }
  if(!length(prior.tau$tau_a)==2){

    stop("`prior.tau.mono$tau_a` must have length 2.")
  }
  if(!prior.tau$tau_a[2]>0){

    stop("`prior.tau.mono$tau_a[2]` is the SD and cannot be negative.")
  }
  if(!is.numeric(prior.tau$tau_b)){

    stop("`prior.tau.mono$tau_b` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_b))){

    stop("`prior.tau.mono$tau_b` must not be NA.")
  }
  if(!length(prior.tau$tau_b)==2){

    stop("`prior.tau.mono$tau_b` must have length 2.")
  }
  if(!prior.tau$tau_b[2]>0){

    stop("`prior.tau.mono$tau_b[2]` is the SD and cannot be negative.")
  }

}

#'@keywords internal
test.prior.tau <- function(prior.tau){
  names_ptau <- names(prior.tau)
  if(!(all(c("tau_a1", "tau_a2", "tau_b1", "tau_b2", "tau_eta")%in%names_ptau))){

    stop("`prior.tau` must have entries \"tau_a1\",\"tau_b1\",\"tau_a2\",\"tau_b2\", and \"tau_eta\".")
  }
  if(!is.numeric(prior.tau$tau_a1)){

    stop("`prior.tau$tau_a1` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_a1))){

    stop("`prior.tau$tau_a1` must not be NA.")
  }
  if(!length(prior.tau$tau_a1)==2){

    stop("`prior.tau$tau_a1` must have length 2.")
  }
  if(!prior.tau$tau_a1[2]>0){

    stop("`prior.tau$tau_a1[2]` is the SD and cannot be negative.")
  }
  if(!is.numeric(prior.tau$tau_b1)){

    stop("`prior.tau$tau_b1` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_b1))){

    stop("`prior.tau$tau_b1` must not be NA.")
  }
  if(!length(prior.tau$tau_b1)==2){

    stop("`prior.tau$tau_b1` must have length 2.")
  }
  if(!prior.tau$tau_b1[2]>0){

    stop("`prior.tau$tau_b1[2]` is the SD and cannot be negative.")
  }
  if(!is.numeric(prior.tau$tau_a2)){

    stop("`prior.tau$tau_a2` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_a2))){

    stop("`prior.tau$tau_a2` must not be NA.")
  }
  if(!length(prior.tau$tau_a2)==2){

    stop("`prior.tau$tau_a2` must have length 2.")
  }
  if(!prior.tau$tau_a2[2]>0){

    stop("`prior.tau$tau_a2[2]` is the SD and cannot be negative.")
  }

  if(!is.numeric(prior.tau$tau_b2)){

    stop("`prior.tau$tau_b2` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_b2))){

    stop("`prior.tau$tau_b2` must not be NA.")
  }
  if(!length(prior.tau$tau_b2)==2){

    stop("`prior.tau$tau_b2` must have length 2.")
  }
  if(!prior.tau$tau_b2[2]>0){

    stop("`prior.tau$tau_b2[2]` is the SD and cannot be negative.")
  }
  if(!is.numeric(prior.tau$tau_eta)){

    stop("`prior.tau$tau_eta` must be numeric.")
  }
  if(any(is.na(prior.tau$tau_eta))){

    stop("`prior.tau$tau_eta` must not be NA.")
  }
  if(!length(prior.tau$tau_eta)==2){

    stop("`prior.tau$tau_eta` must have length 2.")
  }
  if(!prior.tau$tau_eta[2]>0){

    stop("`prior.tau$tau_eta[2]` is the SD and cannot be negative.")
  }
}

#'@keywords internal
#cleans historical data from NA informations
clean.na.hist <- function(x){
  res <- x
  res$dose1[which(is.na(x$dose1))] <- rep(0, length(which(is.na(x$dose1))))
  res$dose2[which(is.na(x$dose2))] <- rep(0, length(which(is.na(x$dose2))))
  return(res)
}

#'@keywords internal
# check for inpurt format of historical data
is.historical.cov.data <- function(x){
  na <- names(x)
  if(!all(c("dose1", "dose2", "n.pat", "n.dlt", "trial", "covar") %in% na)){
    return(F)
  }
  x <- clean.na.hist(x)
  len <- length(x$dose1)
  if(!is.num(x$dose1, len=len, low=0, lB=T, uB=F)){
    return(F)
  }
  if(!is.num(x$dose2, len=len, low=0, lB=T, uB=F)){
    return(F)
  }
  if(!is.num(x$n.pat, len=len, low=0, lB=T, uB=F)|
     !is.wholenumber(x$n.pat)){
    return(F)
  }
  if(!is.num(x$n.dlt, len=len, low=0, lB=T, uB=F)|
     !is.wholenumber(x$n.dlt)){
    return(F)
  }
  if(!is.num(x$covar, len=len, low=0, up = 1, lB=T, uB=T)|
     !is.wholenumber(x$covar)){
    return(F)
  }
  if(!is.numeric(x$trial) & !is.character(x$trial)){
    return(F)
  }
  if(any(is.na(x$trial))){
    return(F)
  }
  if(!length(x$trial)==len){
    return(F)
  }
  return(T)
}

#'@keywords internal
# check for inpurt format of historical data
is.historical.data <- function(x){
  na <- names(x)
  if(!all(c("dose1", "dose2", "n.pat", "n.dlt", "trial") %in% na)){
    return(F)
  }
  x <- clean.na.hist(x)
  len <- length(x$dose1)
  if(!is.num(x$dose1, len=len, low=0, lB=T, uB=F)){
    return(F)
  }
  if(!is.num(x$dose2, len=len, low=0, lB=T, uB=F)){
    return(F)
  }
  if(!is.num(x$n.pat, len=len, low=0, lB=T, uB=F)|
     !is.wholenumber(x$n.pat)){
    return(F)
  }
  if(!is.num(x$n.dlt, len=len, low=0, lB=T, uB=F)|
     !is.wholenumber(x$n.dlt)){
    return(F)
  }
  if(!is.numeric(x$trial) & !is.character(x$trial)){
    return(F)
  }
  if(any(is.na(x$trial))){
    return(F)
  }
  if(!length(x$trial)==len){
    return(F)
  }
  return(T)
}

#'@keywords internal
# remove observations from historical data
remove.noninf.cov.obs <- function(x, idx){
  len <- length(x$dose1)-length(idx)
  allid <- 1:length(x$dose1)
  keep <- allid[which(!allid%in%idx)]
  res <- list(
    "dose1" = x$dose1[keep],
    "dose2" = x$dose2[keep],
    "n.pat" = x$n.pat[keep],
    "n.dlt" = x$n.dlt[keep],
    "trial" = x$trial[keep],
    "covar" = x$covar[keep]
  )
  return(res)
}

#'@keywords internal
# remove observations from historical data
remove.noninf.obs <- function(x, idx){
  len <- length(x$dose1)-length(idx)
  allid <- 1:length(x$dose1)
  keep <- allid[which(!allid%in%idx)]
  res <- list(
    "dose1" = x$dose1[keep],
    "dose2" = x$dose2[keep],
    "n.pat" = x$n.pat[keep],
    "n.dlt" = x$n.dlt[keep],
    "trial" = x$trial[keep]
  )
  return(res)
}

#'@keywords internal
# transforms prior input format to the outputted matrix
# mostly to give cleaner excel sheets as output
prior_mat_out <- function(m, t){
  res <- matrix(nrow=10, ncol=2)
  na <- c("a1", "b1", "a2", "b2", "eta")
  rownames(res) <- c(paste0("mu_", na), paste0("tau_", na))
  colnames(res) <- c("mean", "SD")
  res[1, ] <- m$mu_a1
  res[2, ] <- m$mu_b1
  res[3, ] <- m$mu_a2
  res[4, ] <- m$mu_b2
  res[5, ] <- m$mu_eta
  res[6, ] <- t$tau_a1
  res[7, ] <- t$tau_b1
  res[8, ] <- t$tau_a2
  res[9, ] <- t$tau_b2
  res[10, ] <- t$tau_eta
  return(res)
}

#'@keywords internal
prior_mat_out_cov <- function(m, t){
  res <- matrix(nrow=14, ncol=2)
  na <- c("a1", "b1", "a2", "b2", "eta", "g1", "g2")
  rownames(res) <- c(paste0("mu_", na), paste0("tau_", na))
  colnames(res) <- c("mean", "SD")
  res[1, ] <- m$mu_a1
  res[2, ] <- m$mu_b1
  res[3, ] <- m$mu_a2
  res[4, ] <- m$mu_b2
  res[5, ] <- m$mu_eta
  res[6, ] <- m$mu_c1
  res[7, ] <- m$mu_c2
  res[8, ] <- t$tau_a1
  res[9, ] <- t$tau_b1
  res[10, ] <- t$tau_a2
  res[11, ] <- t$tau_b2
  res[12, ] <- t$tau_eta
  res[13, ] <- t$tau_c1
  res[14, ] <- t$tau_c2
  return(res)
}


#'@keywords internal
# construct data matrix that can be outputted
# mostly for cleaner outputs
data_matrix_jointBLRM <- function(data){

  dose1 <- data$dose1
  dose2 <- data$dose2
  N.Pat <- data$n.pat
  N.DLT <- data$n.dlt
  Study <- data$trial
  nobs <- length(dose1)

  if(all(N.Pat==0)){
    return("No Input Data")
  }else{
    res <- matrix(NA, nrow=5, ncol=nobs)
    res[1,]<-dose1
    res[2,]<-dose2
    res[3,]<-N.Pat
    res[4,]<-N.DLT
    res[5,]<-Study
    rownames(res) <- c("Dose 1", "Dose 2", "N Pat.", "N DLT", "Trial")
    colnames(res) <- paste0("Obs", 1:nobs)
    return(res)
  }

}


#'@keywords internal
# construct data matrix that can be outputted
# mostly for cleaner outputs
data_matrix_covariate_jointBLRM <- function(data){

  dose1 <- data$dose1
  dose2 <- data$dose2
  N.Pat <- data$n.pat
  N.DLT <- data$n.dlt
  Study <- data$trial
  cov <- data$covar
  nobs <- length(dose1)

  if(all(N.Pat==0)){
    return("No Input Data")
  }else{
    res <- matrix(NA, nrow=6, ncol=nobs)
    res[1,]<-dose1
    res[2,]<-dose2
    res[3,]<-N.Pat
    res[4,]<-N.DLT
    res[5,]<-Study
    res[6,]<-cov
    rownames(res) <- c("Dose 1", "Dose 2", "N Pat.", "N DLT", "Trial", "Covariate")
    colnames(res) <- paste0("Obs", 1:nobs)
    return(res)
  }

}
#'@keywords internal
# allows to write all input parameters
# to a clean output matrix.
# depending on activated options
# the contents may vary slightly.
input_out_jointBLRM <- function(
  active.mono1.a,
  active.mono1.b,
  active.mono2.a,
  active.mono2.b,
  active.combi.a,
  active.combi.b,

  dose.ref1,
  dose.ref2,
  seed,
  dosing.intervals,
  saturating,
  ewoc,
  loss.weights,
  dynamic.weights,

  esc.rule,
  esc.comp.max,

  cohort.size.mono1.a,
  cohort.prob.mono1.a,
  cohort.size.mono1.b,
  cohort.prob.mono1.b,
  cohort.size.mono2.a,
  cohort.prob.mono2.a,
  cohort.size.mono2.b,
  cohort.prob.mono2.b,
  cohort.size.combi.a,
  cohort.prob.combi.a,
  cohort.size.combi.b,
  cohort.prob.combi.b,


  esc.step.mono1.a,
  esc.step.mono2.a ,
  esc.step.mono1.b,
  esc.step.mono2.b,
  esc.step.combi.a1,
  esc.step.combi.b1,
  esc.step.combi.a2,
  esc.step.combi.b2,

  esc.constrain.mono1.a,
  esc.constrain.mono2.a,
  esc.constrain.mono1.b,
  esc.constrain.mono2.b,
  esc.constrain.combi.a1,
  esc.constrain.combi.b1,
  esc.constrain.combi.a2,
  esc.constrain.combi.b2,

  start.dose.mono1.a,
  start.dose.mono2.a,
  start.dose.mono1.b,
  start.dose.mono2.b,
  start.dose.combi.a1,
  start.dose.combi.a2,
  start.dose.combi.b1,
  start.dose.combi.b2,

  max.n.mono1.a,
  max.n.mono1.b,
  max.n.mono2.a,
  max.n.mono2.b,
  max.n.combi.a,
  max.n.combi.b,
  decision.combi.a,
  decision.combi.b,
  decision.mono1.a,
  decision.mono1.b,
  decision.mono2.a,
  decision.mono2.b,

  mtd.enforce.mono1.a,
  mtd.enforce.mono2.a,
  mtd.enforce.mono1.b,
  mtd.enforce.mono2.b,
  mtd.enforce.combi.a,
  mtd.enforce.combi.b,

  backfill.mono1.a,
  backfill.mono1.b,
  backfill.mono2.a,
  backfill.mono2.b,
  backfill.combi.a,
  backfill.combi.b,

  backfill.size.mono1.a,
  backfill.size.mono1.b,
  backfill.size.mono2.a,
  backfill.size.mono2.b,
  backfill.size.combi.a,
  backfill.size.combi.b,

  backfill.prob.mono1.a,
  backfill.prob.mono1.b,
  backfill.prob.mono2.a,
  backfill.prob.mono2.b,
  backfill.prob.combi.a,
  backfill.prob.combi.b,
  backfill.start.mono1.a = NULL,
  backfill.start.mono1.b = NULL,
  backfill.start.mono2.a = NULL,
  backfill.start.mono2.b = NULL,
  backfill.start.combi.a1 = NULL,
  backfill.start.combi.a2 = NULL,
  backfill.start.combi.b1 = NULL,
  backfill.start.combi.b2 = NULL
){
  nrow <- 7
  n.str <- c("seed", "dosing.intervals", "esc.rule", "esc.comp.max", "dose.ref1", "dose.ref2", "saturating")
  ncol <- 2
  if(esc.rule%in%c("dynamic.loss", "dynamic")){
    nrow <- nrow + 4
    ncol <- 4
    n.str <- c(n.str, "dynamic.weights[1,]", "dynamic.weights[2,]", "dynamic.weights[3,]", "dynamic.weights[4,]")
  }else if(esc.rule%in%c("loss")){
    nrow <- nrow + 1
    ncol <- 4
    n.str <- c(n.str, "loss.weights")
  }else{
    nrow <- nrow + 1
    ncol <- 2
    n.str <- c(n.str, "ewoc.threshold")
  }

  if(active.mono1.a){
    nrow <- nrow + 12 + 4
    n.str <- c(n.str, "start.dose.mono1.a", "esc.step.mono1.a", "esc.constrain.mono1.a", "max.n.mono1.a", "cohort.size.mono1.a", "cohort.prob.mono1.a",
               "mtd.decision.mono1.a$target.prob","mtd.decision.mono1.a$pat.at.mtd","mtd.decision.mono1.a$min.pat",
               "mtd.decision.mono1.a$min.dlt", "mtd.decision.mono1.a$rule", "mtd.enforce.mono1.a",
               "backfill.mono1.a", "backfill.size.mono1.a", "backfill.prob.mono1.a",
               "backfill.start.mono1.a")
    ncol <- max(ncol, length(cohort.size.mono1.a), length(backfill.size.mono1.a))
  }


  if(active.mono1.b){
    nrow <- nrow + 12 + 4
    n.str <- c(n.str, "start.dose.mono1.b", "esc.step.mono1.b","esc.constrain.mono1.b", "max.n.mono1.b", "cohort.size.mono1.b", "cohort.prob.mono1.b",
               "mtd.decision.mono1.b$target.prob","mtd.decision.mono1.b$pat.at.mtd","mtd.decision.mono1.b$min.pat",
               "mtd.decision.mono1.b$min.dlt", "mtd.decision.mono1.b$rule", "mtd.enforce.mono1.b",
               "backfill.mono1.b", "backfill.size.mono1.b", "backfill.prob.mono1.b",
               "backfill.start.mono1.b")
    ncol <- max(ncol, length(cohort.size.mono1.b), length(backfill.size.mono1.b))
  }

  if(active.mono2.a){
    nrow <- nrow + 12 + 4
    n.str <- c(n.str, "start.dose.mono2.a", "esc.step.mono2.a", "esc.constrain.mono2.a","max.n.mono2.a", "cohort.size.mono2.a", "cohort.prob.mono2.a",
               "mtd.decision.mono2.a$target.prob","mtd.decision.mono2.a$pat.at.mtd","mtd.decision.mono2.a$min.pat",
               "mtd.decision.mono2.a$min.dlt", "mtd.decision.mono2.a$rule", "mtd.enforce.mono2.a",
               "backfill.mono2.a", "backfill.size.mono2.a", "backfill.prob.mono2.a",
               "backfill.start.mono2.a")
    ncol <- max(ncol, length(cohort.size.mono2.a), length(backfill.size.mono2.a))
  }

  if(active.mono2.b){
    nrow <- nrow + 12 + 4
    n.str <- c(n.str, "start.dose.mono2.b", "esc.step.mono2.b","esc.constrain.mono2.b", "max.n.mono2.b", "cohort.size.mono2.b", "cohort.prob.mono2.b",
               "mtd.decision.mono2.b$target.prob","mtd.decision.mono2.b$pat.at.mtd","mtd.decision.mono2.b$min.pat",
               "mtd.decision.mono2.b$min.dlt", "mtd.decision.mono2.b$rule", "mtd.enforce.mono2.b",
               "backfill.mono2.b", "backfill.size.mono2.b", "backfill.prob.mono2.b",
               "backfill.start.mono2.b")
    ncol <- max(ncol, length(cohort.size.mono2.b), length(backfill.size.mono2.b))
  }

  if(active.combi.a){
    nrow <- nrow + 15 + 5
    n.str <- c(n.str, "start.dose.combi.a1","start.dose.combi.a2",
               "esc.step.combi.a1","esc.step.combi.a2",
               "esc.constrain.combi.a1","esc.constrain.combi.a2", "max.n.combi.a",
               "cohort.size.combi.a", "cohort.prob.combi.a",
               "mtd.decision.combi.a$target.prob","mtd.decision.combi.a$pat.at.mtd","mtd.decision.combi.a$min.pat",
               "mtd.decision.combi.a$min.dlt", "mtd.decision.combi.a$rule", "mtd.enforce.combi.a",
               "backfill.combi.a", "backfill.size.combi.a", "backfill.prob.combi.a",
               "backfill.start.combi.a1", "backfill.start.combi.a2")
    ncol <- max(ncol, length(cohort.size.combi.a), length(backfill.size.combi.a))
  }

  if(active.combi.b){
    nrow <- nrow + 15 + 5
    n.str <- c(n.str, "start.dose.combi.b1","start.dose.combi.b2",
               "esc.step.combi.b1","esc.step.combi.b2",
               "esc.constrain.combi.b1","esc.constrain.combi.b2","max.n.combi.b",
               "cohort.size.combi.b", "cohort.prob.combi.b",
               "mtd.decision.combi.b$target.prob","mtd.decision.combi.b$pat.at.mtd","mtd.decision.combi.b$min.pat",
               "mtd.decision.combi.b$min.dlt", "mtd.decision.combi.b$rule", "mtd.enforce.combi.b",
               "backfill.combi.b", "backfill.size.combi.b", "backfill.prob.combi.b",
               "backfill.start.combi.b1", "backfill.start.combi.b2")

    ncol <- max(ncol, length(cohort.size.combi.b), length(backfill.size.combi.b))
  }

  res <- matrix("-", nrow=nrow, ncol=ncol)
  rownames(res) <- n.str
  colnames(res) <- c("Value", rep("-", times = ncol-1))

  res[1,1] <- seed
  if(esc.rule%in%c("dynamic.loss", "dynamic")){
    res[2,1:3] <- dosing.intervals[c(1, 2, 3)]
  }else if(esc.rule%in%c("loss")){
    res[2,1:3] <- dosing.intervals[c(1, 2, 3)]
  }else{
    res[2,1:2] <- dosing.intervals[c(1, 2)]
  }
  res[3,1] <- esc.rule
  res[4,1] <- esc.comp.max
  res[5,1] <- dose.ref1
  res[6,1] <- dose.ref2
  res[7,1] <- saturating

  curr <- 8
  if(esc.rule%in%c("dynamic.loss", "dynamic")){
    res[8,1:4] <- dynamic.weights[1, ]
    res[9,1:4] <- dynamic.weights[2, ]
    res[10,1:4] <- dynamic.weights[3, ]
    res[11,1:4] <- dynamic.weights[4, ]
    curr <- 12
  }else if(esc.rule%in%c("loss")){
    res[8,1:4] <- loss.weights
    curr <- 9
  }else{
    res[8, 1] <- ewoc
    curr <- 9
  }

  if(active.mono1.a){
    len <- length(cohort.size.mono1.a)
    len2 <- length(backfill.size.mono1.a)
    res[curr, 1] <- start.dose.mono1.a
    res[curr+1, 1] <- esc.step.mono1.a
    res[curr+2, 1] <- esc.constrain.mono1.a
    res[curr+3, 1] <- max.n.mono1.a
    res[curr+4, 1:len] <- cohort.size.mono1.a
    res[curr+5, 1:len] <- cohort.prob.mono1.a
    res[curr+6, 1] <- decision.mono1.a$TARGET.PROB
    res[curr+7, 1] <- decision.mono1.a$PAT.AT.MTD
    res[curr+8, 1] <- decision.mono1.a$MIN.PAT
    res[curr+9, 1] <- decision.mono1.a$MIN.DLT
    res[curr+10, 1] <- decision.mono1.a$RULE
    res[curr+11, 1] <- mtd.enforce.mono1.a
    res[curr+12, 1] <- backfill.mono1.a
    res[curr+13, 1:len2] <- backfill.size.mono1.a
    res[curr+14, 1:len2] <- backfill.prob.mono1.a
    res[curr+15, 1] <- backfill.start.mono1.a
    curr <- curr + 12 + 4
  }

  if(active.mono1.b){
    len <- length(cohort.size.mono1.b)
    len2 <- length(backfill.size.mono1.b)
    res[curr, 1] <- start.dose.mono1.b
    res[curr+1, 1] <- esc.step.mono1.b
    res[curr+2, 1] <- esc.constrain.mono1.b
    res[curr+3, 1] <- max.n.mono1.b
    res[curr+4, 1:len] <- cohort.size.mono1.b
    res[curr+5, 1:len] <- cohort.prob.mono1.b
    res[curr+6, 1] <- decision.mono1.b$TARGET.PROB
    res[curr+7, 1] <- decision.mono1.b$PAT.AT.MTD
    res[curr+8, 1] <- decision.mono1.b$MIN.PAT
    res[curr+9, 1] <- decision.mono1.b$MIN.DLT
    res[curr+10, 1] <- decision.mono1.b$RULE
    res[curr+11, 1] <- mtd.enforce.mono1.b
    res[curr+12, 1] <- backfill.mono1.b
    res[curr+13, 1:len2] <- backfill.size.mono1.b
    res[curr+14, 1:len2] <- backfill.prob.mono1.b
    res[curr+15, 1] <- backfill.start.mono1.b
    curr <- curr + 12 + 4
  }

  if(active.mono2.a){
    len <- length(cohort.size.mono2.a)
    len2 <- length(backfill.size.mono2.a)
    res[curr, 1] <- start.dose.mono2.a
    res[curr+1, 1] <- esc.step.mono2.a
    res[curr+2, 1] <- esc.constrain.mono2.a
    res[curr+3, 1] <- max.n.mono2.a
    res[curr+4, 1:len] <- cohort.size.mono2.a
    res[curr+5, 1:len] <- cohort.prob.mono2.a
    res[curr+6, 1] <- decision.mono2.a$TARGET.PROB
    res[curr+7, 1] <- decision.mono2.a$PAT.AT.MTD
    res[curr+8, 1] <- decision.mono2.a$MIN.PAT
    res[curr+9, 1] <- decision.mono2.a$MIN.DLT
    res[curr+10, 1] <- decision.mono2.a$RULE
    res[curr+11, 1] <- mtd.enforce.mono2.a
    res[curr+12, 1] <- backfill.mono2.a
    res[curr+13, 1:len2] <- backfill.size.mono2.a
    res[curr+14, 1:len2] <- backfill.prob.mono2.a
    res[curr+15, 1] <- backfill.start.mono2.a
    curr <- curr + 12 + 4
  }

  if(active.mono2.b){
    len <- length(cohort.size.mono2.b)
    len2 <- length(backfill.size.mono2.b)
    res[curr, 1] <- start.dose.mono2.b
    res[curr+1, 1] <- esc.step.mono2.b
    res[curr+2, 1] <- esc.constrain.mono2.b
    res[curr+3, 1] <- max.n.mono2.b
    res[curr+4, 1:len] <- cohort.size.mono2.b
    res[curr+5, 1:len] <- cohort.prob.mono2.b
    res[curr+6, 1] <- decision.mono2.b$TARGET.PROB
    res[curr+7, 1] <- decision.mono2.b$PAT.AT.MTD
    res[curr+8, 1] <- decision.mono2.b$MIN.PAT
    res[curr+9, 1] <- decision.mono2.b$MIN.DLT
    res[curr+10, 1] <- decision.mono2.b$RULE
    res[curr+11, 1] <- mtd.enforce.mono2.b
    res[curr+12, 1] <- backfill.mono2.b
    res[curr+13, 1:len2] <- backfill.size.mono2.b
    res[curr+14, 1:len2] <- backfill.prob.mono2.b
    res[curr+15, 1] <- backfill.start.mono2.b
    curr <- curr + 12 + 4
  }


  if(active.combi.a){
    len <- length(cohort.size.combi.a)
    len2 <- length(backfill.size.combi.a)
    res[curr, 1] <- start.dose.combi.a1
    res[curr+1, 1] <- start.dose.combi.a2
    res[curr+2, 1] <- esc.step.combi.a1
    res[curr+3, 1] <- esc.step.combi.a2
    res[curr+4, 1] <- esc.constrain.combi.a1
    res[curr+5, 1] <- esc.constrain.combi.a2
    res[curr+6, 1] <- max.n.combi.a
    res[curr+7, 1:len] <- cohort.size.combi.a
    res[curr+8, 1:len] <- cohort.prob.combi.a
    res[curr+9, 1] <- decision.combi.a$TARGET.PROB
    res[curr+10, 1] <- decision.combi.a$PAT.AT.MTD
    res[curr+11, 1] <- decision.combi.a$MIN.PAT
    res[curr+12, 1] <- decision.combi.a$MIN.DLT
    res[curr+13, 1] <- decision.combi.a$RULE
    res[curr+14, 1] <- mtd.enforce.combi.a
    res[curr+15, 1] <- backfill.combi.a
    res[curr+16, 1:len2] <- backfill.size.combi.a
    res[curr+17, 1:len2] <- backfill.prob.combi.a
    res[curr+18, 1] <- backfill.start.combi.a1
    res[curr+19, 1] <- backfill.start.combi.a2
    curr <- curr + 15 + 5
  }

  if(active.combi.b){
    len <- length(cohort.size.combi.b)
    len2 <- length(backfill.size.combi.b)
    res[curr, 1] <- start.dose.combi.b1
    res[curr+1, 1] <- start.dose.combi.b2
    res[curr+2, 1] <- esc.step.combi.b1
    res[curr+3, 1] <- esc.step.combi.b2
    res[curr+4, 1] <- esc.constrain.combi.b1
    res[curr+5, 1] <- esc.constrain.combi.b2
    res[curr+6, 1] <- max.n.combi.b
    res[curr+7, 1:len] <- cohort.size.combi.b
    res[curr+8, 1:len] <- cohort.prob.combi.b
    res[curr+9, 1] <- decision.combi.b$TARGET.PROB
    res[curr+10, 1] <- decision.combi.b$PAT.AT.MTD
    res[curr+11, 1] <- decision.combi.b$MIN.PAT
    res[curr+12, 1] <- decision.combi.b$MIN.DLT
    res[curr+13, 1] <- decision.combi.b$RULE
    res[curr+14, 1] <- mtd.enforce.combi.b
    res[curr+15, 1] <- backfill.combi.b
    res[curr+16, 1:len2] <- backfill.size.combi.b
    res[curr+17, 1:len2] <- backfill.prob.combi.b
    res[curr+18, 1] <- backfill.start.combi.b1
    res[curr+19, 1] <- backfill.start.combi.b2
    curr <- curr + 15 + 5
  }

  return(res)

}

#'@keywords internal
 get_int_probs_BLRM <- function(smpl, ints){
    nint <- length(ints)+1
    ret <- rep(0, nint)
    for(i in 1:nint){
      if(i==1){
        ret[i] <- mean(ifelse(smpl<ints[i], 1, 0))
      }else if(i==nint){
        ret[i] <- mean(ifelse(smpl>=ints[i-1], 1, 0))
      }else{
        ret[i] <- mean(ifelse(smpl>=ints[i-1] & smpl<ints[i],1,0))
      }
    }
    return(ret)
 }


 #'@keywords internal
 # allows to write all input parameters
 # to a clean output matrix.
 # depending on activated options
 # the contents may vary slightly.
 input_out_covjointBLRM <- function(
   active.mono1.a,
   active.mono1.b,
   active.mono2.a,
   active.mono2.b,
   active.combi.a,
   active.combi.b,

   dose.ref1,
   dose.ref2,
   seed,
   dosing.intervals,
   saturating,
   ewoc,
   loss.weights,
   dynamic.weights,

   esc.rule,
   esc.comp.max,

   cohort.size.mono1.a,
   cohort.prob.mono1.a,
   cohort.size.mono1.b,
   cohort.prob.mono1.b,
   cohort.size.mono2.a,
   cohort.prob.mono2.a,
   cohort.size.mono2.b,
   cohort.prob.mono2.b,
   cohort.size.combi.a,
   cohort.prob.combi.a,
   cohort.size.combi.b,
   cohort.prob.combi.b,


   esc.step.mono1.a,
   esc.step.mono2.a ,
   esc.step.mono1.b,
   esc.step.mono2.b,
   esc.step.combi.a1,
   esc.step.combi.b1,
   esc.step.combi.a2,
   esc.step.combi.b2,

   esc.constrain.mono1.a,
   esc.constrain.mono2.a,
   esc.constrain.mono1.b,
   esc.constrain.mono2.b,
   esc.constrain.combi.a1,
   esc.constrain.combi.b1,
   esc.constrain.combi.a2,
   esc.constrain.combi.b2,

   start.dose.mono1.a,
   start.dose.mono2.a,
   start.dose.mono1.b,
   start.dose.mono2.b,
   start.dose.combi.a1,
   start.dose.combi.a2,
   start.dose.combi.b1,
   start.dose.combi.b2,

   max.n.mono1.a,
   max.n.mono1.b,
   max.n.mono2.a,
   max.n.mono2.b,
   max.n.combi.a,
   max.n.combi.b,
   decision.combi.a,
   decision.combi.b,
   decision.mono1.a,
   decision.mono1.b,
   decision.mono2.a,
   decision.mono2.b,

   mtd.enforce.mono1.a,
   mtd.enforce.mono2.a,
   mtd.enforce.mono1.b,
   mtd.enforce.mono2.b,
   mtd.enforce.combi.a,
   mtd.enforce.combi.b,

   backfill.mono1.a,
   backfill.mono1.b,
   backfill.mono2.a,
   backfill.mono2.b,
   backfill.combi.a,
   backfill.combi.b,

   backfill.size.mono1.a,
   backfill.size.mono1.b,
   backfill.size.mono2.a,
   backfill.size.mono2.b,
   backfill.size.combi.a,
   backfill.size.combi.b,

   backfill.prob.mono1.a,
   backfill.prob.mono1.b,
   backfill.prob.mono2.a,
   backfill.prob.mono2.b,
   backfill.prob.combi.a,
   backfill.prob.combi.b,
   backfill.start.mono1.a = NULL,
   backfill.start.mono1.b = NULL,
   backfill.start.mono2.a = NULL,
   backfill.start.mono2.b = NULL,
   backfill.start.combi.a1 = NULL,
   backfill.start.combi.a2 = NULL,
   backfill.start.combi.b1 = NULL,
   backfill.start.combi.b2 = NULL,
   two_sided1 = FALSE,
   two_sided2 = FALSE,
   covar.mono1.a = FALSE,
   covar.mono1.b = FALSE,
   covar.mono2.a = FALSE,
   covar.mono2.b = FALSE,
   covar.combi.a = FALSE,
   covar.combi.b = FALSE
 ){
   nrow <- 9
   n.str <- c("seed", "dosing.intervals", "esc.rule", "esc.comp.max", "dose.ref1", "dose.ref2", "saturating",
              "two_sided1", "two_sided2")
   ncol <- 2
   if(esc.rule%in%c("dynamic.loss", "dynamic")){
     nrow <- nrow + 4
     ncol <- 4
     n.str <- c(n.str, "dynamic.weights[1,]", "dynamic.weights[2,]", "dynamic.weights[3,]", "dynamic.weights[4,]")
   }else if(esc.rule%in%c("loss")){
     nrow <- nrow + 1
     ncol <- 4
     n.str <- c(n.str, "loss.weights")
   }else{
     nrow <- nrow + 1
     ncol <- 2
     n.str <- c(n.str, "ewoc.threshold")
   }

   if(active.mono1.a){
     nrow <- nrow + 12 + 5
     n.str <- c(n.str, "start.dose.mono1.a", "esc.step.mono1.a", "esc.constrain.mono1.a", "max.n.mono1.a", "cohort.size.mono1.a", "cohort.prob.mono1.a",
                "mtd.decision.mono1.a$target.prob","mtd.decision.mono1.a$pat.at.mtd","mtd.decision.mono1.a$min.pat",
                "mtd.decision.mono1.a$min.dlt", "mtd.decision.mono1.a$rule", "mtd.enforce.mono1.a",
                "backfill.mono1.a", "backfill.size.mono1.a", "backfill.prob.mono1.a",
                "backfill.start.mono1.a", "covar.mono1.a")
     ncol <- max(ncol, length(cohort.size.mono1.a), length(backfill.size.mono1.a))
   }


   if(active.mono1.b){
     nrow <- nrow + 12 + 5
     n.str <- c(n.str, "start.dose.mono1.b", "esc.step.mono1.b","esc.constrain.mono1.b", "max.n.mono1.b", "cohort.size.mono1.b", "cohort.prob.mono1.b",
                "mtd.decision.mono1.b$target.prob","mtd.decision.mono1.b$pat.at.mtd","mtd.decision.mono1.b$min.pat",
                "mtd.decision.mono1.b$min.dlt", "mtd.decision.mono1.b$rule", "mtd.enforce.mono1.b",
                "backfill.mono1.b", "backfill.size.mono1.b", "backfill.prob.mono1.b",
                "backfill.start.mono1.b", "covar.mono1.b")
     ncol <- max(ncol, length(cohort.size.mono1.b), length(backfill.size.mono1.b))
   }

   if(active.mono2.a){
     nrow <- nrow + 12 + 5
     n.str <- c(n.str, "start.dose.mono2.a", "esc.step.mono2.a", "esc.constrain.mono2.a","max.n.mono2.a", "cohort.size.mono2.a", "cohort.prob.mono2.a",
                "mtd.decision.mono2.a$target.prob","mtd.decision.mono2.a$pat.at.mtd","mtd.decision.mono2.a$min.pat",
                "mtd.decision.mono2.a$min.dlt", "mtd.decision.mono2.a$rule", "mtd.enforce.mono2.a",
                "backfill.mono2.a", "backfill.size.mono2.a", "backfill.prob.mono2.a",
                "backfill.start.mono2.a", "covar.mono2.a")
     ncol <- max(ncol, length(cohort.size.mono2.a), length(backfill.size.mono2.a))
   }

   if(active.mono2.b){
     nrow <- nrow + 12 + 5
     n.str <- c(n.str, "start.dose.mono2.b", "esc.step.mono2.b","esc.constrain.mono2.b", "max.n.mono2.b", "cohort.size.mono2.b", "cohort.prob.mono2.b",
                "mtd.decision.mono2.b$target.prob","mtd.decision.mono2.b$pat.at.mtd","mtd.decision.mono2.b$min.pat",
                "mtd.decision.mono2.b$min.dlt", "mtd.decision.mono2.b$rule", "mtd.enforce.mono2.b",
                "backfill.mono2.b", "backfill.size.mono2.b", "backfill.prob.mono2.b",
                "backfill.start.mono2.b", "covar.mono2.b")
     ncol <- max(ncol, length(cohort.size.mono2.b), length(backfill.size.mono2.b))
   }

   if(active.combi.a){
     nrow <- nrow + 15 + 6
     n.str <- c(n.str, "start.dose.combi.a1","start.dose.combi.a2",
                "esc.step.combi.a1","esc.step.combi.a2",
                "esc.constrain.combi.a1","esc.constrain.combi.a2", "max.n.combi.a",
                "cohort.size.combi.a", "cohort.prob.combi.a",
                "mtd.decision.combi.a$target.prob","mtd.decision.combi.a$pat.at.mtd","mtd.decision.combi.a$min.pat",
                "mtd.decision.combi.a$min.dlt", "mtd.decision.combi.a$rule", "mtd.enforce.combi.a",
                "backfill.combi.a", "backfill.size.combi.a", "backfill.prob.combi.a",
                "backfill.start.combi.a1", "backfill.start.combi.a2", "covar.combi.a")
     ncol <- max(ncol, length(cohort.size.combi.a), length(backfill.size.combi.a))
   }

   if(active.combi.b){
     nrow <- nrow + 15 + 6
     n.str <- c(n.str, "start.dose.combi.b1","start.dose.combi.b2",
                "esc.step.combi.b1","esc.step.combi.b2",
                "esc.constrain.combi.b1","esc.constrain.combi.b2","max.n.combi.b",
                "cohort.size.combi.b", "cohort.prob.combi.b",
                "mtd.decision.combi.b$target.prob","mtd.decision.combi.b$pat.at.mtd","mtd.decision.combi.b$min.pat",
                "mtd.decision.combi.b$min.dlt", "mtd.decision.combi.b$rule", "mtd.enforce.combi.b",
                "backfill.combi.b", "backfill.size.combi.b", "backfill.prob.combi.b",
                "backfill.start.combi.b1", "backfill.start.combi.b2", "covar.combi.b")

     ncol <- max(ncol, length(cohort.size.combi.b), length(backfill.size.combi.b))
   }

   res <- matrix("-", nrow=nrow, ncol=ncol)
   rownames(res) <- n.str
   colnames(res) <- c("Value", rep("-", times = ncol-1))

   res[1,1] <- seed
   if(esc.rule%in%c("dynamic.loss", "dynamic")){
     res[2,1:3] <- dosing.intervals[c(1, 2, 3)]
   }else if(esc.rule%in%c("loss")){
     res[2,1:3] <- dosing.intervals[c(1, 2, 3)]
   }else{
     res[2,1:2] <- dosing.intervals[c(1, 2)]
   }
   res[3,1] <- esc.rule
   res[4,1] <- esc.comp.max
   res[5,1] <- dose.ref1
   res[6,1] <- dose.ref2
   res[7,1] <- saturating

   curr <- 8
   if(esc.rule%in%c("dynamic.loss", "dynamic")){
     res[8,1:4] <- dynamic.weights[1, ]
     res[9,1:4] <- dynamic.weights[2, ]
     res[10,1:4] <- dynamic.weights[3, ]
     res[11,1:4] <- dynamic.weights[4, ]
     curr <- 12
   }else if(esc.rule%in%c("loss")){
     res[8,1:4] <- loss.weights
     curr <- 9
   }else{
     res[8, 1] <- ewoc
     curr <- 9
   }

   if(active.mono1.a){
     len <- length(cohort.size.mono1.a)
     len2 <- length(backfill.size.mono1.a)
     res[curr, 1] <- start.dose.mono1.a
     res[curr+1, 1] <- esc.step.mono1.a
     res[curr+2, 1] <- esc.constrain.mono1.a
     res[curr+3, 1] <- max.n.mono1.a
     res[curr+4, 1:len] <- cohort.size.mono1.a
     res[curr+5, 1:len] <- cohort.prob.mono1.a
     res[curr+6, 1] <- decision.mono1.a$TARGET.PROB
     res[curr+7, 1] <- decision.mono1.a$PAT.AT.MTD
     res[curr+8, 1] <- decision.mono1.a$MIN.PAT
     res[curr+9, 1] <- decision.mono1.a$MIN.DLT
     res[curr+10, 1] <- decision.mono1.a$RULE
     res[curr+11, 1] <- mtd.enforce.mono1.a
     res[curr+12, 1] <- backfill.mono1.a
     res[curr+13, 1:len2] <- backfill.size.mono1.a
     res[curr+14, 1:len2] <- backfill.prob.mono1.a
     res[curr+15, 1] <- backfill.start.mono1.a
     res[curr+16, 1] <- covar.mono1.a
     curr <- curr + 12 + 5
   }

   if(active.mono1.b){
     len <- length(cohort.size.mono1.b)
     len2 <- length(backfill.size.mono1.b)
     res[curr, 1] <- start.dose.mono1.b
     res[curr+1, 1] <- esc.step.mono1.b
     res[curr+2, 1] <- esc.constrain.mono1.b
     res[curr+3, 1] <- max.n.mono1.b
     res[curr+4, 1:len] <- cohort.size.mono1.b
     res[curr+5, 1:len] <- cohort.prob.mono1.b
     res[curr+6, 1] <- decision.mono1.b$TARGET.PROB
     res[curr+7, 1] <- decision.mono1.b$PAT.AT.MTD
     res[curr+8, 1] <- decision.mono1.b$MIN.PAT
     res[curr+9, 1] <- decision.mono1.b$MIN.DLT
     res[curr+10, 1] <- decision.mono1.b$RULE
     res[curr+11, 1] <- mtd.enforce.mono1.b
     res[curr+12, 1] <- backfill.mono1.b
     res[curr+13, 1:len2] <- backfill.size.mono1.b
     res[curr+14, 1:len2] <- backfill.prob.mono1.b
     res[curr+15, 1] <- backfill.start.mono1.b
     res[curr+16, 1] <- covar.mono1.b
     curr <- curr + 12 + 5
   }

   if(active.mono2.a){
     len <- length(cohort.size.mono2.a)
     len2 <- length(backfill.size.mono2.a)
     res[curr, 1] <- start.dose.mono2.a
     res[curr+1, 1] <- esc.step.mono2.a
     res[curr+2, 1] <- esc.constrain.mono2.a
     res[curr+3, 1] <- max.n.mono2.a
     res[curr+4, 1:len] <- cohort.size.mono2.a
     res[curr+5, 1:len] <- cohort.prob.mono2.a
     res[curr+6, 1] <- decision.mono2.a$TARGET.PROB
     res[curr+7, 1] <- decision.mono2.a$PAT.AT.MTD
     res[curr+8, 1] <- decision.mono2.a$MIN.PAT
     res[curr+9, 1] <- decision.mono2.a$MIN.DLT
     res[curr+10, 1] <- decision.mono2.a$RULE
     res[curr+11, 1] <- mtd.enforce.mono2.a
     res[curr+12, 1] <- backfill.mono2.a
     res[curr+13, 1:len2] <- backfill.size.mono2.a
     res[curr+14, 1:len2] <- backfill.prob.mono2.a
     res[curr+15, 1] <- backfill.start.mono2.a
     res[curr+16, 1] <- covar.mono2.a
     curr <- curr + 12 + 5
   }

   if(active.mono2.b){
     len <- length(cohort.size.mono2.b)
     len2 <- length(backfill.size.mono2.b)
     res[curr, 1] <- start.dose.mono2.b
     res[curr+1, 1] <- esc.step.mono2.b
     res[curr+2, 1] <- esc.constrain.mono2.b
     res[curr+3, 1] <- max.n.mono2.b
     res[curr+4, 1:len] <- cohort.size.mono2.b
     res[curr+5, 1:len] <- cohort.prob.mono2.b
     res[curr+6, 1] <- decision.mono2.b$TARGET.PROB
     res[curr+7, 1] <- decision.mono2.b$PAT.AT.MTD
     res[curr+8, 1] <- decision.mono2.b$MIN.PAT
     res[curr+9, 1] <- decision.mono2.b$MIN.DLT
     res[curr+10, 1] <- decision.mono2.b$RULE
     res[curr+11, 1] <- mtd.enforce.mono2.b
     res[curr+12, 1] <- backfill.mono2.b
     res[curr+13, 1:len2] <- backfill.size.mono2.b
     res[curr+14, 1:len2] <- backfill.prob.mono2.b
     res[curr+15, 1] <- backfill.start.mono2.b
     res[curr+16, 1] <- covar.mono2.b
     curr <- curr + 12 + 5
   }


   if(active.combi.a){
     len <- length(cohort.size.combi.a)
     len2 <- length(backfill.size.combi.a)
     res[curr, 1] <- start.dose.combi.a1
     res[curr+1, 1] <- start.dose.combi.a2
     res[curr+2, 1] <- esc.step.combi.a1
     res[curr+3, 1] <- esc.step.combi.a2
     res[curr+4, 1] <- esc.constrain.combi.a1
     res[curr+5, 1] <- esc.constrain.combi.a2
     res[curr+6, 1] <- max.n.combi.a
     res[curr+7, 1:len] <- cohort.size.combi.a
     res[curr+8, 1:len] <- cohort.prob.combi.a
     res[curr+9, 1] <- decision.combi.a$TARGET.PROB
     res[curr+10, 1] <- decision.combi.a$PAT.AT.MTD
     res[curr+11, 1] <- decision.combi.a$MIN.PAT
     res[curr+12, 1] <- decision.combi.a$MIN.DLT
     res[curr+13, 1] <- decision.combi.a$RULE
     res[curr+14, 1] <- mtd.enforce.combi.a
     res[curr+15, 1] <- backfill.combi.a
     res[curr+16, 1:len2] <- backfill.size.combi.a
     res[curr+17, 1:len2] <- backfill.prob.combi.a
     res[curr+18, 1] <- backfill.start.combi.a1
     res[curr+19, 1] <- backfill.start.combi.a2
     res[curr+20, 1] <- covar.combi.a
     curr <- curr + 15 + 6
   }

   if(active.combi.b){
     len <- length(cohort.size.combi.b)
     len2 <- length(backfill.size.combi.b)
     res[curr, 1] <- start.dose.combi.b1
     res[curr+1, 1] <- start.dose.combi.b2
     res[curr+2, 1] <- esc.step.combi.b1
     res[curr+3, 1] <- esc.step.combi.b2
     res[curr+4, 1] <- esc.constrain.combi.b1
     res[curr+5, 1] <- esc.constrain.combi.b2
     res[curr+6, 1] <- max.n.combi.b
     res[curr+7, 1:len] <- cohort.size.combi.b
     res[curr+8, 1:len] <- cohort.prob.combi.b
     res[curr+9, 1] <- decision.combi.b$TARGET.PROB
     res[curr+10, 1] <- decision.combi.b$PAT.AT.MTD
     res[curr+11, 1] <- decision.combi.b$MIN.PAT
     res[curr+12, 1] <- decision.combi.b$MIN.DLT
     res[curr+13, 1] <- decision.combi.b$RULE
     res[curr+14, 1] <- mtd.enforce.combi.b
     res[curr+15, 1] <- backfill.combi.b
     res[curr+16, 1:len2] <- backfill.size.combi.b
     res[curr+17, 1:len2] <- backfill.prob.combi.b
     res[curr+18, 1] <- backfill.start.combi.b1
     res[curr+19, 1] <- backfill.start.combi.b2
     res[curr+20, 1] <- covar.combi.b
     curr <- curr + 15 + 6
   }

   return(res)

 }

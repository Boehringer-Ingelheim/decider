#---------------------------------------------------------------------------
#Utility function for deciding on the next dose level
#---------------------------------------------------------------------------
#'@keywords internal

get_ewoc_dec_jointBLRM <- function(
  tox,
  curr.dose1,
  curr.dose2,
  type,
  mode,
  ewoc,
  comp.max,
  esc.step1,
  esc.step2,
  max.next1=NULL,
  max.next2=NULL
){

  #read out available dose levels
  ndose <- length(tox[,1])
  dnames <- strsplit(rownames(tox), split="+", fixed = TRUE)
  d1 <- vector("numeric", ndose)
  d2 <- vector("numeric", ndose)
  for(d in 1:ndose){
    d1[d] <- as.numeric(dnames[[d]][1])
    d2[d] <- as.numeric(dnames[[d]][2])
  }

  #check whether at least some doses satisfy ewoc
  idx.d.ewoc <- which(tox[, 3]<ewoc)
  if(length(idx.d.ewoc)==0){
    return(c(0,0))
  }

  #doses that fulfill ewoc principle
  d1.ew <- d1[idx.d.ewoc]
  d2.ew <- d2[idx.d.ewoc]
  if(mode=="opt"){
    ptar.ew <- tox[idx.d.ewoc,2]
  }
  #to avoid that the escalation step does not work as expected
  eps.tol <- .Machine$double.eps

  if(type=="mono1"){

    #perform escalation for mono component 1
    if(is.null(max.next1)){
      idx.d1.poss <- which((d1.ew-eps.tol)<=(curr.dose1*esc.step1))
    }else{
      idx.d1.poss <- which((d1.ew-eps.tol)<=(max.next1))
    }
    if(length(idx.d1.poss)==0){
      return(c(0,0))
    }else if(length(idx.d1.poss)==1){
      return(c(d1.ew[idx.d1.poss], 0))
    }else{
      if(mode=="max"){
        return(c(max(d1.ew[idx.d1.poss]), 0))
      }else{
        d1.rem <- d1.ew[idx.d1.poss]
        maxptar <- max(ptar.ew[idx.d1.poss])
        idx.ptar <- which(ptar.ew[idx.d1.poss]==maxptar)
        return(c(max(d1.rem[idx.ptar]), 0))
      }
    }

  }else if(type=="mono2"){

    #perform escalation for mono component 2
    if(is.null(max.next2)){
      idx.d2.poss <- which((d2.ew-eps.tol)<=(curr.dose2*esc.step2))
    }else{
      idx.d2.poss <- which((d2.ew-eps.tol)<=(max.next2))
    }
    if(length(idx.d2.poss)==0){
      return(c(0,0))
    }else if(length(idx.d2.poss)==1){
      return(c(0, d2.ew[idx.d2.poss]))
    }else{
      if(mode=="max"){
        return(c(0, max(d2.ew[idx.d2.poss])))
      }else{
        d2.rem <- d2.ew[idx.d2.poss]
        maxptar <- max(ptar.ew[idx.d2.poss])
        idx.ptar <- which(ptar.ew[idx.d2.poss]==maxptar)
        return(c(0, max(d2.rem[idx.ptar])))
      }
    }

  }else{

    #perform combination therapy escalation
    if(is.null(max.next1) & is.null(max.next2)){
      idx.poss <- which( ((d1.ew-eps.tol)<=(curr.dose1*esc.step1))
                        & ((d2.ew-eps.tol)<=(curr.dose2*esc.step2)))
    }else if (is.null(max.next1) & !is.null(max.next2)){
      idx.poss <- which( ((d1.ew-eps.tol)<=(curr.dose1*esc.step1))
                         & ((d2.ew-eps.tol)<=(max.next2)))
    }else if (!is.null(max.next1) & is.null(max.next2)){
      idx.poss <- which( ((d1.ew-eps.tol)<=(max.next1))
                         & ((d2.ew-eps.tol)<=(curr.dose2*esc.step2)))
    }else{
      idx.poss <- which( ((d1.ew-eps.tol)<=(max.next1))
                         & ((d2.ew-eps.tol)<=(max.next2)))
    }

    if(length(idx.poss)==0){
      return(c(0,0))
    }else if(length(idx.poss)==1){
      return(c(d1.ew[idx.poss],
               d2.ew[idx.poss]))
    }else{
      #more than 1 dose combi remains, so, additional rules are applied.
      d1.rem <- d1.ew[idx.poss]
      d2.rem <- d2.ew[idx.poss]

      #decision rule is depending on mode:  opt or max
      if(mode=="max"){
        #for max, maximize both components in the given order
        if(comp.max == 1){
          md1 <- max(d1.rem)
          md2 <- max(d2.rem[which(d1.rem==md1)])
          return(c(md1, md2))
        }else{
          md2 <- max(d2.rem)
          md1 <- max(d1.rem[which(d2.rem==md2)])
          return(c(md1, md2))
        }
      }else{
        #mode is opt
        ptar.rem <- ptar.ew[idx.poss]
        mtar <- max(ptar.rem)
        idx_mtar <- which(ptar.rem==mtar)
        if(length(idx_mtar)==1){
          return(c(d1.rem[idx_mtar],
                   d2.rem[idx_mtar]))
        }else{
          d1mtar <- d1.rem[idx_mtar]
          d2mtar <- d2.rem[idx_mtar]
          if(comp.max==1){
            md1 <- max(d1mtar)
            idx_md1 <- which(d1mtar==md1)
            if(length(idx_md1)==1){
              return(c(d1mtar[idx_md1],
                       d2mtar[idx_md1]))
            }else{
              md2 <- max(d2mtar[idx_md1])
              return(c(md1, md2))
            }
          }else{
            md2 <- max(d2mtar)
            idx_md2 <- which(d2mtar==md2)
            if(length(idx_md2)==1){
              return(c(d1mtar[idx_md2],
                       d2mtar[idx_md2]))
            }else{
              md1 <- max(d1mtar[idx_md2])
              return(c(md1, md2))
            }
          }

        }
      }

    }
  }

}


#-------------------------------------------------------------------------------
#function for loss escalation decisions
#-------------------------------------------------------------------------------

get_loss_dec_jointBLRM <- function(
  tox,
  curr.dose1,
  curr.dose2,
  type,
  loss.weights,
  comp.max,
  esc.step1,
  esc.step2,
  max.next1=NULL,
  max.next2=NULL
){

  #read out available dose levels
  ndose <- length(tox[,1])
  dnames <- strsplit(rownames(tox), split="+", fixed = TRUE)
  d1 <- vector("numeric", ndose)
  d2 <- vector("numeric", ndose)
  for(d in 1:ndose){
    d1[d] <- as.numeric(dnames[[d]][1])
    d2[d] <- as.numeric(dnames[[d]][2])
  }
  #to avoid that the escalation step does not work as expected
  eps.tol <- .Machine$double.eps

  if(type=="mono1"){
    #-------------------
    #mono 1 escalation
    #-------------------
    if(is.null(max.next1)){
      idx_poss <- which((d1-eps.tol)<=(curr.dose1*esc.step1))
    }else{
      idx_poss <- which((d1-eps.tol)<=(max.next1))
    }
    d1.poss <- d1[idx_poss]
    bayes.risk <- vector("numeric", length = length(idx_poss))
    for(level in 1:length(idx_poss)){
      bayes.risk[level] <- sum(tox[idx_poss[level], ]*loss.weights)
    }
    minbr <- min(bayes.risk)
    d1min <- d1.poss[which(bayes.risk==minbr)]
    return(c(max(d1min), 0))

  }else if(type=="mono2"){
    #-------------------
    #mono 2 escalation
    #-------------------
    if(is.null(max.next2)){
      idx_poss <- which((d2-eps.tol)<=(curr.dose2*esc.step2))
    }else{
      idx_poss <- which((d2-eps.tol)<=(max.next2))
    }
    d2.poss <- d2[idx_poss]
    bayes.risk <- vector("numeric", length = length(idx_poss))
    for(level in 1:length(idx_poss)){
      bayes.risk[level] <- sum(tox[idx_poss[level], ]*loss.weights)
    }
    minbr <- min(bayes.risk)
    d2min <- d2.poss[which(bayes.risk==minbr)]
    return(c(0, max(d2min)))

  }else{
    #-------------------
    #combination escalation
    #-------------------
    if(is.null(max.next1) & is.null(max.next2)){
      idx_poss <- which( ((d1-eps.tol)<=(curr.dose1*esc.step1))
                         & ((d2-eps.tol)<=(curr.dose2*esc.step2)))
    }else if (is.null(max.next1) & !is.null(max.next2)){
      idx_poss <- which( ((d1-eps.tol)<=(curr.dose1*esc.step1))
                         & ((d2-eps.tol)<=(max.next2)))
    }else if (!is.null(max.next1) & is.null(max.next2)){
      idx_poss <- which( ((d1-eps.tol)<=(max.next1))
                         & ((d2-eps.tol)<=(curr.dose2*esc.step2)))
    }else{
      idx_poss <- which( ((d1-eps.tol)<=(max.next1))
                         & ((d2-eps.tol)<=(max.next2)))
    }

    d1.poss <- d1[idx_poss]
    d2.poss <- d2[idx_poss]

    bayes.risk <- vector("numeric", length = length(idx_poss))
    for(level in 1:length(idx_poss)){
      bayes.risk[level] <- sum(tox[idx_poss[level], ]*loss.weights)
    }
    minbr <- min(bayes.risk)
    idx_min <-which(bayes.risk==minbr)
    if(length(idx_min)==1){
      return(c(d1.poss[idx_min],
               d2.poss[idx_min]))
    }else{
      d1.minbr <- d1.poss[idx_min]
      d2.minbr <- d2.poss[idx_min]
      if(comp.max==1){
        md1 <- max(d1.minbr)
        d2.md1 <- d2.minbr[which(d1.minbr==md1)]
        return(c(md1, max(d2.md1)))
      }else{
        md2 <- max(d2.minbr)
        d1.md2 <- d1.minbr[which(d2.minbr==md2)]
        return(c(max(d1.md2), md2))
      }
    }
  }

}

#-------------------------------------------------------------------------------
#function for dynamic loss escalation decisions
#-------------------------------------------------------------------------------

get_dloss_dec_jointBLRM <- function(
  tox,
  curr.dose1,
  curr.dose2,
  type,
  dynamic.weights,
  refprobs,
  comp.max,
  esc.step1,
  esc.step2,
  max.next1=NULL,
  max.next2=NULL
){
  #only compute dynamic weights and call usual loss escalation.

  #get interval probabilities of the reference dose
  if(type=="mono1"){
    ip_ref <- refprobs[1, ]
  }else if(type=="mono2"){
    ip_ref <- refprobs[2, ]
  }else{
    ip_ref <- refprobs[3, ]
  }

  l_ast <- ip_ref[1]*dynamic.weights[1,] + ip_ref[2]*dynamic.weights[2,] +
    ip_ref[3]*dynamic.weights[3,] + ip_ref[4]*dynamic.weights[4,]

  return(get_loss_dec_jointBLRM(
    tox=tox,
    curr.dose1=curr.dose1,
    curr.dose2=curr.dose2,
    type=type,
    loss.weights=l_ast,
    comp.max=comp.max,
    esc.step1=esc.step1,
    esc.step2=esc.step2,
    max.next1=max.next1,
    max.next2=max.next2
  ))

}



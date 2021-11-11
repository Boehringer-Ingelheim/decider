# --------------------------
# function: main_jointBLRM
# --------------------------
#This is the main routine for joint BLRM trials:
#The function is responsible for sampling from the model and
#determines the escalation decision.

#To avoid copying large objects, sampling itself is performed via the
#function post_tox_jointBLRM_sim, which only returns the interval probabilities.
#Afterwards, these are used to obtain the decision according to the specified
#escalation rule, which is then returned by the main function.


#'@keywords internal
main_jointBLRM <- function(dose1,
                           dose2,
                           n.pat,
                           n.dlt,
                           study,
                           prior.tau,
                           prior.mu,
                           saturating = FALSE,
                           seed = as.numeric(Sys.time()),
                           file.name,
                           working.path,

                           dose.ref1 = NULL,
                           dose.ref2 = NULL,
                           study.interest,
                           type.interest = NULL,
                           loss.weights = c(1, 0, 1, 2),
                           dynamic.weights = NULL,
                           esc.rule,
                           esc.comp.max = 1,
                           esc.step1,
                           esc.step2,
                           max.next1 = NULL,
                           max.next2 = NULL,
                           dose1.interest,
                           dose2.interest,
                           curr.dose1,
                           curr.dose2,

                           #argument to check whether current dose
                           #can be back-filled in next iteration
                           check.prev.dose = FALSE,

                           dosing.intervals = c(0.16, 0.33, 0.6),
                           ewoc = 0.25,
                           chains,
                           iter,
                           refresh,
                           warmup,
                           adapt_delta,
                           max_treedepth
                           ) {

  #checks whether the posterior tox of the reference dose need to be
  #computed by post_tox_jointBLRM
  if(tolower(esc.rule)=="dynamic.loss"){
    d.loss <- TRUE
  }else{
    d.loss <- FALSE
  }

  #compute interval probabilities of doses of interest
  post_tox <- post_tox_jointBLRM_sim(
                                 study.interest=study.interest,
                                 type.interest=type.interest,
                                 dose1.interest=dose1.interest,
                                 dose2.interest=dose2.interest,
                                 dose1=dose1,
                                 dose2=dose2,

                                 file.name = file.name,
                                 working.path = working.path,
                                 dose.ref1=dose.ref1,
                                 dose.ref2=dose.ref2,
                                 n.pat=n.pat,
                                 n.dlt=n.dlt,
                                 n.study=study,
                                 d.loss=d.loss,
                                 dosing.intervals = dosing.intervals,
                                 prior.mu = prior.mu,
                                 prior.tau = prior.tau,
                                 saturating = saturating,
                                 iter = iter,
                                 adapt_delta = adapt_delta,
                                 max_treedepth = max_treedepth,
                                 warmup = warmup,
                                 refresh = refresh,
                                 chains = chains,
                                 seed=sample.int(.Machine$integer.max, 1))

  #to extract target probability of current dose (for MTD decision)
  curr.dose.name <- paste( curr.dose1, curr.dose2, sep="+")

  #check given escalation rule and call function to make the decision.
  if(tolower(esc.rule)=="ewoc"|
     tolower(esc.rule)=="ewoc.max"|
     tolower(esc.rule)=="ewoc.opt"){
    #---------------------
    #EWOC escalation
    #---------------------

    #EWOC can be used in modes opt and max.
    #opt is the default, as it is better when handling combination escalation
    if(tolower(esc.rule)=="ewoc"|tolower(esc.rule)=="ewoc.opt"){
      mode <- "opt"
    }else{
      mode <- "max"
    }

    #current target prob
    curr.ptar <- post_tox[curr.dose.name, 2]
    #get next dose
    next.d <- get_ewoc_dec_jointBLRM(
      tox = post_tox,
      curr.dose1=curr.dose1,
      curr.dose2=curr.dose2,
      type = type.interest,
      mode = mode,
      ewoc=ewoc,
      comp.max=esc.comp.max,
      esc.step1 = esc.step1,
      esc.step2 = esc.step2,
      max.next1 = max.next1,
      max.next2 = max.next2)

    #check whether current dose is still okay
    #and could in theory be back-filled in
    #next iteration
    if(check.prev.dose){
      bf.check <- (post_tox[curr.dose.name, 3] < ewoc)
    }else{
      bf.check <- FALSE
    }

    #return decision and current target probability
    return(list(
      "next.d"= next.d,
      "curr.ptar" = curr.ptar,
      "bf.check" = bf.check)
    )
  }else if(tolower(esc.rule)=="loss"){
    #---------------------
    #Loss escalation
    #---------------------

    #current target prob
    curr.ptar <- post_tox[curr.dose.name, 2]
    #get next dose
    next.d <- get_loss_dec_jointBLRM(
      tox=post_tox,
      curr.dose1=curr.dose1,
      curr.dose2=curr.dose2,
      type = type.interest,
      loss.weights=loss.weights,
      comp.max=esc.comp.max,
      esc.step1=esc.step1,
      esc.step2=esc.step2,
      max.next1 = max.next1,
      max.next2 = max.next2)

    #check whether current dose is still okay
    #and could in theory be back-filled in
    #next iteration
    if(check.prev.dose){
      bf.check <- (curr.dose1<next.d[1] & curr.dose2<=next.d[2]+2*.Machine$double.eps) |
        (curr.dose1<=next.d[1]+2*.Machine$double.eps & curr.dose2<next.d[2])
    }else{
      bf.check <- FALSE
    }

    #return decision and current target probability
    return(list(
      "next.d" = next.d,
      "curr.ptar"=curr.ptar,
      "bf.check" = bf.check))

  }else{
    #------------------------
    #dynamic loss escalation
    #------------------------
    #current target prob
    curr.ptar <- post_tox$result[curr.dose.name, 2]

    #get next dose
    next.d <- get_dloss_dec_jointBLRM(
      tox=post_tox$result,
      curr.dose1=curr.dose1,
      curr.dose2=curr.dose2,
      type = type.interest,
      dynamic.weights=dynamic.weights,
      refprobs=post_tox$refprobs,
      comp.max=esc.comp.max,
      esc.step1=esc.step1,
      esc.step2=esc.step2,
      max.next1 = max.next1,
      max.next2 = max.next2)

    #check whether current dose is still okay
    #and could in theory be back-filled in
    #next iteration
    if(check.prev.dose){
      bf.check <- (curr.dose1<next.d[1] & curr.dose2<=next.d[2]+2*.Machine$double.eps) |
        (curr.dose1<=next.d[1]+2*.Machine$double.eps & curr.dose2<next.d[2])
    }else{
      bf.check <- FALSE
    }

    #return decision and current target probability
    return(list(
      "next.d" = next.d,
      "curr.ptar" = curr.ptar,
      "bf.check" = bf.check))

  }

}


#----------------------------------------------------------------------------------------
#END OF FILE        END OF FILE           END OF FILE           END OF FILE     END OF FILE
#----------------------------------------------------------------------------------------

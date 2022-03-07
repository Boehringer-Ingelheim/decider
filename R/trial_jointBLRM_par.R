#----------------------------------------
#function: trial_jointBLRM
#----------------------------------------
#'@keywords internal
trial_jointBLRM_par<- function(          doses.mono1.a = c(0),
                                     doses.mono2.a = c(1),
                                     doses.combi.a = rbind(c(1), c(1)),
                                     doses.mono1.b = c(0),
                                     doses.mono2.b = c(1),
                                     doses.combi.b = rbind(c(1), c(1)),
                                     start.dose.mono1.a= NULL,
                                     start.dose.mono2.a= NULL,
                                     start.dose.mono1.b= NULL,
                                     start.dose.mono2.b= NULL,
                                     start.dose.combi.a1= NULL,
                                     start.dose.combi.a2= NULL,
                                     start.dose.combi.b1 = NULL,
                                     start.dose.combi.b2 = NULL,
                                     #seed = as.numeric(Sys.time()),

                                     historical.data= NULL,

                                     active.mono1.a = FALSE,
                                     active.mono1.b = FALSE,
                                     active.mono2.a = FALSE,
                                     active.mono2.b = FALSE,
                                     active.combi.a = FALSE,
                                     active.combi.b = FALSE,

                                     cohort.queue = rep(c(1,2,3,4,5,6), times = 2000),

                                     tox.combi.a = NULL,
                                     tox.mono1.a = NULL,
                                     tox.mono2.a = NULL,
                                     tox.combi.b = NULL,
                                     tox.mono1.b = NULL,
                                     tox.mono2.b = NULL,

                                     cohort.size.mono1.a = c(3),
                                     cohort.prob.mono1.a = c(1),
                                     cohort.size.mono2.a = c(3),
                                     cohort.prob.mono2.a = c(1),
                                     cohort.size.mono1.b = c(3),
                                     cohort.prob.mono1.b = c(1),
                                     cohort.size.mono2.b = c(3),
                                     cohort.prob.mono2.b = c(1),
                                     cohort.size.combi.a = c(3),
                                     cohort.prob.combi.a = c(1),
                                     cohort.size.combi.b = c(3),
                                     cohort.prob.combi.b = c(1),

                                    esc.rule,
                                    esc.comp.max,

                                     esc.step.mono1.a,
                                     esc.step.mono2.a,
                                     esc.step.mono1.b,
                                     esc.step.mono2.b,
                                     esc.step.combi.a1,
                                     esc.step.combi.b1,
                                     esc.step.combi.a2,
                                     esc.step.combi.b2,

                                    esc.constrain.mono1.a=FALSE,
                                    esc.constrain.mono2.a=FALSE,
                                    esc.constrain.mono1.b=FALSE,
                                    esc.constrain.mono2.b=FALSE,
                                    esc.constrain.combi.a1=FALSE,
                                    esc.constrain.combi.b1=FALSE,
                                    esc.constrain.combi.a2=FALSE,
                                    esc.constrain.combi.b2=FALSE,


                                    prior.tau = list( "tau_a1" = c( 0.25,log(4)/1.96),
                                                   "tau_b1" = c( 0.125,log(4)/1.96),
                                                   "tau_a2" = c( 0.25,log(4)/1.96),
                                                   "tau_b2" = c( 0.125,log(4)/1.96),
                                                   "tau_eta" = c( 0.125,log(4)/1.96)
                                     ),

                                     prior.mu = list(
                                         "mu_a1" = c(-0.6931472, 3),
                                         "mu_b1" = c(0, 1),
                                         "mu_a2" = c(-0.6931472, 3),
                                         "mu_b2" = c(0, 1),
                                         "mu_eta" = c(0, 1.121)
                                     ),
                                    saturating = FALSE,

                                     dose.ref1 = NULL,
                                     dose.ref2 = NULL,

                                     dosing.intervals = c(0.16, 0.33, 0.6),
                                     loss.weights = c(1, 0, 1, 2),
                                     dynamic.weights = NULL,
                                     ewoc = 0.25,

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

                                    working.path,
                                    file.name,
                                     chains = 3,
                                     iter = 5000,
                                     refresh = 0,
                                     adapt_delta = 0.8,
                                     warmup = 2500,
                                     max_treedepth = 14
){


    #set.seed(seed)

    #--------------------------------------
    #preparations to start the actual trial
    #--------------------------------------
    #reset the booleans for the stopping conditons
    combi.is.MTD.1 <- FALSE
    mono.1.is.MTD <- FALSE
    mono.2.is.MTD <- FALSE
    combi.stopped.1 <- FALSE
    mono.1.stopped <- FALSE
    mono.2.stopped <- FALSE
    combi.is.MTD.2 <- FALSE
    mono.4.is.MTD <- FALSE
    mono.5.is.MTD <- FALSE
    combi.stopped.2 <- FALSE
    mono.4.stopped <- FALSE
    mono.5.stopped <- FALSE

    #ensure named tox vectors
    if(active.mono1.a){
        names(tox.mono1.a) <- paste0(doses.mono1.a)
    }

    if(active.mono2.a){
        names(tox.mono2.a) <- paste0(doses.mono2.a)
    }

    if(active.mono1.b){
        names(tox.mono1.b) <- paste0(doses.mono1.b)
    }

    if(active.mono2.b){
        names(tox.mono2.b) <- paste0(doses.mono2.b)
    }

    if(active.combi.a){
        names(tox.combi.a) <- paste(doses.combi.a[1,], doses.combi.a[2,], sep = "+")
    }

    if(active.combi.b){
        names(tox.combi.b) <- paste(doses.combi.b[1,], doses.combi.b[2,], sep = "+")
    }

    if(esc.rule%in% c("ewoc", "ewoc.opt", "ewoc.max")){
        dosing.intervals <- dosing.intervals[1:2]
    }

    #reset various counters for DLTS and patients at each trial/dose
    n.pat <- c(0,0,0,0,0,0) #patients observed in each trial
    n.dlt <- c(0,0,0,0,0,0) #total number of DLTs per trial

    #counters for how much patients were treated with each dose. (for the MTD decision)
    patients.at.dose.combi.1 <- rep(0, times = length(tox.combi.a))
    names(patients.at.dose.combi.1) <- paste(doses.combi.a[1,],doses.combi.a[2,], sep= "+")
    patients.at.dose.mono.1 <- rep(0, times = length(doses.mono1.a))
    names(patients.at.dose.mono.1) <- doses.mono1.a
    patients.at.dose.mono.2 <- rep(0, times = length(doses.mono2.a))
    names(patients.at.dose.mono.2) <- doses.mono2.a
    patients.at.dose.combi.2 <- rep(0, times = length(tox.combi.b))
    names(patients.at.dose.combi.2) <- paste(doses.combi.b[1,],doses.combi.b[2,], sep= "+")
    patients.at.dose.mono.4 <- rep(0, times = length(doses.mono1.b))
    names(patients.at.dose.mono.4) <- doses.mono1.b
    patients.at.dose.mono.5 <- rep(0, times = length(doses.mono2.b))
    names(patients.at.dose.mono.5) <- doses.mono2.b


    #reset the vars for current doses
    current.dose.mono.2 <- start.dose.mono2.a
    current.dose.mono.1 <- start.dose.mono1.a
    current.dose.1.combi.1 <- start.dose.combi.a1
    current.dose.2.combi.1 <- start.dose.combi.a2
    current.dose.mono.5 <- start.dose.mono2.b
    current.dose.mono.4 <- start.dose.mono1.b
    current.dose.1.combi.2 <- start.dose.combi.b1
    current.dose.2.combi.2 <- start.dose.combi.b2

    #initialize available dose levels and constrains if activated.
    if(esc.constrain.mono1.a & active.mono1.a){
        lv.mono1.a <- sort(as.numeric(unique(doses.mono1.a)))
        nlv.mono1.a <- length(lv.mono1.a)
        lvid.mono1.a <- which(lv.mono1.a == current.dose.mono.1)
        if(lvid.mono1.a==nlv.mono1.a){
            max.next.mono1.a <- lv.mono1.a[lvid.mono1.a]
        }else{
            max.next.mono1.a <- lv.mono1.a[lvid.mono1.a+1]
        }
    }else{
        max.next.mono1.a <- NULL
    }

    if(esc.constrain.mono1.b & active.mono1.b){
        lv.mono1.b <- sort(as.numeric(unique(doses.mono1.b)))
        nlv.mono1.b <- length(lv.mono1.b)
        lvid.mono1.b <- which(lv.mono1.b == current.dose.mono.4)
        if(lvid.mono1.b==nlv.mono1.b){
            max.next.mono1.b <- lv.mono1.b[lvid.mono1.b]
        }else{
            max.next.mono1.b <- lv.mono1.b[lvid.mono1.b+1]
        }
    }else{
        max.next.mono1.b <- NULL
    }

    #initialize available dose levels and constrains if activated.
    if(esc.constrain.mono2.a & active.mono2.a){
        lv.mono2.a <- sort(as.numeric(unique(doses.mono2.a)))
        nlv.mono2.a <- length(lv.mono2.a)
        lvid.mono2.a <- which(lv.mono2.a == current.dose.mono.2)
        if(lvid.mono2.a==nlv.mono2.a){
            max.next.mono2.a <- lv.mono2.a[lvid.mono2.a]
        }else{
            max.next.mono2.a <- lv.mono2.a[lvid.mono2.a+1]
        }
    }else{
        max.next.mono2.a <- NULL
    }

    if(esc.constrain.mono2.b & active.mono2.b){
        lv.mono2.b <- sort(as.numeric(unique(doses.mono2.b)))
        nlv.mono2.b <- length(lv.mono2.b)
        lvid.mono2.b <- which(lv.mono2.b == current.dose.mono.5)
        if(lvid.mono2.b==nlv.mono2.b){
            max.next.mono2.b <- lv.mono2.b[lvid.mono2.b]
        }else{
            max.next.mono2.b <- lv.mono2.b[lvid.mono2.b+1]
        }
    }else{
        max.next.mono2.b <- NULL
    }

    #initialize available dose levels and constrains if activated.
    if(esc.constrain.combi.a1 & active.combi.a){
        lv.combi.a1 <- sort(as.numeric(unique(doses.combi.a[1, ])))
        nlv.combi.a1 <- length(lv.combi.a1)
        lvid.combi.a1 <- which(lv.combi.a1 == current.dose.1.combi.1)
        if(lvid.combi.a1==nlv.combi.a1){
            max.next.combi.a1 <- lv.combi.a1[lvid.combi.a1]
        }else{
            max.next.combi.a1 <- lv.combi.a1[lvid.combi.a1+1]
        }
    }else{
        max.next.combi.a1 <- NULL
    }
    if(esc.constrain.combi.a2 & active.combi.a){
        lv.combi.a2 <- sort(as.numeric(unique(doses.combi.a[2, ])))
        nlv.combi.a2 <- length(lv.combi.a2)
        lvid.combi.a2 <- which(lv.combi.a2 == current.dose.2.combi.1)
        if(lvid.combi.a2==nlv.combi.a2){
            max.next.combi.a2 <- lv.combi.a2[lvid.combi.a2]
        }else{
            max.next.combi.a2 <- lv.combi.a2[lvid.combi.a2+1]
        }
    }else{
        max.next.combi.a2 <- NULL
    }

    #initialize available dose levels and constrains if activated.
    if(esc.constrain.combi.b1 & active.combi.b){
        lv.combi.b1 <- sort(as.numeric(unique(doses.combi.b[1, ])))
        nlv.combi.b1 <- length(lv.combi.b1)
        lvid.combi.b1 <- which(lv.combi.b1 == current.dose.1.combi.2)
        if(lvid.combi.b1==nlv.combi.b1){
            max.next.combi.b1 <- lv.combi.b1[lvid.combi.b1]
        }else{
            max.next.combi.b1 <- lv.combi.b1[lvid.combi.b1+1]
        }
    }else{
        max.next.combi.b1 <- NULL
    }
    if(esc.constrain.combi.b2 & active.combi.b){
        lv.combi.b2 <- sort(as.numeric(unique(doses.combi.b[2, ])))
        nlv.combi.b2 <- length(lv.combi.b2)
        lvid.combi.b2 <- which(lv.combi.b2 == current.dose.2.combi.2)
        if(lvid.combi.b2==nlv.combi.b2){
            max.next.combi.b2 <- lv.combi.b2[lvid.combi.b2]
        }else{
            max.next.combi.b2 <- lv.combi.b2[lvid.combi.b2+1]
        }
    }else{
        max.next.combi.b2 <- NULL
    }



    #to evaluate operating characteristics
    trials.n.under <- c(0,0,0,0,0,0)
    trials.n.over <- c(0,0,0,0,0,0)
    trials.n.target <- c(0,0,0,0,0,0)
    trials.dlt.under <- c(0,0,0,0,0,0)
    trials.dlt.over <- c(0,0,0,0,0,0)
    trials.dlt.target <- c(0,0,0,0,0,0)


    #read out historical data
    #NOTE: - Also renames the studies to ensure internal naming conventions work.
    #      - Order of studies and their numbers are matched with cohort.queue
    if(is.list(historical.data)){
        #extract data from lists
        data.doses.1 <- historical.data$dose1
        tmp1 <- which(data.doses.1==0)
        data.doses.1[tmp1] <- rep(NA, length = length(tmp1))
        data.doses.2 <- historical.data$dose2
        tmp2 <- which(data.doses.2==0)
        data.doses.2[tmp2] <- rep(NA, length = length(tmp2))
        data.DLT <- historical.data$n.dlt
        data.n.patients <- historical.data$n.pat

        #changes all study names to numbers, so that
        #actually simulated trials have number greater than 6.
        data.study.raw <- harmonize_vecnames_jointBLRM(historical.data$trial)
        data.study <- rep(NA, length=length(data.study.raw))

        #studies not in the simulation get numbers 1 to (amount of these studies)
        new_st <- as.numeric(rownames(table(data.study.raw[which(data.study.raw>6)])))
        indent_std_num <- length(new_st)
        for(k in 1:indent_std_num){
            data.study[which(data.study.raw == new_st[k])] <-
                rep(k, length = length(which(data.study.raw %in% new_st[k])))
        }

        #check which of studies 1-6 appear in historical data, and process their data accordingly
        #studies without historical data are afterwards sorted according to the cohort queue
        count <- 1+indent_std_num
        temp <- c(active.mono1.a, active.mono2.a, active.combi.a,
                  active.mono1.b, active.mono2.b, active.combi.b)

        for(i in 1:6){

            #check for studies and extract given historical data

            #study1
            if(i %in% data.study.raw & i ==1){
                data.n.study.mono.1 <- count

                #include historical data to the data vectors of the studies
                indices_curr <- which(data.study.raw == 1)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[1] <- sum(data.n.patients[indices_curr])
                n.dlt[1] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h])%in%doses.mono1.a){
                        patients.at.dose.mono.1[paste0(data.doses.1[h])] <-
                            patients.at.dose.mono.1[paste0(data.doses.1[h])]+data.n.patients[h]
                        if(tox.mono1.a[paste0(data.doses.1[h])]<dosing.intervals[1]){
                            trials.n.under[1] <- trials.n.under[1] + data.n.patients[h]
                            trials.dlt.under[1] <- trials.dlt.under[1] + data.DLT[h]
                        }else if(tox.mono1.a[paste0(data.doses.1[h])]>=dosing.intervals[2]){
                            trials.n.over[1] <- trials.n.over[1] + data.n.patients[h]
                            trials.dlt.over[1] <- trials.dlt.over[1] + data.DLT[h]
                        }else{
                            trials.n.target[1] <- trials.n.target[1] + data.n.patients[h]
                            trials.dlt.target[1] <- trials.dlt.target[1] + data.DLT[h]
                        }
                    }
                }
                temp[1] <- FALSE
                count <- count+1
            }

            #study 2
            if(i %in% data.study.raw & i ==2){
                data.n.study.mono.2 <- count
                indices_curr <- which(data.study.raw == 2)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[2] <- sum(data.n.patients[indices_curr])
                n.dlt[2] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.2[h])%in%doses.mono2.a){
                        patients.at.dose.mono.2[paste0(data.doses.2[h])] <-
                            patients.at.dose.mono.2[paste0(data.doses.2[h])]+data.n.patients[h]
                        if(tox.mono2.a[paste0(data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[2] <- trials.n.under[2] + data.n.patients[h]
                            trials.dlt.under[2] <- trials.dlt.under[2] + data.DLT[h]
                        }else if(tox.mono2.a[paste0(data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[2] <- trials.n.over[2] + data.n.patients[h]
                            trials.dlt.over[2] <- trials.dlt.over[2] + data.DLT[h]
                        }else{
                            trials.n.target[2] <- trials.n.target[2] + data.n.patients[h]
                            trials.dlt.target[2] <- trials.dlt.target[2] + data.DLT[h]
                        }
                    }
                }
                temp[2] <- FALSE
                count <- count+1

            }


            #study 3
            if(i %in% data.study.raw & i ==3){
                data.n.study.combi.1 <- count
                indices_curr <- which(data.study.raw == 3)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[3] <- sum(data.n.patients[indices_curr])
                n.dlt[3] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h], "+", data.doses.2[h])
                       %in% paste0(doses.combi.a[1,], "+", doses.combi.a[2,])){
                        index_curr_d <- paste0(data.doses.1[h], "+", data.doses.2[h])
                        patients.at.dose.combi.1[index_curr_d] <- patients.at.dose.combi.1[index_curr_d]+data.n.patients[h]
                        if(tox.combi.a[paste0(data.doses.1[h], "+", data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[3] <- trials.n.under[3] + data.n.patients[h]
                            trials.dlt.under[3] <- trials.dlt.under[3] + data.DLT[h]
                        }else if(tox.combi.a[paste0(data.doses.1[h], "+", data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[3] <- trials.n.over[3] + data.n.patients[h]
                            trials.dlt.over[3] <- trials.dlt.over[3] + data.DLT[h]
                        }else{
                            trials.n.target[3] <- trials.n.target[3] + data.n.patients[h]
                            trials.dlt.target[3] <- trials.dlt.target[3] + data.DLT[h]
                        }
                    }
                }
                temp[3] <- FALSE
                count <- count+1

            }


            #study 4
            if(i %in% data.study.raw & i ==4){
                data.n.study.mono.4 <- count
                #include historical data to the data vectors of the studies
                indices_curr <- which(data.study.raw == 4)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[4] <- sum(data.n.patients[indices_curr])
                n.dlt[4] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h])%in%doses.mono1.b){
                        patients.at.dose.mono.4[paste0(data.doses.1[h])] <-
                            patients.at.dose.mono.4[paste0(data.doses.1[h])]+data.n.patients[h]
                        if(tox.mono1.b[paste0(data.doses.1[h])]<dosing.intervals[1]){
                            trials.n.under[4] <- trials.n.under[4] + data.n.patients[h]
                            trials.dlt.under[4] <- trials.dlt.under[4] + data.DLT[h]
                        }else if(tox.mono1.b[paste0(data.doses.1[h])]>=dosing.intervals[2]){
                            trials.n.over[4] <- trials.n.over[4] + data.n.patients[h]
                            trials.dlt.over[4] <- trials.dlt.over[4] + data.DLT[h]
                        }else{
                            trials.n.target[4] <- trials.n.target[4] + data.n.patients[h]
                            trials.dlt.target[4] <- trials.dlt.target[4] + data.DLT[h]
                        }
                    }
                }
                temp[4] <- FALSE
                count <- count+1

            }

            #study 5
            if(i %in% data.study.raw & i ==5){
                data.n.study.mono.5 <- count
                indices_curr <- which(data.study.raw == 5)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[5] <- sum(data.n.patients[indices_curr])
                n.dlt[5] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.2[h])%in%doses.mono2.b){
                        patients.at.dose.mono.5[paste0(data.doses.2[h])] <-
                            patients.at.dose.mono.5[paste0(data.doses.2[h])]+data.n.patients[h]
                        if(tox.mono2.b[paste0(data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[5] <- trials.n.under[5] + data.n.patients[h]
                            trials.dlt.under[5] <- trials.dlt.under[5] + data.DLT[h]
                        }else if(tox.mono2.b[paste0(data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[5] <- trials.n.over[5] + data.n.patients[h]
                            trials.dlt.over[5] <- trials.dlt.over[5] + data.DLT[h]
                        }else{
                            trials.n.target[5] <- trials.n.target[5] + data.n.patients[h]
                            trials.dlt.target[5] <- trials.dlt.target[5] + data.DLT[h]
                        }
                    }
                }
                temp[5] <- FALSE
                count <- count+1

            }


            #study 6
            if(i %in% data.study.raw & i ==6){
                data.n.study.combi.2 <- count
                indices_curr <- which(data.study.raw == 6)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[6] <- sum(data.n.patients[indices_curr])
                n.dlt[6] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h], "+", data.doses.2[h])%in%
                       paste0(doses.combi.b[1,], "+", doses.combi.b[2,])){
                        index_curr_d <- paste0(data.doses.1[h], "+", data.doses.2[h])
                        patients.at.dose.combi.2[index_curr_d] <- patients.at.dose.combi.2[index_curr_d] +
                            data.n.patients[h]
                        if(tox.combi.b[paste0(data.doses.1[h], "+", data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[6] <- trials.n.under[6] + data.n.patients[h]
                            trials.dlt.under[6] <- trials.dlt.under[6] + data.DLT[h]
                        }else if(tox.combi.b[paste0(data.doses.1[h], "+", data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[6] <- trials.n.over[6] + data.n.patients[h]
                            trials.dlt.over[6] <- trials.dlt.over[6] + data.DLT[h]
                        }else{
                            trials.n.target[6] <- trials.n.target[6] + data.n.patients[h]
                            trials.dlt.target[6] <- trials.dlt.target[6] + data.DLT[h]
                        }

                    }
                }
                temp[6] <- FALSE
                count <- count+1

            }

        }

        #sort remaining studies according to the cohort queue
        index <-1
        #loop over queue to see the order. Trials detected in the hist data
        #were already deactivated before (temp[i] is FALSE)
        while(temp[1]||temp[2]||temp[3]||temp[4]||temp[5]||temp[6]){
            if(temp[1] && cohort.queue[index]==1){
                data.n.study.mono.1 <- count

                #include historical data to the data vectors of the studies
                indices_curr <- which(data.study.raw == 1)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[1] <- sum(data.n.patients[indices_curr])
                n.dlt[1] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h])%in%doses.mono1.a){
                        patients.at.dose.mono.1[paste0(data.doses.1[h])] <-
                            patients.at.dose.mono.1[paste0(data.doses.1[h])]+data.n.patients[h]
                        if(tox.mono1.a[paste0(data.doses.1[h])]<dosing.intervals[1]){
                            trials.n.under[1] <- trials.n.under[1] + data.n.patients[h]
                            trials.dlt.under[1] <- trials.dlt.under[1] + data.DLT[h]
                        }else if(tox.mono1.a[paste0(data.doses.1[h])]>=dosing.intervals[2]){
                            trials.n.over[1] <- trials.n.over[1] + data.n.patients[h]
                            trials.dlt.over[1] <- trials.dlt.over[1] + data.DLT[h]
                        }else{
                            trials.n.target[1] <- trials.n.target[1] + data.n.patients[h]
                            trials.dlt.target[1] <- trials.dlt.target[1] + data.DLT[h]
                        }
                    }
                }
                temp[1] <- FALSE
                count <- count+1
            }
            if(temp[4] && cohort.queue[index]==4){
                data.n.study.mono.4 <- count
                #include historical data to the data vectors of the studies
                indices_curr <- which(data.study.raw == 4)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[4] <- sum(data.n.patients[indices_curr])
                n.dlt[4] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h])%in%doses.mono1.b){
                        patients.at.dose.mono.4[paste0(data.doses.1[h])] <-
                            patients.at.dose.mono.4[paste0(data.doses.1[h])]+data.n.patients[h]
                        if(tox.mono1.b[paste0(data.doses.1[h])]<dosing.intervals[1]){
                            trials.n.under[4] <- trials.n.under[4] + data.n.patients[h]
                            trials.dlt.under[4] <- trials.dlt.under[4] + data.DLT[h]
                        }else if(tox.mono1.b[paste0(data.doses.1[h])]>=dosing.intervals[2]){
                            trials.n.over[4] <- trials.n.over[4] + data.n.patients[h]
                            trials.dlt.over[4] <- trials.dlt.over[4] + data.DLT[h]
                        }else{
                            trials.n.target[4] <- trials.n.target[4] + data.n.patients[h]
                            trials.dlt.target[4] <- trials.dlt.target[4] + data.DLT[h]
                        }
                    }
                }
                temp[4] <- FALSE
                count <- count+1
            }
            if(temp[2] && cohort.queue[index]==2){
                data.n.study.mono.2 <- count
                indices_curr <- which(data.study.raw == 2)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[2] <- sum(data.n.patients[indices_curr])
                n.dlt[2] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.2[h])%in%doses.mono2.a){
                        patients.at.dose.mono.2[paste0(data.doses.2[h])] <-
                            patients.at.dose.mono.2[paste0(data.doses.2[h])]+data.n.patients[h]
                        if(tox.mono2.a[paste0(data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[2] <- trials.n.under[2] + data.n.patients[h]
                            trials.dlt.under[2] <- trials.dlt.under[2] + data.DLT[h]
                        }else if(tox.mono2.a[paste0(data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[2] <- trials.n.over[2] + data.n.patients[h]
                            trials.dlt.over[2] <- trials.dlt.over[2] + data.DLT[h]
                        }else{
                            trials.n.target[2] <- trials.n.target[2] + data.n.patients[h]
                            trials.dlt.target[2] <- trials.dlt.target[2] + data.DLT[h]
                        }
                    }
                }
                temp[2] <- FALSE
                count <- count+1
            }
            if(temp[5] && cohort.queue[index]==5){
                data.n.study.mono.5 <- count
                indices_curr <- which(data.study.raw == 5)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[5] <- sum(data.n.patients[indices_curr])
                n.dlt[5] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.2[h])%in%doses.mono2.b){
                        patients.at.dose.mono.5[paste0(data.doses.2[h])] <-
                            patients.at.dose.mono.5[paste0(data.doses.2[h])]+data.n.patients[h]
                        if(tox.mono2.b[paste0(data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[5] <- trials.n.under[5] + data.n.patients[h]
                            trials.dlt.under[5] <- trials.dlt.under[5] + data.DLT[h]
                        }else if(tox.mono2.b[paste0(data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[5] <- trials.n.over[5] + data.n.patients[h]
                            trials.dlt.over[5] <- trials.dlt.over[5] + data.DLT[h]
                        }else{
                            trials.n.target[5] <- trials.n.target[5] + data.n.patients[h]
                            trials.dlt.target[5] <- trials.dlt.target[5] + data.DLT[h]
                        }
                    }
                }
                temp[5] <- FALSE
                count <- count+1
            }
            if(temp[3] && cohort.queue[index]==3){
                data.n.study.combi.1 <- count
                indices_curr <- which(data.study.raw == 3)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[3] <- sum(data.n.patients[indices_curr])
                n.dlt[3] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h], "+", data.doses.2[h])%in%
                       paste0(doses.combi.a[1,], "+", doses.combi.a[2,])){
                        index_curr_d <-paste0(data.doses.1[h], "+", data.doses.2[h])
                        patients.at.dose.combi.1[index_curr_d] <- patients.at.dose.combi.1[index_curr_d] +
                            data.n.patients[h]
                        if(tox.combi.a[paste0(data.doses.1[h], "+", data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[3] <- trials.n.under[3] + data.n.patients[h]
                            trials.dlt.under[3] <- trials.dlt.under[3] + data.DLT[h]
                        }else if(tox.combi.a[paste0(data.doses.1[h], "+", data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[3] <- trials.n.over[3] + data.n.patients[h]
                            trials.dlt.over[3] <- trials.dlt.over[3] + data.DLT[h]
                        }else{
                            trials.n.target[3] <- trials.n.target[3] + data.n.patients[h]
                            trials.dlt.target[3] <- trials.dlt.target[3] + data.DLT[h]
                        }
                    }
                }
                temp[3] <- FALSE
                count <- count+1
            }
            if(temp[6] && cohort.queue[index]==6){
                data.n.study.combi.2 <- count
                indices_curr <- which(data.study.raw == 6)
                data.study[indices_curr] <- rep(count, length = length(indices_curr))
                n.pat[6] <- sum(data.n.patients[indices_curr])
                n.dlt[6] <- sum(data.DLT[indices_curr])

                for(h in indices_curr){
                    if(paste0(data.doses.1[h], "+", data.doses.2[h])%in%
                       paste0(doses.combi.b[1,], "+", doses.combi.b[2,])){
                        index_curr_d <-paste0(data.doses.1[h], "+", data.doses.2[h])
                        patients.at.dose.combi.2[index_curr_d] <- patients.at.dose.combi.2[index_curr_d] +
                            data.n.patients[h]
                        if(tox.combi.b[paste0(data.doses.1[h], "+", data.doses.2[h])]<dosing.intervals[1]){
                            trials.n.under[6] <- trials.n.under[6] + data.n.patients[h]
                            trials.dlt.under[6] <- trials.dlt.under[6] + data.DLT[h]
                        }else if(tox.combi.b[paste0(data.doses.1[h], "+", data.doses.2[h])]>=dosing.intervals[2]){
                            trials.n.over[6] <- trials.n.over[6] + data.n.patients[h]
                            trials.dlt.over[6] <- trials.dlt.over[6] + data.DLT[h]
                        }else{
                            trials.n.target[6] <- trials.n.target[6] + data.n.patients[h]
                            trials.dlt.target[6] <- trials.dlt.target[6] + data.DLT[h]
                        }

                    }
                }
                temp[6] <- FALSE
                count <- count+1
            }
            index <- index +1
        }

        #save number of previously included cohort from historical data
        curr.data.entry <- length(data.doses.1)

    } else{

        #No historical data was given,
        #just initialize the data vectors
        data.doses.1 <- c(0)
        data.doses.2 <- c(0)
        data.DLT <- c(0)
        data.n.patients <- c(0)
        data.study <- c(0)
        temp <- c(active.mono1.a, active.mono2.a, active.combi.a,
                  active.mono1.b, active.mono2.b, active.combi.b)
        index <-1
        count <- 1

        #Same while loop as above: Adjust numbering to order of
        #cohort queue.
        while(temp[1]||temp[2]||temp[3]||temp[4]||temp[5]||temp[6]){
            if(temp[1] && cohort.queue[index]==1){
                data.n.study.mono.1 <- count
                temp[1] <- FALSE
                count <- count+1
            }
            if(temp[4] && cohort.queue[index]==4){
                data.n.study.mono.4 <- count
                temp[4] <- FALSE
                count <- count+1
            }
            if(temp[2] && cohort.queue[index]==2){
                data.n.study.mono.2 <- count
                temp[2] <- FALSE
                count <- count+1
            }
            if(temp[5] && cohort.queue[index]==5){
                data.n.study.mono.5 <- count
                temp[5] <- FALSE
                count <- count+1
            }
            if(temp[3] && cohort.queue[index]==3){
                data.n.study.combi.1 <- count
                temp[3] <- FALSE
                count <- count+1
            }
            if(temp[6] && cohort.queue[index]==6){
                data.n.study.combi.2 <- count
                temp[6] <- FALSE
                count <- count+1
            }
            index <- index +1
        }
        #no previous data
        curr.data.entry <- c(0)
    }


    #pre-sort the vectors, and add sufficiently many zeros to ensure
    #that all data that might occur in case of 1-patient cohorts can be saved.
    prep.data.doses.1 <- c(data.doses.1[which(is.na(data.doses.2))],
                           rep(0, times = length(which(is.na(data.doses.1)))),
                           data.doses.1[which(!is.na(data.doses.1) & !is.na(data.doses.2))],
                           rep(0, times = 2*(max.n.mono1.a+max.n.mono1.b+max.n.mono2.a+max.n.mono2.b+
                                   max.n.combi.a+max.n.combi.b)+1) )
    prep.data.doses.2 <- c(rep(0, times = length(which(is.na(data.doses.2)))),
                           data.doses.2[which(is.na(data.doses.1))],
                           data.doses.2[which(!is.na(data.doses.1) & !is.na(data.doses.2))],
                           rep(0, times = 2*(max.n.mono1.a+max.n.mono1.b+max.n.mono2.a+max.n.mono2.b+
                                   max.n.combi.a+max.n.combi.b)+1) )
    prep.data.DLT <- c(data.DLT[which(is.na(data.doses.2))],
                       data.DLT[which(is.na(data.doses.1))],
                       data.DLT[which(!is.na(data.doses.1) & !is.na(data.doses.2))],
                       rep(0, times = 2*(max.n.mono1.a+max.n.mono1.b+max.n.mono2.a+max.n.mono2.b+
                               max.n.combi.a+max.n.combi.b)+1) )
    prep.data.n.patients <- c(data.n.patients[which(is.na(data.doses.2))],
                              data.n.patients[which(is.na(data.doses.1))],
                              data.n.patients[which(!is.na(data.doses.1) & !is.na(data.doses.2))],
                              rep(0, times = 2*(max.n.mono1.a+max.n.mono1.b+max.n.mono2.a+max.n.mono2.b+
                                      max.n.combi.a+max.n.combi.b)+1) )
    prep.data.study <- c(data.study[which(is.na(data.doses.2))],
                         data.study[which(is.na(data.doses.1))],
                         data.study[which(!is.na(data.doses.1) & !is.na(data.doses.2))],
                         rep(0, times = 2*(max.n.mono1.a+max.n.mono1.b+max.n.mono2.a+max.n.mono2.b+
                                 max.n.combi.a+max.n.combi.b)+1) )
    length.of.data <- length(prep.data.doses.1)

    #variables for back-fill cohorts
    if(active.mono1.a & backfill.mono1.a){
        #stores previous dose
        prev.dose.mono1.a <- 0
        #stores whether previous dose satisfies EWOC
        prev.dose.ewoc.mono1.a <- FALSE
        #stores whether previous dose is lower
        prev.dose.esc.mono1.a <- FALSE
        #stores which doses were backfilled
        backfilled.mono1.a <- rep(FALSE, times = length(doses.mono1.a))
    }
    if(active.mono1.b & backfill.mono1.b){
        #stores previous dose
        prev.dose.mono1.b <- 0
        #stores whether previous dose satisfies EWOC
        prev.dose.ewoc.mono1.b <- FALSE
        #stores whether previous dose is lower
        prev.dose.esc.mono1.b <- FALSE
        #stores which doses were backfilled
        backfilled.mono1.b <- rep(FALSE, times = length(doses.mono1.b))
    }

    if(active.mono2.a & backfill.mono2.a){
        #stores previous dose
        prev.dose.mono2.a <- 0
        #stores whether previous dose satisfies EWOC
        prev.dose.ewoc.mono2.a <- FALSE
        #stores whether previous dose is lower
        prev.dose.esc.mono2.a <- FALSE
        #stores which doses were backfilled
        backfilled.mono2.a <- rep(FALSE, times = length(doses.mono2.a))
    }
    if(active.mono2.b & backfill.mono2.b){
        #stores previous dose
        prev.dose.mono2.b <- 0
        #stores whether previous dose satisfies EWOC
        prev.dose.ewoc.mono2.b <- FALSE
        #stores whether previous dose is lower
        prev.dose.esc.mono2.b <- FALSE
        #stores which doses were backfilled
        backfilled.mono2.b <- rep(FALSE, times = length(doses.mono2.b))
    }

    if(active.combi.a & backfill.combi.a){
        #stores previous dose
        prev.dose.combi.a1 <- 0
        prev.dose.combi.a2 <- 0
        #stores whether previous dose satisfies EWOC
        prev.dose.ewoc.combi.a <- FALSE
        #stores whether previous dose is lower
        prev.dose.esc.combi.a <- FALSE
        #stores which doses were backfilled
        backfilled.combi.a <- rep(FALSE, times = length(doses.combi.a[1,]))
    }

    if(active.combi.b & backfill.combi.b){
        #stores previous dose
        prev.dose.combi.b1 <- 0
        prev.dose.combi.b2 <- 0
        #stores whether previous dose satisfies EWOC
        prev.dose.ewoc.combi.b <- FALSE
        #stores whether previous dose is lower
        prev.dose.esc.combi.b <- FALSE
        #stores which doses were backfilled
        backfilled.combi.b <- rep(FALSE, times = length(doses.combi.b[1,]))
    }

    #counter variable that reads out the type of the arriving cohorts
    cohort.counter <- 0

    #empty list for results.
    result_simulation <- list()

    #sample seed for MCMC call for each potential cohort
    mcmc_seeds <- sample.int(.Machine$integer.max, length(cohort.queue))


    #as long as one of the trials is active, has not found an MTD, has not
    #stopped due to EWOC, and has not run out of patients, simulate
    #the next cohort using the while loop.
    #The loop cannot run infinitely, as sim_jointBLRM has previously checked
    #whether sufficiently many cohorts can be enrolled in each active trial
    #to ensure that the maximal patient number can be reached.
    while( ((!combi.is.MTD.1[1]) & n.pat[3] < max.n.combi.a & !combi.stopped.1 & active.combi.a) ||
           ((!combi.is.MTD.2[1]) & n.pat[6] < max.n.combi.b & !combi.stopped.2 & active.combi.b) ||
           (!(mono.1.is.MTD[1]) & n.pat[1] < max.n.mono1.a & !mono.1.stopped & active.mono1.a)||
           (!(mono.4.is.MTD[1]) & n.pat[4] < max.n.mono1.b & !mono.4.stopped & active.mono1.b)||
           (!(mono.2.is.MTD[1]) & n.pat[2] < max.n.mono2.a & !mono.2.stopped & active.mono2.a)||
           (!(mono.5.is.MTD[1]) & n.pat[5] < max.n.mono2.b & !mono.5.stopped & active.mono2.b) ){

        #increment counter that determines the next cohort type from the queue
        cohort.counter <- cohort.counter+1

        #check which trial type the cohort belongs to and whether this trial is active
        if((cohort.queue[cohort.counter]==3) &
           (!(combi.is.MTD.1[1]) &
            n.pat[3] < max.n.combi.a &
            !combi.stopped.1 &
            active.combi.a)){
            #--------------------------------------------------------------
            #case: combi
            #--------------------------------------------------------------
            curr.data.entry <- curr.data.entry+1
            curr.dose.name <- paste( current.dose.1.combi.1, current.dose.2.combi.1, sep="+")

            #read out the current true dlt-probability
            curr.tox.combi.1 <- tox.combi.a[curr.dose.name]

            #sample the number of patients
            if(length(cohort.size.combi.a) == 1){
                curr.n.pat <- cohort.size.combi.a
            }else{
                curr.n.pat <- sample(cohort.size.combi.a, size=1, prob=cohort.prob.combi.a)
            }
            #save the number
            n.pat[3] <- n.pat[3] + curr.n.pat

            #sample the number of DLTs that are observed
            curr.dlt <- rbinom(1, curr.n.pat,  prob = curr.tox.combi.1)
            n.dlt[3] <- n.dlt[3] + curr.dlt

            #save the simulated treatment and result
            patients.at.dose.combi.1[curr.dose.name] <- curr.n.pat + patients.at.dose.combi.1[curr.dose.name]
            if(curr.tox.combi.1 < dosing.intervals[1]){  #the cohort was treated with underdose
                trials.n.under[3] <- as.numeric(trials.n.under[3]) + curr.n.pat
                trials.dlt.under[3] <- as.numeric(trials.dlt.under[3]) + curr.dlt
            } else if (dosing.intervals[1] <= curr.tox.combi.1 & curr.tox.combi.1 < dosing.intervals[2]){
                #the cohort was treated with dose in target toxicicty interval
                trials.n.target[3] <- as.numeric(trials.n.target[3]) + curr.n.pat
                trials.dlt.target[3] <- as.numeric(trials.dlt.target[3]) + curr.dlt
            } else { # the cohort was given an overdose
                trials.n.over[3] <- as.numeric(trials.n.over[3]) + curr.n.pat
                trials.dlt.over[3] <- as.numeric(trials.dlt.over[3]) + curr.dlt
            }

            #----------------------------------------------------------------------------
            # prepare data and execute the main computation step
            #----------------------------------------------------------------------------
            #update current trial data
            prep.data.doses.1[curr.data.entry] <- current.dose.1.combi.1
            prep.data.doses.2[curr.data.entry] <- current.dose.2.combi.1
            prep.data.DLT[curr.data.entry] <- curr.dlt
            prep.data.n.patients[curr.data.entry] <- curr.n.pat
            prep.data.study[curr.data.entry] <- data.n.study.combi.1

            #check for back-fill cohort and simulate it if allowed
            if(backfill.combi.a){
                #if there is a previous dose level
                #and if it is at least equal to the back-fill start dose
                if(prev.dose.combi.a1 > 0 & prev.dose.combi.a2 > 0 &
                   prev.dose.combi.a1 + 2*.Machine$double.eps>=backfill.start.combi.a1 &
                   prev.dose.combi.a2 + 2*.Machine$double.eps>=backfill.start.combi.a2){
                    #if it is actually contained in doses
                    idx.prev.dose.combi.a <- which(doses.combi.a[1, ]==prev.dose.combi.a1 &
                                                   doses.combi.a[2, ]==prev.dose.combi.a2)
                    if(!length(idx.prev.dose.combi.a)==0){
                        #if there was an escalation to reach the current dose,
                        #and if the previous dose was not backfilled.
                        if(prev.dose.esc.combi.a & prev.dose.ewoc.combi.a&
                           !backfilled.combi.a[idx.prev.dose.combi.a]){
                            #generate backfill cohort for dose
                            prev.dose.name.combi.a <- paste0(prev.dose.combi.a1, "+",
                                                             prev.dose.combi.a2)
                            prev.tox.combi.a <- tox.combi.a[prev.dose.name.combi.a]
                            if(length(backfill.size.combi.a)==1){
                                bf.pat.combi.a <- backfill.size.combi.a
                            }else{
                                bf.pat.combi.a <- sample(backfill.size.combi.a,
                                                         size=1, prob=backfill.prob.combi.a)
                            }
                            backfilled.combi.a[idx.prev.dose.combi.a] <- TRUE
                            #when there is at least 1 back-fill patient
                            if(bf.pat.combi.a>0){
                                #generate and save cohort data
                                n.pat[3] <- n.pat[3] + bf.pat.combi.a
                                bf.dlt.combi.a <- rbinom(1, bf.pat.combi.a,
                                                         prob = prev.tox.combi.a)
                                n.dlt[3] <- n.dlt[3] + bf.dlt.combi.a

                                #add to data frame
                                curr.data.entry <- curr.data.entry + 1
                                prep.data.doses.1[curr.data.entry] <- prev.dose.combi.a1
                                prep.data.doses.2[curr.data.entry] <- prev.dose.combi.a2
                                prep.data.DLT[curr.data.entry] <- bf.dlt.combi.a
                                prep.data.n.patients[curr.data.entry] <- bf.pat.combi.a
                                prep.data.study[curr.data.entry] <- data.n.study.combi.1

                                #save data
                                patients.at.dose.combi.1[prev.dose.name.combi.a] <- bf.pat.combi.a +
                                    patients.at.dose.combi.1[prev.dose.name.combi.a]
                                if(prev.tox.combi.a < dosing.intervals[1]){  #the cohort was treated with underdose
                                    trials.n.under[3] <- as.numeric(trials.n.under[3]) + bf.pat.combi.a
                                    trials.dlt.under[3] <- as.numeric(trials.dlt.under[3]) + bf.dlt.combi.a
                                } else if (dosing.intervals[1] <= prev.tox.combi.a &
                                           prev.tox.combi.a < dosing.intervals[2]){
                                    #the cohort was treated with dose in target toxicicty interval
                                    trials.n.target[3] <- as.numeric(trials.n.target[3]) + bf.pat.combi.a
                                    trials.dlt.target[3] <- as.numeric(trials.dlt.target[3]) + bf.dlt.combi.a
                                } else { # the cohort was given an overdose
                                    trials.n.over[3] <- as.numeric(trials.n.over[3]) + bf.pat.combi.a
                                    trials.dlt.over[3] <- as.numeric(trials.dlt.over[3]) + bf.dlt.combi.a
                                }

                            }
                        }
                    }
                }
            }

             #Call the main function
            current.results <- main_jointBLRM(dose1 = prep.data.doses.1[1:curr.data.entry],
                                              dose2 = prep.data.doses.2[1:curr.data.entry],
                                              n.pat = prep.data.n.patients[1:curr.data.entry],
                                              n.dlt = prep.data.DLT[1:curr.data.entry],
                                              study = prep.data.study[1:curr.data.entry],
                                              prior.tau=prior.tau,
                                              prior.mu = prior.mu,
                                              saturating = saturating,
                                              file.name = file.name,
                                              working.path = working.path,
                                              seed = mcmc_seeds[cohort.counter],

                                              dose.ref1 = dose.ref1,
                                              dose.ref2 = dose.ref2,
                                              study.interest = data.n.study.combi.1,
                                              type.interest = "combi",
                                              loss.weights = loss.weights,
                                              dynamic.weights = dynamic.weights,
                                              esc.rule = esc.rule,
                                              esc.comp.max = esc.comp.max,
                                              esc.step1 = esc.step.combi.a1,
                                              esc.step2 = esc.step.combi.a2,
                                              max.next1 = max.next.combi.a1,
                                              max.next2 = max.next.combi.a2,

                                              dose1.interest = doses.combi.a[1, ],
                                              dose2.interest = doses.combi.a[2, ],

                                              #ToDO
                                              curr.dose1=current.dose.1.combi.1,
                                              curr.dose2=current.dose.2.combi.1,

                                              check.prev.dose = backfill.combi.a,
                                              dosing.intervals = dosing.intervals,
                                              ewoc = ewoc,
                                              chains = chains,
                                              iter=iter,
                                              refresh=refresh,
                                              warmup=warmup,
                                              adapt_delta=adapt_delta,
                                              max_treedepth=max_treedepth
            )

            #obtain recommended next dose
            pot.nd <- current.results$next.d
            next.combi.1.dose.1 <- pot.nd[1]
            next.combi.1.dose.2 <- pot.nd[2]

            #check whether the trial is stopped by escalation rule
            if(all(pot.nd==c(0,0))){
                combi.stopped.1 <- TRUE
            }else{

                #otherwise, check whether MTD decision rules are fulfilled
                target.prob <- current.results$curr.ptar
                curr.n.treated <- patients.at.dose.combi.1[curr.dose.name]

                #rule=1
                if(decision.combi.a$RULE == 1){
                    #result is true or false, i.e. either the dose is declared MTD or not.
                    combi.is.MTD.1 <- (next.combi.1.dose.1==current.dose.1.combi.1) &
                        (next.combi.1.dose.2==current.dose.2.combi.1) &
                        (n.dlt[3]  >= decision.combi.a$MIN.DLT) &
                        (curr.n.treated >= decision.combi.a$PAT.AT.MTD) &
                        ((n.pat[3] >= decision.combi.a$MIN.PAT) | (target.prob >= decision.combi.a$TARGET.PROB ))

                } else { #rule = 2
                    #same as above
                    combi.is.MTD.1 <- (next.combi.1.dose.1==current.dose.1.combi.1) &
                        (next.combi.1.dose.2==current.dose.2.combi.1) &
                        (n.dlt[3]  >= decision.combi.a$MIN.DLT) &
                        (n.pat[3] >= decision.combi.a$MIN.PAT) &
                        ((curr.n.treated>= decision.combi.a$PAT.AT.MTD) | (target.prob >= decision.combi.a$TARGET.PROB ))

                }
                #update data that controls back-fill cohorts
                if(backfill.combi.a){
                    prev.dose.combi.a1 <- current.dose.1.combi.1
                    prev.dose.combi.a2 <- current.dose.2.combi.1
                    prev.dose.esc.combi.a <- ((current.dose.1.combi.1 < next.combi.1.dose.1 &
                                               current.dose.2.combi.1 <= next.combi.1.dose.2 + 2*.Machine$double.eps)|
                                              (current.dose.2.combi.1 < next.combi.1.dose.2 &
                                               current.dose.1.combi.1 <= next.combi.1.dose.1 + 2*.Machine$double.eps))
                    prev.dose.ewoc.combi.a <- current.results$bf.check
                }
                #continue with the new doses
                current.dose.1.combi.1 <- next.combi.1.dose.1
                current.dose.2.combi.1 <- next.combi.1.dose.2

                #reset levels if constrain active
                if(esc.constrain.combi.a1){
                    lvid.combi.a1 <- which(lv.combi.a1 == current.dose.1.combi.1)
                    if(lvid.combi.a1==nlv.combi.a1){
                        max.next.combi.a1 <- lv.combi.a1[lvid.combi.a1]
                    }else{
                        max.next.combi.a1 <- lv.combi.a1[lvid.combi.a1+1]
                    }
                }

                if(esc.constrain.combi.a2){
                    lvid.combi.a2 <- which(lv.combi.a2 == current.dose.2.combi.1)
                    if(lvid.combi.a2==nlv.combi.a2){
                        max.next.combi.a2 <- lv.combi.a2[lvid.combi.a2]
                    }else{
                        max.next.combi.a2 <- lv.combi.a2[lvid.combi.a2+1]
                    }
                }

                if(n.pat[3]>=max.n.combi.a & mtd.enforce.combi.a==TRUE){
                    combi.is.MTD.1 <- TRUE
                }

            }
    #---------------------------------------------------------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------------------------------------------------------

            }else if((cohort.queue[cohort.counter]==6) &
                     (!(combi.is.MTD.2[1]) &
                      n.pat[6] < max.n.combi.b &
                      !combi.stopped.2 &
                      active.combi.b)){
                #--------------------------------------------------------------
                #case: combi 2
                #--------------------------------------------------------------
                #sample the number of patients and DLTs
                curr.data.entry <- curr.data.entry+1
                curr.dose.name <- paste( current.dose.1.combi.2, current.dose.2.combi.2, sep="+")
                curr.tox.combi.2 <- tox.combi.b[curr.dose.name]
                if(length(cohort.size.combi.b) == 1){
                    curr.n.pat <- cohort.size.combi.b
                }else{
                    curr.n.pat <- sample(cohort.size.combi.b, size=1, prob=cohort.prob.combi.b)
                }
                n.pat[6] <- n.pat[6] + curr.n.pat
                curr.dlt <- rbinom(1, curr.n.pat,  prob = curr.tox.combi.2)
                n.dlt[6] <- n.dlt[6] + curr.dlt

                #save relevant data
                patients.at.dose.combi.2[curr.dose.name] <- curr.n.pat + patients.at.dose.combi.2[curr.dose.name]

                if(curr.tox.combi.2 < dosing.intervals[1]){
                    #the cohort was treated with underdose
                    trials.n.under[6] <- as.numeric(trials.n.under[6]) + curr.n.pat
                    trials.dlt.under[6] <- as.numeric(trials.dlt.under[6]) + curr.dlt
                } else if (dosing.intervals[1] <= curr.tox.combi.2 & curr.tox.combi.2 < dosing.intervals[2]){
                    #the cohort was treated with dose in target toxicicty interval
                    trials.n.target[6] <- as.numeric(trials.n.target[6]) + curr.n.pat
                    trials.dlt.target[6] <- as.numeric(trials.dlt.target[6]) + curr.dlt
                } else {
                    # the cohort was given an overdose
                    trials.n.over[6] <- as.numeric(trials.n.over[6]) + curr.n.pat
                    trials.dlt.over[6] <- as.numeric(trials.dlt.over[6]) + curr.dlt
                }

                #----------------------------------------------------------------------------
                # prepare data and execute the main computation step
                #----------------------------------------------------------------------------
                #Prepare current data in vectors
                prep.data.doses.1[curr.data.entry] <- current.dose.1.combi.2
                prep.data.doses.2[curr.data.entry] <- current.dose.2.combi.2
                prep.data.DLT[curr.data.entry] <- curr.dlt
                prep.data.n.patients[curr.data.entry] <- curr.n.pat
                prep.data.study[curr.data.entry] <- data.n.study.combi.2

                #check for back-fill cohort and simulate it if allowed
                if(backfill.combi.b){
                    #if there is a previous dose level
                    #and if it is at least equal to the back-fill start dose
                    if(prev.dose.combi.b1 > 0 & prev.dose.combi.b2 > 0 &
                       prev.dose.combi.b1 + 2*.Machine$double.eps>=backfill.start.combi.b1 &
                       prev.dose.combi.b2 + 2*.Machine$double.eps>=backfill.start.combi.b2){
                        #if it is actually contained in doses
                        idx.prev.dose.combi.b <- which(doses.combi.b[1, ]==prev.dose.combi.b1 &
                                                           doses.combi.b[2, ]==prev.dose.combi.b2)
                        if(!length(idx.prev.dose.combi.b)==0){
                            #if there was an escalation to reach the current dose,
                            #and if the previous dose was not backfilled.
                            if(prev.dose.esc.combi.b & prev.dose.ewoc.combi.b&
                               !backfilled.combi.b[idx.prev.dose.combi.b]){
                                #generate backfill cohort for dose
                                prev.dose.name.combi.b <- paste0(prev.dose.combi.b1, "+",
                                                                 prev.dose.combi.b2)
                                prev.tox.combi.b <- tox.combi.b[prev.dose.name.combi.b]
                                if(length(backfill.size.combi.b)==1){
                                    bf.pat.combi.b <- backfill.size.combi.b
                                }else{
                                    bf.pat.combi.b <- sample(backfill.size.combi.b,
                                                             size=1, prob=backfill.prob.combi.b)
                                }
                                backfilled.combi.b[idx.prev.dose.combi.b] <- TRUE
                                #when there is at least 1 back-fill patient
                                if(bf.pat.combi.b>0){
                                    #generate and save cohort data
                                    n.pat[6] <- n.pat[6] + bf.pat.combi.b
                                    bf.dlt.combi.b <- rbinom(1, bf.pat.combi.b,
                                                             prob = prev.tox.combi.b)
                                    n.dlt[6] <- n.dlt[6] + bf.dlt.combi.b

                                    #add to data frame
                                    curr.data.entry <- curr.data.entry + 1
                                    prep.data.doses.1[curr.data.entry] <- prev.dose.combi.b1
                                    prep.data.doses.2[curr.data.entry] <- prev.dose.combi.b2
                                    prep.data.DLT[curr.data.entry] <- bf.dlt.combi.b
                                    prep.data.n.patients[curr.data.entry] <- bf.pat.combi.b
                                    prep.data.study[curr.data.entry] <- data.n.study.combi.2

                                    #save data
                                    patients.at.dose.combi.2[prev.dose.name.combi.b] <- bf.pat.combi.b +
                                        patients.at.dose.combi.2[prev.dose.name.combi.b]
                                    if(prev.tox.combi.b < dosing.intervals[1]){  #the cohort was treated with underdose
                                        trials.n.under[6] <- as.numeric(trials.n.under[6]) + bf.pat.combi.b
                                        trials.dlt.under[6] <- as.numeric(trials.dlt.under[6]) + bf.dlt.combi.b
                                    } else if (dosing.intervals[1] <= prev.tox.combi.b &
                                               prev.tox.combi.b < dosing.intervals[2]){
                                        #the cohort was treated with dose in target toxicicty interval
                                        trials.n.target[6] <- as.numeric(trials.n.target[6]) + bf.pat.combi.b
                                        trials.dlt.target[6] <- as.numeric(trials.dlt.target[6]) + bf.dlt.combi.b
                                    } else { # the cohort was given an overdose
                                        trials.n.over[6] <- as.numeric(trials.n.over[6]) + bf.pat.combi.b
                                        trials.dlt.over[6] <- as.numeric(trials.dlt.over[6]) + bf.dlt.combi.b
                                    }

                                }
                            }
                        }
                    }
                }

                #Call the main function
                current.results <- main_jointBLRM(dose1 = prep.data.doses.1[1:curr.data.entry],
                                                  dose2 = prep.data.doses.2[1:curr.data.entry],
                                                  n.pat = prep.data.n.patients[1:curr.data.entry],
                                                  n.dlt = prep.data.DLT[1:curr.data.entry],
                                                  study = prep.data.study[1:curr.data.entry],
                                                  prior.tau=prior.tau,
                                                  prior.mu = prior.mu,
                                                  saturating = saturating,
                                                  file.name = file.name,
                                                  working.path = working.path,


                                                  seed = mcmc_seeds[cohort.counter],

                                                  dose.ref1 = dose.ref1,
                                                  dose.ref2 = dose.ref2,
                                                  study.interest = data.n.study.combi.2,
                                                  type.interest = "combi",
                                                  loss.weights = loss.weights,
                                                  dynamic.weights = dynamic.weights,
                                                  esc.rule = esc.rule,
                                                  esc.comp.max = esc.comp.max,
                                                  esc.step1 = esc.step.combi.b1,
                                                  esc.step2 = esc.step.combi.b2,
                                                  max.next1 = max.next.combi.b1,
                                                  max.next2 = max.next.combi.b2,
                                                  dose1.interest = doses.combi.b[1, ],
                                                  dose2.interest = doses.combi.b[2, ],

                                                  curr.dose1=current.dose.1.combi.2,
                                                  curr.dose2=current.dose.2.combi.2,

                                                  check.prev.dose = backfill.combi.b,
                                                  dosing.intervals = dosing.intervals,
                                                  ewoc = ewoc,
                                                  chains = chains,
                                                  iter=iter,
                                                  refresh=refresh,
                                                  warmup=warmup,
                                                  adapt_delta=adapt_delta,
                                                  max_treedepth=max_treedepth
                )

                #obtain recommended next dose
                pot.nd <- current.results$next.d
                next.combi.2.dose.1 <- pot.nd[1]
                next.combi.2.dose.2 <- pot.nd[2]

                #check whether the trial was stopped
                if(all(pot.nd==c(0,0))){
                    combi.stopped.2 <- TRUE
                }else{
                        #Check whether this dose is the MTD
                        target.prob <- current.results$curr.ptar
                        curr.n.treated <- patients.at.dose.combi.2[curr.dose.name]

                        if(decision.combi.b$RULE == 1){

                            combi.is.MTD.2 <- (next.combi.2.dose.1==current.dose.1.combi.2) &
                                (next.combi.2.dose.2==current.dose.2.combi.2) &
                                (n.dlt[6]  >= decision.combi.b$MIN.DLT) &
                                (curr.n.treated >= decision.combi.b$PAT.AT.MTD) &
                                ((n.pat[6] >= decision.combi.b$MIN.PAT) | (target.prob >= decision.combi.b$TARGET.PROB ))


                        } else{

                            combi.is.MTD.2 <- (next.combi.2.dose.1==current.dose.1.combi.2) &
                                (next.combi.2.dose.2==current.dose.2.combi.2) &
                                (n.dlt[6]  >= decision.combi.b$MIN.DLT) &
                                (n.pat[6] >= decision.combi.b$MIN.PAT) &
                                ((curr.n.treated>= decision.combi.b$PAT.AT.MTD) | (target.prob >= decision.combi.b$TARGET.PROB ))

                        }

                        #update data that controls back-fill cohorts
                        if(backfill.combi.b){
                            prev.dose.combi.b1 <- current.dose.1.combi.2
                            prev.dose.combi.b2 <- current.dose.2.combi.2
                            prev.dose.esc.combi.b <- ((current.dose.1.combi.2 < next.combi.2.dose.1 &
                                                           current.dose.2.combi.2 <= next.combi.2.dose.2 + 2*.Machine$double.eps)|
                                                          (current.dose.2.combi.2 < next.combi.2.dose.2 &
                                                               current.dose.1.combi.2 <= next.combi.2.dose.1 + 2*.Machine$double.eps))
                            prev.dose.ewoc.combi.b <- current.results$bf.check
                        }

                        #continue with the new doses
                        current.dose.1.combi.2 <- next.combi.2.dose.1
                        current.dose.2.combi.2 <- next.combi.2.dose.2

                        #reset constrain levels if constrain active
                        if(esc.constrain.combi.b1){
                            lvid.combi.b1 <- which(lv.combi.b1 == current.dose.1.combi.2)
                            if(lvid.combi.b1==nlv.combi.b1){
                                max.next.combi.b1 <- lv.combi.b1[lvid.combi.b1]
                            }else{
                                max.next.combi.b1 <- lv.combi.b1[lvid.combi.b1+1]
                            }
                        }

                        if(esc.constrain.combi.b2){
                            lvid.combi.b2 <- which(lv.combi.b2 == current.dose.2.combi.2)
                            if(lvid.combi.b2==nlv.combi.b2){
                                max.next.combi.b2 <- lv.combi.b2[lvid.combi.b2]
                            }else{
                                max.next.combi.b2 <- lv.combi.b2[lvid.combi.b2+1]
                            }
                        }


                        if(n.pat[6]>=max.n.combi.b & mtd.enforce.combi.b==TRUE){
                            combi.is.MTD.2 <- TRUE
                        }
            }

#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------

        } else if ((cohort.queue[cohort.counter]==1) &
                   (!(mono.1.is.MTD[1]) &
                    n.pat[1] < max.n.mono1.a &
                    !mono.1.stopped & active.mono1.a)) {
            #-----------------------------------------------------------------------------------------------------------
            #CASE: MONO COMPONENT 1
            #-----------------------------------------------------------------------------------------------------------

            #similar to the processing of combination cohorts
            curr.data.entry <- curr.data.entry+1
            curr.dose.name <- toString(current.dose.mono.1)
            curr.tox.mono.1 <- tox.mono1.a[curr.dose.name]

            #sample the number of patients and DLTs
            if(length(cohort.size.mono1.a) == 1){
                curr.n.pat <- cohort.size.mono1.a
            }else{
                curr.n.pat <- sample(cohort.size.mono1.a, size=1, prob=cohort.prob.mono1.a)
            }
            n.pat[1] <- n.pat[1] + curr.n.pat
            curr.dlt <- rbinom(1, curr.n.pat,  prob = curr.tox.mono.1)
            n.dlt[1] <- n.dlt[1] + curr.dlt

            #save data
            patients.at.dose.mono.1[curr.dose.name ] <- curr.n.pat + patients.at.dose.mono.1[curr.dose.name ]
            if(curr.tox.mono.1 < dosing.intervals[1]){  #the cohort was treated with underdose
                trials.n.under[1] <- as.numeric(trials.n.under[1]) + curr.n.pat
                trials.dlt.under[1] <- as.numeric(trials.dlt.under[1]) + curr.dlt
            } else if (dosing.intervals[1] <= curr.tox.mono.1 & curr.tox.mono.1 < dosing.intervals[2]){
                #the cohort was treated with dose in target toxicicty interval
                trials.n.target[1] <- as.numeric(trials.n.target[1]) + curr.n.pat
                trials.dlt.target[1] <- as.numeric(trials.dlt.target[1]) + curr.dlt
            } else { # the cohort was given an overdose
                trials.n.over[1] <- as.numeric(trials.n.over[1]) + curr.n.pat
                trials.dlt.over[1] <- as.numeric(trials.dlt.over[1]) + curr.dlt
            }

            #----------------------------------------------------------------------------
            # prepare data and execute the main computation step
            #----------------------------------------------------------------------------
            #Prepare current data
            prep.data.doses.1[curr.data.entry] <- current.dose.mono.1
            prep.data.doses.2[curr.data.entry] <- 0
            prep.data.DLT[curr.data.entry] <- curr.dlt
            prep.data.n.patients[curr.data.entry] <- curr.n.pat
            prep.data.study[curr.data.entry] <- data.n.study.mono.1

            #check for back-fill cohort and simulate it if allowed
            if(backfill.mono1.a){
                #if there is a previous dose level
                #and if it is at least equal to the back-fill start dose
                if(prev.dose.mono1.a > 0 &
                   prev.dose.mono1.a + 2*.Machine$double.eps>=backfill.start.mono1.a){
                    #if it is actually contained in doses
                    idx.prev.dose.mono1.a <- which(doses.mono1.a==prev.dose.mono1.a)
                    if(!length(idx.prev.dose.mono1.a)==0){
                        #if there was an escalation to reach the current dose,
                        #and if the previous dose was not backfilled.
                        if(prev.dose.esc.mono1.a & prev.dose.ewoc.mono1.a&
                            !backfilled.mono1.a[idx.prev.dose.mono1.a]){
                            #generate backfill cohort for dose
                            prev.tox.mono1.a <- tox.mono1.a[toString(prev.dose.mono1.a)]
                            if(length(backfill.size.mono1.a)==1){
                                bf.pat.mono1.a <- backfill.size.mono1.a
                            }else{
                                bf.pat.mono1.a <- sample(backfill.size.mono1.a,
                                                         size=1, prob=backfill.prob.mono1.a)
                            }
                            backfilled.mono1.a[idx.prev.dose.mono1.a] <- TRUE
                            #when there is at least 1 back-fill patient
                            if(bf.pat.mono1.a>0){
                                #generate and save cohort data
                                n.pat[1] <- n.pat[1] + bf.pat.mono1.a
                                bf.dlt.mono1.a <- rbinom(1, bf.pat.mono1.a,
                                                    prob = prev.tox.mono1.a)
                                n.dlt[1] <- n.dlt[1] + bf.dlt.mono1.a

                                #add to data frame
                                curr.data.entry <- curr.data.entry + 1
                                prep.data.doses.1[curr.data.entry] <- prev.dose.mono1.a
                                prep.data.doses.2[curr.data.entry] <- 0
                                prep.data.DLT[curr.data.entry] <- bf.dlt.mono1.a
                                prep.data.n.patients[curr.data.entry] <- bf.pat.mono1.a
                                prep.data.study[curr.data.entry] <- data.n.study.mono.1

                                #save data
                                patients.at.dose.mono.1[toString(prev.dose.mono1.a) ] <- bf.pat.mono1.a +
                                    patients.at.dose.mono.1[toString(prev.dose.mono1.a)]
                                if(prev.tox.mono1.a < dosing.intervals[1]){  #the cohort was treated with underdose
                                    trials.n.under[1] <- as.numeric(trials.n.under[1]) + bf.pat.mono1.a
                                    trials.dlt.under[1] <- as.numeric(trials.dlt.under[1]) + bf.dlt.mono1.a
                                } else if (dosing.intervals[1] <= prev.tox.mono1.a &
                                           prev.tox.mono1.a < dosing.intervals[2]){
                                    #the cohort was treated with dose in target toxicicty interval
                                    trials.n.target[1] <- as.numeric(trials.n.target[1]) + bf.pat.mono1.a
                                    trials.dlt.target[1] <- as.numeric(trials.dlt.target[1]) + bf.dlt.mono1.a
                                } else { # the cohort was given an overdose
                                    trials.n.over[1] <- as.numeric(trials.n.over[1]) + bf.pat.mono1.a
                                    trials.dlt.over[1] <- as.numeric(trials.dlt.over[1]) + bf.dlt.mono1.a
                                }

                            }
                        }
                    }
                }
            }

            #Call the main function
            current.results <- main_jointBLRM(dose1 = prep.data.doses.1[1:curr.data.entry],
                                              dose2 = prep.data.doses.2[1:curr.data.entry],
                                              n.pat = prep.data.n.patients[1:curr.data.entry],
                                              n.dlt = prep.data.DLT[1:curr.data.entry],
                                              study = prep.data.study[1:curr.data.entry],
                                              prior.tau = prior.tau,
                                              prior.mu = prior.mu,
                                              saturating = saturating,
                                              file.name = file.name,
                                              working.path = working.path,


                                              seed = mcmc_seeds[cohort.counter],

                                              dose.ref1 = dose.ref1,
                                              dose.ref2 = dose.ref2,
                                              study.interest = data.n.study.mono.1, ##
                                              type.interest = "mono1",              ##
                                              loss.weights = loss.weights,
                                              dynamic.weights = dynamic.weights,
                                              esc.rule = esc.rule,
                                              esc.comp.max = esc.comp.max,
                                              esc.step1 = esc.step.mono1.a,        ##
                                              esc.step2 = esc.step.mono1.a,        ##

                                              max.next1 = max.next.mono1.a,
                                              max.next2 = NULL,
                                              dose1.interest = doses.mono1.a,      ##
                                              dose2.interest = rep(0, length(doses.mono1.a)),  ##

                                              curr.dose1=current.dose.mono.1,          ##
                                              curr.dose2=0,                       ##

                                              check.prev.dose = backfill.mono1.a,
                                              dosing.intervals = dosing.intervals,
                                              ewoc = ewoc,
                                              chains = chains,
                                              iter=iter,
                                              refresh=refresh,
                                              warmup=warmup,
                                              adapt_delta=adapt_delta,
                                              max_treedepth=max_treedepth)

            #obtain recommended next dose
            pot.nd <- current.results$next.d
            next.mono.dose.1 <- pot.nd[1]

            #check whether the trial was stopped
            if(all(pot.nd==c(0,0))){
                mono.1.stopped <- TRUE
            }else{

                #decide if dose is MTD
                target.prob <- current.results$curr.ptar
                curr.n.treated <- patients.at.dose.mono.1[curr.dose.name]
                if(decision.mono1.a$RULE == 1){

                    mono.1.is.MTD <- (next.mono.dose.1==current.dose.mono.1) &
                        (n.dlt[1]  >= decision.mono1.a$MIN.DLT) &
                        (curr.n.treated >= decision.mono1.a$PAT.AT.MTD) &
                        ((n.pat[1] >= decision.mono1.a$MIN.PAT) | (target.prob >= decision.mono1.a$TARGET.PROB ))


                } else{

                    mono.1.is.MTD <- (next.mono.dose.1==current.dose.mono.1) &
                        (n.dlt[1]  >= decision.mono1.a$MIN.DLT) &
                        (n.pat[1] >= decision.mono1.a$MIN.PAT) &
                        ((curr.n.treated>= decision.mono1.a$PAT.AT.MTD) | (target.prob >= decision.mono1.a$TARGET.PROB ))

                }

                #update data that controls back-fill cohorts
                if(backfill.mono1.a){
                    prev.dose.mono1.a <- current.dose.mono.1
                    prev.dose.esc.mono1.a <- (current.dose.mono.1 < next.mono.dose.1)
                    prev.dose.ewoc.mono1.a <- current.results$bf.check
                }

                #continue with the new doses
                current.dose.mono.1 <- next.mono.dose.1

                if(esc.constrain.mono1.a){
                    lvid.mono1.a <- which(lv.mono1.a == current.dose.mono.1)
                    if(lvid.mono1.a==nlv.mono1.a){
                        max.next.mono1.a <- lv.mono1.a[lvid.mono1.a]
                    }else{
                        max.next.mono1.a <- lv.mono1.a[lvid.mono1.a+1]
                    }
                }


                if(n.pat[1]>=max.n.mono1.a & mtd.enforce.mono1.a==TRUE){
                    mono.1.is.MTD <- TRUE
                }

            }

#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------

        } else if ((cohort.queue[cohort.counter]==4) &
                   (!(mono.4.is.MTD[1]) &
                    n.pat[4] < max.n.mono1.b &
                    !mono.4.stopped & active.mono1.b)) {
            #-----------------------------------------------------------------------------------------------------------
            #CASE: MONO COMPONENT 4
            #-----------------------------------------------------------------------------------------------------------
            #sample numbers of patients and DLTs
            curr.data.entry <- curr.data.entry+1
            curr.dose.name <- toString(current.dose.mono.4)
            curr.tox.mono.4 <- tox.mono1.b[curr.dose.name]
            if(length(cohort.size.mono1.b) == 1){
                curr.n.pat <- cohort.size.mono1.b
            }else{
                curr.n.pat <- sample(cohort.size.mono1.b, size=1, prob=cohort.prob.mono1.b)
            }
            n.pat[4] <- n.pat[4] + curr.n.pat
            curr.dlt <- rbinom(1, curr.n.pat,  prob = curr.tox.mono.4)
            n.dlt[4] <- n.dlt[4] + curr.dlt

            #save data
            patients.at.dose.mono.4[curr.dose.name] <- curr.n.pat + patients.at.dose.mono.4[curr.dose.name]
            if(curr.tox.mono.4 < dosing.intervals[1]){  #the cohort was treated with underdose
                trials.n.under[4] <- as.numeric(trials.n.under[4]) + curr.n.pat
                trials.dlt.under[4] <- as.numeric(trials.dlt.under[4]) + curr.dlt
            } else if (dosing.intervals[1] <= curr.tox.mono.4 & curr.tox.mono.4 < dosing.intervals[2]){
                #the cohort was treated with dose in target toxicicty interval
                trials.n.target[4] <- as.numeric(trials.n.target[4]) + curr.n.pat
                trials.dlt.target[4] <- as.numeric(trials.dlt.target[4]) + curr.dlt
            } else { # the cohort was given an overdose
                trials.n.over[4] <- as.numeric(trials.n.over[4]) + curr.n.pat
                trials.dlt.over[4] <- as.numeric(trials.dlt.over[4]) + curr.dlt
            }

            #----------------------------------------------------------------------------
            # prepare data and execute the main computation step
            #----------------------------------------------------------------------------
            #Prepare current data
            prep.data.doses.1[curr.data.entry] <- current.dose.mono.4
            prep.data.doses.2[curr.data.entry] <- 0
            prep.data.DLT[curr.data.entry] <- curr.dlt
            prep.data.n.patients[curr.data.entry] <- curr.n.pat
            prep.data.study[curr.data.entry] <- data.n.study.mono.4


            #check for back-fill cohort and simulate it if allowed
            if(backfill.mono1.b){
                #if there is a previous dose level
                #and if it is at least equal to the back-fill start dose
                if(prev.dose.mono1.b > 0 &
                   prev.dose.mono1.b + 2*.Machine$double.eps>=backfill.start.mono1.b){
                    #if it is actually contained in doses
                    idx.prev.dose.mono1.b <- which(doses.mono1.b==prev.dose.mono1.b)
                    if(!length(idx.prev.dose.mono1.b)==0){
                        #if there was an escalation to reach the current dose,
                        #and if the previous dose was not backfilled.
                        if(prev.dose.esc.mono1.b & prev.dose.ewoc.mono1.b&
                           !backfilled.mono1.b[idx.prev.dose.mono1.b]){
                            #generate backfill cohort for dose
                            prev.tox.mono1.b <- tox.mono1.b[toString(prev.dose.mono1.b)]
                            if(length(backfill.size.mono1.b)==1){
                                bf.pat.mono1.b <- backfill.size.mono1.b
                            }else{
                                bf.pat.mono1.b <- sample(backfill.size.mono1.b,
                                                         size=1, prob=backfill.prob.mono1.b)
                            }
                            backfilled.mono1.b[idx.prev.dose.mono1.b] <- TRUE
                            #when there is at least 1 back-fill patient
                            if(bf.pat.mono1.b>0){
                                #generate and save cohort data
                                n.pat[4] <- n.pat[4] + bf.pat.mono1.b
                                bf.dlt.mono1.b <- rbinom(1, bf.pat.mono1.b,
                                                         prob = prev.tox.mono1.b)
                                n.dlt[4] <- n.dlt[4] + bf.dlt.mono1.b

                                #add to data frame
                                curr.data.entry <- curr.data.entry + 1
                                prep.data.doses.1[curr.data.entry] <- prev.dose.mono1.b
                                prep.data.doses.2[curr.data.entry] <- 0
                                prep.data.DLT[curr.data.entry] <- bf.dlt.mono1.b
                                prep.data.n.patients[curr.data.entry] <- bf.pat.mono1.b
                                prep.data.study[curr.data.entry] <- data.n.study.mono.4

                                #save data
                                patients.at.dose.mono.4[toString(prev.dose.mono1.b) ] <- bf.pat.mono1.b +
                                    patients.at.dose.mono.4[toString(prev.dose.mono1.b)]
                                if(prev.tox.mono1.b < dosing.intervals[1]){  #the cohort was treated with underdose
                                    trials.n.under[4] <- as.numeric(trials.n.under[4]) + bf.pat.mono1.b
                                    trials.dlt.under[4] <- as.numeric(trials.dlt.under[4]) + bf.dlt.mono1.b
                                } else if (dosing.intervals[1] <= prev.tox.mono1.b &
                                           prev.tox.mono1.b < dosing.intervals[2]){
                                    #the cohort was treated with dose in target toxicicty interval
                                    trials.n.target[4] <- as.numeric(trials.n.target[4]) + bf.pat.mono1.b
                                    trials.dlt.target[4] <- as.numeric(trials.dlt.target[4]) + bf.dlt.mono1.b
                                } else { # the cohort was given an overdose
                                    trials.n.over[4] <- as.numeric(trials.n.over[4]) + bf.pat.mono1.b
                                    trials.dlt.over[4] <- as.numeric(trials.dlt.over[4]) + bf.dlt.mono1.b
                                }

                            }
                        }
                    }
                }
            }

            #Call the main function
            current.results <- main_jointBLRM(dose1 = prep.data.doses.1[1:curr.data.entry],
                                              dose2 = prep.data.doses.2[1:curr.data.entry],
                                              n.pat = prep.data.n.patients[1:curr.data.entry],
                                              n.dlt = prep.data.DLT[1:curr.data.entry],
                                              study = prep.data.study[1:curr.data.entry],
                                              prior.tau=prior.tau,
                                              prior.mu = prior.mu,
                                              saturating = saturating,
                                              file.name = file.name,
                                              working.path = working.path,


                                              seed = mcmc_seeds[cohort.counter],

                                              dose.ref1 = dose.ref1,
                                              dose.ref2 = dose.ref2,
                                              study.interest = data.n.study.mono.4, ##
                                              type.interest = "mono1",              ##
                                              loss.weights = loss.weights,
                                              dynamic.weights = dynamic.weights,
                                              esc.rule = esc.rule,
                                              esc.comp.max = esc.comp.max,
                                              esc.step1 = esc.step.mono1.b,        ##
                                              esc.step2 = esc.step.mono1.b,        ##

                                              max.next1 = max.next.mono1.b,
                                              max.next2 = NULL,

                                              dose1.interest = doses.mono1.b,      ##
                                              dose2.interest = rep(0, length(doses.mono1.b)),  ##

                                              curr.dose1=current.dose.mono.4,          ##
                                              curr.dose2=0,                       ##

                                              check.prev.dose = backfill.mono1.b,
                                              dosing.intervals = dosing.intervals,
                                              ewoc = ewoc,
                                              chains = chains,
                                              iter=iter,
                                              refresh=refresh,
                                              warmup=warmup,
                                              adapt_delta=adapt_delta,
                                              max_treedepth=max_treedepth)

            #obtain recommended next dose
            pot.nd <- current.results$next.d
            next.mono.dose.4 <- pot.nd[1]

            #check whether the trial was stopped
            if(all(pot.nd==c(0,0))){
                mono.4.stopped <- TRUE
            }else{

                #check if dose is MTD
                target.prob <- current.results$curr.ptar
                curr.n.treated <- patients.at.dose.mono.4[curr.dose.name]

                if(decision.mono1.b$RULE == 1){

                    mono.4.is.MTD <- (next.mono.dose.4==current.dose.mono.4) &
                        (n.dlt[4]  >= decision.mono1.b$MIN.DLT) &
                        (curr.n.treated >= decision.mono1.b$PAT.AT.MTD) &
                        ((n.pat[4] >= decision.mono1.b$MIN.PAT) | (target.prob >= decision.mono1.b$TARGET.PROB ))


                } else {

                    mono.4.is.MTD <- (next.mono.dose.4==current.dose.mono.4) &
                        (n.dlt[4]  >= decision.mono1.b$MIN.DLT) &
                        (n.pat[4] >= decision.mono1.b$MIN.PAT) &
                        ((curr.n.treated>= decision.mono1.b$PAT.AT.MTD) | (target.prob >= decision.mono1.b$TARGET.PROB ))

                }

                #update data that controls back-fill cohorts
                if(backfill.mono1.b){
                    prev.dose.mono1.b <- current.dose.mono.4
                    prev.dose.esc.mono1.b <- (current.dose.mono.4 < next.mono.dose.4)
                    prev.dose.ewoc.mono1.b <- current.results$bf.check
                }

                #continue with the new dose
                current.dose.mono.4 <- next.mono.dose.4

                if(esc.constrain.mono1.b){
                    lvid.mono1.b <- which(lv.mono1.b == current.dose.mono.4)
                    if(lvid.mono1.b==nlv.mono1.b){
                        max.next.mono1.b <- lv.mono1.b[lvid.mono1.b]
                    }else{
                        max.next.mono1.b <- lv.mono1.b[lvid.mono1.b+1]
                    }
                }

                if(n.pat[4]>=max.n.mono1.b & mtd.enforce.mono1.b==TRUE){
                    mono.4.is.MTD <- TRUE
                }

            }

#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------

        } else if ((cohort.queue[cohort.counter]==2) &
                   (!(mono.2.is.MTD[1]) &
                    n.pat[2] < max.n.mono2.a &
                    !mono.2.stopped &
                    active.mono2.a)){
            #------------------------------------------------------------------------------------------------------
            #MONO COMPONENT 2
            #------------------------------------------------------------------------------------------------------
            #sample number of DLTs and patients
            curr.data.entry <- curr.data.entry+1
            curr.dose.name <-toString(current.dose.mono.2)
            curr.tox.mono.2 <- tox.mono2.a[curr.dose.name]
            if(length(cohort.size.mono2.a)<=1){
                curr.n.pat <- cohort.size.mono2.a
            }else{
                curr.n.pat <- sample(cohort.size.mono2.a, size=1, prob=cohort.prob.mono2.a)
            }

            n.pat[2] <- n.pat[2] + curr.n.pat
            curr.dlt <- rbinom(1, curr.n.pat,  prob = curr.tox.mono.2)
            n.dlt[2] <- n.dlt[2] + curr.dlt

            #save the data
            patients.at.dose.mono.2[curr.dose.name] <- curr.n.pat + patients.at.dose.mono.2[curr.dose.name]
            if(curr.tox.mono.2 < dosing.intervals[1]){  #the cohort was treated with underdose
                trials.n.under[2] <- as.numeric(trials.n.under[2]) + curr.n.pat
                trials.dlt.under[2] <- as.numeric(trials.dlt.under[2]) + curr.dlt
            } else if (dosing.intervals[1] <= curr.tox.mono.2 & curr.tox.mono.2 < dosing.intervals[2]){
                #the cohort was treated with dose in target toxicicty interval
                trials.n.target[2] <- as.numeric(trials.n.target[2]) + curr.n.pat
                trials.dlt.target[2] <- as.numeric(trials.dlt.target[2]) + curr.dlt
            } else { # the cohort was given an overdose
                trials.n.over[2] <- as.numeric(trials.n.over[2]) + curr.n.pat
                trials.dlt.over[2] <- as.numeric(trials.dlt.over[2]) + curr.dlt
            }

            #------------------------------------------------------------------------------
            # prepare data and execute the main computation step
            #------------------------------------------------------------------------------
            #Prepare current data
            prep.data.doses.1[curr.data.entry] <- 0
            prep.data.doses.2[curr.data.entry] <- current.dose.mono.2
            prep.data.DLT[curr.data.entry] <- curr.dlt
            prep.data.n.patients[curr.data.entry] <- curr.n.pat
            prep.data.study[curr.data.entry] <- data.n.study.mono.2

            #check for back-fill cohort and simulate it if allowed
            if(backfill.mono2.a){
                #if there is a previous dose level
                #and if it is at least equal to the back-fill start dose
                if(prev.dose.mono2.a > 0 &
                   prev.dose.mono2.a + 2*.Machine$double.eps>=backfill.start.mono2.a){
                    #if it is actually contained in doses
                    idx.prev.dose.mono2.a <- which(doses.mono2.a==prev.dose.mono2.a)
                    if(!length(idx.prev.dose.mono2.a)==0){
                        #if there was an escalation to reach the current dose,
                        #and if the previous dose was not backfilled.
                        if(prev.dose.esc.mono2.a & prev.dose.ewoc.mono2.a&
                           !backfilled.mono2.a[idx.prev.dose.mono2.a]){
                            #generate backfill cohort for dose
                            prev.tox.mono2.a <- tox.mono2.a[toString(prev.dose.mono2.a)]
                            if(length(backfill.size.mono2.a)==1){
                                bf.pat.mono2.a <- backfill.size.mono2.a
                            }else{
                                bf.pat.mono2.a <- sample(backfill.size.mono2.a,
                                                         size=1, prob=backfill.prob.mono2.a)
                            }
                            backfilled.mono2.a[idx.prev.dose.mono2.a] <- TRUE
                            #when there is at least 1 back-fill patient
                            if(bf.pat.mono2.a>0){
                                #generate and save cohort data
                                n.pat[2] <- n.pat[2] + bf.pat.mono2.a
                                bf.dlt.mono2.a <- rbinom(1, bf.pat.mono2.a,
                                                         prob = prev.tox.mono2.a)
                                n.dlt[2] <- n.dlt[2] + bf.dlt.mono2.a

                                #add to data frame
                                curr.data.entry <- curr.data.entry + 1
                                prep.data.doses.1[curr.data.entry] <- 0
                                prep.data.doses.2[curr.data.entry] <- prev.dose.mono2.a
                                prep.data.DLT[curr.data.entry] <- bf.dlt.mono2.a
                                prep.data.n.patients[curr.data.entry] <- bf.pat.mono2.a
                                prep.data.study[curr.data.entry] <- data.n.study.mono.2

                                #save data
                                patients.at.dose.mono.2[toString(prev.dose.mono2.a) ] <- bf.pat.mono2.a +
                                    patients.at.dose.mono.2[toString(prev.dose.mono2.a)]
                                if(prev.tox.mono2.a < dosing.intervals[1]){  #the cohort was treated with underdose
                                    trials.n.under[2] <- as.numeric(trials.n.under[2]) + bf.pat.mono2.a
                                    trials.dlt.under[2] <- as.numeric(trials.dlt.under[2]) + bf.dlt.mono2.a
                                } else if (dosing.intervals[1] <= prev.tox.mono2.a &
                                           prev.tox.mono2.a < dosing.intervals[2]){
                                    #the cohort was treated with dose in target toxicicty interval
                                    trials.n.target[2] <- as.numeric(trials.n.target[2]) + bf.pat.mono2.a
                                    trials.dlt.target[2] <- as.numeric(trials.dlt.target[2]) + bf.dlt.mono2.a
                                } else { # the cohort was given an overdose
                                    trials.n.over[2] <- as.numeric(trials.n.over[2]) + bf.pat.mono2.a
                                    trials.dlt.over[2] <- as.numeric(trials.dlt.over[2]) + bf.dlt.mono2.a
                                }

                            }
                        }
                    }
                }
            }

            #Call the main function
            current.results <- main_jointBLRM(dose1 = prep.data.doses.1[1:curr.data.entry],
                                              dose2 = prep.data.doses.2[1:curr.data.entry],
                                              n.pat = prep.data.n.patients[1:curr.data.entry],
                                              n.dlt = prep.data.DLT[1:curr.data.entry],
                                              study = prep.data.study[1:curr.data.entry],
                                              prior.tau=prior.tau,
                                              prior.mu = prior.mu,
                                              saturating = saturating,
                                              file.name = file.name,
                                              working.path = working.path,


                                              seed = mcmc_seeds[cohort.counter],

                                              dose.ref1 = dose.ref1,
                                              dose.ref2 = dose.ref2,
                                              study.interest = data.n.study.mono.2, ##
                                              type.interest = "mono2",              ##
                                              loss.weights = loss.weights,
                                              dynamic.weights = dynamic.weights,
                                              esc.rule = esc.rule,
                                              esc.comp.max = esc.comp.max,
                                              esc.step1 = esc.step.mono2.a,        ##
                                              esc.step2 = esc.step.mono2.a,        ##
                                              max.next1 = NULL,
                                              max.next2 = max.next.mono2.a,
                                              dose1.interest = rep(0, length(doses.mono2.a)),      ##
                                              dose2.interest = doses.mono2.a,  ##

                                              curr.dose1=0,          ##
                                              curr.dose2=current.dose.mono.2,                       ##

                                              check.prev.dose = backfill.mono2.a,
                                              dosing.intervals = dosing.intervals,
                                              ewoc = ewoc,
                                              chains = chains,
                                              iter=iter,
                                              refresh=refresh,
                                              warmup=warmup,
                                              adapt_delta=adapt_delta,
                                              max_treedepth=max_treedepth)

            #obtain recommended next dose
            pot.nd <- current.results$next.d
            next.mono.dose.2 <- pot.nd[2]

            #check whether the trial was stopped
            if(all(pot.nd==c(0,0))){
                mono.2.stopped <- TRUE
            }else{

                #check if MTD is reached
                target.prob <- current.results$curr.ptar
                curr.n.treated <- patients.at.dose.mono.2[curr.dose.name]
                if(decision.mono2.a$RULE == 1){

                    mono.2.is.MTD <- (next.mono.dose.2==current.dose.mono.2) &
                        (n.dlt[2]  >= decision.mono2.a$MIN.DLT) &
                        (curr.n.treated >= decision.mono2.a$PAT.AT.MTD) &
                        ((n.pat[2] >= decision.mono2.a$MIN.PAT) | (target.prob >= decision.mono2.a$TARGET.PROB ))


                } else {

                    mono.2.is.MTD <- (next.mono.dose.2==current.dose.mono.2) &
                        (n.dlt[2]  >= decision.mono2.a$MIN.DLT) &
                        (n.pat[2] >= decision.mono2.a$MIN.PAT) &
                        ((curr.n.treated>= decision.mono2.a$PAT.AT.MTD) | (target.prob >= decision.mono2.a$TARGET.PROB ))

                }
                #update data that controls back-fill cohorts
                if(backfill.mono2.a){
                    prev.dose.mono2.a <- current.dose.mono.2
                    prev.dose.esc.mono2.a <- (current.dose.mono.2 < next.mono.dose.2)
                    prev.dose.ewoc.mono2.a <- current.results$bf.check
                }
                #continue with the new doses
                current.dose.mono.2 <- next.mono.dose.2
                if(esc.constrain.mono2.a){
                    lvid.mono2.a <- which(lv.mono2.a == current.dose.mono.2)
                    if(lvid.mono2.a==nlv.mono2.a){
                        max.next.mono2.a <- lv.mono2.a[lvid.mono2.a]
                    }else{
                        max.next.mono2.a <- lv.mono2.a[lvid.mono2.a+1]
                    }
                }

                if(n.pat[2]>=max.n.mono2.a & mtd.enforce.mono2.a==TRUE){
                    mono.2.is.MTD <- TRUE
                }

            }
#---------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------

        } else if ((cohort.queue[cohort.counter]==5) &
                   (!(mono.5.is.MTD[1]) &
                    n.pat[5] < max.n.mono2.b &
                    !mono.5.stopped & active.mono2.b)){
            #------------------------------------------------------------------------------------------------------
            #MONO COMPONENT 5
            #------------------------------------------------------------------------------------------------------
            curr.data.entry <- curr.data.entry+1
            curr.dose.name <- toString(current.dose.mono.5)
            curr.tox.mono.5 <- tox.mono2.b[curr.dose.name]
            if(length(cohort.size.mono2.b)<=1){
                curr.n.pat <- cohort.size.mono2.b
            }else{
                curr.n.pat <- sample(cohort.size.mono2.b, size=1, prob=cohort.prob.mono2.b)
            }
            n.pat[5] <- n.pat[5] + curr.n.pat
            curr.dlt <- rbinom(1, curr.n.pat,  prob = curr.tox.mono.5)
            n.dlt[5] <- n.dlt[5] + curr.dlt

            #save data
            patients.at.dose.mono.5[curr.dose.name] <- curr.n.pat + patients.at.dose.mono.5[curr.dose.name]
            if(curr.tox.mono.5 < dosing.intervals[1]){  #the cohort was treated with underdose
                trials.n.under[5] <- as.numeric(trials.n.under[5]) + curr.n.pat
                trials.dlt.under[5] <- as.numeric(trials.dlt.under[5]) + curr.dlt
            } else if (dosing.intervals[1] <= curr.tox.mono.5 & curr.tox.mono.5 < dosing.intervals[2]){
                #the cohort was treated with dose in target toxicicty interval
                trials.n.target[5] <- as.numeric(trials.n.target[5]) + curr.n.pat
                trials.dlt.target[5] <- as.numeric(trials.dlt.target[5]) + curr.dlt
            } else { # the cohort was given an overdose
                trials.n.over[5] <- as.numeric(trials.n.over[5]) + curr.n.pat
                trials.dlt.over[5] <- as.numeric(trials.dlt.over[5]) + curr.dlt
            }

            #------------------------------------------------------------------------------
            # prepare data and execute the main computation step
            #------------------------------------------------------------------------------
            #Prepare current data
            prep.data.doses.1[curr.data.entry] <- 0
            prep.data.doses.2[curr.data.entry] <- current.dose.mono.5
            prep.data.DLT[curr.data.entry] <- curr.dlt
            prep.data.n.patients[curr.data.entry] <- curr.n.pat
            prep.data.study[curr.data.entry] <- data.n.study.mono.5

            #check for back-fill cohort and simulate it if allowed
            if(backfill.mono2.b){
                #if there is a previous dose level
                #and if it is at least equal to the back-fill start dose
                if(prev.dose.mono2.b > 0 &
                   prev.dose.mono2.b + 2*.Machine$double.eps>=backfill.start.mono2.b){
                    #if it is actually contained in doses
                    idx.prev.dose.mono2.b <- which(doses.mono2.b==prev.dose.mono2.b)
                    if(!length(idx.prev.dose.mono2.b)==0){
                        #if there was an escalation to reach the current dose,
                        #and if the previous dose was not backfilled.
                        if(prev.dose.esc.mono2.b & prev.dose.ewoc.mono2.b&
                           !backfilled.mono2.b[idx.prev.dose.mono2.b]){
                            #generate backfill cohort for dose
                            prev.tox.mono2.b <- tox.mono2.b[toString(prev.dose.mono2.b)]
                            if(length(backfill.size.mono2.b)==1){
                                bf.pat.mono2.b <- backfill.size.mono2.b
                            }else{
                                bf.pat.mono2.b <- sample(backfill.size.mono2.b,
                                                         size=1, prob=backfill.prob.mono2.b)
                            }
                            backfilled.mono2.b[idx.prev.dose.mono2.b] <- TRUE
                            #when there is at least 1 back-fill patient
                            if(bf.pat.mono2.b>0){
                                #generate and save cohort data
                                n.pat[5] <- n.pat[5] + bf.pat.mono2.b
                                bf.dlt.mono2.b <- rbinom(1, bf.pat.mono2.b,
                                                         prob = prev.tox.mono2.b)
                                n.dlt[5] <- n.dlt[5] + bf.dlt.mono2.b

                                #add to data frame
                                curr.data.entry <- curr.data.entry + 1
                                prep.data.doses.1[curr.data.entry] <- 0
                                prep.data.doses.2[curr.data.entry] <- prev.dose.mono2.b
                                prep.data.DLT[curr.data.entry] <- bf.dlt.mono2.b
                                prep.data.n.patients[curr.data.entry] <- bf.pat.mono2.b
                                prep.data.study[curr.data.entry] <- data.n.study.mono.5

                                #save data
                                patients.at.dose.mono.5[toString(prev.dose.mono2.b) ] <- bf.pat.mono2.b +
                                    patients.at.dose.mono.5[toString(prev.dose.mono2.b)]
                                if(prev.tox.mono2.b < dosing.intervals[1]){  #the cohort was treated with underdose
                                    trials.n.under[5] <- as.numeric(trials.n.under[5]) + bf.pat.mono2.b
                                    trials.dlt.under[5] <- as.numeric(trials.dlt.under[5]) + bf.dlt.mono2.b
                                } else if (dosing.intervals[1] <= prev.tox.mono2.b &
                                           prev.tox.mono2.b < dosing.intervals[2]){
                                    #the cohort was treated with dose in target toxicicty interval
                                    trials.n.target[5] <- as.numeric(trials.n.target[5]) + bf.pat.mono2.b
                                    trials.dlt.target[5] <- as.numeric(trials.dlt.target[5]) + bf.dlt.mono2.b
                                } else { # the cohort was given an overdose
                                    trials.n.over[5] <- as.numeric(trials.n.over[5]) + bf.pat.mono2.b
                                    trials.dlt.over[5] <- as.numeric(trials.dlt.over[5]) + bf.dlt.mono2.b
                                }

                            }
                        }
                    }
                }
            }

            #Call the main function
            current.results <- main_jointBLRM(dose1 = prep.data.doses.1[1:curr.data.entry],
                                              dose2 = prep.data.doses.2[1:curr.data.entry],
                                              n.pat = prep.data.n.patients[1:curr.data.entry],
                                              n.dlt = prep.data.DLT[1:curr.data.entry],
                                              study = prep.data.study[1:curr.data.entry],
                                              prior.tau=prior.tau,
                                              prior.mu = prior.mu,
                                              saturating = saturating,
                                              file.name = file.name,
                                              working.path = working.path,


                                              seed = mcmc_seeds[cohort.counter],

                                              dose.ref1 = dose.ref1,
                                              dose.ref2 = dose.ref2,
                                              study.interest = data.n.study.mono.5, ##
                                              type.interest = "mono2",              ##
                                              loss.weights = loss.weights,
                                              dynamic.weights = dynamic.weights,
                                              esc.rule = esc.rule,
                                              esc.comp.max = esc.comp.max,
                                              esc.step1 = esc.step.mono2.b,        ##
                                              esc.step2 = esc.step.mono2.b,        ##
                                              max.next1 = NULL,
                                              max.next2 = max.next.mono2.b,
                                              dose1.interest = rep(0, length(doses.mono2.b)),      ##
                                              dose2.interest = doses.mono2.b,  ##

                                              curr.dose1=0,          ##
                                              curr.dose2=current.dose.mono.5,                       ##

                                              check.prev.dose = backfill.mono2.b,
                                              dosing.intervals = dosing.intervals,
                                              ewoc = ewoc,
                                              chains = chains,
                                              iter=iter,
                                              refresh=refresh,
                                              warmup=warmup,
                                              adapt_delta=adapt_delta,
                                              max_treedepth=max_treedepth)

            #obtain recommended next dose
            pot.nd <- current.results$next.d
            next.mono.dose.5 <- pot.nd[2]

            #check whether the trial was stopped
            if(all(pot.nd==c(0,0))){
                mono.5.stopped <- TRUE
            }else{

                #check for MTD
                target.prob <- current.results$curr.ptar
                curr.n.treated <- patients.at.dose.mono.5[curr.dose.name]
                if(decision.mono2.b$RULE == 1){

                    mono.5.is.MTD <- (next.mono.dose.5==current.dose.mono.5) &
                        (n.dlt[5]  >= decision.mono2.b$MIN.DLT) &
                        (curr.n.treated >= decision.mono2.b$PAT.AT.MTD) &
                        ((n.pat[5] >= decision.mono2.b$MIN.PAT) | (target.prob >= decision.mono2.b$TARGET.PROB ))

                } else {

                    mono.5.is.MTD <- (next.mono.dose.5==current.dose.mono.5) &
                        (n.dlt[5]  >= decision.mono2.b$MIN.DLT) &
                        (n.pat[5] >= decision.mono2.b$MIN.PAT) &
                        ((curr.n.treated>= decision.mono2.b$PAT.AT.MTD) | (target.prob >= decision.mono2.b$TARGET.PROB ))

                }

                #update data that controls back-fill cohorts
                if(backfill.mono2.b){
                    prev.dose.mono2.b <- current.dose.mono.5
                    prev.dose.esc.mono2.b <- (current.dose.mono.5 < next.mono.dose.5)
                    prev.dose.ewoc.mono2.b <- current.results$bf.check
                }

                #continue with the new doses
                current.dose.mono.5 <- next.mono.dose.5
                if(esc.constrain.mono2.b){
                    lvid.mono2.b <- which(lv.mono2.b == current.dose.mono.5)
                    if(lvid.mono2.b==nlv.mono2.b){
                        max.next.mono2.b <- lv.mono2.b[lvid.mono2.b]
                    }else{
                        max.next.mono2.b <- lv.mono2.b[lvid.mono2.b+1]
                    }
                }


                if(n.pat[5]>=max.n.mono2.b & mtd.enforce.mono2.b==TRUE){
                    mono.5.is.MTD <- TRUE
                }

            }

        }
        #---------------------
        #END OF IF-Conditions
        #---------------------

    }
    #------------------------------------------------------------------------
    #END OF THE WHILE-LOOP
    #------------------------------------------------------------------------

    #save the simulation results and prepare output
    result_simulation[["n.pat"]] <- n.pat
    result_simulation[["n.dlt"]] <- n.dlt

    result_simulation[["trials.n.under"]] <- trials.n.under
    result_simulation[["trials.dlt.under"]] <- trials.dlt.under

    result_simulation[["trials.n.over"]] <- trials.n.over
    result_simulation[["trials.dlt.over"]] <- trials.dlt.over

    result_simulation[["trials.n.target"]] <- trials.n.target
    result_simulation[["trials.dlt.target"]] <- trials.dlt.target


    #for all possible trials: check whether the trial was activated, and,
    #if yes, save the corresponding variables and current doses
    if(active.mono1.a==TRUE){
        result_simulation[['mono.1.stopped']] <- mono.1.stopped
        result_simulation[['mono.1.is.MTD']] <- mono.1.is.MTD
        result_simulation[['current.dose.mono.1']] <- current.dose.mono.1
        result_simulation[['curr.tox.mono.1']] <- curr.tox.mono.1
        result_simulation[['pat.d.mono.1']] <- patients.at.dose.mono.1
    }
    if(active.mono2.a== TRUE){
        result_simulation[['mono.2.stopped']] <- mono.2.stopped
        result_simulation[['mono.2.is.MTD']] <- mono.2.is.MTD
        result_simulation[['current.dose.mono.2']] <- current.dose.mono.2
        result_simulation[['curr.tox.mono.2']] <- curr.tox.mono.2
        result_simulation[['pat.d.mono.2']] <- patients.at.dose.mono.2
    }
    if(active.combi.a == TRUE){
        result_simulation[['combi.stopped.1']] <- combi.stopped.1
        result_simulation[['combi.is.MTD.1']] <- combi.is.MTD.1
        result_simulation[['current.dose.1.combi.1']] <- current.dose.1.combi.1
        result_simulation[['current.dose.2.combi.1']] <- current.dose.2.combi.1
        result_simulation[['curr.tox.combi.1']] <- curr.tox.combi.1
        result_simulation[['pat.d.combi.1']] <- patients.at.dose.combi.1
    }
    if(active.mono1.b==TRUE){
        result_simulation[['mono.4.stopped']] <- mono.4.stopped
        result_simulation[['mono.4.is.MTD']] <- mono.4.is.MTD
        result_simulation[['current.dose.mono.4']] <- current.dose.mono.4
        result_simulation[['curr.tox.mono.4']] <- curr.tox.mono.4
        result_simulation[['pat.d.mono.4']] <- patients.at.dose.mono.4
    }
    if(active.mono2.b== TRUE){
        result_simulation[['mono.5.stopped']] <- mono.5.stopped
        result_simulation[['mono.5.is.MTD']] <- mono.5.is.MTD
        result_simulation[['current.dose.mono.5']] <- current.dose.mono.5
        result_simulation[['curr.tox.mono.5']] <- curr.tox.mono.5
        result_simulation[['pat.d.mono.5']] <- patients.at.dose.mono.5
    }
    if(active.combi.b == TRUE){
        result_simulation[['combi.stopped.2']] <- combi.stopped.2
        result_simulation[['combi.is.MTD.2']] <- combi.is.MTD.2
        result_simulation[['current.dose.1.combi.2']] <- current.dose.1.combi.2
        result_simulation[['current.dose.2.combi.2']] <- current.dose.2.combi.2
        result_simulation[['curr.tox.combi.2']] <- curr.tox.combi.2
        result_simulation[['pat.d.combi.2']] <- patients.at.dose.combi.2
    }

    #return the results.
    return(result_simulation)

}

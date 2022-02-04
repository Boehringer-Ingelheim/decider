#context("scenario_jointBLRM")


test_that("scenario_jointBLRM() input checks", {

  #doses of interest
  expect_error(
    suppressMessages(scenario_jointBLRM(
      data = NULL,
      historical.data = NULL,
      doses.of.interest = c(1, 2, 3),
      dose.ref1 = 1,
      dose.ref2 = 1,
      trials.of.interest = c(1, 2, 3),
      types.of.interest = NULL,
      esc.rule = c("ewoc", "loss", "dynamic.loss"),
      dosing.intervals = c(0.16, 0.33, 0.6),
      ewoc.threshold = 0.25,
      loss.weights = c(1, 0, 1, 2),
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
  )


  #data tests
  wrong_data <- list(
    123,
    c(1, 2, 3),
    "foo",
    list(c(3), c(1), c(3), c(0), c(1)),
    list(dose1 =c(3), c(1), c(3), c(0), c(1)),
    list(dose1 =c(3), dose2 = c(1), c(3), c(0), c(1)),
    list(dose1 =c(3), dose2 = c(1), n.dlt=c(0), c(2), c(1)),
    list(dose1 =c(3), dose2 = c(1), n.dlt=c(0), n.pat = c(2), c(1)),
    list(dose1 =c(-13), dose2 = c(1), n.dlt=c(0), n.pat = c(2), trial= c(1)),
    list(dose1 =c(3), dose2 = c(-1), n.dlt=c(0), n.pat = c(2), trial= c(1)),
    list(dose1 =c(3), dose2 = c(1), n.dlt=c(12), n.pat = c(1), trial= c(1)),
    list(dose1 =c(13), dose2 = c(1), n.dlt=c(1.3), n.pat = c(4), trial= c(1)),
    list(dose1 =c(13), dose2 = c(1), n.dlt=c(1), n.pat = c(4.5), trial= c(1)),
    list(dose1 =c(3), dose2 = c(1), n.dlt=c(1), n.pat = c(4, 5), trial= c(1)),
    matrix(data = c(1, 2, 3, 3), ncol=2)
  )

  for(i in 1:length(wrong_data)){
   expect_error(
    suppressMessages(scenario_jointBLRM(
      data = wrong_data[[i]],
      historical.data = NULL,
      doses.of.interest = rbind(c(1, 2, 3),
                                c(0, 0, 0)),
      dose.ref1 = 1,
      dose.ref2 = 1,
      trials.of.interest = c(1, 2, 3),
      types.of.interest = NULL,
      esc.rule = c("ewoc", "loss", "dynamic.loss"),
      dosing.intervals = c(0.16, 0.33, 0.6),
      ewoc.threshold = 0.25,
      loss.weights = c(1, 0, 1, 2),
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
    )
  }

  for(i in 1:length(wrong_data)){
    expect_error(
      suppressMessages(scenario_jointBLRM(
        data = NULL,
        historical.data = wrong_data[[i]],
        doses.of.interest = rbind(c(1, 2, 3),
                                  c(0, 0, 0)),
        dose.ref1 = 1,
        dose.ref2 = 1,
        trials.of.interest = c(1, 2, 3),
        types.of.interest = NULL,
        esc.rule = c("ewoc", "loss", "dynamic.loss"),
        dosing.intervals = c(0.16, 0.33, 0.6),
        ewoc.threshold = 0.25,
        loss.weights = c(1, 0, 1, 2),
        dynamic.weights = rbind(
          c(0.32, 0, 0.32, 0.36),
          c(0.29, 0, 0.31, 0.4),
          c(0.27, 0, 0.33, 0.4),
          c(0.2,  0, 0.3,  0.5)
        ),
        prior.mu = list(
          mu_a1 =  c(logit(0.33), 2),
          mu_b1 =  c(0,          1),
          mu_a2 =  c(logit(0.33), 2),
          mu_b2 =  c(0,          1),
          mu_eta = c(0,          1.121)
        ),
        prior.tau = list(
          tau_a1 =  c(log(0.25),  log(2) / 1.96),
          tau_b1 =  c(log(0.125), log(2) /
                        1.96),
          tau_a2 =  c(log(0.25),  log(2) /
                        1.96),
          tau_b2 =  c(log(0.125), log(2) /
                        1.96),
          tau_eta = c(log(0.125), log(2) /
                        1.96)
        ),
        saturating = FALSE,
        probs = c(0.025, 0.5, 0.975),
        iter = 2000,
        warmup = 1000,
        refresh = 0,
        #floor(iter/10),
        adapt_delta = 0.8,
        max_treedepth = 15,
        chains = 4,
        seed = 4832042L,

        path = NULL,
        file.name = NULL,
        plot.decisions = FALSE,
        plot.combi.heatmap = TRUE,
        plot.int.probs.loss = FALSE,
        plot.return = FALSE,
        plot.file.format = "pdf",
        plot.width = 12,
        plot.height = 12,
        plot.unit = "mm",
        output.scen.config = FALSE
      ))
    )
  }

  #doses of interest
  wrong_doses <- list(
    "wrong",
    list(c(1, 2, 3), c(1, 2, 3)),
    c(123, 213, 231),
    rbind(c(12,12), c(-1 , -1)),
    rbind(c(-12,-12), c(1 , 1))
  )
  for(i in 1:length(wrong_data)){
    expect_error(
      suppressMessages(scenario_jointBLRM(
        data = NULL,
        historical.data = NULL,
        doses.of.interest = wrong_doses[[i]],
        dose.ref1 = 1,
        dose.ref2 = 1,
        trials.of.interest = c(1, 2, 3),
        types.of.interest = NULL,
        esc.rule = c("ewoc", "loss", "dynamic.loss"),
        dosing.intervals = c(0.16, 0.33, 0.6),
        ewoc.threshold = 0.25,
        loss.weights = c(1, 0, 1, 2),
        dynamic.weights = rbind(
          c(0.32, 0, 0.32, 0.36),
          c(0.29, 0, 0.31, 0.4),
          c(0.27, 0, 0.33, 0.4),
          c(0.2,  0, 0.3,  0.5)
        ),
        prior.mu = list(
          mu_a1 =  c(logit(0.33), 2),
          mu_b1 =  c(0,          1),
          mu_a2 =  c(logit(0.33), 2),
          mu_b2 =  c(0,          1),
          mu_eta = c(0,          1.121)
        ),
        prior.tau = list(
          tau_a1 =  c(log(0.25),  log(2) / 1.96),
          tau_b1 =  c(log(0.125), log(2) /
                        1.96),
          tau_a2 =  c(log(0.25),  log(2) /
                        1.96),
          tau_b2 =  c(log(0.125), log(2) /
                        1.96),
          tau_eta = c(log(0.125), log(2) /
                        1.96)
        ),
        saturating = FALSE,
        probs = c(0.025, 0.5, 0.975),
        iter = 2000,
        warmup = 1000,
        refresh = 0,
        #floor(iter/10),
        adapt_delta = 0.8,
        max_treedepth = 15,
        chains = 4,
        seed = 4832042L,

        path = NULL,
        file.name = NULL,
        plot.decisions = FALSE,
        plot.combi.heatmap = TRUE,
        plot.int.probs.loss = FALSE,
        plot.return = FALSE,
        plot.file.format = "pdf",
        plot.width = 12,
        plot.height = 12,
        plot.unit = "mm",
        output.scen.config = FALSE
      ))
    )
  }

  #doses of interest
  wrong_doses <- list(
    "wrong",
    list(c(1, 2, 3), c(1, 2, 3)),
    c(123, 213, 231),
    rbind(c(12,12), c(-1 , -1)),
    rbind(c(-12,-12), c(1 , 1))
  )
  for(i in 1:length(wrong_data)){
    expect_error(
      suppressMessages(scenario_jointBLRM(
        data = NULL,
        historical.data = NULL,
        doses.of.interest = rbind(c(12,12), c(1 , 1)),
        dose.ref1 = wrong_doses[[i]],
        dose.ref2 = 1,
        trials.of.interest = c(1, 2, 3),
        types.of.interest = NULL,
        esc.rule = c("ewoc", "loss", "dynamic.loss"),
        dosing.intervals = c(0.16, 0.33, 0.6),
        ewoc.threshold = 0.25,
        loss.weights = c(1, 0, 1, 2),
        dynamic.weights = rbind(
          c(0.32, 0, 0.32, 0.36),
          c(0.29, 0, 0.31, 0.4),
          c(0.27, 0, 0.33, 0.4),
          c(0.2,  0, 0.3,  0.5)
        ),
        prior.mu = list(
          mu_a1 =  c(logit(0.33), 2),
          mu_b1 =  c(0,          1),
          mu_a2 =  c(logit(0.33), 2),
          mu_b2 =  c(0,          1),
          mu_eta = c(0,          1.121)
        ),
        prior.tau = list(
          tau_a1 =  c(log(0.25),  log(2) / 1.96),
          tau_b1 =  c(log(0.125), log(2) /
                        1.96),
          tau_a2 =  c(log(0.25),  log(2) /
                        1.96),
          tau_b2 =  c(log(0.125), log(2) /
                        1.96),
          tau_eta = c(log(0.125), log(2) /
                        1.96)
        ),
        saturating = FALSE,
        probs = c(0.025, 0.5, 0.975),
        iter = 2000,
        warmup = 1000,
        refresh = 0,
        #floor(iter/10),
        adapt_delta = 0.8,
        max_treedepth = 15,
        chains = 4,
        seed = 4832042L,

        path = NULL,
        file.name = NULL,
        plot.decisions = FALSE,
        plot.combi.heatmap = TRUE,
        plot.int.probs.loss = FALSE,
        plot.return = FALSE,
        plot.file.format = "pdf",
        plot.width = 12,
        plot.height = 12,
        plot.unit = "mm",
        output.scen.config = FALSE
      ))
    )

    expect_error(
      suppressMessages(scenario_jointBLRM(
        data = NULL,
        historical.data = NULL,
        doses.of.interest = rbind(c(12,12), c(1 , 1)),
        dose.ref2 = wrong_doses[[i]],
        dose.ref1 = 1,
        trials.of.interest = c(1, 2, 3),
        types.of.interest = NULL,
        esc.rule = c("ewoc", "loss", "dynamic.loss"),
        dosing.intervals = c(0.16, 0.33, 0.6),
        ewoc.threshold = 0.25,
        loss.weights = c(1, 0, 1, 2),
        dynamic.weights = rbind(
          c(0.32, 0, 0.32, 0.36),
          c(0.29, 0, 0.31, 0.4),
          c(0.27, 0, 0.33, 0.4),
          c(0.2,  0, 0.3,  0.5)
        ),
        prior.mu = list(
          mu_a1 =  c(logit(0.33), 2),
          mu_b1 =  c(0,          1),
          mu_a2 =  c(logit(0.33), 2),
          mu_b2 =  c(0,          1),
          mu_eta = c(0,          1.121)
        ),
        prior.tau = list(
          tau_a1 =  c(log(0.25),  log(2) / 1.96),
          tau_b1 =  c(log(0.125), log(2) /
                        1.96),
          tau_a2 =  c(log(0.25),  log(2) /
                        1.96),
          tau_b2 =  c(log(0.125), log(2) /
                        1.96),
          tau_eta = c(log(0.125), log(2) /
                        1.96)
        ),
        saturating = FALSE,
        probs = c(0.025, 0.5, 0.975),
        iter = 2000,
        warmup = 1000,
        refresh = 0,
        #floor(iter/10),
        adapt_delta = 0.8,
        max_treedepth = 15,
        chains = 4,
        seed = 4832042L,

        path = NULL,
        file.name = NULL,
        plot.decisions = FALSE,
        plot.combi.heatmap = TRUE,
        plot.int.probs.loss = FALSE,
        plot.return = FALSE,
        plot.file.format = "pdf",
        plot.width = 12,
        plot.height = 12,
        plot.unit = "mm",
        output.scen.config = FALSE
      ))
    )
  }

  #---------------------------
  #trials of interest
  #---------------------------
  expect_error(
    suppressMessages(scenario_jointBLRM(
      data = NULL,
      historical.data = NULL,
      doses.of.interest = rbind(c(12,12), c(1 , 1)),
      dose.ref1 = 2,
      dose.ref2 = 1,
      trials.of.interest = NULL,
      types.of.interest = NULL,
      esc.rule = c("ewoc", "loss", "dynamic.loss"),
      dosing.intervals = c(0.16, 0.33, 0.6),
      ewoc.threshold = 0.25,
      loss.weights = c(1, 0, 1, 2),
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
  )

  #---------------------------
  #trials of interest
  #---------------------------
  expect_error(
    suppressMessages(scenario_jointBLRM(
      data = NULL,
      historical.data = NULL,
      doses.of.interest = rbind(c(12,12), c(1 , 1)),
      dose.ref1 = 2,
      dose.ref2 = 1,
      trials.of.interest = c(1, 2, 3),
      types.of.interest = "type",
      esc.rule = c("ewoc", "loss", "dynamic.loss"),
      dosing.intervals = c(0.16, 0.33, 0.6),
      ewoc.threshold = 0.25,
      loss.weights = c(1, 0, 1, 2),
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
  )

  expect_error(
    suppressMessages(scenario_jointBLRM(
      data = NULL,
      historical.data = NULL,
      doses.of.interest = rbind(c(12,12), c(1 , 1)),
      dose.ref1 = 2,
      dose.ref2 = 1,
      trials.of.interest = c(1, 2, 3),
      types.of.interest = 123298,
      esc.rule = c("ewoc", "loss", "dynamic.loss"),
      dosing.intervals = c(0.16, 0.33, 0.6),
      ewoc.threshold = 0.25,
      loss.weights = c(1, 0, 1, 2),
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
  )


  #---------------------------
  #dosing intervals
  #---------------------------
  expect_error(
    suppressMessages(scenario_jointBLRM(
      data = NULL,
      historical.data = NULL,
      doses.of.interest = rbind(c(12,12), c(1 , 1)),
      dose.ref1 = 2,
      dose.ref2 = 1,
      trials.of.interest = c(1, 2, 3),
      types.of.interest = NULL,
      esc.rule = c("ewoc"),
      dosing.intervals = c(2, 0.16, 0.33, 0.6),
      ewoc.threshold = 0.25,
      loss.weights = c(1, 0, 1, 2),
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
  )


  #---------------------------
  #ewoc threshold
  #---------------------------
  expect_error(
    suppressMessages(scenario_jointBLRM(
      data = NULL,
      historical.data = NULL,
      doses.of.interest = rbind(c(12,12), c(1 , 1)),
      dose.ref1 = 2,
      dose.ref2 = 1,
      trials.of.interest = c(1, 2, 3),
      types.of.interest = NULL,
      esc.rule = c("ewoc"),
      dosing.intervals = c(0.16, 0.33, 0.6),
      ewoc.threshold = "none",
      loss.weights = c(1, 0, 1, 2),
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
  )


  #---------------------------
  #loss stuff
  #---------------------------
  for(i in 1:length(wrong_data)){
  expect_error(
    suppressMessages(scenario_jointBLRM(
      data = NULL,
      historical.data = NULL,
      doses.of.interest = rbind(c(12,12), c(1 , 1)),
      dose.ref1 = 2,
      dose.ref2 = 1,
      trials.of.interest = c(1, 2, 3),
      types.of.interest = NULL,
      esc.rule = c("ewoc"),
      dosing.intervals = c(0.16, 0.33, 0.6),
      ewoc.threshold = 0.25,
      loss.weights = wrong_data[[i]],
      dynamic.weights = rbind(
        c(0.32, 0, 0.32, 0.36),
        c(0.29, 0, 0.31, 0.4),
        c(0.27, 0, 0.33, 0.4),
        c(0.2,  0, 0.3,  0.5)
      ),
      prior.mu = list(
        mu_a1 =  c(logit(0.33), 2),
        mu_b1 =  c(0,          1),
        mu_a2 =  c(logit(0.33), 2),
        mu_b2 =  c(0,          1),
        mu_eta = c(0,          1.121)
      ),
      prior.tau = list(
        tau_a1 =  c(log(0.25),  log(2) / 1.96),
        tau_b1 =  c(log(0.125), log(2) /
                      1.96),
        tau_a2 =  c(log(0.25),  log(2) /
                      1.96),
        tau_b2 =  c(log(0.125), log(2) /
                      1.96),
        tau_eta = c(log(0.125), log(2) /
                      1.96)
      ),
      saturating = FALSE,
      probs = c(0.025, 0.5, 0.975),
      iter = 2000,
      warmup = 1000,
      refresh = 0,
      #floor(iter/10),
      adapt_delta = 0.8,
      max_treedepth = 15,
      chains = 4,
      seed = 4832042L,

      path = NULL,
      file.name = NULL,
      plot.decisions = FALSE,
      plot.combi.heatmap = TRUE,
      plot.int.probs.loss = FALSE,
      plot.return = FALSE,
      plot.file.format = "pdf",
      plot.width = 12,
      plot.height = 12,
      plot.unit = "mm",
      output.scen.config = FALSE
    ))
  )

    expect_error(
      suppressMessages(scenario_jointBLRM(
        data = NULL,
        historical.data = NULL,
        doses.of.interest = rbind(c(12,12), c(1 , 1)),
        dose.ref1 = 2,
        dose.ref2 = 1,
        trials.of.interest = c(1, 2, 3),
        types.of.interest = NULL,
        esc.rule = c("ewoc"),
        dosing.intervals = c(0.16, 0.33, 0.6),
        ewoc.threshold = 0.25,
        loss.weights = c(1, 0, 1, 2),
        dynamic.weights = wrong_data[[i]],
        prior.mu = list(
          mu_a1 =  c(logit(0.33), 2),
          mu_b1 =  c(0,          1),
          mu_a2 =  c(logit(0.33), 2),
          mu_b2 =  c(0,          1),
          mu_eta = c(0,          1.121)
        ),
        prior.tau = list(
          tau_a1 =  c(log(0.25),  log(2) / 1.96),
          tau_b1 =  c(log(0.125), log(2) /
                        1.96),
          tau_a2 =  c(log(0.25),  log(2) /
                        1.96),
          tau_b2 =  c(log(0.125), log(2) /
                        1.96),
          tau_eta = c(log(0.125), log(2) /
                        1.96)
        ),
        saturating = FALSE,
        probs = c(0.025, 0.5, 0.975),
        iter = 2000,
        warmup = 1000,
        refresh = 0,
        #floor(iter/10),
        adapt_delta = 0.8,
        max_treedepth = 15,
        chains = 4,
        seed = 4832042L,

        path = NULL,
        file.name = NULL,
        plot.decisions = FALSE,
        plot.combi.heatmap = TRUE,
        plot.int.probs.loss = FALSE,
        plot.return = FALSE,
        plot.file.format = "pdf",
        plot.width = 12,
        plot.height = 12,
        plot.unit = "mm",
        output.scen.config = FALSE
      ))
    )
  }


#end of test_that call
})

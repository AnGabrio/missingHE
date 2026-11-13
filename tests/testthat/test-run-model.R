test_that("check output from running model internal functions", {

  testthat::skip_on_cran()

  load(test_path("testdata/sm.RData"))

  n <- sum(sm$data_set$data_raw$n)
  n.trt <- length(names(sm$data_set$data_raw$n))
  names.trt <- names(sm$data_set$data_raw$n)
  data_model <- list("n" = n,
                     "eff" = sm$data_set$data_raw$e, "cost" = sm$data_set$data_raw$c,
                     "m_eff" = sm$data_set$data_raw$me, "m_cost" = sm$data_set$data_raw$mc,
                     "X_e_fixed" = sm$data_set$data_raw$x_e_fixed,
                     "X_c_fixed" = sm$data_set$data_raw$x_c_fixed,
                     "Z_e_fixed" = sm$data_set$data_raw$z_e_fixed,
                     "Z_c_fixed" = sm$data_set$data_raw$z_c_fixed,
                     "mean_cov_e_fixed" = apply(sm$data_set$data_raw$x_e_fixed, 2, mean, na.rm = TRUE),
                     "mean_cov_c_fixed" = apply(sm$data_set$data_raw$x_c_fixed, 2, mean, na.rm = TRUE),
                     "mean_z_e_fixed" = apply(sm$data_set$data_raw$z_e_fixed, 2, mean, na.rm = TRUE),
                     "mean_z_c_fixed" = apply(sm$data_set$data_raw$z_c_fixed, 2, mean, na.rm = TRUE),
                     "pe_fixed" = ncol(sm$data_set$data_raw$x_e_fixed),
                     "pc_fixed" = ncol(sm$data_set$data_raw$x_c_fixed),
                     "ze_fixed" = ncol(sm$data_set$data_raw$z_e_fixed),
                     "zc_fixed" = ncol(sm$data_set$data_raw$z_c_fixed),
                     "X_e_random" = NULL,
                     "X_c_random" = NULL,
                     "mean_cov_e_random" = NULL,
                     "mean_cov_c_random" = NULL,
                     "pe_random" = 0, "pc_random" = 0,
                     "clus_e" = NULL, "clus_c" = NULL,
                     "n_clus_e" = NULL, "n_clus_c" = NULL,
                     "Z_e_random" = NULL, "Z_c_random" = NULL,
                     "mean_z_e_random" = NULL,
                     "mean_z_c_random" = NULL,
                     "ze_random" = 0, "zc_random" = 0,
                     "clus_me" = NULL, "clus_mc" = NULL,
                     "n_clus_me" = NULL, "n_clus_mc" = NULL,
                     "trt_pos_e" = sm$data_set$data_raw$trt_pos_e,
                     "trt_pos_c" = sm$data_set$data_raw$trt_pos_c,
                     "trt_pos_me" = sm$data_set$data_raw$trt_pos_me,
                     "trt_pos_mc" = sm$data_set$data_raw$trt_pos_mc,
                     "n_trt" = sm$data_set$data_raw$n,
                     "trt_index" = sm$data_set$data_raw$trt_index,
                     "trt_lev" = names.trt,
                     "efft" = sm$data_set$data_raw$efft,
                     "costt" = sm$data_set$data_raw$costt,
                     "m_efft" = sm$data_set$data_raw$m_efft,
                     "m_costt" = sm$data_set$data_raw$m_costt,
                     "clus_e_lev" = sm$data_set$data_raw$clus_e_lev,
                     "clus_c_lev" = sm$data_set$data_raw$clus_c_lev,
                     "clus_me_lev" = sm$data_set$data_raw$clus_me_lev,
                     "clus_mc_lev" = sm$data_set$data_raw$clus_mc_lev)
  type <- sm$type
  dist_e <- sm$model_output$dist_e
  dist_c <- sm$model_output$dist_c
  n.iter <- 10
  model_info <- list("is_c_random_c" = FALSE, "is_int_c_random_c" = FALSE,
                     "is_me_random_e" = FALSE, "is_mc_random_c" = FALSE,
                     "is_int_me_random_e" = FALSE, "is_int_mc_random_c" = FALSE,
                     "ind_random" = sm$model_output$ind_random,
                     "ind_fixed" = sm$model_output$ind_fixed,
                     "model_e_fixed" = "trt",
                     "model_c_fixed" = "trt",
                     "model_e_random" = character(0),
                     "model_c_random" = character(0),
                     "model_me_random" = character(0),
                     "model_mc_random" = character(0),
                     "dic" = TRUE,
                     "pd" = FALSE,
                     "n.iter.pd" = 10,
                     "n.adapt" = 100,
                     "inits" = NULL,
                     "n.chains" = 2,
                     "n.iter" = n.iter,
                     "n.burnin" = floor(10 / 2),
                     "n.thin" = 1,
                     "prior" = list("default"),
                     "n.mci" = 10,
                     "prob" = c(0.025, 0.975),
                     "progress.bar" = "text",
                     "refresh" = 10000/50,
                     "RNGname" = "Wichmann-Hill",
                     "quiet" = FALSE,
                     "digits" = 3,
                     "save_model" = FALSE)
  run_sm <- run_selection(data_model = data_model, type = type ,
                            dist_e = dist_e, dist_c = dist_c,
                            model_info = model_info)

  expect_length(run_sm, 26)
  expect_type(run_sm, "list")
  expect_named(run_sm, c("summary", "model", "mean_effects", "mean_costs",
                         "sd_effects", "sd_costs", "effects_fixed", "costs_fixed",
                         "pmis_effects", "mis_effects_fixed", "pmis_costs", "mis_costs_fixed",
                         "effects_random", "costs_random", "mis_effects_random", "mis_costs_random",
                         "imputed", "imputed_trt", "loglik", "type", "method", "ind_fixed",
                         "ind_random", "dist_e", "dist_c", "filein"))

  expect_length(run_sm$filein, 1)
  expect_type(run_sm$filein, "character")
  expect_length(run_sm$dist_c, 1)
  expect_type(run_sm$dist_c, "character")
  expect_equal(run_sm$dist_c, dist_c)
  expect_length(run_sm$dist_e, 1)
  expect_type(run_sm$dist_e, "character")
  expect_equal(run_sm$dist_e, dist_e)
  expect_length(run_sm$ind_random, 1)
  expect_type(run_sm$ind_random, "logical")
  expect_equal(run_sm$ind_random, TRUE)
  expect_length(run_sm$ind_fixed, 1)
  expect_type(run_sm$ind_fixed, "logical")
  expect_equal(run_sm$ind_fixed, TRUE)
  expect_length(run_sm$method, 1)
  expect_type(run_sm$method, "character")
  expect_equal(run_sm$method, "SELECTION")
  expect_length(run_sm$type, 1)
  expect_type(run_sm$type, "character")
  expect_equal(run_sm$type, type)
  expect_length(run_sm$loglik, 4)
  expect_type(run_sm$loglik, "list")
  expect_named(run_sm$loglik, c("effects", "costs",
                                "missing indicators effects", "missing indicators costs" ))
  expect_shape(run_sm$loglik$effects, dim = c(n.iter, n))
  expect_shape(run_sm$loglik$costs, dim = c(n.iter, n))
  expect_shape(run_sm$loglik$`missing indicators effects`, dim = c(n.iter, n))
  expect_shape(run_sm$loglik$`missing indicators costs`, dim = c(n.iter, n))
  expect_length(run_sm$imputed_trt, 2)
  expect_type(run_sm$imputed_trt, "list")
  expect_named(run_sm$imputed_trt, c("effects", "costs"))
  expect_length(run_sm$imputed_trt$effects, 4)
  expect_type(run_sm$imputed_trt$effects, "list")
  expect_named(run_sm$imputed_trt$effects, c("avg", "ql", "qu", "index"))
  expect_length(run_sm$imputed_trt$costs, 4)
  expect_type(run_sm$imputed_trt$costs, "list")
  expect_named(run_sm$imputed_trt$costs, c("avg", "ql", "qu", "index"))
  expect_length(run_sm$imputed_trt$effects$avg, n.trt)
  expect_type(run_sm$imputed_trt$effects$avg, "list")
  expect_named(run_sm$imputed_trt$effects$avg, names.trt)
  expect_length(run_sm$imputed_trt$effects$ql, n.trt)
  expect_type(run_sm$imputed_trt$effects$ql, "list")
  expect_named(run_sm$imputed_trt$effects$ql, names.trt)
  expect_length(run_sm$imputed_trt$effects$qu, n.trt)
  expect_type(run_sm$imputed_trt$effects$qu, "list")
  expect_named(run_sm$imputed_trt$effects$qu, names.trt)
  expect_length(run_sm$imputed_trt$effects$index, n.trt)
  expect_type(run_sm$imputed_trt$effects$index, "list")
  expect_named(run_sm$imputed_trt$effects$index, names.trt)
  expect_vector(run_sm$imputed_trt$effects$avg$`1`)
  expect_vector(run_sm$imputed_trt$effects$avg$`2`)
  expect_vector(run_sm$imputed_trt$effects$ql$`1`)
  expect_vector(run_sm$imputed_trt$effects$ql$`2`)
  expect_vector(run_sm$imputed_trt$effects$qu$`1`)
  expect_vector(run_sm$imputed_trt$effects$qu$`2`)
  expect_vector(run_sm$imputed_trt$effects$index$`1`)
  expect_vector(run_sm$imputed_trt$effects$index$`2`)
  expect_length(run_sm$imputed_trt$costs$avg, n.trt)
  expect_type(run_sm$imputed_trt$costs$avg, "list")
  expect_named(run_sm$imputed_trt$costs$avg, names.trt)
  expect_length(run_sm$imputed_trt$costs$ql, n.trt)
  expect_type(run_sm$imputed_trt$costs$ql, "list")
  expect_named(run_sm$imputed_trt$costs$ql, names.trt)
  expect_length(run_sm$imputed_trt$costs$qu, n.trt)
  expect_type(run_sm$imputed_trt$costs$qu, "list")
  expect_named(run_sm$imputed_trt$costs$qu, names.trt)
  expect_length(run_sm$imputed_trt$costs$index, n.trt)
  expect_type(run_sm$imputed_trt$costs$index, "list")
  expect_named(run_sm$imputed_trt$costs$index, names.trt)
  expect_vector(run_sm$imputed_trt$costs$avg$`1`)
  expect_vector(run_sm$imputed_trt$costs$avg$`2`)
  expect_vector(run_sm$imputed_trt$costs$ql$`1`)
  expect_vector(run_sm$imputed_trt$costs$ql$`2`)
  expect_vector(run_sm$imputed_trt$costs$qu$`1`)
  expect_vector(run_sm$imputed_trt$costs$qu$`2`)
  expect_vector(run_sm$imputed_trt$costs$index$`1`)
  expect_vector(run_sm$imputed_trt$costs$index$`2`)
  expect_length(run_sm$imputed, 2)
  expect_type(run_sm$imputed, "list")
  expect_named(run_sm$imputed, c("effects", "costs"))
  expect_length(run_sm$imputed$effects, 4)
  expect_type(run_sm$imputed$effects, "list")
  expect_named(run_sm$imputed$effects, c("avg", "ql", "qu", "index"))
  expect_vector(run_sm$imputed$effects$avg)
  expect_vector(run_sm$imputed$effects$ql)
  expect_vector(run_sm$imputed$effects$qu)
  expect_vector(run_sm$imputed$effects$index)
  expect_length(run_sm$imputed$costs, 4)
  expect_type(run_sm$imputed$costs, "list")
  expect_named(run_sm$imputed$costs, c("avg", "ql", "qu", "index"))
  expect_vector(run_sm$imputed$costs$avg)
  expect_vector(run_sm$imputed$costs$ql)
  expect_vector(run_sm$imputed$costs$qu)
  expect_vector(run_sm$imputed$costs$index)
  expect_null(run_sm$mis_costs_random)
  expect_null(run_sm$mis_effects_random)
  expect_null(run_sm$costs_random)
  expect_null(run_sm$effects_random)
  expect_shape(run_sm$mis_costs_fixed, nrow = n.iter)
  expect_shape(run_sm$mis_effects_fixed, nrow = n.iter)
  expect_shape(run_sm$pmis_costs, dim = c(n.iter, 1))
  expect_shape(run_sm$pmis_effects, dim = c(n.iter, 1))
  expect_shape(run_sm$costs_fixed, nrow = n.iter)
  expect_shape(run_sm$effects_fixed, nrow = n.iter)
  expect_shape(run_sm$sd_costs, dim = c(n.iter, 1))
  expect_shape(run_sm$sd_effects, dim = c(n.iter, 1))
  expect_shape(run_sm$mean_costs, dim = c(n.iter, n.trt))
  expect_shape(run_sm$mean_effects, dim = c(n.iter, n.trt))
  expect_true(is.matrix(run_sm$summary))
  expect_contains(colnames(run_sm$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                              "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_sm$summary, ncol = 9)
  expect_length(run_sm$model, 6)
  expect_type(run_sm$model, "list")
  expect_named(run_sm$model, c("model", "BUGSoutput", "parameters.to.save",
                                        "model.file", "n.iter", "DIC"))
  expect_true(run_sm$model$DIC)
  expect_equal(run_sm$model$n.iter, n.iter)
  expect_type(run_sm$model$parameters.to.save, "character")
  expect_vector(run_sm$model$parameters.to.save)
  expect_type(run_sm$model$model, "list")
  expect_named(run_sm$model$model, c("ptr", "data", "model", "state", "nchain",
                                              "iter", "sync", "recompile"))
  expect_type(run_sm$model$BUGSoutput, "list")
  expect_named(run_sm$model$BUGSoutput, c("n.chains", "n.iter", "n.burnin", "n.thin", "n.keep",
                                                   "n.sims", "sims.array", "sims.list", "sims.matrix",
                                                   "summary", "mean", "sd", "median", "root.short",
                                                   "long.short", "dimension.short", "indexes.short",
                                                   "last.values", "program", "model.file", "isDIC",
                                                   "DICbyR", "pV", "DIC", "time2run" ))
  expect_gt(run_sm$model$BUGSoutput$time2run, 0)
  expect_length(run_sm$model$BUGSoutput$DIC, 1)
  expect_length(run_sm$model$BUGSoutput$pV, 1)
  expect_equal(run_sm$model$BUGSoutput$DICbyR, run_sm$model$DIC)
  expect_equal(run_sm$model$BUGSoutput$isDIC, run_sm$model$DIC)
  expect_equal(run_sm$model$BUGSoutput$program, "jags")
  expect_type(run_sm$model$BUGSoutput$last.values, "list")
  expect_type(run_sm$model$BUGSoutput$indexes.short, "list")
  expect_vector(run_sm$model$BUGSoutput$dimension.short)
  expect_type(run_sm$model$BUGSoutput$long.short, "list")
  expect_vector(run_sm$model$BUGSoutput$root.short)
  expect_type(run_sm$model$BUGSoutput$median, "list")
  expect_type(run_sm$model$BUGSoutput$sd, "list")
  expect_type(run_sm$model$BUGSoutput$mean, "list")
  expect_true(is.matrix(run_sm$model$BUGSoutput$summary))
  expect_contains(colnames(run_sm$model$BUGSoutput$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                              "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_sm$model$BUGSoutput$summary, ncol = 9)
  expect_shape(run_sm$model$BUGSoutput$sims.matrix, nrow = n.iter)
  expect_type(run_sm$model$BUGSoutput$sims.list, "list")
  expect_shape(run_sm$model$BUGSoutput$sims.array, dim = c(n.iter - run_sm$model$BUGSoutput$n.burnin,
                                                                    run_sm$model$BUGSoutput$n.chains,
                                                                    dim(run_sm$model$BUGSoutput$sims.matrix)[2]))
  expect_equal(run_sm$model$BUGSoutput$n.sims, n.iter)
  expect_equal(run_sm$model$BUGSoutput$n.keep, n.iter - run_sm$model$BUGSoutput$n.burnin)
  expect_equal(run_sm$model$BUGSoutput$n.iter, n.iter)



  load(test_path("testdata/pmm.RData"))

  n <- sum(pmm$data_set$data_raw$n)
  n.trt <- length(names(pmm$data_set$data_raw$n))
  names.trt <- names(pmm$data_set$data_raw$n)
  n.pattern <- max(pmm$data_set$data_raw$n_obs_pat)
  data_model <- list("n" = n,
                     "eff" = pmm$data_set$data_raw$e, "cost" = pmm$data_set$data_raw$c,
                     "m_eff" = pmm$data_set$data_raw$me, "m_cost" = pmm$data_set$data_raw$mc,
                     "d" = pmm$data_set$data_raw$d, "n_patterns" = n.pattern,
                     "X_e_fixed" = pmm$data_set$data_raw$x_e_fixed,
                     "X_c_fixed" = pmm$data_set$data_raw$x_c_fixed,
                     "mean_cov_e_fixed" = apply(pmm$data_set$data_raw$x_e_fixed, 2, mean, na.rm = TRUE),
                     "mean_cov_c_fixed" = apply(pmm$data_set$data_raw$x_c_fixed, 2, mean, na.rm = TRUE),
                     "pe_fixed" = ncol(pmm$data_set$data_raw$x_e_fixed),
                     "pc_fixed" = ncol(pmm$data_set$data_raw$x_c_fixed),
                     "X_e_random" = NULL,
                     "X_c_random" = NULL,
                     "mean_cov_e_random" = NULL,
                     "mean_cov_c_random" = NULL,
                     "pe_random" = 0, "pc_random" = 0,
                     "clus_e" = NULL, "clus_c" = NULL,
                     "n_clus_e" = NULL, "n_clus_c" = NULL,
                     "trt_pos_e" = pmm$data_set$data_raw$trt_pos_e,
                     "trt_pos_c" = pmm$data_set$data_raw$trt_pos_c,
                     "n_trt" = pmm$data_set$data_raw$n,
                     "trt_index" = pmm$data_set$data_raw$trt_index,
                     "trt_lev" = names.trt,
                     "efft" = pmm$data_set$data_raw$efft,
                     "costt" = pmm$data_set$data_raw$costt,
                     "m_efft" = pmm$data_set$data_raw$m_efft,
                     "m_costt" = pmm$data_set$data_raw$m_costt,
                     "dt" = pmm$data_set$data_raw$dt,
                     "n_patternst" = pmm$data_set$data_raw$n_obs_pat,
                     "clus_e_lev" = pmm$data_set$data_raw$clus_e_lev,
                     "clus_c_lev" = pmm$data_set$data_raw$clus_c_lev)
  type <- pmm$type
  dist_e <- pmm$model_output$dist_e
  dist_c <- pmm$model_output$dist_c
  n.iter <- 10
  restriction <- pmm$restriction
  model_info <- list("is_c_random_c" = FALSE, "is_int_c_random_c" = FALSE,
                     "ind_random" = pmm$model_output$ind_random,
                     "ind_fixed" = pmm$model_output$ind_fixed,
                     "model_e_fixed" = "trt",
                     "model_c_fixed" = "trt",
                     "model_e_random" = character(0),
                     "model_c_random" = character(0),
                     "dic" = TRUE,
                     "pd" = FALSE,
                     "n.iter.pd" = 10,
                     "n.adapt" = 100,
                     "inits" = NULL,
                     "n.chains" = 2,
                     "n.iter" = n.iter,
                     "n.burnin" = floor(10 / 2),
                     "n.thin" = 1,
                     "prior" = list("default"),
                     "patterns.prior" = NULL,
                     "n.mci" = 10,
                     "prob" = c(0.025, 0.975),
                     "restriction" = restriction,
                     "progress.bar" = "text",
                     "refresh" = 10000/50,
                     "RNGname" = "Wichmann-Hill",
                     "quiet" = FALSE,
                     "digits" = 3,
                     "save_model" = FALSE)
  run_pmm <- run_pattern(data_model = data_model, type = type ,
                          dist_e = dist_e, dist_c = dist_c,
                          model_info = model_info)

  expect_length(run_pmm, 22)
  expect_type(run_pmm, "list")
  expect_named(run_pmm, c("summary", "model", "mean_effects", "mean_costs",
                         "sd_effects", "sd_costs", "effects_fixed", "costs_fixed",
                         "p_pattern", "effects_random", "costs_random",
                         "imputed", "imputed_trt", "loglik", "type", "method", "ind_fixed",
                         "ind_random", "dist_e", "dist_c", "filein", "d_mod"))

  expect_length(run_pmm$filein, 1)
  expect_type(run_pmm$filein, "character")
  expect_length(run_pmm$dist_c, 1)
  expect_type(run_pmm$dist_c, "character")
  expect_equal(run_pmm$dist_c, dist_c)
  expect_length(run_pmm$dist_e, 1)
  expect_type(run_pmm$dist_e, "character")
  expect_equal(run_pmm$dist_e, dist_e)
  expect_length(run_pmm$ind_random, 1)
  expect_type(run_pmm$ind_random, "logical")
  expect_equal(run_pmm$ind_random, TRUE)
  expect_length(run_pmm$ind_fixed, 1)
  expect_type(run_pmm$ind_fixed, "logical")
  expect_equal(run_pmm$ind_fixed, TRUE)
  expect_length(run_pmm$method, 1)
  expect_type(run_pmm$method, "character")
  expect_equal(run_pmm$method, "PATTERN")
  expect_length(run_pmm$type, 1)
  expect_type(run_pmm$type, "character")
  expect_equal(run_pmm$type, type)
  expect_length(run_pmm$loglik, 3)
  expect_type(run_pmm$loglik, "list")
  expect_named(run_pmm$loglik, c("effects", "costs", "pattern indicators"))
  expect_shape(run_pmm$loglik$effects, dim = c(n.iter, n))
  expect_shape(run_pmm$loglik$costs, dim = c(n.iter, n))
  expect_shape(run_pmm$loglik$`pattern indicators`, dim = c(n.iter, n))
  expect_length(run_pmm$imputed_trt, 2)
  expect_type(run_pmm$imputed_trt, "list")
  expect_named(run_pmm$imputed_trt, c("effects", "costs"))
  expect_length(run_pmm$imputed_trt$effects, 4)
  expect_type(run_pmm$imputed_trt$effects, "list")
  expect_named(run_pmm$imputed_trt$effects, c("avg", "ql", "qu", "index"))
  expect_length(run_pmm$imputed_trt$costs, 4)
  expect_type(run_pmm$imputed_trt$costs, "list")
  expect_named(run_pmm$imputed_trt$costs, c("avg", "ql", "qu", "index"))
  expect_length(run_pmm$imputed_trt$effects$avg, n.trt)
  expect_type(run_pmm$imputed_trt$effects$avg, "list")
  expect_named(run_pmm$imputed_trt$effects$avg, names.trt)
  expect_length(run_pmm$imputed_trt$effects$ql, n.trt)
  expect_type(run_pmm$imputed_trt$effects$ql, "list")
  expect_named(run_pmm$imputed_trt$effects$ql, names.trt)
  expect_length(run_pmm$imputed_trt$effects$qu, n.trt)
  expect_type(run_pmm$imputed_trt$effects$qu, "list")
  expect_named(run_pmm$imputed_trt$effects$qu, names.trt)
  expect_length(run_pmm$imputed_trt$effects$index, n.trt)
  expect_type(run_pmm$imputed_trt$effects$index, "list")
  expect_named(run_pmm$imputed_trt$effects$index, names.trt)
  expect_vector(run_pmm$imputed_trt$effects$avg$`1`)
  expect_vector(run_pmm$imputed_trt$effects$avg$`2`)
  expect_vector(run_pmm$imputed_trt$effects$ql$`1`)
  expect_vector(run_pmm$imputed_trt$effects$ql$`2`)
  expect_vector(run_pmm$imputed_trt$effects$qu$`1`)
  expect_vector(run_pmm$imputed_trt$effects$qu$`2`)
  expect_vector(run_pmm$imputed_trt$effects$index$`1`)
  expect_vector(run_pmm$imputed_trt$effects$index$`2`)
  expect_length(run_pmm$imputed_trt$costs$avg, n.trt)
  expect_type(run_pmm$imputed_trt$costs$avg, "list")
  expect_named(run_pmm$imputed_trt$costs$avg, names.trt)
  expect_length(run_pmm$imputed_trt$costs$ql, n.trt)
  expect_type(run_pmm$imputed_trt$costs$ql, "list")
  expect_named(run_pmm$imputed_trt$costs$ql, names.trt)
  expect_length(run_pmm$imputed_trt$costs$qu, n.trt)
  expect_type(run_pmm$imputed_trt$costs$qu, "list")
  expect_named(run_pmm$imputed_trt$costs$qu, names.trt)
  expect_length(run_pmm$imputed_trt$costs$index, n.trt)
  expect_type(run_pmm$imputed_trt$costs$index, "list")
  expect_named(run_pmm$imputed_trt$costs$index, names.trt)
  expect_vector(run_pmm$imputed_trt$costs$avg$`1`)
  expect_vector(run_pmm$imputed_trt$costs$avg$`2`)
  expect_vector(run_pmm$imputed_trt$costs$ql$`1`)
  expect_vector(run_pmm$imputed_trt$costs$ql$`2`)
  expect_vector(run_pmm$imputed_trt$costs$qu$`1`)
  expect_vector(run_pmm$imputed_trt$costs$qu$`2`)
  expect_vector(run_pmm$imputed_trt$costs$index$`1`)
  expect_vector(run_pmm$imputed_trt$costs$index$`2`)
  expect_length(run_pmm$imputed, 2)
  expect_type(run_pmm$imputed, "list")
  expect_named(run_pmm$imputed, c("effects", "costs"))
  expect_length(run_pmm$imputed$effects, 4)
  expect_type(run_pmm$imputed$effects, "list")
  expect_named(run_pmm$imputed$effects, c("avg", "ql", "qu", "index"))
  expect_vector(run_pmm$imputed$effects$avg)
  expect_vector(run_pmm$imputed$effects$ql)
  expect_vector(run_pmm$imputed$effects$qu)
  expect_vector(run_pmm$imputed$effects$index)
  expect_length(run_pmm$imputed$costs, 4)
  expect_type(run_pmm$imputed$costs, "list")
  expect_named(run_pmm$imputed$costs, c("avg", "ql", "qu", "index"))
  expect_vector(run_pmm$imputed$costs$avg)
  expect_vector(run_pmm$imputed$costs$ql)
  expect_vector(run_pmm$imputed$costs$qu)
  expect_vector(run_pmm$imputed$costs$index)
  expect_null(run_pmm$costs_random)
  expect_null(run_pmm$effects_random)
  expect_shape(run_pmm$costs_fixed, nrow = n.iter)
  expect_shape(run_pmm$effects_fixed, nrow = n.iter)
  expect_shape(run_pmm$sd_costs, dim = c(n.iter, n.pattern))
  expect_shape(run_pmm$sd_effects, dim = c(n.iter, n.pattern))
  expect_shape(run_pmm$mean_costs, dim = c(n.iter, n.trt))
  expect_shape(run_pmm$mean_effects, dim = c(n.iter, n.trt))
  expect_shape(run_pmm$p_pattern, dim = c(n.iter, n.pattern))
  expect_vector(run_pmm$d_mod)
  expect_length(run_pmm$d_mod, n)
  expect_true(is.matrix(run_pmm$summary))
  expect_contains(colnames(run_pmm$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                              "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_pmm$summary, ncol = 9)
  expect_length(run_pmm$model, 6)
  expect_type(run_pmm$model, "list")
  expect_named(run_pmm$model, c("model", "BUGSoutput", "parameters.to.save",
                               "model.file", "n.iter", "DIC"))
  expect_true(run_pmm$model$DIC)
  expect_equal(run_pmm$model$n.iter, n.iter)
  expect_type(run_pmm$model$parameters.to.save, "character")
  expect_vector(run_pmm$model$parameters.to.save)
  expect_type(run_pmm$model$model, "list")
  expect_named(run_pmm$model$model, c("ptr", "data", "model", "state", "nchain",
                                     "iter", "sync", "recompile"))
  expect_type(run_pmm$model$BUGSoutput, "list")
  expect_named(run_pmm$model$BUGSoutput, c("n.chains", "n.iter", "n.burnin", "n.thin", "n.keep",
                                          "n.sims", "sims.array", "sims.list", "sims.matrix",
                                          "summary", "mean", "sd", "median", "root.short",
                                          "long.short", "dimension.short", "indexes.short",
                                          "last.values", "program", "model.file", "isDIC",
                                          "DICbyR", "pV", "DIC", "time2run" ))
  expect_gt(run_pmm$model$BUGSoutput$time2run, 0)
  expect_length(run_pmm$model$BUGSoutput$DIC, 1)
  expect_length(run_pmm$model$BUGSoutput$pV, 1)
  expect_equal(run_pmm$model$BUGSoutput$DICbyR, run_pmm$model$DIC)
  expect_equal(run_pmm$model$BUGSoutput$isDIC, run_pmm$model$DIC)
  expect_equal(run_pmm$model$BUGSoutput$program, "jags")
  expect_type(run_pmm$model$BUGSoutput$last.values, "list")
  expect_type(run_pmm$model$BUGSoutput$indexes.short, "list")
  expect_vector(run_pmm$model$BUGSoutput$dimension.short)
  expect_type(run_pmm$model$BUGSoutput$long.short, "list")
  expect_vector(run_pmm$model$BUGSoutput$root.short)
  expect_type(run_pmm$model$BUGSoutput$median, "list")
  expect_type(run_pmm$model$BUGSoutput$sd, "list")
  expect_type(run_pmm$model$BUGSoutput$mean, "list")
  expect_true(is.matrix(run_pmm$model$BUGSoutput$summary))
  expect_contains(colnames(run_pmm$model$BUGSoutput$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                                               "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_pmm$model$BUGSoutput$summary, ncol = 9)
  expect_shape(run_pmm$model$BUGSoutput$sims.matrix, nrow = n.iter)
  expect_type(run_pmm$model$BUGSoutput$sims.list, "list")
  expect_shape(run_pmm$model$BUGSoutput$sims.array, dim = c(n.iter - run_pmm$model$BUGSoutput$n.burnin,
                                                           run_pmm$model$BUGSoutput$n.chains,
                                                           dim(run_pmm$model$BUGSoutput$sims.matrix)[2]))
  expect_equal(run_pmm$model$BUGSoutput$n.sims, n.iter)
  expect_equal(run_pmm$model$BUGSoutput$n.keep, n.iter - run_pmm$model$BUGSoutput$n.burnin)
  expect_equal(run_pmm$model$BUGSoutput$n.iter, n.iter)



  load(test_path("testdata/hm.RData"))

  n <- sum(hm$data_set$data_raw$n)
  n.trt <- length(names(hm$data_set$data_raw$n))
  names.trt <- names(hm$data_set$data_raw$n)
  se <-1; sc <- 0
  data_model <- list("n" = n,
                     "eff" = hm$data_set$data_raw$e, "cost" = hm$data_set$data_raw$c,
                     "m_eff" = hm$data_set$data_raw$me, "m_cost" = hm$data_set$data_raw$mc,
                     "s_eff" = hm$data_set$data_raw$se, "s_cost" = hm$data_set$data_raw$sc,
                     "se" = se, "sc" = sc,
                     "sde" = hm$se.prior, "sdc" = hm$sc.prior,
                     "X_e_fixed" = hm$data_set$data_raw$x_e_fixed,
                     "X_c_fixed" = hm$data_set$data_raw$x_c_fixed,
                     "Z_e_fixed" = hm$data_set$data_raw$z_e_fixed,
                     "Z_c_fixed" = hm$data_set$data_raw$z_c_fixed,
                     "mean_cov_e_fixed" = apply(hm$data_set$data_raw$x_e_fixed, 2, mean, na.rm = TRUE),
                     "mean_cov_c_fixed" = apply(hm$data_set$data_raw$x_c_fixed, 2, mean, na.rm = TRUE),
                     "mean_z_e_fixed" = apply(hm$data_set$data_raw$z_e_fixed, 2, mean, na.rm = TRUE),
                     "mean_z_c_fixed" = apply(hm$data_set$data_raw$z_c_fixed, 2, mean, na.rm = TRUE),
                     "pe_fixed" = ncol(hm$data_set$data_raw$x_e_fixed),
                     "pc_fixed" = ncol(hm$data_set$data_raw$x_c_fixed),
                     "ze_fixed" = ncol(hm$data_set$data_raw$z_e_fixed),
                     "zc_fixed" = ncol(hm$data_set$data_raw$z_c_fixed),
                     "X_e_random" = NULL,
                     "X_c_random" = NULL,
                     "mean_cov_e_random" = NULL,
                     "mean_cov_c_random" = NULL,
                     "pe_random" = 0, "pc_random" = 0,
                     "clus_e" = NULL, "clus_c" = NULL,
                     "n_clus_e" = NULL, "n_clus_c" = NULL,
                     "Z_e_random" = NULL, "Z_c_random" = NULL,
                     "mean_z_e_random" = NULL,
                     "mean_z_c_random" = NULL,
                     "ze_random" = 0, "zc_random" = 0,
                     "clus_se" = NULL, "clus_sc" = NULL,
                     "n_clus_se" = NULL, "n_clus_sc" = NULL,
                     "trt_pos_e" = hm$data_set$data_raw$trt_pos_e,
                     "trt_pos_c" = hm$data_set$data_raw$trt_pos_c,
                     "trt_pos_se" = hm$data_set$data_raw$trt_pos_se,
                     "trt_pos_sc" = hm$data_set$data_raw$trt_pos_sc,
                     "n_trt" = hm$data_set$data_raw$n,
                     "trt_index" = hm$data_set$data_raw$trt_index,
                     "trt_lev" = names.trt,
                     "efft" = hm$data_set$data_raw$efft,
                     "costt" = hm$data_set$data_raw$costt,
                     "m_efft" = hm$data_set$data_raw$m_efft,
                     "m_costt" = hm$data_set$data_raw$m_costt,
                     "s_efft" = hm$data_set$data_raw$s_efft,
                     "s_costt" = hm$data_set$data_raw$s_costt,
                     "clus_e_lev" = hm$data_set$data_raw$clus_e_lev,
                     "clus_c_lev" = hm$data_set$data_raw$clus_c_lev,
                     "clus_se_lev" = hm$data_set$data_raw$clus_se_lev,
                     "clus_sc_lev" = hm$data_set$data_raw$clus_sc_lev)
  type <- hm$type
  dist_e <- hm$model_output$dist_e
  dist_c <- hm$model_output$dist_c
  n.iter <- 10
  model_info <- list("is_c_random_c" = FALSE, "is_int_c_random_c" = FALSE,
                     "ind_random" = hm$model_output$ind_random,
                     "ind_fixed" = hm$model_output$ind_fixed,
                     "model_e_fixed" = "trt",
                     "model_c_fixed" = "trt",
                     "model_e_random" = character(0),
                     "model_c_random" = character(0),
                     "model_se_random" = character(0),
                     "model_sc_random" = character(0),
                     "dic" = TRUE,
                     "pd" = FALSE,
                     "n.iter.pd" = 10,
                     "n.adapt" = 100,
                     "inits" = NULL,
                     "n.chains" = 2,
                     "n.iter" = n.iter,
                     "n.burnin" = floor(10 / 2),
                     "n.thin" = 1,
                     "prior" = list("default"),
                     "n.mci" = 10,
                     "prob" = c(0.025, 0.975),
                     "progress.bar" = "text",
                     "refresh" = 10000/50,
                     "RNGname" = "Wichmann-Hill",
                     "quiet" = FALSE,
                     "digits" = 3,
                     "save_model" = FALSE)
  run_hm <- run_hurdle(data_model = data_model, type = type ,
                          dist_e = dist_e, dist_c = dist_c,
                          model_info = model_info)

  expect_length(run_hm, 30)
  expect_type(run_hm, "list")
  expect_named(run_hm, c("summary", "model", "mean_effects", "mean_costs",
                         "mean_nostr_effects", "mean_nostr_costs",
                         "sd_effects", "sd_costs", "effects_fixed", "costs_fixed",
                         "pstr_effects", "str_effects_fixed", "pstr_costs", "str_costs_fixed",
                         "effects_random", "costs_random", "str_effects_random", "str_costs_random",
                         "imputed", "imputed_trt", "loglik", "type", "method", "ind_fixed",
                         "ind_random", "dist_e", "dist_c", "se", "sc", "filein"))

  expect_length(run_hm$filein, 1)
  expect_type(run_hm$filein, "character")
  expect_length(run_hm$dist_c, 1)
  expect_type(run_hm$dist_c, "character")
  expect_equal(run_hm$dist_c, dist_c)
  expect_length(run_hm$dist_e, 1)
  expect_type(run_hm$dist_e, "character")
  expect_equal(run_hm$dist_e, dist_e)
  expect_length(run_hm$ind_random, 1)
  expect_type(run_hm$ind_random, "logical")
  expect_equal(run_hm$ind_random, TRUE)
  expect_length(run_hm$ind_fixed, 1)
  expect_type(run_hm$ind_fixed, "logical")
  expect_equal(run_hm$ind_fixed, TRUE)
  expect_length(run_hm$method, 1)
  expect_type(run_hm$method, "character")
  expect_equal(run_hm$method, "HURDLE")
  expect_length(run_hm$type, 1)
  expect_type(run_hm$type, "character")
  expect_equal(run_hm$type, type)
  expect_length(run_hm$loglik, 4)
  expect_type(run_hm$loglik, "list")
  expect_named(run_hm$loglik, c("effects", "costs",
                                "structural indicators effects", "structural indicators costs" ))
  expect_shape(run_hm$loglik$effects, dim = c(n.iter, n))
  expect_shape(run_hm$loglik$costs, dim = c(n.iter, n))
  expect_shape(run_hm$loglik$`structural indicators effects`, dim = c(n.iter, n))
  expect_shape(run_hm$loglik$`structural indicators costs`, dim = c(n.iter, n))
  expect_length(run_hm$imputed_trt, 2)
  expect_type(run_hm$imputed_trt, "list")
  expect_named(run_hm$imputed_trt, c("effects", "costs"))
  expect_length(run_hm$imputed_trt$effects, 4)
  expect_type(run_hm$imputed_trt$effects, "list")
  expect_named(run_hm$imputed_trt$effects, c("avg", "ql", "qu", "index"))
  expect_length(run_hm$imputed_trt$costs, 4)
  expect_type(run_hm$imputed_trt$costs, "list")
  expect_named(run_hm$imputed_trt$costs, c("avg", "ql", "qu", "index"))
  expect_length(run_hm$imputed_trt$effects$avg, n.trt)
  expect_type(run_hm$imputed_trt$effects$avg, "list")
  expect_named(run_hm$imputed_trt$effects$avg, names.trt)
  expect_length(run_hm$imputed_trt$effects$ql, n.trt)
  expect_type(run_hm$imputed_trt$effects$ql, "list")
  expect_named(run_hm$imputed_trt$effects$ql, names.trt)
  expect_length(run_hm$imputed_trt$effects$qu, n.trt)
  expect_type(run_hm$imputed_trt$effects$qu, "list")
  expect_named(run_hm$imputed_trt$effects$qu, names.trt)
  expect_length(run_hm$imputed_trt$effects$index, n.trt)
  expect_type(run_hm$imputed_trt$effects$index, "list")
  expect_named(run_hm$imputed_trt$effects$index, names.trt)
  expect_vector(run_hm$imputed_trt$effects$avg$`1`)
  expect_vector(run_hm$imputed_trt$effects$avg$`2`)
  expect_vector(run_hm$imputed_trt$effects$ql$`1`)
  expect_vector(run_hm$imputed_trt$effects$ql$`2`)
  expect_vector(run_hm$imputed_trt$effects$qu$`1`)
  expect_vector(run_hm$imputed_trt$effects$qu$`2`)
  expect_vector(run_hm$imputed_trt$effects$index$`1`)
  expect_vector(run_hm$imputed_trt$effects$index$`2`)
  expect_length(run_hm$imputed_trt$costs$avg, n.trt)
  expect_type(run_hm$imputed_trt$costs$avg, "list")
  expect_named(run_hm$imputed_trt$costs$avg, names.trt)
  expect_length(run_hm$imputed_trt$costs$ql, n.trt)
  expect_type(run_hm$imputed_trt$costs$ql, "list")
  expect_named(run_hm$imputed_trt$costs$ql, names.trt)
  expect_length(run_hm$imputed_trt$costs$qu, n.trt)
  expect_type(run_hm$imputed_trt$costs$qu, "list")
  expect_named(run_hm$imputed_trt$costs$qu, names.trt)
  expect_length(run_hm$imputed_trt$costs$index, n.trt)
  expect_type(run_hm$imputed_trt$costs$index, "list")
  expect_named(run_hm$imputed_trt$costs$index, names.trt)
  expect_vector(run_hm$imputed_trt$costs$avg$`1`)
  expect_vector(run_hm$imputed_trt$costs$avg$`2`)
  expect_vector(run_hm$imputed_trt$costs$ql$`1`)
  expect_vector(run_hm$imputed_trt$costs$ql$`2`)
  expect_vector(run_hm$imputed_trt$costs$qu$`1`)
  expect_vector(run_hm$imputed_trt$costs$qu$`2`)
  expect_vector(run_hm$imputed_trt$costs$index$`1`)
  expect_vector(run_hm$imputed_trt$costs$index$`2`)
  expect_length(run_hm$imputed, 2)
  expect_type(run_hm$imputed, "list")
  expect_named(run_hm$imputed, c("effects", "costs"))
  expect_length(run_hm$imputed$effects, 4)
  expect_type(run_hm$imputed$effects, "list")
  expect_named(run_hm$imputed$effects, c("avg", "ql", "qu", "index"))
  expect_vector(run_hm$imputed$effects$avg)
  expect_vector(run_hm$imputed$effects$ql)
  expect_vector(run_hm$imputed$effects$qu)
  expect_vector(run_hm$imputed$effects$index)
  expect_length(run_hm$imputed$costs, 4)
  expect_type(run_hm$imputed$costs, "list")
  expect_named(run_hm$imputed$costs, c("avg", "ql", "qu", "index"))
  expect_vector(run_hm$imputed$costs$avg)
  expect_vector(run_hm$imputed$costs$ql)
  expect_vector(run_hm$imputed$costs$qu)
  expect_vector(run_hm$imputed$costs$index)
  expect_null(run_hm$mis_costs_random)
  expect_null(run_hm$mis_effects_random)
  expect_null(run_hm$costs_random)
  expect_null(run_hm$effects_random)
  expect_shape(run_hm$str_costs_fixed, nrow = n.iter)
  expect_shape(run_hm$str_effects_fixed, nrow = n.iter)
  expect_shape(run_hm$pstr_costs, dim = c(n.iter, 1))
  expect_shape(run_hm$pstr_effects, dim = c(n.iter, 1))
  expect_shape(run_hm$costs_fixed, nrow = n.iter)
  expect_shape(run_hm$effects_fixed, nrow = n.iter)
  expect_shape(run_hm$sd_costs, dim = c(n.iter, 1))
  expect_shape(run_hm$sd_effects, dim = c(n.iter, 1))
  expect_shape(run_hm$mean_costs, dim = c(n.iter, n.trt))
  expect_shape(run_hm$mean_effects, dim = c(n.iter, n.trt))
  expect_true(is.matrix(run_hm$summary))
  expect_contains(colnames(run_hm$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                              "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_hm$summary, ncol = 9)
  expect_length(run_hm$model, 6)
  expect_type(run_hm$model, "list")
  expect_named(run_hm$model, c("model", "BUGSoutput", "parameters.to.save",
                               "model.file", "n.iter", "DIC"))
  expect_true(run_hm$model$DIC)
  expect_equal(run_hm$model$n.iter, n.iter)
  expect_type(run_hm$model$parameters.to.save, "character")
  expect_vector(run_hm$model$parameters.to.save)
  expect_type(run_hm$model$model, "list")
  expect_named(run_hm$model$model, c("ptr", "data", "model", "state", "nchain",
                                     "iter", "sync", "recompile"))
  expect_type(run_hm$model$BUGSoutput, "list")
  expect_named(run_hm$model$BUGSoutput, c("n.chains", "n.iter", "n.burnin", "n.thin", "n.keep",
                                          "n.sims", "sims.array", "sims.list", "sims.matrix",
                                          "summary", "mean", "sd", "median", "root.short",
                                          "long.short", "dimension.short", "indexes.short",
                                          "last.values", "program", "model.file", "isDIC",
                                          "DICbyR", "pV", "DIC", "time2run" ))
  expect_gt(run_hm$model$BUGSoutput$time2run, 0)
  expect_length(run_hm$model$BUGSoutput$DIC, 1)
  expect_length(run_hm$model$BUGSoutput$pV, 1)
  expect_equal(run_hm$model$BUGSoutput$DICbyR, run_hm$model$DIC)
  expect_equal(run_hm$model$BUGSoutput$isDIC, run_hm$model$DIC)
  expect_equal(run_hm$model$BUGSoutput$program, "jags")
  expect_type(run_hm$model$BUGSoutput$last.values, "list")
  expect_type(run_hm$model$BUGSoutput$indexes.short, "list")
  expect_vector(run_hm$model$BUGSoutput$dimension.short)
  expect_type(run_hm$model$BUGSoutput$long.short, "list")
  expect_vector(run_hm$model$BUGSoutput$root.short)
  expect_type(run_hm$model$BUGSoutput$median, "list")
  expect_type(run_hm$model$BUGSoutput$sd, "list")
  expect_type(run_hm$model$BUGSoutput$mean, "list")
  expect_true(is.matrix(run_hm$model$BUGSoutput$summary))
  expect_contains(colnames(run_hm$model$BUGSoutput$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                                               "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_hm$model$BUGSoutput$summary, ncol = 9)
  expect_shape(run_hm$model$BUGSoutput$sims.matrix, nrow = n.iter)
  expect_type(run_hm$model$BUGSoutput$sims.list, "list")
  expect_shape(run_hm$model$BUGSoutput$sims.array, dim = c(n.iter - run_hm$model$BUGSoutput$n.burnin,
                                                           run_hm$model$BUGSoutput$n.chains,
                                                           dim(run_hm$model$BUGSoutput$sims.matrix)[2]))
  expect_equal(run_hm$model$BUGSoutput$n.sims, n.iter)
  expect_equal(run_hm$model$BUGSoutput$n.keep, n.iter - run_hm$model$BUGSoutput$n.burnin)
  expect_equal(run_hm$model$BUGSoutput$n.iter, n.iter)


  load(test_path("testdata/lmm.RData"))

  n <- sum(lmm$data_set$data_raw$n)
  n.trt <- length(names(lmm$data_set$data_raw$n))
  names.trt <- names(lmm$data_set$data_raw$n)
  data_model <- list("n" = n, "n_long" = lmm$data_set$data_raw$n_long,
                     "eff" = lmm$data_set$data_raw$e, "cost" = lmm$data_set$data_raw$c,
                     "eff_long" = lmm$data_set$data_raw$e_long,
                     "cost_long" = lmm$data_set$data_raw$c_long,
                     "time" = lmm$data_set$data_raw$time_long,
                     "m_eff" = lmm$data_set$data_raw$count_me,
                     "m_cost" = lmm$data_set$data_raw$count_mc,
                     "m_eff_long" = lmm$data_set$data_raw$count_me_long,
                     "m_cost_long" = lmm$data_set$data_raw$count_mc_long,
                     "X_e_fixed" = lmm$data_set$data_raw$x_e_fixed,
                     "X_c_fixed" = lmm$data_set$data_raw$x_c_fixed,
                     "Z_e_fixed" = lmm$data_set$data_raw$z_e_fixed,
                     "Z_c_fixed" = lmm$data_set$data_raw$z_c_fixed,
                     "mean_cov_e_fixed" = apply(lmm$data_set$data_raw$x_e_fixed, 2, mean, na.rm = TRUE),
                     "mean_cov_c_fixed" = apply(lmm$data_set$data_raw$x_c_fixed, 2, mean, na.rm = TRUE),
                     "mean_z_e_fixed" = apply(lmm$data_set$data_raw$z_e_fixed, 2, mean, na.rm = TRUE),
                     "mean_z_c_fixed" = apply(lmm$data_set$data_raw$z_c_fixed, 2, mean, na.rm = TRUE),
                     "pe_fixed" = ncol(lmm$data_set$data_raw$x_e_fixed),
                     "pc_fixed" = ncol(lmm$data_set$data_raw$x_c_fixed),
                     "ze_fixed" = ncol(lmm$data_set$data_raw$z_e_fixed),
                     "zc_fixed" = ncol(lmm$data_set$data_raw$z_c_fixed),
                     "X_e_random" = NULL,
                     "X_c_random" = NULL,
                     "mean_cov_e_random" = NULL,
                     "mean_cov_c_random" = NULL,
                     "pe_random" = 0, "pc_random" = 0,
                     "clus_e" = NULL, "clus_c" = NULL,
                     "n_clus_e" = NULL, "n_clus_c" = NULL,
                     "Z_e_random" = NULL, "Z_c_random" = NULL,
                     "mean_z_e_random" = NULL,
                     "mean_z_c_random" = NULL,
                     "ze_random" = 0, "zc_random" = 0,
                     "clus_me" = NULL, "clus_mc" = NULL,
                     "n_clus_me" = NULL, "n_clus_mc" = NULL,
                     "trt_pos_e" = lmm$data_set$data_raw$trt_pos_e,
                     "trt_pos_c" = lmm$data_set$data_raw$trt_pos_c,
                     "trt_pos_me" = lmm$data_set$data_raw$trt_pos_me,
                     "trt_pos_mc" = lmm$data_set$data_raw$trt_pos_mc,
                     "n_trt" = lmm$data_set$data_raw$n,
                     "n_trt_long" = lmm$data_set$data_raw$n_long,
                     "trt_index" = lmm$data_set$data_raw$trt_index,
                     "trt_index_long"= lmm$data_set$data_raw$trt_index_long,
                     "trt_lev" = names.trt,
                     "timet" = lmm$data_set$data_raw$timet,
                     "efft" = lmm$data_set$data_raw$efft,
                     "costt" = lmm$data_set$data_raw$costt,
                     "m_efft" = lmm$data_set$data_raw$m_efft,
                     "m_costt" = lmm$data_set$data_raw$m_costt,
                     "count_met" = lmm$data_set$data_raw$count_met,
                     "count_mct" = lmm$data_set$data_raw$count_mct,
                     "clus_e_lev" = lmm$data_set$data_raw$clus_e_lev,
                     "clus_c_lev" = lmm$data_set$data_raw$clus_c_lev,
                     "clus_me_lev" = lmm$data_set$data_raw$clus_me_lev,
                     "clus_mc_lev" = lmm$data_set$data_raw$clus_mc_lev)
  type <- lmm$type
  dist_e <- lmm$model_output$dist_e
  dist_c <- lmm$model_output$dist_c
  n.iter <- 10
  time_dep <- lmm$time_dep
  n.time <- max(lmm$data_set$data_raw$time_long)
  model_info <- list("is_c_random_c" = FALSE, "is_int_c_random_c" = FALSE,
                     "is_me_random_e" = FALSE, "is_mc_random_c" = FALSE,
                     "is_int_me_random_e" = FALSE, "is_int_mc_random_c" = FALSE,
                     "ind_random" = lmm$model_output$ind_random,
                     "ind_fixed" = lmm$model_output$ind_fixed,
                     "ind_time_fixed" = FALSE,
                     "model_e_fixed" = "trt",
                     "model_c_fixed" = "trt",
                     "time_dep" = time_dep,
                     "model_e_random" = character(0),
                     "model_c_random" = character(0),
                     "model_me_random" = character(0),
                     "model_mc_random" = character(0),
                     "dic" = TRUE,
                     "pd" = FALSE,
                     "n.iter.pd" = 10,
                     "n.adapt" = 100,
                     "inits" = NULL,
                     "n.chains" = 2,
                     "n.iter" = n.iter,
                     "n.burnin" = floor(10 / 2),
                     "n.thin" = 1,
                     "prior" = list("default"),
                     "n.mci" = 10,
                     "prob" = c(0.025, 0.975),
                     "progress.bar" = "text",
                     "refresh" = 10000/50,
                     "RNGname" = "Wichmann-Hill",
                     "quiet" = FALSE,
                     "digits" = 3,
                     "save_model" = FALSE)
  run_lmm <- run_lmdm(data_model = data_model, type = type ,
                          dist_e = dist_e, dist_c = dist_c,
                          model_info = model_info)
  unlink("lmdm.txt")

  expect_length(run_lmm, 26)
  expect_type(run_lmm, "list")
  expect_named(run_lmm, c("summary", "model", "mean_effects", "mean_costs",
                         "sd_effects", "sd_costs", "effects_fixed", "costs_fixed",
                         "pmis_effects", "mis_effects_fixed", "pmis_costs", "mis_costs_fixed",
                         "effects_random", "costs_random", "mis_effects_random", "mis_costs_random",
                         "imputed", "imputed_trt", "loglik", "type", "method", "ind_fixed",
                         "ind_random", "dist_e", "dist_c", "filein"))

  expect_length(run_lmm$filein, 1)
  expect_type(run_lmm$filein, "character")
  expect_length(run_lmm$dist_c, 1)
  expect_type(run_lmm$dist_c, "character")
  expect_equal(run_lmm$dist_c, dist_c)
  expect_length(run_lmm$dist_e, 1)
  expect_type(run_lmm$dist_e, "character")
  expect_equal(run_lmm$dist_e, dist_e)
  expect_length(run_lmm$ind_random, 1)
  expect_type(run_lmm$ind_random, "logical")
  expect_equal(run_lmm$ind_random, TRUE)
  expect_length(run_lmm$ind_fixed, 1)
  expect_type(run_lmm$ind_fixed, "logical")
  expect_equal(run_lmm$ind_fixed, FALSE)
  expect_length(run_lmm$method, 1)
  expect_type(run_lmm$method, "character")
  expect_equal(run_lmm$method, "LMDM")
  expect_length(run_lmm$type, 1)
  expect_type(run_lmm$type, "character")
  expect_equal(run_lmm$type, type)
  expect_length(run_lmm$loglik, 4)
  expect_type(run_lmm$loglik, "list")
  expect_named(run_lmm$loglik, c("effects", "costs",
                                "missing indicators effects", "missing indicators costs" ))
  expect_shape(run_lmm$loglik$effects, dim = c(n.iter, n, n.time))
  expect_shape(run_lmm$loglik$costs, dim = c(n.iter, n, n.time))
  expect_shape(run_lmm$loglik$`missing indicators effects`, dim = c(n.iter, n, n.time))
  expect_shape(run_lmm$loglik$`missing indicators costs`, dim = c(n.iter, n, n.time))
  expect_length(run_lmm$imputed_trt, 2)
  expect_type(run_lmm$imputed_trt, "list")
  expect_named(run_lmm$imputed_trt, c("effects", "costs"))
  expect_length(run_lmm$imputed_trt$effects, 4)
  expect_type(run_lmm$imputed_trt$effects, "list")
  expect_named(run_lmm$imputed_trt$effects, c("avg", "ql", "qu", "index"))
  expect_length(run_lmm$imputed_trt$costs, 4)
  expect_type(run_lmm$imputed_trt$costs, "list")
  expect_named(run_lmm$imputed_trt$costs, c("avg", "ql", "qu", "index"))
  expect_length(run_lmm$imputed_trt$effects$avg, n.trt)
  expect_type(run_lmm$imputed_trt$effects$avg, "list")
  expect_named(run_lmm$imputed_trt$effects$avg, names.trt)
  expect_length(run_lmm$imputed_trt$effects$ql, n.trt)
  expect_type(run_lmm$imputed_trt$effects$ql, "list")
  expect_named(run_lmm$imputed_trt$effects$ql, names.trt)
  expect_length(run_lmm$imputed_trt$effects$qu, n.trt)
  expect_type(run_lmm$imputed_trt$effects$qu, "list")
  expect_named(run_lmm$imputed_trt$effects$qu, names.trt)
  expect_length(run_lmm$imputed_trt$effects$index, n.trt)
  expect_type(run_lmm$imputed_trt$effects$index, "list")
  expect_named(run_lmm$imputed_trt$effects$index, names.trt)
  expect_vector(run_lmm$imputed_trt$effects$avg$`1`)
  expect_vector(run_lmm$imputed_trt$effects$avg$`2`)
  expect_vector(run_lmm$imputed_trt$effects$ql$`1`)
  expect_vector(run_lmm$imputed_trt$effects$ql$`2`)
  expect_vector(run_lmm$imputed_trt$effects$qu$`1`)
  expect_vector(run_lmm$imputed_trt$effects$qu$`2`)
  expect_vector(run_lmm$imputed_trt$effects$index$`1`)
  expect_vector(run_lmm$imputed_trt$effects$index$`2`)
  expect_length(run_lmm$imputed_trt$costs$avg, n.trt)
  expect_type(run_lmm$imputed_trt$costs$avg, "list")
  expect_named(run_lmm$imputed_trt$costs$avg, names.trt)
  expect_length(run_lmm$imputed_trt$costs$ql, n.trt)
  expect_type(run_lmm$imputed_trt$costs$ql, "list")
  expect_named(run_lmm$imputed_trt$costs$ql, names.trt)
  expect_length(run_lmm$imputed_trt$costs$qu, n.trt)
  expect_type(run_lmm$imputed_trt$costs$qu, "list")
  expect_named(run_lmm$imputed_trt$costs$qu, names.trt)
  expect_length(run_lmm$imputed_trt$costs$index, n.trt)
  expect_type(run_lmm$imputed_trt$costs$index, "list")
  expect_named(run_lmm$imputed_trt$costs$index, names.trt)
  expect_vector(run_lmm$imputed_trt$costs$avg$`1`)
  expect_vector(run_lmm$imputed_trt$costs$avg$`2`)
  expect_vector(run_lmm$imputed_trt$costs$ql$`1`)
  expect_vector(run_lmm$imputed_trt$costs$ql$`2`)
  expect_vector(run_lmm$imputed_trt$costs$qu$`1`)
  expect_vector(run_lmm$imputed_trt$costs$qu$`2`)
  expect_vector(run_lmm$imputed_trt$costs$index$`1`)
  expect_vector(run_lmm$imputed_trt$costs$index$`2`)
  expect_length(run_lmm$imputed, 2)
  expect_type(run_lmm$imputed, "list")
  expect_named(run_lmm$imputed, c("effects", "costs"))
  expect_length(run_lmm$imputed$effects, 4)
  expect_type(run_lmm$imputed$effects, "list")
  expect_named(run_lmm$imputed$effects, c("avg", "ql", "qu", "index"))
  expect_vector(run_lmm$imputed$effects$avg)
  expect_vector(run_lmm$imputed$effects$ql)
  expect_vector(run_lmm$imputed$effects$qu)
  expect_vector(run_lmm$imputed$effects$index)
  expect_length(run_lmm$imputed$costs, 4)
  expect_type(run_lmm$imputed$costs, "list")
  expect_named(run_lmm$imputed$costs, c("avg", "ql", "qu", "index"))
  expect_vector(run_lmm$imputed$costs$avg)
  expect_vector(run_lmm$imputed$costs$ql)
  expect_vector(run_lmm$imputed$costs$qu)
  expect_vector(run_lmm$imputed$costs$index)
  expect_null(run_lmm$mis_costs_random)
  expect_null(run_lmm$mis_effects_random)
  expect_null(run_lmm$costs_random)
  expect_null(run_lmm$effects_random)
  expect_shape(run_lmm$mis_costs_fixed, nrow = n.iter)
  expect_shape(run_lmm$mis_effects_fixed, nrow = n.iter)
  expect_shape(run_lmm$pmis_costs, dim = c(n.iter, 3, n.time))
  expect_shape(run_lmm$pmis_effects, dim = c(n.iter, 3, n.time))
  expect_type(run_lmm$costs_fixed, "list")
  expect_shape(run_lmm$costs_fixed$beta, dim = c(n.iter, 2, n.time))
  expect_shape(run_lmm$costs_fixed$beta_f, dim = c(n.iter, n.time))
  expect_shape(run_lmm$costs_fixed$beta_te, dim = c(n.iter, n.time))
  expect_shape(run_lmm$costs_fixed$beta_tc, dim = c(n.iter, n.time))
  expect_type(run_lmm$effects_fixed, "list")
  expect_shape(run_lmm$effects_fixed$alpha, dim = c(n.iter, 2, n.time))
  expect_shape(run_lmm$effects_fixed$alpha_te, dim = c(n.iter, n.time))
  expect_shape(run_lmm$effects_fixed$alpha_tc, dim = c(n.iter, n.time))
  expect_shape(run_lmm$sd_costs, dim = c(n.iter, n.time))
  expect_shape(run_lmm$sd_effects, dim = c(n.iter, n.time))
  expect_shape(run_lmm$mean_costs, dim = c(n.iter, n.trt, n.time))
  expect_shape(run_lmm$mean_effects, dim = c(n.iter, n.trt, n.time))
  expect_true(is.matrix(run_lmm$summary))
  expect_contains(colnames(run_lmm$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                              "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_lmm$summary, ncol = 9)
  expect_length(run_lmm$model, 6)
  expect_type(run_lmm$model, "list")
  expect_named(run_lmm$model, c("model", "BUGSoutput", "parameters.to.save",
                               "model.file", "n.iter", "DIC"))
  expect_true(run_lmm$model$DIC)
  expect_equal(run_lmm$model$n.iter, n.iter)
  expect_type(run_lmm$model$parameters.to.save, "character")
  expect_vector(run_lmm$model$parameters.to.save)
  expect_type(run_lmm$model$model, "list")
  expect_named(run_lmm$model$model, c("ptr", "data", "model", "state", "nchain",
                                     "iter", "sync", "recompile"))
  expect_type(run_lmm$model$BUGSoutput, "list")
  expect_named(run_lmm$model$BUGSoutput, c("n.chains", "n.iter", "n.burnin", "n.thin", "n.keep",
                                          "n.sims", "sims.array", "sims.list", "sims.matrix",
                                          "summary", "mean", "sd", "median", "root.short",
                                          "long.short", "dimension.short", "indexes.short",
                                          "last.values", "program", "model.file", "isDIC",
                                          "DICbyR", "pV", "DIC", "time2run" ))
  expect_gt(run_lmm$model$BUGSoutput$time2run, 0)
  expect_length(run_lmm$model$BUGSoutput$DIC, 1)
  expect_length(run_lmm$model$BUGSoutput$pV, 1)
  expect_equal(run_lmm$model$BUGSoutput$DICbyR, run_lmm$model$DIC)
  expect_equal(run_lmm$model$BUGSoutput$isDIC, run_lmm$model$DIC)
  expect_equal(run_lmm$model$BUGSoutput$program, "jags")
  expect_type(run_lmm$model$BUGSoutput$last.values, "list")
  expect_type(run_lmm$model$BUGSoutput$indexes.short, "list")
  expect_vector(run_lmm$model$BUGSoutput$dimension.short)
  expect_type(run_lmm$model$BUGSoutput$long.short, "list")
  expect_vector(run_lmm$model$BUGSoutput$root.short)
  expect_type(run_lmm$model$BUGSoutput$median, "list")
  expect_type(run_lmm$model$BUGSoutput$sd, "list")
  expect_type(run_lmm$model$BUGSoutput$mean, "list")
  expect_true(is.matrix(run_lmm$model$BUGSoutput$summary))
  expect_contains(colnames(run_lmm$model$BUGSoutput$summary), c("mean", "sd", "2.5%", "25%", "50%",
                                                               "75%", "97.5%", "Rhat", "n.eff"))
  expect_shape(run_lmm$model$BUGSoutput$summary, ncol = 9)
  expect_shape(run_lmm$model$BUGSoutput$sims.matrix, nrow = n.iter)
  expect_type(run_lmm$model$BUGSoutput$sims.list, "list")
  expect_shape(run_lmm$model$BUGSoutput$sims.array, dim = c(n.iter - run_lmm$model$BUGSoutput$n.burnin,
                                                           run_lmm$model$BUGSoutput$n.chains,
                                                           dim(run_lmm$model$BUGSoutput$sims.matrix)[2]))
  expect_equal(run_lmm$model$BUGSoutput$n.sims, n.iter)
  expect_equal(run_lmm$model$BUGSoutput$n.keep, n.iter - run_lmm$model$BUGSoutput$n.burnin)
  expect_equal(run_lmm$model$BUGSoutput$n.iter, n.iter)


})


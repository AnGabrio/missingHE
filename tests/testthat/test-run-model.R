test_that("check output from running model internal functions", {

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
                     "quiet" = FALSE)
  run_sm <- run_selection(data_model = data_model, type = type ,
                            dist_e = dist_e, dist_c = dist_c,
                            model_info = model_info)
  unlink("selection.txt")

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
  expect_equal(run_sm$filein, "selection.txt")
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
  expect_equal(run_sm$model$model.file, "selection.txt")
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
  expect_equal(run_sm$model$BUGSoutput$model.file, sm$model_output$filein)
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


})


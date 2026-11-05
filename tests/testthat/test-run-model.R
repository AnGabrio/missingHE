test_that("check output from running model internal functions", {

  load(test_path("testdata/sm.RData"))

  data_model <- list("n" = sum(sm$data_set$data_raw$n),
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
                     "trt_lev" = names(sm$data_set$data_raw$n),
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
                     "n.iter" = 10,
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



})








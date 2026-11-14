test_that("check output from writing model file and prior updating internal functions", {

  testthat::skip_on_cran()
  testthat::skip_if_not(
    isTRUE(tryCatch(requireNamespace("rjags", quietly = TRUE), error = function(e) FALSE)),
    "JAGS / rjags is not available"
  )

  load(test_path("testdata/sm.RData"))

  type <- sm$type
  dist_e <- sm$model_output$dist_e
  dist_c <- sm$model_output$dist_c
  model_txt_info <- list("pe_fixed" = ncol(sm$model_output$effects_fixed),
                         "pc_fixed" = ncol(sm$model_output$costs_fixed),
                         "ze_fixed" = ncol(sm$model_output$mis_effects_fixed),
                         "zc_fixed" = ncol(sm$model_output$mis_costs_fixed),
                         "ind_fixed" = sm$model_output$ind_fixed,
                         "pe_random" = 0,
                         "pc_random" = 0,
                         "ze_random" = 0,
                         "zc_random" = 0,
                         "ind_random" = sm$model_output$ind_random,
                         "model_e_random" = character(0),
                         "model_c_random" = character(0),
                         "model_me_random" = character(0),
                         "model_mc_random" = character(0),
                         "is_c_random_c" = FALSE,
                         "is_int_c_random_c" = FALSE,
                         "is_me_random_e" = FALSE,
                         "is_mc_random_c" = FALSE,
                         "is_int_me_random_e" = FALSE,
                         "is_int_mc_random_c" = FALSE,
                         "prior" = sm$prior,
                         "save_model" = TRUE)
  write_sm <- write_selection(type = type , dist_e = dist_e, dist_c = dist_c,
                              model_txt_info = model_txt_info)
  sm_txt <- paste(readLines("selection.txt"), collapse = "\n")
  prior_sm <- prior_selection(type = type, dist_e = dist_e, dist_c = dist_c,
                              model_txt_info = model_txt_info,
                              model_string_jags = sm_txt)
  unlink(write_sm)

  expect_length(write_sm, 1)
  expect_vector(write_sm)
  expect_equal(write_sm, "selection.txt")
  expect_length(prior_sm, 1)
  expect_type(prior_sm, "character")


  load(test_path("testdata/pmm.RData"))

  type <- pmm$type
  dist_e <- pmm$model_output$dist_e
  dist_c <- pmm$model_output$dist_c

  model_txt_info <- list("pe_fixed" = ncol(pmm$model_output$effects_fixed),
                         "pc_fixed" = ncol(pmm$model_output$costs_fixed),
                         "n_patterns" = length(pmm$data_set$data_raw$n_mis_pat),
                         "ind_fixed" = pmm$model_output$ind_fixed,
                         "pe_random" = 0,
                         "pc_random" = 0,
                         "ind_random" = pmm$model_output$ind_random,
                         "model_e_random" = character(0),
                         "model_c_random" = character(0),
                         "is_c_random_c" = FALSE,
                         "is_int_c_random_c" = FALSE,
                         "prior" = pmm$prior,
                         "restriction" = pmm$restriction,
                         "d_or" = unique(pmm$data_set$data_raw$d),
                         "save_model" = TRUE)
  write_pmm <- write_pattern(type = type , dist_e = dist_e, dist_c = dist_c,
                              model_txt_info = model_txt_info)
  pmm_txt <- paste(readLines("pattern.txt"), collapse = "\n")
  prior_pmm <- prior_pattern(type = type, dist_e = dist_e, dist_c = dist_c,
                              model_txt_info = model_txt_info,
                              model_string_jags = pmm_txt)
  unlink(write_pmm)

  expect_length(write_pmm, 1)
  expect_vector(write_pmm)
  expect_equal(write_pmm, "pattern.txt")
  expect_length(prior_pmm, 1)
  expect_type(prior_pmm, "character")



  load(test_path("testdata/hm.RData"))

  type <- hm$type
  dist_e <- hm$model_output$dist_e
  dist_c <- hm$model_output$dist_c
  model_txt_info <- list("pe_fixed" = ncol(hm$model_output$effects_fixed),
                         "pc_fixed" = ncol(hm$model_output$costs_fixed),
                         "ze_fixed" = ncol(hm$model_output$str_effects_fixed),
                         "zc_fixed" = ncol(hm$model_output$str_effects_fixed),
                         "ind_fixed" = hm$model_output$ind_fixed,
                         "pe_random" = 0,
                         "pc_random" = 0,
                         "ze_random" = 0,
                         "zc_random" = 0,
                         "sde" = hm$se.prior, "sdc" = hm$sc.prior,
                         "se" = 1, "sc" = 0,
                         "ind_random" = hm$model_output$ind_random,
                         "model_e_random" = character(0),
                         "model_c_random" = character(0),
                         "model_se_random" = character(0),
                         "model_sc_random" = character(0),
                         "is_c_random_c" = FALSE,
                         "is_int_c_random_c" = FALSE,
                         "is_me_random_e" = FALSE,
                         "is_mc_random_c" = FALSE,
                         "is_int_me_random_e" = FALSE,
                         "is_int_mc_random_c" = FALSE,
                         "prior" = hm$prior,
                         "save_model" = TRUE)
  write_hm <- write_hurdle(type = type , dist_e = dist_e, dist_c = dist_c,
                              model_txt_info = model_txt_info)
  hm_txt <- paste(readLines("hurdle.txt"), collapse = "\n")
  prior_hm <- prior_hurdle(type = type, dist_e = dist_e, dist_c = dist_c,
                             model_txt_info = model_txt_info,
                             model_string_jags = hm_txt)
  unlink(write_hm)

  expect_length(write_hm, 1)
  expect_vector(write_hm)
  expect_equal(write_hm, "hurdle.txt")
  expect_length(prior_hm, 1)
  expect_type(prior_hm, "character")



  load(test_path("testdata/lmm.RData"))

  type <- lmm$type
  dist_e <- lmm$model_output$dist_e
  dist_c <- lmm$model_output$dist_c
  model_txt_info <- list("pe_fixed" = dim(lmm$model_output$effects_fixed$alpha)[2],
                         "pc_fixed" = dim(lmm$model_output$costs_fixed$beta)[2],
                         "ze_fixed" = 1,
                         "zc_fixed" = 1,
                         "ind_fixed" = lmm$model_output$ind_fixed,
                         "pe_random" = 0,
                         "pc_random" = 0,
                         "ze_random" = 0,
                         "zc_random" = 0,
                         "ind_time_fixed" = FALSE,
                         "ind_random" = lmm$model_output$ind_random,
                         "model_e_random" = character(0),
                         "model_c_random" = character(0),
                         "model_me_random" = character(0),
                         "model_mc_random" = character(0),
                         "is_c_random_c" = FALSE,
                         "is_int_c_random_c" = FALSE,
                         "is_me_random_e" = FALSE,
                         "is_mc_random_c" = FALSE,
                         "is_int_me_random_e" = FALSE,
                         "is_int_mc_random_c" = FALSE,
                         "prior" = lmm$prior,
                         "time_dep" = lmm$time_dep,
                         "save_model" = TRUE)
  write_lmdm <- write_lmdm(type = type , dist_e = dist_e, dist_c = dist_c,
                              model_txt_info = model_txt_info)
  lmdm_txt <- paste(readLines("lmdm.txt"), collapse = "\n")
  prior_lmdm <- prior_hurdle(type = type, dist_e = dist_e, dist_c = dist_c,
                           model_txt_info = model_txt_info,
                           model_string_jags = lmdm_txt)
  unlink(write_lmdm)

  expect_length(write_lmdm, 1)
  expect_vector(write_lmdm)
  expect_equal(write_lmdm, "lmdm.txt")
  expect_length(prior_lmdm, 1)
  expect_type(prior_lmdm, "character")


})








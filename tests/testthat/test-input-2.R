# test errors in input arguments

test_that("input errors check within inner functions", {

  load(test_path("testdata/MenSS.RData"))

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = 1, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "'model.me' and/or 'model.mc' must be formula objects")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = 1, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "'model.me' and/or 'model.mc' must be formula objects")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = me ~ trt , model.mc = mc ~ trt + (1 | site) + (1 | site), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Random effects can be included in the formula only through a single expression within brackets")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = me ~ trt + sti, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Missing covariates cannot be included in the model")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = me ~ trt + ager, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide names in the formula that correspond to those in the data")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = mes ~ trt, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please set 'me' as the response in `model.me`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mcs ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please set 'mc' as the response in `model.mc`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = me ~ trt + e*trt, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'e' are allowed")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + c*trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'c' are allowed")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = 1, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.eff`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = 1,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.cost`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt + (0 + 1 | site), model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please add a random intercept to `model.me`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt + (1 | site + marital), model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a single clustering variable for each model")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt + (1 | ager), model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.me`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt + (age | site), model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Only covariates defined as fixed effects can be included as random effects")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt + (me | site), model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please remove 'me' from the random effects of 'model.me'")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt + e + (e*trt | site), model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'e' are allowed")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt + (1 | age), model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please define clustering variables as factors")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + (0 + 1 | site), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please add a random intercept to `model.mc`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + (1 | site + marital), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a single clustering variable for each model")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + (1 | ager), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.mc`")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + (age | site), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Only covariates defined as fixed effects can be included as random effects")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + (mc | site), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please remove 'mc' from the random effects of 'model.mc'")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + c + (c*trt | site), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'c' are allowed")

  expect_error(
    data_read_selection(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        model.me = me ~ trt, model.mc = mc ~ trt + (1 | age), cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please define clustering variables as factors")



})

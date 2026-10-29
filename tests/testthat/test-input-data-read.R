# test errors in input arguments of data_read internal functions

test_that("input errors check within data_read functions", {

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



  expect_error(
    data_read_pattern(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = 1, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.eff`")

  expect_error(
    data_read_pattern(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt,
                        cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "MAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = 1,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.cost`")


  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                        model.se = 1, model.sc = sc ~ trt, se = 1, sc = 0,
                        cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                        type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                        fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                        trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "'model.se' and/or 'model.sc' must be formula objects")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = 1, model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "'model.se' and/or 'model.sc' must be formula objects")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + (1 | site) + (1 | site), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Random effects can be included in the formula only through a single expression within brackets")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + sti, model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Missing covariates cannot be included in the model")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + ager, model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide names in the formula that correspond to those in the data")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = ses ~ trt, model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please set 'se' as the response in `model.se`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = scs ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please set 'sc' as the response in `model.sc`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + e*trt, model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'e' are allowed")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + c*trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'c' are allowed")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = 1, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.eff`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = 1,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.cost`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + (0 + 1 | site), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please add a random intercept to `model.se`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + (1 | site + marital), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a single clustering variable for each model")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + (1 | ager), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.se`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + (age | site), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Only covariates defined as fixed effects can be included as random effects")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + (se | site), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please remove 'se' from the random effects of 'model.se'")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + e + (e*trt | site), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'e' are allowed")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt + (1 | age), model.sc = sc ~ trt, se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please define clustering variables as factors")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + (0 + 1 | site), se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please add a random intercept to `model.sc`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + (1 | site + marital), se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a single clustering variable for each model")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + (1 | ager), se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.sc`")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + (age | site), se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Only covariates defined as fixed effects can be included as random effects")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + (sc | site), se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please remove 'sc' from the random effects of 'model.sc'")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + c + (c*trt | site), se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'c' are allowed")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt + (1 | age), se = 1, sc = 0,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please define clustering variables as factors")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt, se = NULL, sc = NULL,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Structural values in at least one outcome variable are required,
         please provide the structural value to be used in 'se' and/or in 'sc'")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt, se = 100, sc = NULL,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please select a structural value in 'se' that is present in the data")

  expect_error(
    data_read_hurdle(data = MenSS, model.eff = e ~ trt + u.0, model.cost = c ~ trt,
                     model.se = se ~ trt, model.sc = sc ~ trt, se = NULL, sc = 0.05,
                     cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                     type = "SAR", center = FALSE, fixed_e = e ~ trt + u.0,
                     fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                     trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please select a structural value in 'sc' that is present in the data")



  load(test_path("testdata/PBS.RData"))

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = 1, model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "'model.me' and/or 'model.mc' must be formula objects")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = ~ trt, model.mc = 1, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "'model.me' and/or 'model.mc' must be formula objects")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt , model.mc = mc ~ trt + (1 | site) + (1 | site), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Random effects can be included in the formula only through a single expression within brackets")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt , model.cost = c ~ trt,
                   model.me = me ~ trt + c, model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Missing covariates cannot be included in the model")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt , model.cost = c ~ trt,
                   model.me = me ~ trt + ager, model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide names in the formula that correspond to those in the data")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt , model.cost = c ~ trt,
                   model.me = mes ~ trt, model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please set 'me' as the response in `model.me`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt , model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mcs ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please set 'mc' as the response in `model.mc`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt , model.cost = c ~ trt,
                   model.me = me ~ trt + e*trt, model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'e' are allowed")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt , model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + c*trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'c' are allowed")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = 1, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.eff`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = 1,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.cost`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt + (0 + 1 | site), model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please add a random intercept to `model.me`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt + (1 | site + marital), model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a single clustering variable for each model")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt + (1 | ager), model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.me`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt + (age | site), model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Only covariates defined as fixed effects can be included as random effects")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt + (me | site), model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please remove 'me' from the random effects of 'model.me'")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt + e + (e*trt | site), model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'e' are allowed")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt + (1 | age), model.mc = mc ~ trt, cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please define clustering variables as factors")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + (0 + 1 | site), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please add a random intercept to `model.mc`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + (1 | site + marital), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a single clustering variable for each model")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + (1 | ager), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please provide a valid clustering variable in `model.mc`")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + (age | site), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Only covariates defined as fixed effects can be included as random effects")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + (mc | site), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please remove 'mc' from the random effects of 'model.mc'")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + c + (c*trt | site), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "No interaction effects for 'c' are allowed")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = me ~ trt, model.mc = mc ~ trt + (1 | age), cov_matrix = PBS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt ,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "Please define clustering variables as factors")

})


# check output generated vai data_read internal functions

test_that("check output of data_read functions", {

  load(test_path("testdata/MenSS.RData"))

  model.eff <- e ~ trt + u.0
  model.cost <- c ~ trt + e
  model.me <- me ~ 1
  model.mc <- mc ~ 1
  type <- "MAR"
  cov_matrix <- MenSS[,c(1,2,5:8,10,12,13)]
  fixed_e <- e ~ trt + u.0
  fixed_c <- c ~ trt + e
  trt_lev <- c(1, 2)
  trt_pos <- c(1, 2)
  n.trt <- length(unique(MenSS$trt))
  n.trt.names <- unique(MenSS$trt)
  n.ind <- nrow(MenSS)

  data_sm <- data_read_selection(data = MenSS, model.eff = model.eff, model.cost = c ~ trt,
                                 model.me = model.me, model.mc = model.mc,
                                 cov_matrix = cov_matrix, type = type, center = FALSE,
                                 fixed_e = e ~ trt + u.0, fixed_c = c ~ trt + e,
                                 random_e = NULL, random_c = NULL,
                                 trt_lev = trt_lev, trt_pos = trt_pos)

  expect_type(data_sm, "list")
  expect_length(data_sm, 2)
  expect_named(data_sm, c("data_raw", "model_formula"))

  expect_length(data_sm$data_raw, 63)
  expect_type(data_sm$data_raw, "list")
  expect_length(data_sm$model_formula, 8)
  expect_type(data_sm$model_formula, "list")
  expect_type(data_sm$data_raw$data, "list")
  expect_vector(data_sm$data_raw$e)
  expect_vector(data_sm$data_raw$c)
  expect_equal(data_sm$data_raw$e, MenSS$e, ignore_attr = TRUE)
  expect_equal(data_sm$data_raw$c, MenSS$c, ignore_attr = TRUE)
  expect_vector(data_sm$data_raw$me)
  expect_vector(data_sm$data_raw$mc)
  expect_length(data_sm$data_raw$me, n.ind)
  expect_length(data_sm$data_raw$mc, n.ind)
  expect_vector(data_sm$data_raw$n_obs_e)
  expect_vector(data_sm$data_raw$n_obs_c)
  expect_length(data_sm$data_raw$n_obs_e, n.trt)
  expect_length(data_sm$data_raw$n_obs_c, n.trt)
  expect_named(data_sm$data_raw$n_obs_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_sm$data_raw$n_obs_c, n.trt.names, ignore.case = TRUE)
  expect_shape(data_sm$data_raw$n_mis_e, dim = n.trt)
  expect_shape(data_sm$data_raw$n_mis_c, dim = n.trt)
  expect_named(data_sm$data_raw$n_mis_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_sm$data_raw$n_mis_c, n.trt.names, ignore.case = TRUE)
  expect_shape(data_sm$data_raw$n, dim = n.trt)
  expect_named(data_sm$data_raw$n, n.trt.names, ignore.case = TRUE)
  expect_type(data_sm$data_raw$cov_fixed_e, "list")
  expect_type(data_sm$data_raw$cov_fixed_c, "list")
  expect_length(data_sm$data_raw$cov_fixed_e, n.trt)
  expect_length(data_sm$data_raw$cov_fixed_c, n.trt)
  expect_type(data_sm$data_raw$avg_cov_fixed_e, "list")
  expect_type(data_sm$data_raw$avg_cov_fixed_c, "list")
  expect_length(data_sm$data_raw$cov_fixed_me, n.trt)
  expect_length(data_sm$data_raw$cov_fixed_mc, n.trt)
  expect_length(data_sm$data_raw$avg_cov_fixed_me, n.trt)
  expect_length(data_sm$data_raw$avg_cov_fixed_mc, n.trt)
  expect_null(data_sm$data_raw$cov_random_e)
  expect_null(data_sm$data_raw$cov_random_c)
  expect_null(data_sm$data_raw$avg_cov_random_e)
  expect_null(data_sm$data_raw$avg_cov_random_c)
  expect_null(data_sm$data_raw$cov_random_me)
  expect_null(data_sm$data_raw$cov_random_mc)
  expect_null(data_sm$data_raw$avg_cov_random_me)
  expect_null(data_sm$data_raw$avg_cov_random_mc)
  expect_null(data_sm$data_raw$clus_e)
  expect_null(data_sm$data_raw$clus_c)
  expect_null(data_sm$data_raw$clus_me)
  expect_null(data_sm$data_raw$clus_mc)
  expect_null(data_sm$data_raw$n_clus_e)
  expect_null(data_sm$data_raw$n_clus_e)
  expect_null(data_sm$data_raw$n_clus_me)
  expect_null(data_sm$data_raw$n_clus_mc)
  expect_length(data_sm$data_raw$trt_pos_e, 1)
  expect_length(data_sm$data_raw$trt_pos_c, 1)
  expect_length(data_sm$data_raw$trt_pos_me, 1)
  expect_length(data_sm$data_raw$trt_pos_mc, 1)
  expect_type(data_sm$data_raw$trt_index, "list")
  expect_length(data_sm$data_raw$trt_index, n.trt)
  expect_named(data_sm$data_raw$trt_index, n.trt.names, ignore.case = TRUE)
  expect_type(data_sm$data_raw$efft, "list")
  expect_type(data_sm$data_raw$costt, "list")
  expect_type(data_sm$data_raw$m_efft, "list")
  expect_type(data_sm$data_raw$m_costt, "list")
  expect_null(data_sm$data_raw$name_clus_e)
  expect_null(data_sm$data_raw$name_clus_c)
  expect_null(data_sm$data_raw$name_clus_me)
  expect_null(data_sm$data_raw$name_clus_mc)
  expect_shape(data_sm$data_raw$x_e_fixed, nrow = n.ind)
  expect_shape(data_sm$data_raw$x_c_fixed, nrow = n.ind)
  expect_null(data_sm$data_raw$x_e_random)
  expect_null(data_sm$data_raw$x_c_random)
  expect_shape(data_sm$data_raw$z_e_fixed, nrow = n.ind)
  expect_shape(data_sm$data_raw$z_c_fixed, nrow = n.ind)
  expect_null(data_sm$data_raw$z_e_random)
  expect_null(data_sm$data_raw$z_c_random)
  expect_null(data_sm$data_raw$clusn_e)
  expect_null(data_sm$data_raw$clusn_c)
  expect_null(data_sm$data_raw$clusn_me)
  expect_null(data_sm$data_raw$clusn_mc)
  expect_null(data_sm$data_raw$clus_e_lev)
  expect_null(data_sm$data_raw$clus_c_lev)
  expect_null(data_sm$data_raw$clus_me_lev)
  expect_null(data_sm$data_raw$clus_mc_lev)
  expect_equal(data_sm$model_formula$mf_model.e_fixed, model.eff, ignore_attr = TRUE)
  expect_equal(data_sm$model_formula$mf_model.c_fixed, model.cost, ignore_attr = TRUE)
  expect_equal(data_sm$model_formula$mf_model.me_fixed, model.me, ignore_attr = TRUE)
  expect_equal(data_sm$model_formula$mf_model.mc_fixed, model.mc, ignore_attr = TRUE)
  expect_equal(data_sm$model_formula$mf_model.e_random, e ~ 0, ignore_attr = TRUE)
  expect_equal(data_sm$model_formula$mf_model.c_random, c ~ 0, ignore_attr = TRUE)
  expect_equal(data_sm$model_formula$mf_model.me_random, me ~ 0, ignore_attr = TRUE)
  expect_equal(data_sm$model_formula$mf_model.mc_random, mc ~ 0, ignore_attr = TRUE)


  MenSSv2 <- MenSS
  MenSSv2$e <- MenSS$e - 0.01
  MenSSv2$c <- MenSS$c + 0.1
  cov_matrix <- MenSSv2[,c(1,2,5:8,10,12,13)]
  n.trt <- length(unique(MenSSv2$trt))
  n.trt.names <- unique(MenSSv2$trt)
  n.ind <- nrow(MenSSv2)

  data_pmm <- data_read_pattern(data = MenSSv2, model.eff = model.eff, model.cost = c ~ trt,
                                 cov_matrix = cov_matrix, type = type, center = FALSE,
                                 fixed_e = e ~ trt + u.0, fixed_c = c ~ trt + e,
                                 random_e = NULL, random_c = NULL,
                                 trt_lev = trt_lev, trt_pos = trt_pos)

  expect_length(data_pmm$data_raw, 45)
  expect_type(data_pmm$data_raw, "list")
  expect_length(data_pmm$model_formula, 4)
  expect_type(data_pmm$model_formula, "list")
  expect_type(data_pmm$data_raw$data, "list")
  expect_vector(data_pmm$data_raw$e)
  expect_vector(data_pmm$data_raw$c)
  expect_equal(data_pmm$data_raw$e, MenSSv2$e, ignore_attr = TRUE)
  expect_equal(data_pmm$data_raw$c, MenSSv2$c, ignore_attr = TRUE)
  expect_vector(data_pmm$data_raw$me)
  expect_vector(data_pmm$data_raw$mc)
  expect_length(data_pmm$data_raw$me, n.ind)
  expect_length(data_pmm$data_raw$mc, n.ind)
  expect_vector(data_pmm$data_raw$n_obs_e)
  expect_vector(data_pmm$data_raw$n_obs_c)
  expect_vector(data_pmm$data_raw$n_obs_pat)
  expect_length(data_pmm$data_raw$n_obs_e, n.trt)
  expect_length(data_pmm$data_raw$n_obs_c, n.trt)
  expect_length(data_pmm$data_raw$n_obs_pat, n.trt)
  expect_named(data_pmm$data_raw$n_obs_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_pmm$data_raw$n_obs_c, n.trt.names, ignore.case = TRUE)
  expect_named(data_pmm$data_raw$n_obs_pat, n.trt.names, ignore.case = TRUE)
  expect_shape(data_pmm$data_raw$n_mis_e, dim = n.trt)
  expect_shape(data_pmm$data_raw$n_mis_c, dim = n.trt)
  expect_length(data_pmm$data_raw$n_mis_pat, n.trt)
  expect_named(data_pmm$data_raw$n_mis_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_pmm$data_raw$n_mis_c, n.trt.names, ignore.case = TRUE)
  expect_named(data_pmm$data_raw$n_mis_pat, n.trt.names, ignore.case = TRUE)
  expect_shape(data_pmm$data_raw$n, dim = n.trt)
  expect_named(data_pmm$data_raw$n, n.trt.names, ignore.case = TRUE)
  expect_type(data_pmm$data_raw$cov_fixed_e, "list")
  expect_type(data_pmm$data_raw$cov_fixed_c, "list")
  expect_length(data_pmm$data_raw$cov_fixed_e, n.trt)
  expect_length(data_pmm$data_raw$cov_fixed_c, n.trt)
  expect_type(data_pmm$data_raw$avg_cov_fixed_e, "list")
  expect_type(data_pmm$data_raw$avg_cov_fixed_c, "list")
  expect_null(data_pmm$data_raw$cov_random_e)
  expect_null(data_pmm$data_raw$cov_random_c)
  expect_null(data_pmm$data_raw$avg_cov_random_e)
  expect_null(data_pmm$data_raw$avg_cov_random_c)
  expect_null(data_pmm$data_raw$clus_e)
  expect_null(data_pmm$data_raw$clus_c)
  expect_null(data_pmm$data_raw$n_clus_e)
  expect_null(data_pmm$data_raw$n_clus_e)
  expect_length(data_pmm$data_raw$trt_pos_e, 1)
  expect_length(data_pmm$data_raw$trt_pos_c, 1)
  expect_length(data_pmm$data_raw$trt_pos_me, 1)
  expect_length(data_pmm$data_raw$trt_pos_mc, 1)
  expect_type(data_pmm$data_raw$trt_index, "list")
  expect_length(data_pmm$data_raw$trt_index, n.trt)
  expect_named(data_pmm$data_raw$trt_index, n.trt.names, ignore.case = TRUE)
  expect_type(data_pmm$data_raw$efft, "list")
  expect_type(data_pmm$data_raw$costt, "list")
  expect_type(data_pmm$data_raw$m_efft, "list")
  expect_type(data_pmm$data_raw$m_costt, "list")
  expect_type(data_pmm$data_raw$dt, "list")
  expect_null(data_pmm$data_raw$name_clus_e)
  expect_null(data_pmm$data_raw$name_clus_c)
  expect_shape(data_pmm$data_raw$x_e_fixed, nrow = n.ind)
  expect_shape(data_pmm$data_raw$x_c_fixed, nrow = n.ind)
  expect_null(data_pmm$data_raw$x_e_random)
  expect_null(data_pmm$data_raw$x_c_random)
  expect_null(data_pmm$data_raw$clusn_e)
  expect_null(data_pmm$data_raw$clusn_c)
  expect_null(data_pmm$data_raw$clus_e_lev)
  expect_null(data_pmm$data_raw$clus_c_lev)
  expect_equal(data_pmm$model_formula$mf_model.e_fixed, model.eff, ignore_attr = TRUE)
  expect_equal(data_pmm$model_formula$mf_model.c_fixed, model.cost, ignore_attr = TRUE)
  expect_equal(data_pmm$model_formula$mf_model.e_random, e ~ 0, ignore_attr = TRUE)
  expect_equal(data_pmm$model_formula$mf_model.c_random, c ~ 0, ignore_attr = TRUE)



  model.se <- se ~ 1
  model.sc <- sc ~ 1
  se <- 1
  sc <- 0
  type <- "SCAR"

  data_hm <- data_read_hurdle(data = MenSS, model.eff = model.eff, model.cost = c ~ trt,
                                 model.se = model.se, model.sc = model.sc, se = se, sc = sc,
                                 cov_matrix = cov_matrix, type = type, center = FALSE,
                                 fixed_e = e ~ trt + u.0, fixed_c = c ~ trt + e,
                                 random_e = NULL, random_c = NULL,
                                 trt_lev = trt_lev, trt_pos = trt_pos)


  expect_length(data_hm$data_raw, 71)
  expect_type(data_hm$data_raw, "list")
  expect_length(data_hm$model_formula, 8)
  expect_type(data_hm$model_formula, "list")
  expect_type(data_hm$data_raw$data, "list")

  expect_vector(data_hm$data_raw$e)
  expect_vector(data_hm$data_raw$c)
  expect_equal(data_hm$data_raw$e, MenSS$e, ignore_attr = TRUE)
  expect_equal(data_hm$data_raw$c, MenSS$c, ignore_attr = TRUE)
  expect_vector(data_hm$data_raw$me)
  expect_vector(data_hm$data_raw$mc)
  expect_length(data_hm$data_raw$me, n.ind)
  expect_length(data_hm$data_raw$mc, n.ind)
  expect_vector(data_hm$data_raw$se)
  expect_vector(data_hm$data_raw$sc)
  expect_length(data_hm$data_raw$se, n.ind)
  expect_length(data_hm$data_raw$sc, n.ind)
  expect_vector(data_hm$data_raw$n_obs_e)
  expect_vector(data_hm$data_raw$n_obs_c)
  expect_length(data_hm$data_raw$n_obs_e, n.trt)
  expect_length(data_hm$data_raw$n_obs_c, n.trt)
  expect_named(data_hm$data_raw$n_obs_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_hm$data_raw$n_obs_c, n.trt.names, ignore.case = TRUE)
  expect_shape(data_hm$data_raw$n_mis_e, dim = n.trt)
  expect_shape(data_hm$data_raw$n_mis_c, dim = n.trt)
  expect_named(data_hm$data_raw$n_mis_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_hm$data_raw$n_mis_c, n.trt.names, ignore.case = TRUE)
  expect_vector(data_hm$data_raw$n_ns_eff)
  expect_vector(data_hm$data_raw$n_ns_cost)
  expect_length(data_hm$data_raw$n_ns_eff, n.trt)
  expect_length(data_hm$data_raw$n_ns_cost, n.trt)
  expect_named(data_hm$data_raw$n_ns_eff, n.trt.names, ignore.case = TRUE)
  expect_named(data_hm$data_raw$n_ns_cost, n.trt.names, ignore.case = TRUE)
  expect_vector(data_hm$data_raw$n_s_eff)
  expect_vector(data_hm$data_raw$n_s_cost)
  expect_length(data_hm$data_raw$n_s_eff, n.trt)
  expect_length(data_hm$data_raw$n_s_cost, n.trt)
  expect_named(data_hm$data_raw$n_s_eff, n.trt.names, ignore.case = TRUE)
  expect_named(data_hm$data_raw$n_s_cost, n.trt.names, ignore.case = TRUE)
  expect_shape(data_hm$data_raw$n, dim = n.trt)
  expect_named(data_hm$data_raw$n, n.trt.names, ignore.case = TRUE)
  expect_type(data_hm$data_raw$cov_fixed_e, "list")
  expect_type(data_hm$data_raw$cov_fixed_c, "list")
  expect_length(data_hm$data_raw$cov_fixed_e, n.trt)
  expect_length(data_hm$data_raw$cov_fixed_c, n.trt)
  expect_type(data_hm$data_raw$avg_cov_fixed_e, "list")
  expect_type(data_hm$data_raw$avg_cov_fixed_c, "list")
  expect_length(data_hm$data_raw$cov_fixed_se, n.trt)
  expect_length(data_hm$data_raw$cov_fixed_sc, n.trt)
  expect_length(data_hm$data_raw$avg_cov_fixed_se, n.trt)
  expect_length(data_hm$data_raw$avg_cov_fixed_sc, n.trt)
  expect_null(data_hm$data_raw$cov_random_e)
  expect_null(data_hm$data_raw$cov_random_c)
  expect_null(data_hm$data_raw$avg_cov_random_e)
  expect_null(data_hm$data_raw$avg_cov_random_c)
  expect_null(data_hm$data_raw$cov_random_se)
  expect_null(data_hm$data_raw$cov_random_sc)
  expect_null(data_hm$data_raw$avg_cov_random_se)
  expect_null(data_hm$data_raw$avg_cov_random_sc)
  expect_null(data_hm$data_raw$clus_e)
  expect_null(data_hm$data_raw$clus_c)
  expect_null(data_hm$data_raw$clus_se)
  expect_null(data_hm$data_raw$clus_sc)
  expect_null(data_hm$data_raw$n_clus_e)
  expect_null(data_hm$data_raw$n_clus_e)
  expect_null(data_hm$data_raw$n_clus_se)
  expect_null(data_hm$data_raw$n_clus_sc)
  expect_length(data_hm$data_raw$trt_pos_e, 1)
  expect_length(data_hm$data_raw$trt_pos_c, 1)
  expect_length(data_hm$data_raw$trt_pos_se, 1)
  expect_length(data_hm$data_raw$trt_pos_sc, 1)
  expect_type(data_hm$data_raw$trt_index, "list")
  expect_length(data_hm$data_raw$trt_index, n.trt)
  expect_named(data_hm$data_raw$trt_index, n.trt.names, ignore.case = TRUE)
  expect_type(data_hm$data_raw$efft, "list")
  expect_type(data_hm$data_raw$costt, "list")
  expect_type(data_hm$data_raw$m_efft, "list")
  expect_type(data_hm$data_raw$m_costt, "list")
  expect_type(data_hm$data_raw$s_efft, "list")
  expect_type(data_hm$data_raw$s_costt, "list")
  expect_null(data_hm$data_raw$name_clus_e)
  expect_null(data_hm$data_raw$name_clus_c)
  expect_null(data_hm$data_raw$name_clus_se)
  expect_null(data_hm$data_raw$name_clus_sc)
  expect_shape(data_hm$data_raw$x_e_fixed, nrow = n.ind)
  expect_shape(data_hm$data_raw$x_c_fixed, nrow = n.ind)
  expect_null(data_hm$data_raw$x_e_random)
  expect_null(data_hm$data_raw$x_c_random)
  expect_shape(data_hm$data_raw$z_e_fixed, nrow = n.ind)
  expect_shape(data_hm$data_raw$z_c_fixed, nrow = n.ind)
  expect_null(data_hm$data_raw$z_e_random)
  expect_null(data_hm$data_raw$z_c_random)
  expect_null(data_hm$data_raw$clusn_e)
  expect_null(data_hm$data_raw$clusn_c)
  expect_null(data_hm$data_raw$clusn_se)
  expect_null(data_hm$data_raw$clusn_sc)
  expect_null(data_hm$data_raw$clus_e_lev)
  expect_null(data_hm$data_raw$clus_c_lev)
  expect_null(data_hm$data_raw$clus_se_lev)
  expect_null(data_hm$data_raw$clus_sc_lev)
  expect_equal(data_hm$model_formula$mf_model.e_fixed, model.eff, ignore_attr = TRUE)
  expect_equal(data_hm$model_formula$mf_model.c_fixed, model.cost, ignore_attr = TRUE)
  expect_equal(data_hm$model_formula$mf_model.se_fixed, model.se, ignore_attr = TRUE)
  expect_equal(data_hm$model_formula$mf_model.sc_fixed, model.sc, ignore_attr = TRUE)
  expect_equal(data_hm$model_formula$mf_model.e_random, e ~ 0, ignore_attr = TRUE)
  expect_equal(data_hm$model_formula$mf_model.c_random, c ~ 0, ignore_attr = TRUE)
  expect_equal(data_hm$model_formula$mf_model.se_random, se ~ 0, ignore_attr = TRUE)
  expect_equal(data_hm$model_formula$mf_model.sc_random, sc ~ 0, ignore_attr = TRUE)



  load(test_path("testdata/PBS.RData"))

  model.eff <- fixed_e <- e ~ trt
  model.cost <- fixed_c <- c ~ trt
  cov_matrix <- PBS[,c(1,2,5:8,10,12,13)]
  type <- "MAR"
  n.time <- max(PBS$time)
  n.ind <- length(unique(PBS$id))
  n.rows <- nrow(PBS)

  data_lmdm <- data_read_lmdm(data = PBS, model.eff = model.eff, model.cost = model.cost,
                   model.me = model.me, model.mc = model.mc, cov_matrix = cov_matrix,
                   type = type, center = FALSE, fixed_e = fixed_e,
                   fixed_c = fixed_c, random_e = NULL, random_c = NULL,
                   trt_lev = trt_lev, trt_pos = trt_pos)

  expect_length(data_lmdm$data_raw, 93)
  expect_type(data_lmdm$data_raw, "list")
  expect_length(data_lmdm$model_formula, 8)
  expect_type(data_lmdm$model_formula, "list")
  expect_type(data_lmdm$data_raw$data, "list")

  expect_vector(data_lmdm$data_raw$e_long)
  expect_vector(data_lmdm$data_raw$c_long)
  expect_equal(data_lmdm$data_raw$e_long, PBS$e, ignore_attr = TRUE)
  expect_equal(data_lmdm$data_raw$c_long, PBS$c, ignore_attr = TRUE)
  expect_vector(data_lmdm$data_raw$e)
  expect_vector(data_lmdm$data_raw$c)
  expect_shape(data_lmdm$data_raw$e, dim = c(n.ind, n.time))
  expect_shape(data_lmdm$data_raw$c, dim = c(n.ind, n.time))
  expect_vector(data_lmdm$data_raw$trt)
  expect_length(data_lmdm$data_raw$trt, n.ind)
  expect_vector(data_lmdm$data_raw$me_long)
  expect_vector(data_lmdm$data_raw$mc_long)
  expect_length(data_lmdm$data_raw$me_long, n.rows)
  expect_length(data_lmdm$data_raw$mc_long, n.rows)
  expect_vector(data_lmdm$data_raw$time_long)
  expect_length(data_lmdm$data_raw$time_long, n.rows)
  expect_vector(data_lmdm$data_raw$n_obs_e)
  expect_vector(data_lmdm$data_raw$n_obs_c)
  expect_length(data_lmdm$data_raw$n_obs_e, n.trt)
  expect_length(data_lmdm$data_raw$n_obs_c, n.trt)
  expect_named(data_lmdm$data_raw$n_obs_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_lmdm$data_raw$n_obs_c, n.trt.names, ignore.case = TRUE)
  expect_vector(data_lmdm$data_raw$n_obs_time_e)
  expect_vector(data_lmdm$data_raw$n_obs_time_c)
  expect_length(data_lmdm$data_raw$n_obs_time_e, n.time)
  expect_length(data_lmdm$data_raw$n_obs_time_c, n.time)
  expect_named(data_lmdm$data_raw$n_obs_time_e, unique(PBS$time), ignore.case = TRUE)
  expect_named(data_lmdm$data_raw$n_obs_time_c, unique(PBS$time), ignore.case = TRUE)
  expect_vector(data_lmdm$data_raw$n_obs_time_trt_e)
  expect_vector(data_lmdm$data_raw$n_obs_time_trt_c)
  expect_length(data_lmdm$data_raw$n_obs_time_trt_e, n.time*n.trt)
  expect_length(data_lmdm$data_raw$n_obs_time_trt_c, n.time*n.trt)
  expect_shape(data_lmdm$data_raw$n_mis_e, dim = n.trt)
  expect_shape(data_lmdm$data_raw$n_mis_c, dim = n.trt)
  expect_named(data_lmdm$data_raw$n_mis_e, n.trt.names, ignore.case = TRUE)
  expect_named(data_lmdm$data_raw$n_mis_c, n.trt.names, ignore.case = TRUE)
  expect_vector(data_lmdm$data_raw$n_mis_time_e)
  expect_vector(data_lmdm$data_raw$n_mis_time_c)
  expect_length(data_lmdm$data_raw$n_mis_time_e, n.time)
  expect_length(data_lmdm$data_raw$n_mis_time_c, n.time)
  expect_named(data_lmdm$data_raw$n_mis_time_e, unique(PBS$time), ignore.case = TRUE)
  expect_named(data_lmdm$data_raw$n_mis_time_c, unique(PBS$time), ignore.case = TRUE)
  expect_vector(data_lmdm$data_raw$n_mis_time_trt_e)
  expect_vector(data_lmdm$data_raw$n_mis_time_trt_c)
  expect_length(data_lmdm$data_raw$n_mis_time_trt_e, n.time*n.trt)
  expect_length(data_lmdm$data_raw$n_mis_time_trt_c, n.time*n.trt)
  expect_length(data_lmdm$data_raw$n_long, n.trt)
  expect_named(data_lmdm$data_raw$n_long, n.trt.names, ignore.case = TRUE)
  expect_length(data_lmdm$data_raw$n_long_time, n.time)
  expect_named(data_lmdm$data_raw$n_long_time, unique(PBS$time), ignore.case = TRUE)
  expect_vector(data_lmdm$data_raw$n_long_time_trt)
  expect_length(data_lmdm$data_raw$n_long_time_trt, n.time*n.trt)
  expect_shape(data_lmdm$data_raw$n, dim = n.trt)
  expect_named(data_lmdm$data_raw$n, n.trt.names, ignore.case = TRUE)
  expect_shape(data_lmdm$data_raw$count_me, dim = c(n.ind, n.time, 3))
  expect_shape(data_lmdm$data_raw$count_mc, dim = c(n.ind, n.time, 3))
  expect_vector(data_lmdm$data_raw$count_me_long)
  expect_vector(data_lmdm$data_raw$count_mc_long)
  expect_length(data_lmdm$data_raw$count_me_long, n.rows)
  expect_length(data_lmdm$data_raw$count_mc_long, n.rows)
  expect_type(data_lmdm$data_raw$cov_fixed_e, "list")
  expect_type(data_lmdm$data_raw$cov_fixed_c, "list")
  expect_length(data_lmdm$data_raw$cov_fixed_e, n.trt)
  expect_length(data_lmdm$data_raw$cov_fixed_c, n.trt)
  expect_type(data_lmdm$data_raw$avg_cov_fixed_e, "list")
  expect_type(data_lmdm$data_raw$avg_cov_fixed_c, "list")
  expect_length(data_lmdm$data_raw$cov_fixed_me, n.trt)
  expect_length(data_lmdm$data_raw$cov_fixed_mc, n.trt)
  expect_length(data_lmdm$data_raw$avg_cov_fixed_me, n.trt)
  expect_length(data_lmdm$data_raw$avg_cov_fixed_mc, n.trt)
  expect_null(data_lmdm$data_raw$cov_random_e)
  expect_null(data_lmdm$data_raw$cov_random_c)
  expect_null(data_lmdm$data_raw$avg_cov_random_e)
  expect_null(data_lmdm$data_raw$avg_cov_random_c)
  expect_null(data_lmdm$data_raw$cov_random_me)
  expect_null(data_lmdm$data_raw$cov_random_mc)
  expect_null(data_lmdm$data_raw$avg_cov_random_me)
  expect_null(data_lmdm$data_raw$avg_cov_random_mc)
  expect_null(data_lmdm$data_raw$clus_e)
  expect_null(data_lmdm$data_raw$clus_c)
  expect_null(data_lmdm$data_raw$clus_me)
  expect_null(data_lmdm$data_raw$clus_mc)
  expect_null(data_lmdm$data_raw$n_clus_e)
  expect_null(data_lmdm$data_raw$n_clus_e)
  expect_null(data_lmdm$data_raw$n_clus_me)
  expect_null(data_lmdm$data_raw$n_clus_mc)
  expect_length(data_lmdm$data_raw$trt_pos_e, 1)
  expect_length(data_lmdm$data_raw$trt_pos_c, 1)
  expect_length(data_lmdm$data_raw$trt_pos_me, 1)
  expect_length(data_lmdm$data_raw$trt_pos_mc, 1)
  expect_type(data_lmdm$data_raw$trt_index_long, "list")
  expect_length(data_lmdm$data_raw$trt_index_long, n.trt)
  expect_named(data_lmdm$data_raw$trt_index_long, n.trt.names, ignore.case = TRUE)
  expect_type(data_lmdm$data_raw$trt_index, "list")
  expect_length(data_lmdm$data_raw$trt_index, n.trt)
  expect_named(data_lmdm$data_raw$trt_index, n.trt.names, ignore.case = TRUE)
  expect_type(data_lmdm$data_raw$efft, "list")
  expect_type(data_lmdm$data_raw$costt, "list")
  expect_type(data_lmdm$data_raw$m_efft, "list")
  expect_type(data_lmdm$data_raw$m_costt, "list")
  expect_type(data_lmdm$data_raw$count_met, "list")
  expect_type(data_lmdm$data_raw$count_mct, "list")
  expect_length(data_lmdm$data_raw$count_met, n.trt)
  expect_length(data_lmdm$data_raw$count_mct, n.trt)
  expect_null(data_lmdm$data_raw$name_clus_e)
  expect_null(data_lmdm$data_raw$name_clus_c)
  expect_null(data_lmdm$data_raw$name_clus_me)
  expect_null(data_lmdm$data_raw$name_clus_mc)
  expect_shape(data_lmdm$data_raw$x_e_fixed, nrow = n.ind)
  expect_shape(data_lmdm$data_raw$x_c_fixed, nrow = n.ind)
  expect_null(data_lmdm$data_raw$x_e_random)
  expect_null(data_lmdm$data_raw$x_c_random)
  expect_shape(data_lmdm$data_raw$z_e_fixed, nrow = n.ind)
  expect_shape(data_lmdm$data_raw$z_c_fixed, nrow = n.ind)
  expect_null(data_lmdm$data_raw$z_e_random)
  expect_null(data_lmdm$data_raw$z_c_random)
  expect_null(data_lmdm$data_raw$clusn_e)
  expect_null(data_lmdm$data_raw$clusn_c)
  expect_null(data_lmdm$data_raw$clusn_me)
  expect_null(data_lmdm$data_raw$clusn_mc)
  expect_null(data_lmdm$data_raw$clus_e_lev)
  expect_null(data_lmdm$data_raw$clus_c_lev)
  expect_null(data_lmdm$data_raw$clus_me_lev)
  expect_null(data_lmdm$data_raw$clus_mc_lev)
  expect_type(data_lmdm$data_raw$eff_time_list, "list")
  expect_type(data_lmdm$data_raw$cost_time_list, "list")
  expect_length(data_lmdm$data_raw$eff_time_list, n.time)
  expect_length(data_lmdm$data_raw$cost_time_list, n.time)
  expect_named(data_lmdm$data_raw$eff_time_list, unique(PBS$time), ignore.case = TRUE)
  expect_named(data_lmdm$data_raw$cost_time_list, unique(PBS$time), ignore.case = TRUE)
  expect_type(data_lmdm$data_raw$time_index_long, "list")
  expect_length(data_lmdm$data_raw$time_index_long, n.time)
  expect_named(data_lmdm$data_raw$time_index_long, unique(PBS$time), ignore.case = TRUE)
  expect_type(data_lmdm$data_raw$eff_time_trt_list, "list")
  expect_type(data_lmdm$data_raw$cost_time_trt_list, "list")
  expect_length(data_lmdm$data_raw$eff_time_trt_list, n.time*n.trt)
  expect_length(data_lmdm$data_raw$cost_time_trt_list, n.time*n.trt)
  expect_named(data_lmdm$data_raw$eff_time_trt_list,
               names(data_lmdm$data_raw$n_long_time_trt), ignore.case = TRUE)
  expect_named(data_lmdm$data_raw$cost_time_trt_list,
               names(data_lmdm$data_raw$n_long_time_trt), ignore.case = TRUE)
  expect_type(data_lmdm$data_raw$time_trt_index_long, "list")
  expect_length(data_lmdm$data_raw$time_trt_index_long, n.time*n.trt)
  expect_named(data_lmdm$data_raw$time_trt_index_long,
               names(data_lmdm$data_raw$n_long_time_trt), ignore.case = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.e_fixed, model.eff, ignore_attr = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.c_fixed, model.cost, ignore_attr = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.me_fixed, model.me, ignore_attr = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.mc_fixed, model.mc, ignore_attr = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.e_random, e ~ 0, ignore_attr = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.c_random, c ~ 0, ignore_attr = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.me_random, me ~ 0, ignore_attr = TRUE)
  expect_equal(data_lmdm$model_formula$mf_model.mc_random, mc ~ 0, ignore_attr = TRUE)


})




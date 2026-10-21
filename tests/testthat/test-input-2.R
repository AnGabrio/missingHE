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
                   model.me = 1, model.mc = mc ~ trt, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
                   type = "MAR", center = FALSE, fixed_e = e ~ trt,
                   fixed_c = c ~ trt + e, random_e = NULL, random_c = NULL,
                   trt_lev = c(1, 2), trt_pos = c(1, 2)),
    regexp = "'model.me' and/or 'model.mc' must be formula objects")

  expect_error(
    data_read_lmdm(data = PBS, model.eff = e ~ trt, model.cost = c ~ trt,
                   model.me = ~ trt, model.mc = 1, cov_matrix = MenSS[,c(1,2,5:8,10,12,13)],
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

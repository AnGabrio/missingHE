#' Diagnostic checks for assessing MCMC convergence of Bayesian models fitted in \code{JAGS} using the function \code{\link{selection}}, \code{\link{pattern}} or \code{\link{hurdle}}
#'
#' The focus is restricted to full Bayesian models in cost-effectiveness analyses based on the function \code{\link{selection}}, \code{\link{pattern}} and \code{\link{hurdle}}, 
#' with convergence of the MCMC chains that is assessed through graphical checks of the posterior distribution of the parameters of interest,
#' Examples are density plots, trace plots, autocorrelation plots, etc. Other types of posterior checks are related to some summary MCMC statistics 
#' that are able to detect possible issues in the convergence of the algorithm, such as the potential scale reduction factor or the effective sample size.
#' Different types of diagnostic tools and statistics are used to assess model convergence using functions contained in the package \strong{ggmcmc} and \strong{mcmcplots}. 
#' Graphics and plots are managed using functions contained in the package \strong{ggplot2} and \strong{ggthemes}.
#' @keywords diagnostics MCMC convergence checks
#' @param x An object of class "missingHE" containing the posterior results of a full Bayesian model implemented using the function \code{\link{selection}}, 
#' \code{\link{pattern}} or \code{\link{hurdle}}
#' @param type Type of diagnostic check to be plotted for the model parameter selected Available choices include: 'histogram' for histogram plots,
#' 'denplot' for density plots, 'traceplot' for trace plots, 'acf' for autocorrelation plots, 'running' for running mean plots,
#' 'compare' for comparing the distribution of the whole chain with only its last part, 'cross' for crosscorrelation plots, 'Rhat' for the potential scale reduction factor, 'geweke' for the geweke diagnostic,
#' 'pairs' for posterior correlation among the parameters,'caterpillar' for caterpillar plots. In addition the class 'summary' provides an overview of some of the most popular
#' diagnostic checks for each parameter selected.
#' @param param Name of the family of parameters to process, as given by a regular expression. For example the mean parameters 
#' for the effect and cost variables can be specified using 'mu.e' and 'mu.c', respectively. Different types of
#' models may have different parameters depending on the assumed distributions and missing data assumptions. 
#' To see a complete list of all possible parameters by types of models assumed see details.
#' @param theme Type of ggplot theme among some pre-defined themes. For a full list of available themes see details.
#' @param ... Additional parameters that can be provided to manage the graphical output of \code{diagnostic}.
#' @return A \strong{ggplot} object containing the plots specified in the argument \code{type}
#' @seealso \code{\link[ggmcmc]{ggs}} \code{\link{selection}} \code{\link{selection}} \code{\link{hurdle}}.
#' @details Depending on the types of plots specified in the argument \code{type}, the output of \code{diagnostic} can produce
#' different combinations of MCMC visual posterior checks for the family of parameters indicated in the argument \code{param}.
#' For a full list of the available plots see the description of the argument \code{type} or see the corresponding plots in the package \strong{ggmcmc}.
#' 
#' The parameters that can be assessed through \code{diagnostic} are only those included in the object \code{x} (see Arguments). Specific character names
#' must be specified in the argument \code{param} according to the specific model implemented. The available names and the parameters associated with them are:
#' \itemize{
#' \item "mu.e" the mean parameters of the effect variables in the two treatment arms.
#' \item "mu.c" the mean parameters of the cost variables in the two treatment arms.
#' \item "mu.e.p" the pattern-specific mean parameters of the effect variables in the two treatment arms (only with the function \code{pattern}).
#' \item "mu.c.p" the pattern-specific mean parameters of the cost variables in the two treatment arms (only with the function \code{pattern}).
#' \item "sd.e" the standard deviation parameters of the effect variables in the two treatment arms.
#' \item "sd.c" the standard deviation parameters of the cost variables in the two treatment arms.
#' \item "alpha" the regression intercept and covariate coefficient parameters for the effect variables in the two treatment arms.
#' \item "beta" the regression intercept and covariate coefficient parameters for the cost variables in the two treatment arms.
#' \item "p.e" the probability parameters of the missingness or structural values mechanism for the effect variables in the two treatment arms 
#' (only with the function \code{selection} or \code{hurdle}).
#' \item "p.c" the probability parameters of the missingness or structural values mechanism for the cost variables in the two treatment arms 
#' (only with the function \code{selection} or \code{hurdle}).
#' \item "gamma.e" the regression intercept and covariate coefficient parameters of the missingness or structural values mechanism
#'  for the effect variables in the two treatment arms (only with the function \code{selection} or \code{hurdle}).
#' \item "gamma.c" the regression intercept and covariate coefficient parameters of the missingness or structural values mechanism 
#'  for the cost variables in the two treatment arms (only with the function \code{selection} or \code{hurdle}).
#' \item "pattern" the probabilities associated with the missingness patterns in the data (only with the function \code{pattern}).
#' \item "delta.e" the mnar parameters of the missingness mechanism for the effect variables in the two treatment arms 
#' (only with the function \code{selection} or \code{pattern}).
#' \item "delta.c" the mnar parameters of the missingness mechanism for the cost variables in the two treatment arms 
#' (only with the function \code{selection} or \code{pattern}).
#' \item "all" all available parameters stored in the object \code{x}.
#' }
#' When the object \code{x} is created using the function \code{pattern}, pattern-specific standard deviation ("sd.e", "sd.c") and regression coefficient 
#' parameters ("alpha", "beta") for both outcomes can be visualised. The parameters associated with a missingness mechanism can be accessed only when \code{x}
#' is created using the function \code{selection} or \code{pattern}, while the parameters associated with the model for the structural values mechanism
#' can be accessed only when \code{x} is created using the function \code{hurdle}.
#' 
#' The argument \code{theme} allows to customise the graphical output of the plots generated by \code{diagnostic} and
#' allows to choose among a set of possible pre-defined themes taken form the package \strong{ggtheme}. Those available can
#' be indicated using one of the following character names: "base","calc","economist","excel","few","538","gdocs","hc","par","pander","solarized","stata","tufte","wsj".
#' 
#' @author Andrea Gabrio
#' @references 
#' Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). \emph{Bayesian Data Analysis, 2nd edition}, CRC Press.
#'
#' Brooks, S. Gelman, A. Jones, JL. Meng, XL. (2011). \emph{Handbook of Markov Chain Monte Carlo}, CRC/Chapman and Hall.
#' @import ggplot2 coda
#' @importFrom stats quantile
#' @export 
#' @examples 
#' #For examples see the function selection, pattern or hurdle
#' #
#' #

diagnostic <- function(x, type = "histogram", param = "all", theme = NULL, ...) {
  exArgs <- list(...)
  if(class(x) != "missingHE") {
    stop("Only objects of class 'missingHE' can be used")
  }
  if(!isTRUE(requireNamespace("ggmcmc")) | !isTRUE(requireNamespace("mcmcplots")) | !isTRUE(requireNamespace("coda")) | !isTRUE(requireNamespace("ggthemes"))) {
    stop("You need to install the R packages 'ggmcmc', 'mcmcplots' and 'coda'. Please run in your R terminal:\n install.packages('ggmcmc', 'mcmcplots', 'coda', 'ggthemes')")
  }
  if(length(theme) != 0) {
    theme_names = c("base", "calc", "economist", "excel", "few", "538", "gdocs", "hc", "par", "pander", "solarized", "stata", "tufte", "wsj")
    if(!theme %in% theme_names) {
      stop("You must provide one of the available theme styles")
    } 
  }
  par_hurdle_sar_e <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "p.e", "gamma.e")
  par_hurdle_sar_c <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "p.c", "gamma.c")
  par_hurdle_sar_ec <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "p.e", "p.c", "gamma.e", "gamma.c")
  par_selection_e <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "p.e", "p.c", "gamma.e", "gamma.c", "delta.e")
  par_selection_c <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "p.e", "p.c", "gamma.e", "gamma.c", "delta.c")
  par_selection_ec <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "p.e", "p.c", "gamma.e", "gamma.c", "delta.e", "delta.c")
  par_pattern <- c("all", "mu.e", "mu.c", "mu.e.p", "mu.c.p","sd.e", "sd.c", "alpha", "beta", "pattern")
  par_pattern_e <- c("all", "mu.e", "mu.c", "mu.e.p", "mu.c.p","sd.e", "sd.c", "alpha", "beta", "pattern", "delta.e")
  par_pattern_c <- c("all", "mu.e", "mu.c", "mu.e.p", "mu.c.p","sd.e", "sd.c", "alpha", "beta", "pattern", "delta.c")
  par_pattern_ec <- c("all", "mu.e", "mu.c", "mu.e.p", "mu.c.p","sd.e", "sd.c", "alpha", "beta", "pattern", "delta.e", "delta.c")
  if(x$model_output$type == "SELECTION_e") {
    if(!param %in% par_selection_e) {stop("You must provide valid parameter names contained in the output of selection") }
  } else if(x$model_output$type == "SELECTION_c") {
    if(!param %in% par_selection_c) {stop("You must provide valid parameter names contained in the output of selection") }
  } else if(x$model_output$type == "SELECTION_ec") {
    if(!param %in% par_selection_ec) {stop("You must provide valid parameter names contained in the output of selection") }
  }
  if(x$model_output$type == "HURDLE_e") {
    if(!param %in% par_hurdle_sar_e) {stop("You must provide valid parameter names contained in the output of hurdle") }
  } else if(x$model_output$type == "HURDLE_c") {
    if(!param %in% par_hurdle_sar_c) {stop("You must provide valid parameter names contained in the output of hurdle") }
  } else if(x$model_output$type == "HURDLE_ec") {
    if(!param %in% par_hurdle_sar_ec) {stop("You must provide valid parameter names contained in the output of hurdle") }
  }
  if(x$model_output$type == "PATTERN") {
    if(!param %in% par_pattern) {stop("You must provide valid parameter names contained in the output of pattern") }
  } else if(x$model_output$type == "PATTERN_e") {
    if(!param %in% par_pattern_e) {stop("You must provide valid parameter names contained in the output of pattern") }
  } else if(x$model_output$type == "PATTERN_c") {
    if(!param %in% par_pattern_c) {stop("You must provide valid parameter names contained in the output of pattern") }
  } else if(x$model_output$type == "PATTERN_ec") {
    if(!param %in% par_pattern_ec) {stop("You must provide valid parameter names contained in the output of pattern") }
  }
  if(length(param) != 1) {
    stop("You can only visualise diagnostic checks for one family of parameters at a time 
         or for all parameters together setting the default value param='all'")
  }
  if(!type %in% c("summary", "histogram", "running", "denplot", "compare", "traceplot", "acf", "cross", "Rhat", "geweke", "caterpillar", "pairs")) {
    stop("Types of diagnostics available for use are 'summary', 'histogram', 'running', 'denplot', 'compare', 'traceplot', 'acf', 'cross', 
         'Rhat', 'geweke', 'caterpillar', 'pairs'")
  }
  labs <- param
  if(length(grep("^SELECTION",x$model_output$type)) == 1 | length(grep("^HURDLE",x$model_output$type)) == 1) {
    labs[pmatch("mu.e",labs)] <- "mu_e"
    labs[pmatch("mu.c",labs)] <- "mu_c"
    labs[pmatch("sd.e",labs)] <- "s_e"
    labs[pmatch("sd.c",labs)] <- "s_c"
    labs[pmatch("alpha",labs)] <- "alpha"
    labs[pmatch("beta",labs)] <- "beta"
    labs[pmatch("p.e",labs)] <- "p_e"
    labs[pmatch("p.c",labs)] <- "p_c"
    labs[pmatch("gamma.e",labs)] <- "gamma_e"
    labs[pmatch("gamma.c",labs)] <- "gamma_c"
  }
  if(length(grep("^SELECTION",x$model_output$type)) == 1) {
    labs[pmatch("delta.e",labs)] <- "delta_e"
    labs[pmatch("delta.c",labs)] <- "delta_c"
  }
  if(length(grep("^PATTERN",x$model_output$type)) == 1) {
    labs[match("mu.e",labs)] <- "mu_e\\[.\\]"
    labs[match("mu.c",labs)] <- "mu_c\\[.\\]"
    labs[match("mu.e.p",labs)] <- "mu_e_p"
    labs[match("mu.c.p",labs)] <- "mu_c_p"
    labs[match("sd.e",labs)] <- "s_e_p"
    labs[match("sd.c",labs)] <- "s_c_p"
    labs[match("alpha",labs)] <- "alpha_p"
    labs[match("beta",labs)] <- "beta_p"
    labs[match("pattern",labs)] <- "p_prob"
    labs[match("delta.e",labs)] <- "Delta_e"
    labs[match("delta.c",labs)] <- "Delta_c"
  }
  mcmc_object <- coda::as.mcmc(x$model_output$`model summary`)
  v_name <- coda::varnames(mcmc_object[, , drop = FALSE])
  check_name_eff <- grepl("eff", v_name)
  check_index_eff <- which(check_name_eff, TRUE)
  check_name_cost <- grepl("cost", v_name)
  check_index_cost <- which(check_name_cost, TRUE)
  check_name_loglik <- grepl("loglik", v_name)
  check_index_loglik <- which(check_name_loglik, TRUE)
  check_name_deviance <- grepl("dev", v_name)
  check_index_deviance <- which(check_name_deviance, TRUE)
  check_index <- c(check_index_cost, check_index_eff, check_index_loglik, check_index_deviance)
  parameters <- v_name[-check_index]
  parameters <- gsub("\\[|\\]", "", parameters)
  parameters <- gsub('[[:digit:]]+', '', parameters)
  parameters <- gsub(",", '', parameters)
  parameters <- paste(unique(parameters))
  if(x$model_output$type == "PATTERN" | x$model_output$type == "PATTERN_e" | x$model_output$type == "PATTERN_c" | x$model_output$type == "PATTERN_ec") {
    parameters <- c(paste(parameters,'1', sep = ""), paste(parameters,'2', sep = ""))
    parameters <- gsub("mu_c1", "mu_c", parameters)
    parameters <- gsub("mu_c2", "mu_c", parameters)
    parameters <- gsub("mu_e1", "mu_e", parameters)
    parameters <- gsub("mu_e2", "mu_e", parameters)
    parameters <- paste(unique(parameters))
    if(x$model_output$type == "PATTERN_e" | x$model_output$type == "PATTERN_ec") {
      parameters <- gsub("Delta_e1", "Delta_e", parameters) 
      parameters <- gsub("Delta_e2", "Delta_e", parameters) 
      parameters <- paste(unique(parameters))
    } 
    if(x$model_output$type == "PATTERN_c" | x$model_output$type == "PATTERN_ec") {
      parameters <- gsub("Delta_c1", "Delta_c", parameters) 
      parameters <- gsub("Delta_c2", "Delta_c", parameters) 
      parameters <- paste(unique(parameters))
    } 
  }
  mcmc_object_subset <- subset(mcmc_object, parameters = parameters)
  if(type == "summary") {
    if(param == "all") {
    parameters.mcmc <- parameters
    } else if(param != "all") {
       parameters.mcmc <- labs
      if(length(grep("^SELECTION",x$model_output$type)) == 1 | length(grep("^HURDLE",x$model_output$type)) == 1) {
        parameters = labs }
      if(length(grep("^PATTERN",x$model_output$type)) == 1) {
        if("mu_e\\[.\\]" %in% labs == TRUE) {
          parameters <- "mu_e" 
          exArgs$leaf.marker = "\\["
        }
        if("mu_c\\[.\\]" %in% labs == TRUE) {
          parameters <- "mu_c" 
          exArgs$leaf.marker = "\\["
        }
        if("mu_e_p" %in% labs == TRUE) {parameters <- c("mu_e_p1", "mu_e_p2") }
        if("mu_c_p" %in% labs == TRUE) {parameters <- c("mu_c_p1", "mu_c_p2") }
        if("s_e_p" %in% labs == TRUE) {parameters <- c("s_e_p1", "s_e_p2") }
        if("s_c_p" %in% labs == TRUE) {parameters <- c("s_c_p1", "s_c_p2") }
        if("alpha_p" %in% labs == TRUE) {parameters <- c("alpha_p1", "alpha_p2") }
        if("beta_p" %in% labs == TRUE) {parameters <- c("beta_p1", "beta_p2") }
        if("p_prob" %in% labs == TRUE) {parameters <- c("p_prob1", "p_prob2") }
        if("Delta_e" %in% labs == TRUE | "Delta_c" %in% labs == TRUE) {parameters <- labs }
        parameters.mcmc <- parameters
      }
    } 
    if(exists("regex", where = exArgs)) {regex = exArgs$regex} else {regex = NULL }
    if(exists("leaf.marker", where = exArgs)) {leaf.marker = exArgs$leaf.marker} else {leaf.marker = "[\\[_]" }
    if(exists("random", where = exArgs)) {random = exArgs$random} else {random = NULL }
    if(exists("dir", where = exArgs)) {dir = exArgs$dir} else {dir = tempdir() }
    if(exists("filename", where = exArgs)) {filename = exArgs$filename} else {filename = "MCMCoutput" }
    if(exists("extension", where = exArgs)) {extension = exArgs$extension} else {extension = "html" }
    if(exists("title", where = exArgs)) {title = exArgs$title} else {title = NULL }
    if(exists("col", where = exArgs)) {col = exArgs$col} else {col = NULL }
    if(exists("lty", where = exArgs)) {lty = exArgs$lty} else {lty = 1 }
    if(exists("xlim", where = exArgs)) {xlim = exArgs$xlim} else {xlim = NULL }
    if(exists("ylim", where = exArgs)) {ylim = exArgs$ylim} else {ylim = NULL }
    if(exists("style", where = exArgs)) {style = exArgs$style} else {style = c("gray", "plain") }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = TRUE }
    ggmcmc_out <- mcmcplots::mcmcplot(mcmc_object_subset, parms = parameters.mcmc, regex = regex, leaf.marker = leaf.marker, random = random, dir = dir, 
                                      filename = filename, extension = extension, title = title, col = col, lty = lty, xlim = xlim, ylim = ylim, 
                                      style = style, greek = greek)
    } else {
      if(param == "all") {
        coda::varnames(mcmc_object_subset) <- paste("model", coda::varnames(mcmc_object_subset), sep=".")
        family <- "model"
      } else {
        family <- labs
      }
    ggmcmc_object <- ggmcmc::ggs(mcmc_object_subset)
  }
  if(type == "histogram") {
    if(exists("bins", where = exArgs)) {bins = exArgs$bins} else {bins = 30 }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_histogram(ggmcmc_object, family = family, bins = bins, greek = greek)
  } else if(type == "denplot") {
    if(exists("rug", where = exArgs)) {rug = exArgs$rug} else {rug = FALSE }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_density(ggmcmc_object, family = family, rug = rug, greek = greek)
  } else if(type == "running") {
    if(exists("original_burnin", where=exArgs)) {original_burnin = exArgs$original_burnin} else {original_burnin = TRUE }
    if(exists("original_thin", where=exArgs)) {original_thin = exArgs$original_thin} else {original_thin = TRUE }
    if(exists("greek", where=exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_running(ggmcmc_object, family = family, original_burnin = original_burnin, original_thin = original_thin, greek = greek)
  } else if(type == "compare") {
    if(exists("partial", where = exArgs)) {partial = exArgs$partial} else {partial = 0.1 }
    if(exists("rug", where = exArgs)) {rug = exArgs$rug} else {rug = FALSE }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_compare_partial(ggmcmc_object, family = family, partial = partial, rug = rug, greek = greek)
  } else if(type == "traceplot") {
    if(exists("original_burnin", where = exArgs)) {original_burnin = exArgs$original_burnin} else {original_burnin = TRUE }
    if(exists("original_thin", where = exArgs)) {original_thin = exArgs$original_thin} else {original_thin = TRUE }
    if(exists("simplify", where = exArgs)) {simplify = exArgs$simplify} else {simplify = NULL }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_traceplot(ggmcmc_object, family = family, original_burnin = original_burnin, original_thin = original_thin, 
                                        simplify = simplify, greek = greek)
  } else if(type == "acf") {
    if(exists("nLags", where = exArgs)) {nLags = exArgs$nLags} else {nLags = 50 }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_autocorrelation(ggmcmc_object, family = family, nLags = nLags, greek = greek)
  } else if(type=="cross") {
    if(exists("absolute_scale", where = exArgs)) {absolute_scale = exArgs$absolute_scale} else {absolute_scale = TRUE }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_crosscorrelation(ggmcmc_object, family = family, absolute_scale = absolute_scale, greek = greek)
  } else if(type == "Rhat") {
    if(exists("scaling", where = exArgs)) {scaling = exArgs$scaling} else {scaling = 1.5 }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_Rhat(ggmcmc_object, family = family, scaling = scaling, greek = greek) + ggplot2::xlab("R_hat")
  } else if(type == "geweke") {
    if(exists("frac1", where = exArgs)) {frac1 = exArgs$frac1} else {frac1 = 0.1 }
    if(exists("frac2", where = exArgs)) {frac2 = exArgs$frac2} else {frac2 = 0.5 }
    if(exists("shadow_limit", where = exArgs)) {shadow_limit = exArgs$shadow_limit} else {shadow_limit = TRUE }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_geweke(ggmcmc_object, family = family, frac1 = frac1, frac2 = frac2, shadow_limit = shadow_limit, greek = greek)
  } else if(type == "caterpillar") {
    if(exists("X", where = exArgs)) {X = exArgs$X} else {X = NA }
    if(exists("thick_ci", where = exArgs)) {thick_ci = exArgs$thick_ci} else {thick_ci = c(0.05, 0.95) }
    if(exists("thin_ci", where = exArgs)) {thin_ci = exArgs$thin_ci} else {thin_ci = c(0.025, 0.975) }
    if(exists("line", where = exArgs)) {line = exArgs$line} else {line = NA }
    if(exists("horizontal", where = exArgs)) {horizontal = exArgs$horizontal} else {horizontal = TRUE }
    if(exists("model_labels", where = exArgs)) {model_labels = exArgs$model_labels} else {model_labels = NULL }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out<-ggmcmc::ggs_caterpillar(ggmcmc_object,family=family, X=X,thick_ci = thick_ci,thin_ci = thin_ci,line = line,horizontal = horizontal,model_labels = model_labels,greek = greek)
  } else if(type == "pairs") {
    if(exists("title", where = exArgs)) {title = exArgs$title} else {title = NULL }
    if(exists("upper", where = exArgs)) {upper = exArgs$upper} else {upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na")}
    if(exists("lower", where = exArgs)) {lower = exArgs$lower} else {lower = list(continuous = "points", combo = "facethist", discrete = "facetbar", na = "na")}
    if(exists("diag", where = exArgs)) {diag = exArgs$diag} else {diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag")}
    if(exists("xlab", where = exArgs)) {xlab = exArgs$xlab} else {xlab = NULL }
    if(exists("ylab", where = exArgs)) {ylab = exArgs$ylab} else {ylab = NULL }
    if(exists("axisLabels", where = exArgs)) {axisLabels = exArgs$axisLabels} else {axisLabels = c("show", "internal", "none") }
    if(exists("labeller", where = exArgs)) {labeller = exArgs$labeller} else {labeller = "label_value" }
    if(exists("showStrips", where = exArgs)) {showStrips = exArgs$showStrips} else {showStrips = NULL }
    if(exists("legend", where = exArgs)) {legend = exArgs$legend} else {legend = NULL }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_pairs(ggmcmc_object, family = family,greek = greek, title = title, upper = upper, lower = lower, diag = diag, 
                                    xlab = xlab, ylab = ylab, axisLabels = axisLabels, labeller = labeller, showStrips = showStrips, legend = legend)
  }
  if(length(theme) != 0 & type != "summary") {
    if(theme == "base") ggmcmc_out <- ggmcmc_out + ggthemes::theme_base()
    if(theme == "calc") ggmcmc_out <- ggmcmc_out + ggthemes::theme_calc()
    if(theme == "economist") ggmcmc_out <- ggmcmc_out + ggthemes::theme_economist()
    if(theme == "excel") ggmcmc_out <- ggmcmc_out + ggthemes::theme_excel()
    if(theme == "few") ggmcmc_out <- ggmcmc_out + ggthemes::theme_few()
    if(theme == "538") ggmcmc_out <- ggmcmc_out + ggthemes::theme_fivethirtyeight()
    if(theme == "gdocs") ggmcmc_out <- ggmcmc_out + ggthemes::theme_gdocs()
    if(theme == "hc") ggmcmc_out <- ggmcmc_out + ggthemes::theme_hc()
    if(theme == "par") ggmcmc_out <- ggmcmc_out + ggthemes::theme_par()
    if(theme == "solarized") ggmcmc_out <- ggmcmc_out + ggthemes::theme_solarized()
    if(theme == "pander") ggmcmc_out <- ggmcmc_out + ggthemes::theme_pander()
    if(theme == "stata") ggmcmc_out <- ggmcmc_out + ggthemes::theme_stata()
    if(theme == "tufte") ggmcmc_out <- ggmcmc_out + ggthemes::theme_tufte()
    if(theme == "wsj") ggmcmc_out <- ggmcmc_out + ggthemes::theme_wsj()
  }
  return(print(ggmcmc_out))
}





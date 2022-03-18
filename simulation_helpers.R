volout0 <- function(vol) {
  0
}
volout1 <- function(vol) {
  0.015 * (100 - vol)
}
volout2 <- function(vol) {
  0.001 * (vol - 100)^2
}
# volume outcome effect:
volouts <- list(v0 = volout0, v1 = volout1, v2 = volout2)

# average hospital volume (in Poisson distribution):
Vol_avg <- 100

# fit model with mgcv::gam:
volout_fit <- function(formula = "O ~ 1",
                       data,
                       provider_var = "provider",
                       volume_var = "Volume",
                       gam_method = "REML") {
  if (inherits(formula, "formula")) {
    formula <- format(formula)
  }
  gam_formula <- as.formula(paste(
    formula,
    # use thin plate regression splines for volume effect:
    " + s(", volume_var, ", bs = 'tp', m = 2)",
    # include random effects:
    " + s(", provider_var, ", bs = 're')"))
  mgcv::gam(formula = gam_formula,
            data    = data,
            family  = binomial(),
            method  = gam_method)
}


# plot estimated smooth effect on probability scale:
smooth_prob_plot <- function(fit, var, ref_prob = NULL,
                             xlim = c(NA, NA), ylim = c(0, NA)) {
  stopifnot(inherits(fit, "gam"))
  stopifnot(is.character(var), var %in% names(fit$var.summary))
  svar <- paste0("s(", var, ")")
  stopifnot(svar %in% names(fit$sp))
  if (is.null(ref_prob)) {
    ref_prob <- mean(mgcv::predict.gam(fit, type = "response"))
  } else {
    stopifnot(is.numeric(ref_prob), length(ref_prob) == 1,
              ! is.na(ref_prob))
  }
  stopifnot(length(xlim) == 2, is.na(xlim) | is.numeric(xlim))
  if (is.na(xlim[1])) {
    xlim[1] <- fit$var.summary[[var]][1]
  }
  if (is.na(xlim[2])) {
    xlim[2] <- fit$var.summary[[var]][3]
  }
  stopifnot(length(ylim) == 2, is.na(ylim) | is.numeric(ylim))
  x <- seq(xlim[1], xlim[2], length = 100)

  plot_data <- data.frame(x = x) # create a data.frame with the right number of rows
  for (v in names(fit$var.summary)) {
    plot_data[[v]] <- fit$var.summary[[v]][1] # This may break for parametric variables that are matrices?
  }
  plot_data[[var]] <- x
  y_clip2 <- ifelse(is.na(ylim[2]), 1, ylim[2])
  A <- mgcv::predict.gam(fit, newdata = plot_data, type = "iterms", unconditional = TRUE, se.fit = TRUE)
  plot_data$y <- plogis(qlogis(ref_prob) + A$fit[, svar])
  plot_data$conf_radius <- A$se.fit[, svar]
  plot_data$CI_lower <-
    plogis(qlogis(ref_prob) + A$fit[, svar] - 2 * plot_data$conf_radius)
  plot_data$CI_upper <-
    pmin(y_clip2,
         plogis(qlogis(ref_prob) + A$fit[, svar] + 2 * plot_data$conf_radius))

  # Colors/denominations for plot
  curve_name <- "estimated association\nwith 95%-confidence-band"
  curve_colour <- "black"
  fill_colour <- rgb(0, 1, 0, 0.2)
  ref_name <- "reference line\n(no effect)"
  ref_size <- 1

  ggplot(data = plot_data, aes(x = !! sym(var))) +
    geom_line(data = plot_data, aes(x = x, y = y, colour = curve_name), size = 2.5) +
    scale_colour_manual(name = " ", values = curve_colour) +
    geom_ribbon(data = plot_data, aes(x = x,
                                      ymin = CI_lower,
                                      ymax = CI_upper, fill = curve_name)) +
    scale_fill_manual(name = " ", values = fill_colour) +
    geom_hline(aes(yintercept = ref_prob, size = ref_name), colour = "blue") +
    scale_size_manual(name = "  ", values = ref_size) +
    scale_y_continuous(
      name = "Probability",
      limits = ylim,
      breaks = scales::pretty_breaks(n = 9),
      labels = scales::label_percent(decimal.mark = ",", suffix = " %", accuracy = 0.01))
}


# extract random intercept standard deviation (with confidence interval) 
# from fitted model:
ri_sd <- function(fit, var = NULL, conf.lev=.95) {
  stopifnot(inherits(fit, "gam"))
  if (is.null(var)) {
    # We look for "s(<varname>, <other args?> bs = 're')" and the like:
    var <- regmatches(format(fit$formula),
                      regexec(
                        "s\\(([^,\\)]*),[^\\)]*bs[[:blank:]]*=[[:blank:]]*(\"re\"|'re')[[:blank:]]*\\)",
                        format(fit$formula)))[[1]][2]
  }
  stopifnot(is.character(var), length(var) == 1, nchar(var) > 0L)
  capture.output(vcomp <- mgcv::gam.vcomp(fit, conf.lev = conf.lev),
                 file = 'NUL')
  if (is.matrix(vcomp)) {
    i <- which(rownames(vcomp) == paste0("s(", var, ")"))
    stopifnot(length(i) == 1L)
    as.data.frame(as.list(vcomp[i, ]))
  } else {
    i <- which(names(vcomp) == paste0("s(", var, ")"))
    stopifnot(length(i) == 1L)
    data.frame(std.dev = vcomp[[i]],
               lower = NA, upper = NA)
  }
}


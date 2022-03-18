library(dplyr)

# directory where simulation results are stored
# and where results of fits will be stored
data_path <- "simdata"

source("./simulation_helpers.R")

pred_plot_data <- function(fit, var, xlim = c(NA, NA)) {
  stopifnot(inherits(fit, "gam"))
  stopifnot(is.character(var), var %in% names(fit$var.summary))
  svar <- paste0("s(", var, ")")
  stopifnot(svar %in% names(fit$sp))
  stopifnot(length(xlim) == 2, is.na(xlim) | is.numeric(xlim))
  if (is.na(xlim[1])) {
    xlim[1] <- fit$var.summary[[var]][1]
  }
  if (is.na(xlim[2])) {
    xlim[2] <- fit$var.summary[[var]][3]
  }
  x <- seq(xlim[1], xlim[2]) # , length = 100)
  # In theory, mgcv::predict.gam can "predict" terms when not all variables are
  # specified (using newdata.guaranteed = TRUE and exclude).  However, this is
  # difficult to set up correctly.  Thus, we use a workaround and use
  # values from var.summary
  plot_data <- data.frame(x = x) # create a data.frame with the right number of rows
  for (v in names(fit$var.summary)) {
    plot_data[[v]] <- fit$var.summary[[v]][1] # This may break for parametric variables that are matrices?
  }
  plot_data[[var]] <- x
  A <- mgcv::predict.gam(fit, newdata = plot_data, type = "iterms", unconditional = TRUE, se.fit = TRUE)
  plot_data$y <- A$fit[, svar]
  plot_data$conf_radius <- A$se.fit[, svar]
  plot_data
}

odds_ratio_CI <- function(fit, dataset, volume_var, vol_pair = c(90, 100)) {
  stopifnot(is.numeric(vol_pair), length(vol_pair) == 2L)
  # build matrix with two observations as described in the paper's appendix
  sub <- dataset[c(1, 1), ]            # ... double the first row of the dataset
  sub[[volume_var]] <- vol_pair        # ... and inject vol_pair

  # predict individual volume effects (logit)
  tmp        <- mgcv::predict.gam(fit, newdata = sub, type = "iterms",
                                  unconditional = FALSE, se.fit = TRUE)
  pred       <- unname(tmp$fit[, paste0("s(",volume_var,")")])

  # compute standard error via delta method
  sub_lp      <- mgcv::predict.gam(fit, newdata = sub, type = "lpmatrix")
  sub_lp[2, ] <- -sub_lp[2,]
  OR          <- exp(pred[1] - pred[2])
  nabla       <- as.vector(colSums(sub_lp)) * OR
  sd          <- as.numeric(sqrt(t(nabla) %*% fit$Vp %*% nabla))

  return(data.frame(lower = OR - 2*sd,
                    OR    = OR,
                    upper = OR + 2*sd,
                    sd    = sd))
}

all_data <- readRDS(file.path(data_path, "alldata.RDS"))
N_hosps <- unique(all_data$N_hosp)
N_hosp_breaks <- if (length(N_hosps) < 8) N_hosps else {
  rev(seq(max(N_hosps), min(N_hosps), -200))
}

all_res <- bind_rows(lapply(
  names(volouts),
  function(volout_name) { # For each volume-outcome relationship
    volout <- volouts[[volout_name]]
    message("> volume outcome relationship ", volout_name)
    v_data <- dplyr::filter(all_data, volout == volout_name)
    v_res <- bind_rows(lapply(
      N_hosps, function(N_hosp) {  # For each number of hospitals
        message(">> number of hospitals: ", N_hosp)
        on.exit(message("\r     \r", appendLF = FALSE))
        N_data <- v_data[v_data$N_hosp == N_hosp, ]
        i_s <- sort(unique(N_data$i))
        pp_data <- data.frame()
        N_res <- bind_rows(lapply(
          i_s,
          function(i) {
            message("\r", i, "/", max(i_s), appendLF = FALSE)
            data <- N_data[N_data$i == i, ]
            # Fit the model
            gam <- volout_fit(y ~ x1 + x2, data, "provider", "Volume")
            s <- summary(gam)
            # When depicting the "true" volume outcome effect, we normalize the predicted
            # probabilities such that the average probability over all of the fit data
            # equals the average probability over the plot (computed for those values of
            # volume from the data).
            ref_prob <- mean(mgcv::predict.gam(gam, type = "response"))
            v_tmp <- volout(data$Volume)
            ref_int <- uniroot(function(x) {mean(plogis(x + v_tmp)) - ref_prob},
                               interval = c(-5, 5))$root

            or_info <- odds_ratio_CI(gam, data, "Volume", vol_pair = c(90, 100))
            pp_data_i <- mutate(pred_plot_data(gam, "Volume"),
                                Vol_eff = volout(Volume), # - mean(volout(data$Volume)),
                                y0 = Vol_eff - mean(v_tmp))
            OR_90_100 = exp(pp_data_i$y[pp_data_i$Volume == 90] -
                              pp_data_i$y[pp_data_i$Volume == 100])
            stopifnot(all.equal(OR_90_100, or_info$OR))
            pp_data <<- bind_rows(
              pp_data,
              mutate(pp_data_i, i = i, p_vol = s$s.pv[1]))
            data.frame(
              p_vol = s$s.pv[1],
              p_re = s$s.pv[2],
              tau = ri_sd(gam)$std.dev,
              beta_0 = s$p.coeff['(Intercept)'],
              beta_1 = s$p.coeff['x1'],
              beta_2 = s$p.coeff['x2'],
              # OR_90_100 = exp(pp_data_i$y[pp_data_i$Volume == 90] -
              #   pp_data_i$y[pp_data_i$Volume == 100]),
              OR_90_100       = or_info$OR,
              OR_90_100_sd    = or_info$sd,
              i = i
            )
          }))
        saveRDS(pp_data, file.path(
          data_path,
          paste0(paste("pp_data", volout_name, N_hosp, "all", sep = "_"),
                 ".RDS")))
        ### To save intermediate results for each N:
        # saveRDS(N_res, file.path(
        #   data_path,
        #   paste0(paste("N_res", volout_name, N_hosp, sep = "_"), ".RDS")))
        mutate(N_res, N_hosp = N_hosp)
      }))
    pvalue_plot_data <- tidyr::pivot_longer(v_res, c("p_vol", "p_re"),
                                            names_to = "statistic",
                                            values_to = "pvalue")
    mutate(v_res, volout = volout_name)
  }))

saveRDS(all_res, file.path(data_path, "allres.RDS"))

capture.output(sessioninfo::session_info(),
               file = file.path(data_path, "sessionInfo_analyze"))

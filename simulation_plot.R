library(dplyr)
library(ggplot2)

options(warn = 1)

# To use system fonts (e.g. Calibri), the following needs to be run once:
# extrafont::font_import()
extrafont::loadfonts(device = "win")

# directory where results are stored
data_path <- "simdata"
plot_path <- "simplots"

stopifnot(dir.exists(data_path))
if (! dir.exists(plot_path)) {
  message("Creating directory '", plot_path, "' for plots.")
  stopifnot(dir.create(plot_path))
}

source("./simulation_helpers.R")

# Set the default theme:
theme_set(
  theme_bw(base_size = 20, base_family = "Calibri") +
    theme(panel.border = element_blank(),
          axis.text = element_text(size = rel(1)),
          axis.ticks = element_blank(),
          axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
          axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
          legend.position = "bottom",
          plot.margin = margin(10, 20, 10, 10),
          strip.background = element_blank())
)
# color-blind-friendly colours
# (from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette)
friendly_red <- "#D55E00"
friendly_green <- "#009E73"
friendly_blue <- "#0072B2"

# Generate a plot of the effect of variable `var`
# on the (logistic) predictor scale
pred_plot_from_data <- function(plot_data, var) {
  # Colors/denominations for plot
  curve_name <- "estimated relationship\nwith confidence band (95%)"
  curve_colour <- "black"
  fill_colour <- rgb(0, 1, 0, 0.2)
  ref_name <- "reference\n(no effect)"
  ref_size <- 1

  ggplot(data = plot_data, aes(x = !! sym(var))) +
    geom_line(data = plot_data, aes(x = x, y = y, colour = curve_name), size = 2.5) +
    scale_colour_manual(name = " ", values = curve_colour) +
    geom_ribbon(data = plot_data,
                aes(x = x,
                    ymin = y - 2 * conf_radius, ymax = y + 2 * conf_radius,
                    fill = curve_name)) +
    scale_fill_manual(name = " ", values = fill_colour) +
    geom_hline(aes(yintercept = 0, size = ref_name), colour = friendly_blue) +
    scale_size_manual(name = "  ", values = ref_size) +
    ylab("volume effect (logistic scale)")
}

all_res <- readRDS(file.path(data_path, "allres.RDS"))
N_hosps <- unique(all_res$N_hosp)
N_hosp_breaks <- if (length(N_hosps) < 8) N_hosps else {
  rev(seq(max(N_hosps), min(N_hosps), -200))
}
(lapply(
  names(volouts),
  function(volout_name) { # For each volume-outcome relationship
    volout <- volouts[[volout_name]]
    message("> volume outcome relationship ", volout_name)
    v_res <- dplyr::filter(all_res, volout == volout_name)
    (lapply(
      N_hosps, function(N_hosp) {  # For each number of hospitals
        message(">> number of hospitals: ", N_hosp)
        on.exit(message("\r     \r", appendLF = FALSE))
        pp_data <- readRDS(file.path(
          data_path,
          paste0(paste("pp_data", volout_name, N_hosp, "all", sep = "_"),
                 ".RDS")))
        N_res <- v_res[v_res$N_hosp == N_hosp, ]
        i_s <- sort(unique(N_res$i))
        (lapply(
          i_s,
          function(i) {
            message("\r", i, "/", max(i_s), appendLF = FALSE)
            pp_data_i <- pp_data[pp_data$i == i, ]
            ###### Plot a single curve
            ggsave(file.path(plot_path, paste0(paste0(
              paste("predplot", volout_name, N_hosp, i, sep = "_"), ".png"))),
                   width = 5.25, height = 5.25,
              pred_plot_from_data(pp_data_i, "Volume") +
                annotate("text", x = mean(pp_data_i$Volume), y = Inf,
                         hjust = 0.5, vjust = 2, size = 5,
                         label = paste0(
                           "p = ", format(pp_data_i$p_vol[1], digits = 3))) +
                (if (volout_name != "v0")
                  geom_line(aes(Volume, y0), colour = friendly_red)) +
                theme(legend.position="none")                       # No legend.
            )
            NULL
          }))
        ###### Plot four curves
        ggsave(file.path(plot_path, paste0(paste0(
          paste("predplot", volout_name, N_hosp, sep = "_"), ".png"))),
          width = 10.5, height = 8,
          pred_plot_from_data(dplyr::filter(pp_data, i <= 4), "Volume") +
            facet_wrap("i", 2, scales = "free",
                       labeller = as_labeller(
                         setNames(
                           paste0("p = ",
                                  sapply(N_res$p_vol, format, digits = 3)),
                           N_res$i)
                         )) +
            (if (volout_name != "v0") geom_line(aes(Volume, y0), colour = friendly_red)) +
            theme(text = element_text(size = 15, family = "Calibri"),
                  legend.text = element_text(size = rel(1)),
                  strip.text = element_text(size = rel(1)))
           # to remove facet labels: strip.text.x = element_blank())
        )
        ###### Plot all curves
        true_data <- mutate(distinct(select(pp_data, Volume, Vol_eff),
                                     Volume, .keep_all = TRUE),
                            y0 = Vol_eff - sum(dpois(Volume, Vol_avg) * Vol_eff))
        ggsave(file.path(plot_path, paste0(paste0(
          paste("predplot", volout_name, N_hosp, "all", sep = "_"), ".png"))),
          width = 5.25, height = 5.25,
          ggplot(pp_data) +
            geom_hline(aes(yintercept = 0, size = "reference\n(no effect)"), colour = friendly_blue) +
            scale_size_manual(name = "  ", values = 1) +
            (if (volout_name != "v0") geom_line(aes(Volume, y0), true_data,
                                                colour = friendly_red, size = 1.5)) +
            geom_line(aes(x, y, colour = p_vol <= 0.05, group = i),
                      alpha = 0.7) +
            ylab("volume effect (logistic scale)") + xlab("Volume") +
            scale_colour_manual(values = c("TRUE" = friendly_green, #"#E69F00",
                                           "FALSE" = "gray80")) +
            theme(legend.position="none",
                  legend.text = element_text(size = rel(1)),
                  strip.text = element_text(size = rel(1)))
       )
        NULL
      }))
    ###### Plot estimates of tau
    ggsave(file.path(plot_path,
                     paste0(paste0(
                       paste("tauplot", volout_name, sep = "_"), ".png"))),
           width = 5.25, height = 5.25,
           ggplot(v_res, aes(as.factor(N_hosp), tau)) +
             geom_hline(yintercept = 0.50, colour = friendly_red) +
             geom_boxplot(outlier.shape = 4, width = 0.4) +
             # geom_point(shape = 4) +
             # geom_jitter(height = 0, width = 0.2) +
             scale_x_discrete("number of hospitals", breaks = N_hosp_breaks) +
             ylab(expression(widehat(tau))) +
             theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
    )
    ###### Plot boxplot of estimated odds ratios
    ggsave(file.path(plot_path,
                     paste0(paste0(
                       paste("ORplot", volout_name, sep = "_"), ".png"))),
           width = 5.25, height = 5.25,
           ggplot(v_res, aes(as.factor(N_hosp), OR_90_100)) +
             geom_hline(yintercept = exp(volout(90) - volout(100)),
                        colour = friendly_red) +
             geom_boxplot(outlier.shape = 4, width = 0.4) +
             ## Use `outlier.shape = NA` if `geom_point` or `geom_jitter` is
             ##used to avoid plotting outliers twice
             # geom_point(shape = 4) +
             # geom_jitter(height = 0, width = 0.2) +
             scale_x_discrete("number of hospitals", breaks = N_hosp_breaks) +
             ylab("OR(90, 100)")
    )
    ###### Plot estimated odds ratios for selected values of N_hosp
    ORplot_Nhosps <- c(100, 200, 500, 1000)
    plot_data <- dplyr::filter(
      ungroup(mutate(group_by(v_res, N_hosp),
                     ID = row_number(desc(OR_90_100)),
                     N_hospf = factor(
                       N_hosp, levels = N_hosps,
                       label = paste0("paste(italic(I),' = ',", N_hosps, ")")))),
      N_hosp %in% ORplot_Nhosps)
    ggsave(file.path(plot_path,
                     paste0(paste0(
                       paste("ORplot_sel", volout_name, sep = "_"), ".png"))),
           width = 10.5, height = 7,
           ggplot(plot_data, aes(x = ID, y = OR_90_100)) +
             facet_wrap("N_hospf", labeller = label_parsed) +
             geom_hline(yintercept = exp(volout(90) - volout(100)),
                        colour = friendly_red) +
             geom_errorbar(aes(ymin = OR_90_100 - 2 * OR_90_100_sd,
                               ymax = OR_90_100 + 2 * OR_90_100_sd)) +
             geom_point() +
             scale_x_continuous(name = NULL, breaks = NULL) +
             ylab("OR(90, 100)")
    )
    ###### Plot p-values of volume effect
    ggsave(file.path(plot_path,
                     paste0(paste0(
                       paste("pvalueplot", volout_name, sep = "_"), ".png"))),
           width = 5.25, height = 5.25,
           ggplot(v_res, aes(as.factor(N_hosp), p_vol)) +
             geom_hline(yintercept = c(0.05, 0.005), linetype = "dashed") +
             geom_boxplot(outlier.shape = 4, width = 0.4) + 
             scale_x_discrete("number of hospitals", breaks = N_hosp_breaks) +
             {if (volout_name == "v0") {
               scale_y_continuous("p-value", limits = c(0, 1))
             } else scale_y_log10("p-value")} +
             scale_shape(name = "Legend", breaks = c("p_vol", "p_re"),
                         labels = c(expression(italic(f) [vol]), expression(tau))) +
             scale_colour_discrete(name = "Legend", breaks = c("p_vol", "p_re"),
                                   labels = c(expression(italic(f) [vol]), expression(tau)))
    )
    NULL
  }))

v0_res <- dplyr::filter(all_res, volout == "v0")
message("Counting false-positive results:")
message("p < 0.05: ", sum(v0_res$p_vol <= 0.05), " times (fraction: ",
        sum(v0_res$p_vol <= 0.05) / nrow(v0_res), ")")
message("p < 0.01: ", sum(v0_res$p_vol <= 0.01), " times (fraction: ",
        sum(v0_res$p_vol <= 0.01) / nrow(v0_res), ")")

capture.output(sessioninfo::session_info(),
               file = file.path(plot_path, "sessionInfo_plot"))

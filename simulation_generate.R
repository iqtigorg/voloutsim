library(dplyr)

set.seed(42L)

# directory where to store the results
data_path <- "simdata"
if (! dir.exists(data_path)) {
  message("Creating directory '", data_path, "' for simulation results.")
  stopifnot(dir.create(data_path))
}

##### Parameters of the simulation

# number of iterations per volume:
n_it <- 50

source("./simulation_helpers.R")

### hospital parameters:

# list of values for number of hospitals:
N_hosps <- seq(100, 1000, 100)
# variance/standard deviation of hospital effects:
tau2 <- 0.25
tau <- sqrt(tau2)
# volume outcome effects are defined in simulation_helpers.R.
# Vol_avg is defined in simulation_helpers.R.

# base risk:
beta_0 <- qlogis(0.1)
# There are two risk factors:
# 1. a discrete risk factor
#    prevalence:
prev_x1 <- 0.3
#    coefficient:
beta_x1 <- 0.3
# 2. a continuous risk factor with distribution N(0,1)
#    coefficient:
beta_x2 <- 0.5

##### The simulation function

# The true model is:
# logit(p) = beta_0 + beta_x1 * x1 + beta_x2 * x2 + theta,
#            theta ~ N(0, tau2)
simulate_data <- function(
  N_hosp, Vol_avg, tau, volout,
  bet_0, prev_x1, beta_x1, beta_x2) {
  Vols <- rpois(n = N_hosp, lambda = Vol_avg) + 1
  
  pats <- bind_rows(
    lapply(
      seq_along(Vols),
      function(i) {
        v <- Vols[i]
        data.frame(
          provider = i,
          Volume = v,
          theta = rnorm(n = 1, sd = tau),
          x1 = sample(c(0, 1), size = v, replace = TRUE,
                      prob = c(1 - prev_x1, prev_x1)),
          x2 = rnorm(n = v)
        )
      }
    ))
  mutate(pats,
         pij = plogis(beta_0 + beta_x1 * x1 + beta_x2 * x2 + theta +
                        volout(Volume)),
         y = rbinom(n = nrow(pats), size = 1, prob = pij),
         provider = factor(provider))
}

##### Run the simulation
all_data <- bind_rows(lapply(
  names(volouts), function(volout_name) { # For each volume-outcome relationship
    volout <- volouts[[volout_name]]
    message("> volume outcome relationship ", volout_name)
    v_data <- bind_rows(lapply(
      N_hosps, function(N_hosp) {  # For each number of hospitals
      message(">> number of hospitals: ", N_hosp)
      on.exit(message("\r     \r", appendLF = FALSE))
      N_data <- bind_rows(lapply(
        1:n_it,
        function(i) {
          message("\r", i, "/", n_it, appendLF = FALSE)
          data <- simulate_data(
            N_hosp, Vol_avg, tau, volout,
            bet_0, prev_x1, beta_x1, beta_x2)
          ### To save data separately for each simulation run:
          # saveRDS(data,
          #         file.path(data_path, paste0(
          #           paste(volout_name, N_hosp, i, sep = "_"), ".RDS")))
          mutate(data, i = i)
        }))
      ### To save data separately for each value of N_hosp:
      # saveRDS(N_data,
      #         file.path(data_path, paste0(volout_name, "_", N_hosp, ".RDS")))
      mutate(N_data, N_hosp = N_hosp)
    }))
    ### Save data separately for each volume outcome effect:
    saveRDS(v_data, file.path(data_path, paste0(volout_name, ".RDS")))
    mutate(v_data, volout = volout_name)
  }))
saveRDS(all_data, file.path(data_path, "alldata.RDS"))

capture.output(sessioninfo::session_info(),
               file = file.path(data_path, "sessionInfo_generate"))

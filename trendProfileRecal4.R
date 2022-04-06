load("adj_epi_mif.rda")
load("model1params.rda")

local({r <- getOption("repos")
r["CRAN"] <- "http://repo.miserver.it.umich.edu/cran"
options(repos=r)})

install.packages("tidyverse")
install.packages('doParallel')
install.packages('parallel')
install.packages('tidyverse')
install.packages('devtools')
install.packages('subplex')
install.packages('pomp')
library(tidyverse)
library(devtools)
library(doParallel)
library(tidyverse)
library(subplex)
library(parallel)
library(pomp)

devtools::install_github('cbreto/panelPomp')

devtools::install_github('hetankevin/haitipkg')
library(haitipkg)
library(pomp)

# set run levels
run_level <- 3
Np <-          switch(run_level,100, 1e3, 5e3)
Nmif <-        switch(run_level, 10, 100, 200)
Nreps_eval <-  switch(run_level,  2,  10,  20)
Nreps_local <- switch(run_level, 10,  20,  40)
Nreps_global <-switch(run_level, 10,  20, 100)
Nsim <-        switch(run_level, 50, 100, 500) 
Nprof <-       switch(run_level,  1,  5,  10)


haiti1mod <- function(vacscen = 'id0') {
  ## get data
  dat <- haiti1_agg_data()
  fc_set <- vac_scen(vacscen)
  ## make covariate table
  covar <- covars(tmin = 0,
                  tmax = nrow(dat) + 573, ## for 11 year forecast
                  byt = 1,
                  degree = 6,
                  nbasis = 6,
                  per = 52.14,
                  data = dat,
                  settings = fc_set)
  if (vacscen == 'id0') {
    depts <- 0
  } else {
    depts <- fc_set$nd
  }
  ## make components pomp object building
  ## rinit
  state_names_base <- c("S", "E", "I", "A", "R")
  state_names <- state_names_base
  if (depts > 1) {
    for (i in 1:depts) {
      state_names <- c(state_names, paste0(state_names_base, i))
    }
  }
  ivp_names <- paste0(state_names, "_0")
  denom <- ivp_names[1]
  for (i in 2:length(ivp_names)) {
    denom <- c(denom, " + ", ivp_names[i])
  }
  denom <- paste(denom, collapse = "")
  frac_ivp <- paste0(" = nearbyint(frac * ", ivp_names, "); \n ")
  state_eqs <- c()
  for (i in 1:length(frac_ivp)) {
    state_eqs <- c(state_eqs, state_names[i], frac_ivp[i])
  }
  state_eqs <- paste(state_eqs, collapse = "")
  
  rinit_paste <- paste0("double frac = pop_0 / (", denom, "); \n ",
                        state_eqs,
                        "incid = 0.0; \n ",
                        "foival = 0.0; \n ",
                        "Str0 = 0.0; \n ",
                        "Sout = 0.0; \n ",
                        "Sin = 0.0; \n ")
  if (depts > 1) {
    rinit_paste <- c(rinit_paste, "incidU = 0.0; \n ", "incidV = 0.0; \n",
                     "asymV = 0.0; \n ", "newV = 0.0; \n ")
  }
  
  rinit <- Csnippet(
    paste(rinit_paste, collapse = "")
  )
  
  ## rprocess
  # transition rates and numbers
  rates_base <- c((2 + depts), (3 + depts), (2 + depts), (2 + depts), (2 + depts)) ## for SEIAR
  if (depts > 1) {
    rates_other <- rep(c(2, 3, 2, 2, 2), depts)
    rates <- c(rates_base, rates_other)
  } else {
    rates <- rates_base
  }
  trans_rates <- c()
  trans_numbers <- c()
  for (i in 1:length(rates)) {
    trans_rates <- c(trans_rates, paste0("double ", state_names[i], "rate[", rates[i], "]; \n "))
    trans_numbers <- c(trans_numbers, paste0("double ", state_names[i], "trans[", rates[i], "]; \n "))
  }
  trans_rates <- paste(trans_rates, collapse = "")
  trans_numbers <- paste(trans_numbers, collapse = "")
  trans_numbers <- paste0(trans_numbers, "\n double dgamma; \n ")
  
  # vaccination rates
  if (depts > 1) {
    vac_rates <- paste0("double eta", 1:depts, " = 0.0; \n ")
  }
  
  # demonitors and time checks
  if (depts > 1) {
    demons <- c("int pop_nv = S + E + I + A + R; \n ",
                paste0("int pop_", 1:depts, " = S", 1:depts, " + E",
                       1:depts, " + I", 1:depts, " + A", 1:depts,
                       " + R", 1:depts, "; \n "),
                "int pop = pop_nv", paste0(" + pop_", 1:depts),
                "; \n int births = rpois(mu*pop*dt); \n ") %>%
      paste(collapse = "")
    tchecks <- paste0("if (vac_tcheck == ", 1:depts, ") { \n eta", 1:depts,
                      " = num_vacc/pop_nv/dt; \n } \n ") %>%
      paste(collapse = "")
  } else {
    demons <- c("int pop = S + E + I + A + R; \n ",
                "int births = rpois(mu*pop*dt); \n ") %>%
      paste(collapse = "")
  }
  
  # seasonal beta and foi #1/10*(-betat)
  beta <- "double mybeta = exp(log(beta1*seas1 + beta2*seas2 + beta3*seas3 + beta4*seas4 + beta5*seas5 + beta6*seas6) - (betat)*((t-215)/(430-215))); \n "
  if (depts > 1) {
    foi_i <- c("I", paste0("+I", 1:depts)) %>%
      paste(collapse = "")
    foi_a <- c("A", paste0("+A", 1:depts)) %>%
      paste(collapse = "")
    foi <- paste0("double foi = pow(", foi_i, "+(1-kappa)*(", foi_a, "), nu)*mybeta/pop; \n ") %>%
      paste(collapse = "")
  } else {
    foi <- "double foi = pow(I, nu) * mybeta / pop; \n "
  }
  
  foi <- paste0(foi,
                "double sig_sq; \n ",
                "if (t < 233) { \n ",
                "sig_sq = sig_sq_epi; \n ",
                "} else { \n ",
                "sig_sq = sig_sq_end; \n ",
                "} \n ",
                "\n dgamma = rgammawn(sig_sq, dt); \n foi = foi * dgamma/dt; \n ")
  
  # theta_k
  if (vacscen != "id0") {
    thetas <- paste0("double theta", 1:depts, " = ve_d", 1:depts, "; \n ")
  }
  
  # compute rates
  s_rates <- c("Srate[0] = foi; \n ", "Srate[1] = delta; \n ")
  e_rates <- c("Erate[0] = sigma*(1-theta0); \n ", "Erate[1] = sigma*theta0; \n ", "Erate[2] = delta; \n ")
  i_rates <- c("Irate[0] = gamma; \n ", "Irate[1] = delta; \n ")
  a_rates <- c("Arate[0] = gamma; \n ", "Arate[1] = delta; \n ")
  r_rates <- c("Rrate[0] = alpha; \n ", "Rrate[1] = delta; \n ")
  if (depts > 1) {
    s_rates <- c(s_rates, paste0("Srate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("S", 1:depts, "rate[0] = foi; \n ",
                        "S", 1:depts, "rate[1] = delta; \n "))
    e_rates <- c(e_rates, paste0("Erate[", 3:(depts + 2), "] = eta", 1:depts, "; \n "),
                 paste0("E", 1:depts, "rate[0] = sigma*(1-theta", 1:depts, "); \n ",
                        "E", 1:depts, "rate[1] = sigma*theta", 1:depts, "; \n ",
                        "E", 1:depts, "rate[2] = delta; \n "))
    i_rates <- c(i_rates, paste0("Irate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("I", 1:depts, "rate[0] = gamma; \n ",
                        "I", 1:depts, "rate[1] = delta; \n "))
    a_rates <- c(a_rates, paste0("Arate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("A", 1:depts, "rate[0] = gamma; \n ",
                        "A", 1:depts, "rate[1] = delta; \n "))
    r_rates <- c(r_rates, paste0("Rrate[", 2:(depts + 1), "] = eta", 1:depts, "; \n "),
                 paste0("R", 1:depts, "rate[0] = alpha; \n ",
                        "R", 1:depts, "rate[1] = delta; \n "))
  }
  rates <- c(s_rates, e_rates, i_rates, a_rates, r_rates) %>%
    paste(collapse = "")
  
  # transition numbers
  numbers <- paste0("reulermultinom(", depts + 2, ",S,&Srate[0],dt,&Strans[0]); \n ",
                    "reulermultinom(", depts + 3, ",E,&Erate[0],dt,&Etrans[0]); \n ",
                    "reulermultinom(", depts + 2, ",I,&Irate[0],dt,&Itrans[0]); \n ",
                    "reulermultinom(", depts + 2, ",A,&Arate[0],dt,&Atrans[0]); \n ",
                    "reulermultinom(", depts + 2, ",R,&Rrate[0],dt,&Rtrans[0]); \n ")
  if (depts > 1) {
    numbers <- c(numbers,
                 paste0("reulermultinom(2,S", 1:depts, ",&S", 1:depts, "rate[0],dt,&S", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(3,E", 1:depts, ",&E", 1:depts, "rate[0],dt,&E", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(2,I", 1:depts, ",&I", 1:depts, "rate[0],dt,&I", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(2,A", 1:depts, ",&A", 1:depts, "rate[0],dt,&A", 1:depts, "trans[0]); \n "),
                 paste0("reulermultinom(2,R", 1:depts, ",&R", 1:depts, "rate[0],dt,&R", 1:depts, "trans[0]); \n ")) %>%
      paste(collapse = "")
  }
  
  # cohorts
  if (depts > 1) {
    coh_unvac <- c("S += -Strans[0]", paste0(" - Strans[", 1:(depts + 1), "]"), " + Rtrans[0] + births; \n ",
                   "E += -Etrans[0]", paste0(" - Etrans[", 1:(depts + 2), "]"), " + Strans[0]; \n ",
                   "I += -Itrans[0]", paste0(" - Itrans[", 1:(depts + 1), "]"), " + Etrans[0]; \n ",
                   "A += -Atrans[0]", paste0(" - Atrans[", 1:(depts + 1), "]"), " + Etrans[1]; \n ",
                   "R += -Rtrans[0]", paste0(" - Rtrans[", 1:(depts + 1), "]"), " + Itrans[0] + Atrans[0]; \n ")
    coh_vacs <- c(paste0("S", 1:depts, " += -S", 1:depts, "trans[0] - S", 1:depts,
                         "trans[1] + R", 1:depts, "trans[0] + Strans[", 2:(depts + 1), "]; \n "),
                  paste0("E", 1:depts, " += -E", 1:depts, "trans[0] - E", 1:depts, "trans[1] - E",
                         1:depts, "trans[2] + S", 1:depts, "trans[0] + Etrans[", 3:(depts + 2), "]; \n "),
                  paste0("I", 1:depts, " += -I", 1:depts, "trans[0] - I", 1:depts,
                         "trans[1] + E", 1:depts, "trans[0] + Itrans[", 2:(depts + 1), "]; \n "),
                  paste0("A", 1:depts, " += -A", 1:depts, "trans[0] - A", 1:depts,
                         "trans[1] + E", 1:depts, "trans[0] + Atrans[", 2:(depts + 1), "]; \n "),
                  paste0("R", 1:depts, " += -R", 1:depts, "trans[0] - R", 1:depts,
                         "trans[1] + A", 1:depts, "trans[0] + I", 1:depts, "trans[0] + Rtrans[", 2:(depts + 1), "]; \n "))
    coh <- c(coh_unvac, coh_vacs) %>%
      paste(collapse = "")
  } else {
    coh <- c("S += -Strans[0]", paste0(" - Strans[", 1:(depts + 1), "]"), " + Rtrans[0] + births; \n ",
             "E += -Etrans[0]", paste0(" - Etrans[", 1:(depts + 2), "]"), " + Strans[0]; \n ",
             "I += -Itrans[0]", paste0(" - Itrans[", 1:(depts + 1), "]"), " + Etrans[0]; \n ",
             "A += -Atrans[0]", paste0(" - Atrans[", 1:(depts + 1), "]"), " + Etrans[1]; \n ",
             "R += -Rtrans[0]", paste0(" - Rtrans[", 1:(depts + 1), "]"), " + Itrans[0] + Atrans[0]; \n ") %>%
      paste(collapse = "")
  }
  
  # incidence
  incids <- c("incid += Etrans[0]; \n ")
  foi_val <- "foival += foi; \n "
  str0 <- "Str0 += Strans[0]; \n "
  sin <- c("Sin += Rtrans[0] + births; \n ")
  sout <- c("Sout += Strans[0] + Strans[1]; \n ")
  last <- c(foi_val, str0, sin)
  
  if (depts > 1) {
    incids <- c("incid += Etrans[0]", paste0(" + E", 1:depts, "trans[0]"), "; \n ")
    incid_u <- "incidU += Etrans[0]; \n "
    incid_v <- c("incidV += E1trans[0]", paste0(" + E", 2:depts, "trans[0]"), "; \n ")
    asym_v <- c("asymV += E1trans[1]", paste0(" + E", 2:depts, "trans[1]"), "; \n ")
    new_v <- c("newV += Strans[2]", paste0(" + Strans[", 3:(depts + 1), "]"),
               paste0(" + Etrans[", 3:(depts + 2), "]"),
               paste0(" + Itrans[", 2:(depts + 1), "]"),
               paste0(" + Atrans[", 2:(depts + 1), "]"),
               paste0(" + Rtrans[", 2:(depts + 1), "]"), "; \n ")
    sout <- c("Sout += Strans[0]", paste0(" + Strans[", 1:(depts + 1), "]"), "; \n ")
    last <- c(last, incid_u, incid_v, asym_v, new_v)
  }
  last <- c(last, incids, sout) %>%
    paste(collapse = "")
  
  if (depts > 1) {
    rproc_paste <- c(trans_rates, trans_numbers, vac_rates, demons, tchecks, beta,
                     foi, thetas, rates, numbers, coh, last) %>%
      paste(collapse = "")
  } else {
    rproc_paste <- c(trans_rates, trans_numbers, demons, beta, foi, rates, numbers,
                     coh, last) %>%
      paste(collapse = "")
  }
  
  rproc <- Csnippet(
    rproc_paste
  )
  
  ## dmeasure
  dmeas <- Csnippet("
    if (ISNA(cases)) {
      lik = (give_log) ? 0 : 1;
    } else {
      double rho = rho_epi;
      double tau = tau_epi;
      if (t > 232) {
        rho = rho_end;
        tau = tau_end;
      }
      lik = dnbinom_mu(cases, tau, rho*incid, give_log);
    }
  ")
  
  ## rmeasure
  rmeas <- Csnippet("
    double rho = rho_epi;
    double tau = tau_epi;
    if (t > 232) {
      rho = rho_end;
      tau = tau_end;
    }
    cases = rnbinom_mu(tau, rho*incid);
    if (cases > 0.0) {
      cases = nearbyint(cases);
    } else {
      cases = 0.0;
    }
  ")
  
  
  ## names
  if (depts > 1) {
    ## state names
    state_names <- c("S", paste0("S", 1:depts),
                     "E", paste0("E", 1:depts),
                     "I", paste0("I", 1:depts),
                     "A", paste0("A", 1:depts),
                     "R", paste0("R", 1:depts),
                     "incid", "incidU", "incidV", "asymV", "newV",
                     "foival", "Str0", "Sout", "Sin")
    
    ## parameter names
    param_names <- c("rho_epi", "sig_sq_epi", "tau_epi", # epidemic
                     "rho_end", "sig_sq_end", "tau_end", # endemic
                     "beta1", "beta2", "beta3", "beta4", "beta5",
                     "beta6", "betat", "nu", "gamma", "sigma", "theta0", "alpha", "mu", "delta",
                     "S_0","E_0","I_0","A_0","R_0", "pop_0", "kappa",
                     paste0("S", 1:depts, "_0"),
                     paste0("E", 1:depts, "_0"),
                     paste0("I", 1:depts, "_0"),
                     paste0("A", 1:depts, "_0"),
                     paste0("R", 1:depts, "_0"))
    
    
    
    ## accum vars
    accum_names <- c("incid","incidU","incidV","asymV","newV",
                     "foival","Str0","Sout","Sin")
  } else {
    ## state names
    state_names <- c("S", "E", "I", "A", "R", "incid", "foival", "Str0", "Sout", "Sin")
    
    ## parameter names
    param_names <- c("rho_epi", "sig_sq_epi", "tau_epi", # epidemic
                     "rho_end", "sig_sq_end", "tau_end", # endemic
                     "beta1", "beta2", "beta3", "beta4", "beta5",
                     "beta6", "betat", "gamma", "sigma", "theta0", "alpha", "mu", "delta",
                     "nu", "pop_0", 
                     "S_0","E_0","I_0","A_0","R_0")
    
    ## accum vars
    accum_names <- c("incid", "foival","Str0","Sout","Sin")
  }
  
  ## partrans
  param_trans <- pomp::parameter_trans(
    log = c("sigma", "gamma", "mu", "delta", "alpha",
            "sig_sq_end", "sig_sq_epi", "tau_end", "tau_epi"),
    logit = c("rho_epi", "rho_end", "nu", "theta0"), 
    barycentric = c("S_0", "E_0", "I_0", "A_0", "R_0")
  )
  
  ## build pomp model
  model1 <- pomp::pomp(
    data = dat,
    times = "week",
    t0 = 0,
    rmeasure = rmeas,
    dmeasure = dmeas,
    rprocess = pomp::euler(step.fun = rproc, delta.t = 1/7),
    covar = pomp::covariate_table(covar, times = "time"),
    partrans = param_trans,
    statenames = state_names,
    paramnames = param_names,
    accumvars = accum_names,
    rinit = rinit
  )
  
  
  
  model1 <- pomp::window(model1, start = 1, end = 430) ## full period
  model1@t0 <- 0
  
  return(model1)
}

mod_1 <- haiti1mod()


prof_if[,-c(1,2, 19,27, 28)] %>%
  filter(loglik>max(loglik)-20,loglik.se<2) %>%
  sapply(range) -> box


unfixed_param_names1 <- c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
                          "nu", "rho", "tau", "E_0", "I_0", "sig_sq")

fixed_params <- unlist(MODEL1_INPUT_PARAMETERS$adj_pars_epi)[c("pop_0","mu", "delta", "sigma", "gamma", "alpha", "theta0", "A_0", "R_0")]

betatSamps <- freeze(seed=876,
                     profile_design(
                       betat=seq(0.151,0.2,length=50),
                       lower=box[1,unfixed_param_names1],
                       upper=box[2,unfixed_param_names1],
                       nprof=Nprof, type="runif"))  



epi_pars <- MODEL1_INPUT_PARAMETERS$adj_pars_epi

starts <- betatSamps %>%
  mutate(parid = seq_along(betat)) %>%
  mutate(theta0 = fixed_params["theta0"], mu = fixed_params["mu"], delta = fixed_params["delta"], 
         sigma = fixed_params["sigma"], alpha = fixed_params["alpha"], 
         gamma = fixed_params["gamma"], A_0 =  fixed_params["A_0"], R_0 = fixed_params["R_0"], 
         S_0 = 1- E_0-I_0-A_0-R_0,  pop_0 = fixed_params["pop_0"]) %>% 
  select(parid, rho, tau, beta1, beta2, beta3, beta4, beta5,
         beta6, betat, nu, gamma, sigma, theta0, alpha, mu, delta, sig_sq,
         S_0, E_0, I_0, A_0, R_0, pop_0)


add_starts <- starts %>%
  dplyr::select(rho, tau, sig_sq) %>%
  dplyr::rename(rho_end = rho, sig_sq_end = sig_sq, tau_end = tau)

starts <- cbind(starts, add_starts) %>%
  dplyr::rename(rho_epi = rho, tau_epi = tau, sig_sq_epi = sig_sq)


rw_sds <- rw.sd(beta1 = 0.02, beta2 = 0.02, beta3 = 0.02,
                beta4 = 0.02, beta5 = 0.02, beta6 = 0.02,
                rho_epi = ifelse(time > 232, 0.0, 0.02),
                rho_end = ifelse(time <= 232, 0.0, 0.02),
                tau_epi = ifelse(time > 232, 0.0, 0.02),
                tau_end = ifelse(time <= 232, 0.0, 0.02),
                sig_sq_epi = ifelse(time > 232, 0.0, 0.02),
                sig_sq_end = ifelse(time <= 232, 0.0, 0.02),
                nu = 0.02,
                E_0 = ifelse(time > 232, 0.0, ivp(0.2)),
                I_0 = ifelse(time > 232, 0.0, ivp(0.2)))


library(doParallel)
library(doRNG)
library(parallel)
no_cores <- detectCores(logical=FALSE) 
cl<- parallel::makeCluster(no_cores, outfile = "")

# m2 is from the betatrendexplorations file
# needed the topline because of the str of
# mif2 requiring a special kind of object
# cl defined in betatrend exploration.Rmd
doParallel::registerDoParallel(cl)
registerDoRNG(12)

starts_if <- starts
print('starting fitting')

foreach(start = iter(starts_if, by = 'row'), 
        .combine=rbind, 
        .inorder = FALSE,
        .packages = c("pomp", "tidyverse", "magrittr")) %dopar% {
          
          joint_mod <- mod_1
          timezero(joint_mod) <- 0
          
          pars <- unlist(start[,-1])
          
          names(pars)[c(11:16, 18, 21:23)] <- c("gamma", "sigma", "theta0", "alpha",
                                                "mu", "delta", "S_0", "A_0", "R_0", "pop_0")
          parOrder <- c("rho_epi", "sig_sq_epi", "tau_epi", # epidemic
                        "rho_end", "sig_sq_end", "tau_end", # endemic
                        "beta1", "beta2", "beta3", "beta4", "beta5",
                        "beta6", "betat", "gamma", "sigma", "theta0", "alpha", "mu", "delta",
                        "nu", "pop_0", 
                        "S_0","E_0","I_0","A_0","R_0")
          
          pars <- pars[order(factor(names(pars), levels = parOrder))]
          coef(joint_mod) <- pars
          mf.mod <- mif2(joint_mod,
                         Nmif = Nmif,
                         rw.sd = rw_sds,
                         Np = Np,
                         cooling.fraction.50 = 0.5,
                         verbose = FALSE)
          
          full.lik <- replicate(Nreps_eval, pfilter(mf.mod, Np = Np))
          full_ll <- sapply(full.lik, logLik)
          full_ll <- logmeanexp(full_ll, se = TRUE)
          
          dummy <- data.frame(model = paste0("if"),
                              parid = start$parid,
                              as.list(coef(mf.mod)),
                              full_loglik = full_ll[1],
                              full_loglik.se = full_ll[2])
          rm(mf.mod, full.lik)
          gc()
          dummy
        }-> resultsRecal4

save(resultsRecal4, file = "trendProfileRecalResults4.rda")

load("adj_epi_mif.rda")

local({r <- getOption("repos")
r["CRAN"] <- "http://repo.miserver.it.umich.edu/cran"
options(repos=r)})

install.packages('doParallel')
install.packages('parallel')
install.packages('tidyverse')
install.packages('devtools')
install.packages('subplex')
install.packages('pomp')
library(devtools)
devtools::install_github('cbreto/panelPomp')
devtools::install_github('hetankevin/haitipkg')
library(haitipkg)
mod_1 = haiti1mod()

prof_if[,-c(1,2, 19,27, 28)] %>%
  filter(loglik>max(loglik)-20,loglik.se<2) %>%
  sapply(range) -> box

# S_0 ??
unfixed_param_names <- c("beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
                         "nu", "rho", "tau", "E_0", "I_0", "sig_sq")

fixed_params <- coef(mod_1)[c("pop_0","mu", "delta", "sigma", "gamma", "alpha", "theta0", "A_0", "R_0")]

rw_sd1 <- rw.sd(beta1 = .02, beta2 = .02, beta3 = .02,
                beta4 = .02, beta5 = .02, beta6 = .02,
                tau = 0.02, rho = 0.02, nu = 0.02, sig_sq = 0.02,
                E_0 = ivp(0.2), I_0 = ivp(0.2))

freeze(seed=1196696958,
       profile_design(
         betat=seq(-0.07,0.04,length=100),
         lower=box[1,unfixed_param_names],
         upper=box[2,unfixed_param_names],
         nprof=10, type="runif"
       )) -> guesses

library(doParallel)
library(doRNG)
library(parallel)
no_cores <- detectCores(logical=FALSE) 
cl<- parallel::makeCluster(no_cores)

# m2 is from the betatrendexplorations file
# needed the topline because of the str of
# mif2 requiring a special kind of object
# cl defined in betatrend exploration.Rmd
doParallel::registerDoParallel(cl)
registerDoRNG(123)
#foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
foreach(i = 1:nrow(guesses), .combine=rbind) %dopar% {
  pars = c(unlist(guesses[i,]),fixed_params, S_0 = 1 - guesses[i,"E_0"] - guesses[i,"I_0"])
  parOrder = c("rho","tau","beta1","beta2","beta3","beta4","beta5","beta6","nu",
               "gamma","sigma","theta0","alpha","mu","delta","sig_sq","S_0","E_0",
               "I_0","A_0","R_0","pop_0","betat")
  pars = pars[order(factor(names(pars), levels = parOrder))]
  mod_1@params = pars
  library(pomp)
  library(tidyverse)
  #mod_1@params = c(unlist(guesses[i,]),fixed_params, S_0 = 1 - guesses[i,"E_0"] - guesses[i,"I_0"])
  mod_1 %>%
    # problem is with params argument
    mif2(
      rw.sd=rw_sd1, Nmif=200,cooling.fraction.50=0.5, Np=5000) -> mf
  # replicate 10 runs of pfilter so we can get se's
  replicate(
    100,
    mf %>% pfilter(Np=5000) %>% logLik()) %>%
    logmeanexp(se=TRUE) -> ll
  mf %>% coef() %>% bind_rows() %>%
    bind_cols(loglik=ll[1],loglik.se=ll[2])
}-> results

save(results, file = "trendResults.rda")
#############################################################################################################
############# Simulation study for 'A Note on Parallel Sampling in Markov Graphs ############################
#############################################################################################################

# clear workspace
rm(list=ls())

### maybe install package before
# library(devtools)
# install_github('VerenaMaier/pergm')

### load libraries
library(microbenchmark)
library(devtools)
library(pergm)
library(ergm)

sessionInfo()

?simulate_networks


### Timing of parallel algorithm
# N = 200, N = 2200, N = 4200
# n_cores = 1, 2, 4, 6, 8, 10, 12, 14
# m = 1, 2, 3
# n_updates for m = 1: N*(N-1)/2*10
# n_updates for m > 1: N*(N-1)/2
# n_rep = 50

theta <- c(-0.1, -0.01, 0.02) # approx. density 0.06

n_cores_seq <- c(1, seq(2, 14, by = 2))
dim_seq <- seq(200, 4200, by = 2000)
nsim <- 3
n_func <- 1
n_rep <- 50

set.seed(190986)
timing_sim_pergm <- data.frame(sim = rep(1:nsim, each = n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  dim = rep(rep(dim_seq, each = n_rep*n_func*length(n_cores_seq)), times = nsim),
  n_cores = rep(rep(n_cores_seq, each = n_rep*n_func), times = nsim*length(dim_seq)),
  func = rep(0, each = nsim*n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  time = rep(0, each = nsim*n_rep*n_func*length(n_cores_seq)*length(dim_seq)))

i=1
for(nsim in 1:nsim){
  print(paste("nsim = ", nsim))
  for(dim in dim_seq){
    print(paste("dim = ", dim))
    for (n_cores in n_cores_seq){
      print(paste("n_cores = ", n_cores))
      
      n_update <- dim*(dim-1)/2
      burnin <- n_update*9
      
      microbench <-
        microbenchmark(
          pergm_parallel = simulate_networks(dim = dim, theta = theta, n_update = n_update,  burnin = burnin, 
            n_cores = n_cores, n_sim = nsim, log_change = FALSE, return_nw = FALSE),
          unit = 'ns', times = n_rep
        )
      
      timing_sim_pergm$time[i:(i+(n_rep*n_func-1))] <- microbench$time
      timing_sim_pergm$func[i:(i+(n_rep*n_func-1))] <- microbench$expr
      
      i <- i + n_rep*n_func
      
      saveRDS(timing_sim_pergm, file = 'timing_sim_pergm.RData')
    }
  }
}





rm(list=ls())

sessionInfo()
rm(list=ls())

####################################################################################################################
####################################################################################################################
####################################################################################################################
### Comparison of pergm, serial and ergm (serial)
# N = 200, N = 2200, N = 4200
# n_cores = 14 (if used)
# m = 1, 2, 3
# n_updates for m = 1: N*N*10
# n_updates for m > 1: N*(N-1)/2
# n_rep = 50

theta <- c(-0.1, -0.01, 0.02) # approx. density 0.06

n_cores_seq <- n_cores<- 14
dim_seq <- seq(200, 4200, by = 2000)
nsim <- 3
n_func <- 3
n_rep <- 50
set.seed(190986)

timing_sim_comp_serial <- data.frame(sim = rep(1:nsim, each = n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  dim = rep(rep(dim_seq, each = n_rep*n_func*length(n_cores_seq)), times = nsim),
  n_cores = rep(rep(n_cores_seq, each = n_rep*n_func), times = nsim*length(dim_seq)),
  func = rep(0, each = nsim*n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  time = rep(0, each = nsim*n_rep*n_func*length(n_cores_seq)*length(dim_seq)))

i=1
for(nsim in 1:nsim){
  print(paste("nsim = ", nsim))
  for(dim in dim_seq){
    print(paste("dim = ", dim))
    adj_mat <- matrix(0, nrow = dim, ncol = dim)
    n_update <- dim*(dim-1)/2
    burnin <- n_update*9
    
    microbench <- microbenchmark(
      pergm_parallel = simulate_networks(dim = dim, theta = theta, n_update = n_update,  burnin = burnin, 
        n_cores = n_cores, n_sim = nsim, log_change = FALSE, return_nw = FALSE),
      pergm_serial = simulate_networks(dim = dim, theta = theta, n_update = n_update,  burnin = burnin, 
        n_cores = 1, n_sim = nsim, log_change = FALSE, return_nw = FALSE),
      ergm_serial = simulate(network(adj_mat, directed=FALSE) ~ edges + kstar(2) + triangle,
        coef=c(theta), statsonly=TRUE, nsim = nsim, seed = 190986,
        control=control.simulate(MCMC.burnin = (burnin + n_update), MCMC.interval = n_update,
          MCMC.prop.weights = 'random')),
      unit = 'ns', times = n_rep
    )
    
    timing_sim_comp_serial$time[i:(i+(n_rep*n_func-1))] <- microbench$time
    timing_sim_comp_serial$func[i:(i+(n_rep*n_func-1))] <- microbench$expr
    
    i <- i + n_rep*n_func
    
    saveRDS(timing_sim_comp_serial, file = 'timing_sim_comp_serial.RData')
  }
}

sessionInfo()
rm(list=ls())

###################################################################################################################
###################################################################################################################
###################################################################################################################
## Comparison of pergm with ergm (PSOCK, MPI)
N = 2200
n_cores = 14
m = 100
n_updates for m = 1: N*N*10
n_updates for m > 1: N*(N-1)/2
n_rep = 50


theta <- c(-0.1, -0.01, 0.02) # approx. density 0.06

n_cores <- n_cores_seq <- 14
dim <- dim_seq <- 2200
nsim <- nsim_seq <- 100
n_func <- 5
n_rep <- 50


set.seed(190986)
timing_sim_comp_parallel <- data.frame(sim = rep(nsim_seq, each = n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  dim = rep(rep(dim_seq, each = n_rep*n_func*length(n_cores_seq)), times = length(nsim_seq)),
  n_cores = rep(rep(n_cores_seq, each = n_rep*n_func), times = length(nsim_seq)*length(dim_seq)),
  func = rep(0, each = length(nsim_seq)*n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  time = rep(0, each = length(nsim_seq)*n_rep*n_func*length(n_cores_seq)*length(dim_seq)))


adj_mat <- matrix(0, nrow = dim, ncol = dim)
n_update <- dim*(dim-1)/2
burnin <- n_update*9


microbench <- microbenchmark(
  ergm_paralell_PSOCK = simulate(network(adj_mat, directed=FALSE) ~ edges + kstar(2) + triangle,
    coef = c(theta), statsonly = TRUE,
    nsim = nsim, seed = 190986,
    control = control.simulate(MCMC.burnin = (burnin+n_update),
      MCMC.interval = n_update, MCMC.prop.weights = 'random',
      parallel=n_cores, parallel.type="PSOCK")),
  ergm_paralell_MPI = simulate(network(adj_mat, directed=FALSE) ~ edges + kstar(2) + triangle,
    coef=c(theta), statsonly=TRUE,
    nsim = nsim, seed = 190986,
    control = control.simulate(MCMC.burnin = (burnin+n_update), MCMC.interval = n_update,
      MCMC.prop.weights = 'random',
      parallel = n_cores, parallel.type="MPI")),
  pergm_parallel = simulate_networks(dim = dim, theta = theta, n_update = n_update,  burnin = burnin, 
    n_cores = n_cores, n_sim = nsim, log_change = FALSE),
  pergm_serial = simulate_networks(dim = dim, theta = theta, n_update = n_update,  burnin = burnin, 
    n_cores = 1, n_sim = nsim, log_change = FALSE, return_nw = FALSE),
  ergm_serial = simulate(network(adj_mat, directed=FALSE) ~ edges + kstar(2) + triangle,
    coef=c(theta), statsonly=TRUE, nsim = nsim, seed = 190986,
    control=control.simulate(MCMC.burnin = (burnin + n_update), MCMC.interval = n_update,
      MCMC.prop.weights = 'random')),
  unit = 'ns', times = n_rep)

timing_sim_comp_parallel$time <- microbench$time
timing_sim_comp_parallel$func <- microbench$expr


saveRDS(timing_sim_comp_parallel, file = 'timing_sim_comp_parallel.RData')

sessionInfo()


rm(list=ls())


###################################################################################################################
###################################################################################################################
###################################################################################################################
## Comparison of pergm with ergm (PSOCK, MPI)
N = 2200
n_cores = 14
m = 100
n_updates for m = 1: N*N*10
n_updates for m > 1: N*(N-1)/2
n_rep = 50


theta <- c(-4, -0.01, 0.2) # approx. density 0.01

n_cores <- n_cores_seq <- 14
dim <- dim_seq <- 2200
nsim <- nsim_seq <- 100
n_func <- 5
n_rep <- 50


set.seed(190986)
timing_sim_comp_parallel <- data.frame(sim = rep(nsim_seq, each = n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  dim = rep(rep(dim_seq, each = n_rep*n_func*length(n_cores_seq)), times = length(nsim_seq)),
  n_cores = rep(rep(n_cores_seq, each = n_rep*n_func), times = length(nsim_seq)*length(dim_seq)),
  func = rep(0, each = length(nsim_seq)*n_rep*n_func*length(n_cores_seq)*length(dim_seq)),
  time = rep(0, each = length(nsim_seq)*n_rep*n_func*length(n_cores_seq)*length(dim_seq)))


adj_mat <- matrix(0, nrow = dim, ncol = dim)
n_update <- dim*(dim-1)/2
burnin <- n_update*9


microbench <- microbenchmark(
  ergm_paralell_PSOCK = simulate(network(adj_mat, directed=FALSE) ~ edges + kstar(2) + triangle,
    coef = c(theta), statsonly = TRUE,
    nsim = nsim, seed = 190986,
    control = control.simulate(MCMC.burnin = (burnin+n_update),
      MCMC.interval = n_update, MCMC.prop.weights = 'random',
      parallel=n_cores, parallel.type="PSOCK")),
  ergm_paralell_MPI = simulate(network(adj_mat, directed=FALSE) ~ edges + kstar(2) + triangle,
    coef=c(theta), statsonly=TRUE,
    nsim = nsim, seed = 190986,
    control = control.simulate(MCMC.burnin = (burnin+n_update), MCMC.interval = n_update,
      MCMC.prop.weights = 'random',
      parallel = n_cores, parallel.type="MPI")),
  pergm_parallel = simulate_networks(dim = dim, theta = theta, n_update = n_update,  burnin = burnin, 
    n_cores = n_cores, n_sim = nsim, log_change = FALSE),
  pergm_serial = simulate_networks(dim = dim, theta = theta, n_update = n_update,  burnin = burnin, 
    n_cores = 1, n_sim = nsim, log_change = FALSE, return_nw = FALSE),
  ergm_serial = simulate(network(adj_mat, directed=FALSE) ~ edges + kstar(2) + triangle,
    coef=c(theta), statsonly=TRUE, nsim = nsim, seed = 190986,
    control=control.simulate(MCMC.burnin = (burnin + n_update), MCMC.interval = n_update,
      MCMC.prop.weights = 'random')),
  unit = 'ns', times = n_rep)

timing_sim_comp_parallel$time <- microbench$time
timing_sim_comp_parallel$func <- microbench$expr


saveRDS(timing_sim_comp_parallel, file = 'timing_sim_comp_parallel_setting2.RData')


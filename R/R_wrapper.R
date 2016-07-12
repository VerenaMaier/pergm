#' @title Simulation of networks in parallel based on log/non-log ERGMs
#'
#' @description Simulates networks based on log/non-log exponential random graph models (ERGMs).
#' Simulation can be done in parallel. Used network statistics (number of edges, 2-stars and triangles)
#' induce a conditional independence structure amongst the edges of the network.
#' @param adjacency A start matrix for the simulation (default: empty matrix)
#' @param dim The dimension of the adjacency matrix as a number (N); either \code{dim} or \code{adjacency} has to be defined
#' @param theta A vector with the estimated coefficients of the edge, 2-star and triangle statistic (given in this order)
#' @param n_update The number of updates for the Markov chain
#' @param burnin The number of burnin which is added to \code{n_update} for the first simulated network
#' @param n_cores The number of cores used for parallelization
#' @param n_sim The number of simulated networks
#' @param log_change True (default) if the log changes of the statistics should be used
#' @param return_nw False (default) if only statistics instead of complete simulated networks should be returned
#' @return A matrix with the number of edges, 2-stars and triangles of each simulated network
#' @importFrom Rcpp, evalCpp
#' @details
#' For simulating networks in parallel we start a single Markov chain simulation but apply parallel computing for each single step of the Markov chain. The central idea is to take advantage of the conditional independence structure and simulate multiple conditionally independent edges in networks simultaneously in parallel.
#' The returned statistics consist of the number of edges, 2-stars and triangles of each simulated network. If \code{log_change} is true, the log transformed statistics are returned. The log transformation provides a concave function such that linearity and degeneracy problems diminishes.
#' The algorithm is implemented in C++ and the parallelization works in C++ via OpenMP (Open Multi-Processing).
#' @references
#' Maier, V., Fürlinger, K., Kauermann G. (2016). A Note on Parallel Sampling in Exponential Random Graph Models (ERGM) \cr (to appear)
#' @examples
#' # number of nodes
#' N <- 10
#' # the given parameter vector of edges, 2-star and triangle
#' theta <- c(-2, 0.01, -0.2)
#' # simulate one network with an empty adjacency matrix as start
#' network_stats <- simulate_networks(dim = N, theta = theta, log_change = TRUE)
#'
simulate_networks <- function(adjacency = NULL,
                              dim = NULL,
                              theta = c(0, 0, 0),
                              n_update = 1000000,
                              burnin = 1000,
                              n_cores = 4,
                              n_sim = 1,
                              log_change = TRUE,
                              return_nw = FALSE){


  # default values
  if(is.null(adjacency) && is.null(dim)) stop("either adjacency or dim has to be specified")
  if(is.null(dim)) dim <- nrow(adjacency)
  if(is.null(adjacency)) adjacency <- matrix(0, nrow = dim, ncol = dim)

  if(length(theta) == 1) theta <- c(theta, 0, 0)
  if(length(theta) == 2) theta <- c(theta, 0)
  if(length(theta) > 3 ) stop("theta represents the parameter of edge, 2-star and triangle (#<4)")

  # check data types
  if(!is.matrix(adjacency)) stop("adjacency has to be a matrix")
  if(!is.numeric(dim)) stop("dim has to be numeric")
  if(any(is.na(theta))) stop("theta cannot be NA")
  if(!is.numeric(n_update)) stop("n_update has to be numeric")
  if(!is.numeric(n_cores)) stop("n_cores has to be numeric")
  if(!is.numeric(n_sim)) stop("nsim has to be numeric")
  if(!is.numeric(burnin)) stop("burnin has to be numeric")
  if(!is.logical(log_change)) stop("log_change has to be boolean")
  if(!is.logical(return_nw)) stop("return_nw has to be boolean")

  ### call function
  if(return_nw == FALSE){
    simulate_networks_fit(adjacency = adjacency, dim = dim, theta = theta, n_update = n_update,
      n_cores = n_cores, nsim = n_sim, burnin = burnin, log_change = log_change)

  }else{
    if(n_sim > 6) stop("not more than 5 network adjacencies can be returned, set return_nw to FALSE")
    simulate_networks_fit_nw(adjacency = adjacency, dim = dim, theta = theta, n_update = n_update,
      n_cores = n_cores, nsim = n_sim, burnin = burnin, log_change = log_change)
  }

}




#' @title Calculation of the log transformed statistics of edges, 2-star and triangle for an observed adjacency matrix
#'
#' @description Calculates the log transformed statistics of edges, 2-star and triangle for an observed adjacency matrix. Calculation is combined in one loop in C++ and can be done in parallel.
#'
#' @param adjacency The adjacency matrix of a network (symmetric with only zeros or ones)
#' @param n_cores The number of cores used for parallelization
#' @return A vector of 3 with the number of edges, log transformed number of 2-stars and log transformed number of triangles
#'
#' @examples
#' # generate symmetric matrix with diag = 0
#' N <- 10
#' set.seed(191919)
#' adjacency <- matrix(rbinom(n = N * N, size = 1, prob = 0.5), N, N)
#' diag(adjacency) <- 0
#' adjacency[lower.tri(adjacency)] <- t(adjacency)[lower.tri(adjacency)]
#'
#' summary_stats_log(adjacency, n_cores = 4)
#'

summary_stats_log <- function(adjacency = NULL, n_cores = 4){

  if(is.null(adjacency)) stop("adjacency matrix has to be specified")
  if(isSymmetric(adjacency) == 0) stop("adjacency matrix has to be symmetric")
  if(!is.numeric(n_cores)) stop("n_cores has to be numeric")
  dim = nrow(adjacency)

  summary_stats_log_c(mat = adjacency , dim = dim, n_cores = n_cores)

}

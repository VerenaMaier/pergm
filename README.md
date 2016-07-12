# pergm
The R package allows to simulate Exponential Random Graph Models (ERGM)
in parallel. For simulating networks in parallel we start a single Markov 
chain simulation but apply parallel computing for each single step of the 
Markov chain. The central idea is to take advantage of the conditional 
independence structure and simulate multiple conditionally independent 
edges in networks simultaneously in parallel.
A log transformation can be used as well. The log transformation provides 
a concave function such that linearity and degeneracy problems diminishes.
The algorithm is implemented in C++ and the parallelization works in C++ 
via OpenMP (Open Multi-Processing).

For details on the methodology see:
Maier, V., Fürlinger, K., Kauermann G. (2016). A Note on Parallel Sampling 
in Exponential Random Graph Models (ERGM) (to appear)
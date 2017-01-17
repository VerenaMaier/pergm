#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <iostream>

// [[Rcpp::plugins(openmp)]]

#include "types_20150828.h"
#include "timing.h"

//' @title Return log-statistics of edges, 2-star, triangle as summary statistics
//'
//' @description Return log-statistics of edges, 2-star, triangle as summary statistics.
//' Calculation is combined in one loop in C++ and can be done in parallel
//'
//' @param mat The adjacency as IntegerMatrix
//' @param dim The dimension of the adjacency as a number
//' @param n_cores The number of used cores for parallelization as number
//' @return A vector of 3 with the number of edges, log number of 2-stars and the log number of triangles
//'
//' @examples
//' # generate symmetric matrix with diag = 0
//' N <- 10
//' set.seed(191919)
//' adjacency <- matrix(rbinom(n = N * N, size = 1, prob = 0.02), N, N)
//' diag(adjacency) <- 0
//' adjacency[lower.tri(adjacency)] <- t(adjacency)[lower.tri(adjacency)]
//'
//' summary_stats_log(adjacency, dim = N, n_cores = 4)
//'
// [[Rcpp::export]]
Rcpp::NumericVector summary_stats_log(Rcpp::IntegerMatrix mat,
                                          int dim,
                                          int n_cores){
  Rcpp::NumericVector stats_log(3);
  int edge = 0;
  double star2 = 0;
  double triangle = 0;
  double sum_star = 0;
  double sum_tria = 0;
  int k, l, p;

  // initialize number of different threads
  omp_set_num_threads(n_cores);

  for(k=0; k<dim; k++)
  {
    #pragma omp parallel for reduction(+:star2, sum_star, triangle, sum_tria, edge)
    for(l=0; l<dim; l++)
    {
      if(l > k)
      {
        edge = edge + mat(k,l);

        for(p=0; p<dim; p++)
        {
          if(p!=k & p!=l)
          {
            sum_star = sum_star + mat(k,p) + mat(l,p);
            sum_tria = sum_tria + mat(k,p) * mat(l,p);
          }
        }
        star2 = star2 + mat(k,l)*log(sum_star + 0.1); // log(change + 1)
        sum_star = 0;
        triangle = triangle + mat(k,l)*log(sum_tria + 0.1); // log(change + 1)
        sum_tria = 0;
      }
    }
  }
  star2 = star2/2;
  triangle = triangle/3;
  stats_log(0) = edge;
  stats_log(1) = star2;
  stats_log(2) = triangle;

  return stats_log;
}

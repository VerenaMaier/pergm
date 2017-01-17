#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <iostream>


// [[Rcpp::plugins(openmp)]]

#include "types_20150828.h"
#include "timing.h"

// main simulate function parallel
template<typename T>
void sim_nw_c_par(Matrix<T>& mat2, size_t dim, size_t n_update, size_t n_cores, int k, bool log_change);
// serial
template<typename T>
void sim_nw_c_ser(Matrix<T>& mat2, size_t dim, size_t n_update, int k, bool log_change);

// simulate thread function
template<typename T>
void sim_nw_thread(Matrix<T>& mat, int dim, int i, int j, unsigned *seed, bool log_change);

// helper functions parallel
template<typename T>
int edge_number(Matrix<T>& mat, int dim, int n_cores);
template<typename T>
int star2_number(Matrix<T>& mat, int dim, int n_cores);
template<typename T>
int triangle_number3(Matrix<T>& mat, int dim, int n_cores);

// helper functions serial
template<typename T>
int edge_number_ser(Matrix<T>& mat, int dim);
template<typename T>
int star2_number_ser(Matrix<T>& mat, int dim);
template<typename T>
int triangle_number_ser(Matrix<T>& mat, int dim);

// helper function for log
template<typename T>
Rcpp::NumericVector stats_log_number(Matrix<T>& mat, int dim, int n_cores);
template<typename T>
Rcpp::NumericVector stats_log_number_ser(Matrix<T>& mat, int dim, int n_cores);

// other helper functions
template<typename T>
T inner_prod(T *v1, T *v2, size_t size);
// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
// to get the same in random_suffle after set.seed(1) in R
inline int randWrapper(const int n) { return floor(unif_rand()*n); }



// global variables
double theta0, theta1, theta2;


/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// MAIN FUNCTION (returns statistics only) ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
//main_20151117

// [[Rcpp::export]]
Rcpp::NumericMatrix simulate_networks_fit(Rcpp::IntegerMatrix adjacency,
                                          size_t dim,
                                          Rcpp::NumericVector theta,
                                          size_t n_update = 1000000,
                                          int n_cores = 4,
                                          int nsim = 1,
                                          size_t burnin = 1000,
                                          bool log_change = true)
{

  // initialize parameters
  Matrix<char> adj_mat(dim,dim);

  Rcpp::NumericMatrix stats(nsim, 3);
  Rcpp::NumericVector log_stats_temp(3);

  size_t n_update_eff;

  // copy adjacency matrix entries
  adj_mat.init_adj(0, dim, adjacency);

  int k;
  theta0 = theta[0];
  theta1 = theta[1];
  theta2 = theta[2];

  // simulate k networks (start for next network from the previous)
  for(k=0; k<nsim; k++){
    // number of updates for the first network
    if(k==0) n_update_eff = burnin + n_update;
    // number of updates for the other networks
    if(k!=0) n_update_eff = n_update;
    if(n_cores < 2)
    {
      sim_nw_c_ser(adj_mat, dim, n_update_eff, k, log_change);
    }else{
      sim_nw_c_par(adj_mat, dim, n_update_eff, n_cores, k, log_change);
    }


    // calculate statistics for the final networks
    if(n_cores < 2)
    {
      if(log_change){
        log_stats_temp = stats_log_number_ser(adj_mat, dim, n_cores);
        stats(k,0) = log_stats_temp(0);
        stats(k,1) = log_stats_temp(1);
        stats(k,2) = log_stats_temp(2);
      }else{
        stats(k,0) = edge_number_ser(adj_mat, dim);
        stats(k,1) = star2_number_ser(adj_mat, dim);
        stats(k,2) = triangle_number_ser(adj_mat, dim);
      }

    }else{
      if(log_change){
        log_stats_temp = stats_log_number(adj_mat, dim, n_cores);
        stats(k,0) = log_stats_temp(0);
        stats(k,1) = log_stats_temp(1);
        stats(k,2) = log_stats_temp(2);
      }else{
        stats(k,0) = edge_number(adj_mat, dim, n_cores);
        stats(k,1) = star2_number(adj_mat, dim, n_cores);
        stats(k,2) = triangle_number3(adj_mat, dim, n_cores);
      }

    }
  }
  return stats;
}



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////// MAIN FUNCTION (returns at most 5 network adjacencies) ///////////////
/////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::List simulate_networks_fit_nw(Rcpp::IntegerMatrix adjacency,
                                    size_t dim,
                                    Rcpp::NumericVector theta,
                                    size_t n_update = 1000000,
                                    int n_cores = 4,
                                    int nsim = 1,
                                    size_t burnin = 1000,
                                    bool log_change = true)
{

  Rcpp::IntegerMatrix adj1(dim, dim);
  Rcpp::IntegerMatrix adj2(dim, dim);
  Rcpp::IntegerMatrix adj3(dim, dim);
  Rcpp::IntegerMatrix adj4(dim, dim);
  Rcpp::IntegerMatrix adj5(dim, dim);


  // initialize parameters
  Matrix<char> adj_mat(dim,dim);

  Rcpp::NumericMatrix stats(nsim, 3);
  Rcpp::NumericVector log_stats_temp(3);

  size_t n_update_eff;

  // copy adjacency matrix entries
  adj_mat.init_adj(0, dim, adjacency);


  int k;
  theta0 = theta[0];
  theta1 = theta[1];
  theta2 = theta[2];

  // simulate k networks (start for next network from the previous)
  for(k=0; k<nsim; k++){
    // number of updates for the first network
    if(k==0) n_update_eff = burnin + n_update;
    // number of updates for the other networks
    if(k!=0) n_update_eff = n_update;
    if(n_cores < 2)
    {
      sim_nw_c_ser(adj_mat, dim, n_update_eff, k, log_change);
    }else{
      sim_nw_c_par(adj_mat, dim, n_update_eff, n_cores, k, log_change);
    }

    if(k == 0) adj_mat.copy_adj(0, dim, adj1);
    if(k == 1) adj_mat.copy_adj(0, dim, adj2);
    if(k == 2) adj_mat.copy_adj(0, dim, adj3);
    if(k == 3) adj_mat.copy_adj(0, dim, adj4);
    if(k == 4) adj_mat.copy_adj(0, dim, adj5);

  }
  return Rcpp::List::create(Rcpp::Named("adj1") = adj1,
    Rcpp::Named("adj2") = adj2,
    Rcpp::Named("adj3") = adj3,
    Rcpp::Named("adj4") = adj4,
    Rcpp::Named("adj5") = adj5);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main simulate function parallel
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
void sim_nw_c_par(Matrix<T>& mat2, size_t dim, size_t n_update, size_t n_cores, int k, bool log_change)
{

  // correct number of updates
  size_t nrounds = (size_t) n_update/dim*2;

    // for small number of updates to at least one round
  if(nrounds == 0) nrounds = 1;
  int dim_ind;

  // ignore one row of ind_mat if dim uneven
  if(dim%2!=0){
    dim_ind = dim-1;
  }else{
    dim_ind = dim;
  }

  Matrix<size_t> ind_mat(nrounds, dim);

  for(int i=0; i<nrounds; i++ ) {
    for(int j=0; j<dim; j++ ) {
      ind_mat(i,j)=j;
    }
    // shuffle row i with fixed seed (in R)
    std::random_shuffle(&(ind_mat(i,0)), &(ind_mat(i,dim-1)) + 1, randWrapper);
  }


  // initialize number of different threads
  omp_set_num_threads(n_cores);

  #pragma omp parallel
  {
    auto tid = omp_get_thread_num();
    auto nth = omp_get_num_threads();

    // seed value for rand_r()
    unsigned seed = 31+31337*tid;

    // this thread is responsible for
    // rows [row_beg, row_end)
    size_t row_beg = (double)tid*dim/(double)nth;
    size_t row_end = (double)(tid+1)*dim/(double)nth;

    if( tid==0 )     { assert(row_beg==0); }
    if( tid==nth-1 ) { assert(row_end==dim); }


  #pragma omp barrier

    for(int round=0; round<nrounds; round++){
      for(auto i=0; i<dim_ind; i+=2){
	      // we'll be updating adj_mat(row,col)
	      int row=ind_mat(round,i);
	      int col=ind_mat(round,i+1);

	      if( row_beg<=row && row<row_end ){
	        sim_nw_thread(mat2, dim, row, col, &seed, log_change);
	      }
	    }
      // wait after each round (update adj_mat)
      #pragma omp barrier
    }
  }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// simulate serial
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void sim_nw_c_ser(Matrix<T>& mat2, size_t dim, size_t n_update, int k, bool log_change)
{
  // seed value for rand_r()
  unsigned seed = 31+31337*k;
  int i, j;

  for(int m=0; m<n_update; m++ )
  {
    // choose one dyad randomly
    do
    {
      i = rand_r(&seed)%dim;
      j = rand_r(&seed)%dim;
    }
    while(i==j);
    sim_nw_thread(mat2, dim, i, j, &seed, log_change);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper thread function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
void sim_nw_thread(Matrix<T>& mat, int dim, int i, int j, unsigned *seed, bool log_change)
{
  int k, deg_i, deg_j, mat_ij;
  double change0, change1, change2;
  double hr, accept_prob;
  int decision;

  mat_ij = mat(i,j);
  mat_ij ^= 1;

  deg_i = mat.rowsum(i);
  deg_j = mat.rowsum(j);



  if(mat_ij) {
    ++deg_i; ++deg_j;
  } else {
    --deg_i; --deg_j;
  }

  // calculate change statistics
  if(mat_ij == 1){
    // change of number of edges
    change0 = 1.0;
    // 2-stars: sum of coloum i and row j
    change1 = deg_i-1 + deg_j-1;
    // 3-stars: sum of coloum i and row j
    //  change3 = 0.5*((deg_i-1)*((deg_i-1)-1) + (deg_j-1)*((deg_j-1)-1));
    // triangle: number of common neighbors
    change2 = 0.0;
    change2 = inner_prod( &(mat(i,0)), &(mat(j,0)), dim);

  }else{
    // edges: change of number of edges
    change0 = -1.0 ;
    // 2-stars: sum of coloum i and row j
    change1 = -(deg_i + deg_j);
    // 3-stars: sum of coloum i and row j
    //  change3 = -(0.5*((deg_i)*((deg_i)-1) + (deg_j)*((deg_j)-1)));
    // triangle: number of common neighbors
    change2 = 0.0;
    change2 = -inner_prod( &(mat(i,0)), &(mat(j,0)), dim);
  }

  // hastings ratio (acceptance probability)
  if(log_change){
    hr = exp((theta0*abs(change0) + theta1*log(abs(change1) + 0.1) + theta2*log(abs(change2) + 0.1))*change0);
  }else{
    hr = exp(theta0*change0 + theta1*change1 + theta2*change2);
  }


  // accept_prob = min{1, hr}
  if(hr < 1) accept_prob = hr;
  else accept_prob = 1;

  // sample decision with acceptance probability
  // decision = 1 means toggle (change of entry 0 to 1 or 1 to 0) is accepted
  // decision = 0 means toggle is rejected
  double r = ((double)rand_r(seed) / (double)(RAND_MAX));
  if(r <= accept_prob) decision = 1;
  else decision = 0;

  // keep old adj. matrix if decision is 0 otherwise keep new adj. matrix
  if(decision == 1)
    {
      mat(i,j) = mat_ij;
      mat(j,i) = mat_ij;

      mat.rowsum(i) = deg_i;
      mat.rowsum(j) = deg_j;
    }
}






/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// NUMBER OF EDGES ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// calculate number of edges via sum of row_sums
template<typename T>
int edge_number(Matrix<T>& mat, int dim, int n_cores){
  int edge = 0;
  int k;
  // initialize number of different threads
  omp_set_num_threads(n_cores);

  #pragma omp parallel for reduction(+:edge)
  for(k=0; k<dim; k++){
      edge = edge + mat.rowsum(k);
    }
  edge = edge/2;
  return edge;
}


/////////////////////////// serial //////////////////////////////////////////////////////
// calculate number of edges via sum of row_sums
template<typename T>
int edge_number_ser(Matrix<T>& mat, int dim){
  int edge = 0;
  int k;

  for(k=0; k<dim; k++){
      edge = edge + mat.rowsum(k);
    }
  edge = edge/2;
  return edge;
}

/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// NUMBER OF 2-STARS /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// calculate number of 2-stars via summation of evaluated binomial coefficient
// (\sum^dim_i \binom(row_sum[i], 2) = \sum^dim_i row_sum[l]*(row_sum[l]-1)/2 )
template<typename T>
int star2_number(Matrix<T>& mat, int dim, int n_cores){
  // 2-stars
  int star2 = 0;
  int l;
  // initialize number of different threads
  omp_set_num_threads(n_cores);
  #pragma omp parallel for reduction (+:star2)
  for(l=0; l<dim; l++)
  {
    star2 = star2 + mat.rowsum(l)*(mat.rowsum(l)-1)/2;
  }
  return star2;
}

/////////////////////////// serial //////////////////////////////////////////////////////
template<typename T>
int star2_number_ser(Matrix<T>& mat, int dim){
  // 2-stars
  int star2 = 0;
  int l;

  for(l=0; l<dim; l++)
  {
    star2 = star2 + mat.rowsum(l)*(mat.rowsum(l)-1)/2;
  }
  return star2;
}


/////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// NUMBER OF TRIANGLES ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// calculate number of triangles via upper matrix multiplication and blocking
template<typename T>
int triangle_number3(Matrix<T>& mat, int dim, int n_cores){

  int triangle = 0;
  int tria;

  int k_b, k, p_b, p, l_b,l;
  int blocksize=100;

  // initialize number of different threads
  omp_set_num_threads(n_cores);

  int last_block = dim%blocksize;
  int last_block_start = dim-last_block;
  int block;

  for(k_b=0; k_b<last_block_start; k_b+=blocksize){
    for(l_b=0; l_b<last_block_start; l_b+=blocksize){


    #pragma omp parallel for reduction(+:triangle)
      for(p_b=0; p_b<last_block_start; p_b+=blocksize){
          block = blocksize;
          for(k=k_b; k<(k_b+block); k++){
            for(l=l_b; l<(l_b+block); l++){
              for(p=p_b; p<(p_b+block); p++){
                tria=mat(k,p)*mat(p,l)*mat(l,k);
                triangle=tria+triangle;
              }
            }
          }
        }
      }
    }

   //if last block is reached do rest, if there are uncomplete blocks
   if(last_block_start < dim){
    for(k=last_block_start; k<dim; k++){
      for(l=0; l<dim; l++){
        for(p=0; p<dim; p++){
         tria=mat(k,p)*mat(p,l)*mat(l,k);
         triangle=tria+triangle;
        }
      }
    }
    for(k=0; k<last_block_start; k++){
      for(l=last_block_start; l<dim; l++){
        for(p=0; p<last_block_start; p++){
          tria=mat(k,p)*mat(p,l)*mat(l,k);
          triangle=tria+triangle;
        }
      }
    }
    for(k=0; k<last_block_start; k++){
      for(l=0; l<dim; l++){
        for(p=last_block_start; p<dim; p++){
          tria=mat(k,p)*mat(p,l)*mat(l,k);
          triangle=tria+triangle;
        }
      }
    }
    triangle = triangle/6;

    return triangle;
  }

  triangle = triangle/6;

  return triangle;
}


/////////////////////////// serial //////////////////////////////////////////////////////
template<typename T>
int triangle_number_ser(Matrix<T>& mat, int dim){
  int triangle = 0;
  int tria, k, l, p;

  for(k=0; k<dim; k++)
  {
    for(l=0; l<dim; l++)
    {
      if(l > k)
      {
        for(p=0; p<dim; p++)
        {
          if(p > l)
          {
           tria = mat(k,p)*mat(p,l)*mat(l,k);
           triangle = tria + triangle;
          }
        }
      }
    }
  }

  return triangle;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper thread function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns log-statistics of edges, 2-star, triangle (combined in one loop)
// parallel
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
Rcpp::NumericVector stats_log_number(Matrix<T>& mat, int dim, int n_cores){
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns log-statistics of edges, 2-star, triangle (combined in one loop)
// serial
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
Rcpp::NumericVector stats_log_number_ser(Matrix<T>& mat, int dim, int n_cores){
  Rcpp::NumericVector stats_log(3);
  int edge = 0;
  double star2 = 0;
  double triangle = 0;
  double sum_star = 0;
  double sum_tria = 0;
  int k, l, p;

  for(k=0; k<dim; k++)
  {
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




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// other helper function
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TYPE>
TYPE inner_prod(TYPE *v1, TYPE *v2, size_t size)
{
  TYPE sum=0;
  for(size_t i=0; i<size; ++i ) {
    sum+=(*v1++)*(*v2++);
  }
  return sum;
}




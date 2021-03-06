#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

#include <unistd.h>
#include <iostream>
using namespace std;

template<typename TYPE=int>
class Matrix
{
private:
  TYPE   *m_data=0;
//  size_t  m_dim=0;
  size_t  m_dim_row=0;
  size_t  m_dim_col=0;
  size_t *m_rowsum=0;

public:
  Matrix(size_t dim_row, size_t dim_col) {
    m_dim_row = dim_row;
    m_dim_col = dim_col;
    m_data = new TYPE[dim_row*dim_col];
    m_rowsum = new size_t[dim_row];
  }

  ~Matrix() {
    delete[] m_data;
    delete[] m_rowsum;
  }

  TYPE& operator()(size_t row, size_t col) {
    return m_data[row*m_dim_col+col];
  }

  size_t& rowsum(size_t row) {
    return m_rowsum[row];
  }

  void init_rows(size_t row_beg, size_t row_end, TYPE val=0) {

    for( auto i=row_beg; i<row_end; ++i ) {
      for( auto j=0; j<m_dim_col; j++ ) {
        m_data[i*m_dim_col+j]=val;
      }
      m_rowsum[i]=(size_t)val*m_dim_col;
    }
  }

  void init_adj(size_t row_beg, size_t row_end, Rcpp::IntegerMatrix& adjacency2) {
    for( auto i=row_beg; i<row_end; ++i ) {
      size_t sum = 0;
      for( auto j=0; j<m_dim_col; j++ ) {
        m_data[i*m_dim_col+j] = adjacency2(i,j);
        sum = sum + adjacency2(i,j);
      }
      m_rowsum[i] = sum;
    }
  }

  void copy_adj(size_t row_beg, size_t row_end, Rcpp::IntegerMatrix& adjacency2) {
    for( auto i=row_beg; i<row_end; ++i ) {
      size_t sum = 0;
      for( auto j=0; j<m_dim_col; j++ ) {
         adjacency2(i,j) = m_data[i*m_dim_col+j];
        sum = sum + adjacency2(i,j);
      }
      m_rowsum[i] = sum;
    }
  }


  void print() {

    for( auto i=0; i<m_dim_row; i++ ) {
      for( auto j=0; j<m_dim_col; j++ ) {
	if( operator()(i,j)!=operator()(j,i) ) {
	  cerr<<"Not symmetric!"<<endl;
	}
      }
    }

   if( m_dim_row<=3000 ) {
      for( auto i=0; i<m_dim_row; i++ ) {
	cout<<rowsum(i)<<": ";
	for( auto j=0; j<m_dim_col; j++ ) {
	  TYPE val = operator()(i,j);
	  if( val )
	    cout<<val<<" 1 ";
	  else
	    cout<<" 0 ";
	}
	cout<<endl;
      }
    }
}

};


#endif /* TYPES_H_INCLUDED*/

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;


class CellAuto {
public:
    CellAuto(Rcpp::IntegerMatrix mat_sp, Rcpp::IntegerMatrix mat_suc, 
             Rcpp::NumericVector c_e, double prob_distur);
    Rcpp::IntegerMatrix returnSp();
    Rcpp::IntegerMatrix returnSucc();
    void update();
    void updateCell(int i, int j);
    void iterate(unsigned int iterations);
private:
    Rcpp::IntegerMatrix m_mat_sp;
    Rcpp::IntegerMatrix m_mat_suc;
    Rcpp::IntegerMatrix m_otherLandscape;
    Rcpp:: NumericVector m_c_e;
    int m_nrow;
    int m_ncol;
    double m_prob_distur;
};


CellAuto::CellAuto(Rcpp::IntegerMatrix mat_sp, Rcpp::IntegerMatrix mat_suc, 
		   Rcpp::NumericVector c_e, double prob_distur):
                                            m_mat_sp(mat_sp), 
					    m_mat_suc(mat_suc), 
					    m_c_e(c_e), 
					    m_prob_distur(prob_distur)
{
  m_nrow = mat_sp.nrow();
  m_ncol = mat_sp.nrow();
}

Rcpp::IntegerMatrix CellAuto::returnSp() {
  return m_mat_sp;
}

Rcpp::IntegerMatrix CellAuto::returnSucc() {
  return m_mat_sp;
}



void CellAuto::iterate( unsigned int iterations ) {
    for ( int i = 0; i < iterations; i++ ) {
        update();
    }
}



void CellAuto::update() {
        for ( int ii = 0; ii < m_nrow; ii++ ) {
            for ( int jj = 0; jj < m_ncol; jj++ ) {
	      CellAuto::updateCell(ii , jj);
            }
        }
}
 

void CellAuto::updateCell(int ii, int jj) {
   bool test;
    for( int i = -1; i < 2; i++ ) {
        for( int j = -1; j < 2; j++ ) {
          int i_n = ( m_nrow + ii + i ) % m_nrow;
   	  int j_n = jj + j;
 	  if(j_n < 0) j_n = 1;
 	  if(j_n > m_ncol-1) j_n = m_ncol -2;
          test = m_c_e[m_mat_sp(i_n,j_n)] > m_c_e[m_mat_sp(ii,jj)];
	  m_mat_sp = (test) ? m_mat_sp(i_n,j_n) : m_mat_sp(ii,jj);
	}
    }
m_mat_sp(ii , jj) = (R::runif(0,1) < m_prob_distur) ? 0 : m_mat_sp(ii , jj);
}


// [[Rcpp::export]]
Rcpp::List UpdateIterR(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
		       NumericVector c_e, double prob_distur, int n){

    CellAuto cells(mat_sp, mat_suc, 
                   c_e,prob_distur);
    cells.iterate(n);
    return Rcpp::List::create(Rcpp::Named("sp") = cells.returnSp(),
                              Rcpp::Named("suc") = cells.returnSucc());
}


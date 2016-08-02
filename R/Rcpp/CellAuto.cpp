#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

int Torus(int i, int ii, int i_max);

int Reflexion(int i, int ii, int i_max);


class CellAuto {
public:
    CellAuto(Rcpp::IntegerMatrix mat_sp, Rcpp::IntegerMatrix mat_suc,
                   int nrow, int ncol,
             Rcpp::NumericVector c_e, double prob_distur);
    Rcpp::IntegerMatrix returnSp() const;
    Rcpp::IntegerMatrix returnSucc() const;
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



int Torus(int i, int ii, int i_max){
 int i_n = ( i_max + ii + i ) % i_max;

 return i_n;
}
int Reflexion(int i, int ii, int i_max){
  int i_n = ii + i;
 if(i_n < 0) i_n = 1;
 if(i_n > i_max-1) i_n = i_max -2;

 return i_n;
}


CellAuto::CellAuto(Rcpp::IntegerMatrix mat_sp, Rcpp::IntegerMatrix mat_suc,
                   int nrow, int ncol,
		   Rcpp::NumericVector c_e, double prob_distur):
                                            m_mat_sp(mat_sp),
					    m_mat_suc(mat_suc),
					    m_c_e(c_e),
					    m_prob_distur(prob_distur),
					    m_nrow(nrow),
					    m_ncol(ncol)
{
}

Rcpp::IntegerMatrix CellAuto::returnSp() const{
  return m_mat_sp;
}

Rcpp::IntegerMatrix CellAuto::returnSucc() const{
  return m_mat_suc;
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
	  int i_n = Torus(i, ii, m_nrow);
	  int j_n = Reflexion(j, jj, m_ncol);
          test = m_c_e[m_mat_sp(i_n,j_n)] > m_c_e[m_mat_sp(ii,jj)];
	  m_mat_sp(ii, jj) = (test) ? m_mat_sp(i_n,j_n) : m_mat_sp(ii,jj);
	}
    }
m_mat_sp(ii , jj) = (R::runif(0,1) < m_prob_distur) ? 0 : m_mat_sp(ii , jj);
}


// [[Rcpp::export]]
Rcpp::List UpdateIterR(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
                       int nrow, int ncol,
		       NumericVector c_e, double prob_distur, int n){

    CellAuto cells(mat_sp, mat_suc,
                   nrow, ncol,
                   c_e,prob_distur);
    cells.iterate(n);
    return Rcpp::List::create(Rcpp::Named("sp") = cells.returnSp(),
                              Rcpp::Named("suc") = cells.returnSucc());
}


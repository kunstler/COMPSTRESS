#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

int Torus(int i, int ii, int i_max);

int Reflexion(int i, int ii, int i_max);

int ColonizeNeigh(int ii, int jj , IntegerMatrix mat_sp,  IntegerMatrix mat_suc,
                  NumericVector c_e, NumericVector c_l, int i_max, int j_max,
		  double K) ;


int KillCellDisturb(int ii, int jj, IntegerMatrix mat_sp,  double prob_distur);


int KillCellStress(int ii, int jj, IntegerMatrix mat_sp,
                   NumericVector c_s, NumericVector ss);

int SucCell(int ii, int jj, IntegerMatrix mat_sp, IntegerMatrix mat_suc,
             double prob_suc);

double Prob_Law(double x, double y, double K);

int Sampel_Law(int ii, int jj, int ni, int nj,
               IntegerMatrix mat_sp, NumericVector c_,
               double K);



class CellAuto {
public:
    CellAuto(Rcpp::IntegerMatrix mat_sp, Rcpp::IntegerMatrix mat_suc,
             int nrow, int ncol,
	     Rcpp::NumericVector c_e, Rcpp::NumericVector c_, Rcpp::NumericVector c_s,
	     Rcpp::NumericVector ss,
	     double prob_distur, double prob_suc, double K);
    Rcpp::IntegerMatrix returnSp() const;
    Rcpp::IntegerMatrix returnSuc() const;
    void update();
    void updateCell(int i, int j);
    void iterate(unsigned int iterations);
private:
    Rcpp::IntegerMatrix m_mat_sp;
    Rcpp::IntegerMatrix m_mat_suc;
    Rcpp::IntegerMatrix m_otherLandscape;
    Rcpp:: NumericVector m_c_e;
    Rcpp:: NumericVector m_c_l;
    Rcpp:: NumericVector m_c_s;
    Rcpp:: NumericVector m_ss;
    int m_nrow;
    int m_ncol;
    double m_prob_distur;
    double m_prob_suc;
    double m_K;
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
		   Rcpp::NumericVector c_e, Rcpp::NumericVector c_l, Rcpp::NumericVector c_s,
		   Rcpp::NumericVector ss,
		   double prob_distur, double prob_suc, double K):
                                            m_mat_sp(mat_sp),
					    m_mat_suc(mat_suc),
					    m_c_e(c_e),
					    m_c_l(c_l),
					    m_c_s(c_s),
					    m_ss(ss),
					    m_nrow(nrow),
					    m_ncol(ncol),
					    m_prob_distur(prob_distur),
					    m_prob_suc(prob_suc),
					    m_K(K)
{
}

Rcpp::IntegerMatrix CellAuto::returnSp() const{
  return m_mat_sp;
}

Rcpp::IntegerMatrix CellAuto::returnSuc() const{
  return m_mat_suc;
}



void CellAuto::iterate( unsigned int iterations ) {
    for ( int i = 0; i < iterations; i++ ) {
        update();
    }
}



void CellAuto::update() {
// Rcpp::Rcout << "m_row and m_ncol = " << m_nrow << " and "<< m_ncol << std::endl;

        for ( int ii = 0; ii < m_nrow; ii++ ) {
            for ( int jj = 0; jj < m_ncol; jj++ ) {
	      // Rcpp::Rcout << "ii and jj = " << ii << " and "<< jj << std::endl;
	      CellAuto::updateCell(ii , jj);
            }
        }
}


void CellAuto::updateCell(int ii, int jj) {
  // Rcpp::Rcout << "start colo" << std::endl;
  m_mat_sp(ii , jj) = ColonizeNeigh(ii, jj,
				    m_mat_sp, m_mat_suc,
	                            m_c_e, m_c_l,
				    m_nrow, m_ncol,
				    m_K);
  // Rcpp::Rcout << "start kill" << std::endl;
  m_mat_sp(ii , jj) = KillCellDisturb(ii, jj, m_mat_sp, m_prob_distur);
  m_mat_sp(ii, jj) = KillCellStress(ii, jj, m_mat_sp,
  				    m_c_s, m_ss);
  m_mat_suc(ii, jj) = SucCell(ii, jj, m_mat_sp, m_mat_suc, m_prob_suc);
}


int ColonizeNeigh(int ii, int jj , IntegerMatrix mat_sp,  IntegerMatrix mat_suc,
                  NumericVector c_e, NumericVector c_l, int i_max, int j_max,
		  double K) {
  int i_n, j_n, res;
    for( int i = -1; i < 2; i++ ) {
        for( int j = -1; j < 2; j++ ) {
            i_n = Torus(i, ii, i_max);
            j_n = Reflexion(j, jj, j_max);

	    if(mat_suc(ii, jj) == 1){
	    	mat_sp(ii, jj) = Sampel_Law(ii, jj, i_n, j_n, mat_sp, c_e, K);
             }else{
	    	mat_sp(ii, jj) = Sampel_Law(ii, jj, i_n, j_n, mat_sp, c_l, K);
             }
        }
    }
    return mat_sp(ii, jj);
}


int KillCellDisturb(int ii, int jj, IntegerMatrix mat_sp,  double prob_distur){
    int res;
    if(R::runif(0,1) < prob_distur){
	res = 0;
    } else { 
        res = mat_sp(ii , jj);
    }
  return res;
}


int KillCellStress(int ii, int jj, IntegerMatrix mat_sp,
                   NumericVector c_s, NumericVector ss){
    bool test;
    int res;
    if(c_s[mat_sp(ii, jj)] != -100){
	if(R::runif(0,1) < (ss[jj] * c_s[mat_sp(ii, jj)])){
	    res = 0;
	} else {
            res =  mat_sp(ii , jj);
	}     
    } else{
	res = 0;
    }
    return res;
}

int SucCell(int ii, int jj, IntegerMatrix mat_sp, IntegerMatrix mat_suc,
             double prob_suc){
  int res;

  res = (R::runif(0,1) < prob_suc) ? 2 : mat_suc(ii , jj);
  if(mat_sp(ii, jj) == 0) res = 1;
  return res;
}

double Prob_Law(double x, double y, double K){

  return  1 / (1 + exp(-K*(x - y)));
}

int Sampel_Law(int ii, int jj, int ni, int nj,
               IntegerMatrix mat_sp, NumericVector c_,
               double K){
  bool test;
  // Rcpp::Rcout << "ni and nj " << ni << " , " << nj << std::endl;
  // Rcpp::Rcout << "val sp " << mat_sp(ni, nj) << std::endl;
  // Rcpp::Rcout << "val c " << c_[mat_sp(ni, nj)]  << " and " << c_[mat_sp(ii, jj)]<< std::endl;
  // Rcpp::Rcout << "p law " << Prob_Law(c_[mat_sp(ni, nj)],c_[mat_sp(ii, jj)], K) << std::endl;
  int res;
  double pval = Prob_Law(c_[mat_sp(ni, nj)],c_[mat_sp(ii, jj)],K);
  if(c_[mat_sp(ii, jj)] != -100){
      if(c_[mat_sp(ni, nj)] != -100 ){
          if(R::runif(0,1) < pval){
	      res = mat_sp(ni, nj);
	  } else {
              res = mat_sp(ii, jj);
	  }
        } else {
          res = mat_sp(ii, jj);
        }
  }else{
      if(c_[mat_sp(ni, nj)] != -100 ){
          res = mat_sp(ni, nj);
        }else{
          res = mat_sp(ii, jj);
        }
  }   
  // Rcpp::Rcout << "val sp col " <<  res << std::endl;

  return res;
}



// [[Rcpp::export]]
Rcpp::List UpdateIterR(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
		       NumericVector c_e, NumericVector c_l,NumericVector c_s,
                       NumericVector ss,
		       double prob_distur, double prob_suc, double K, int n){

    int nrow = mat_sp.nrow();
    int ncol = mat_sp.ncol();

    CellAuto cells(mat_sp, mat_suc,
                   nrow, ncol,
                   c_e, c_l, c_s, ss,
		   prob_distur, prob_suc, K);
    cells.iterate(n);
    return Rcpp::List::create(Rcpp::Named("sp") = cells.returnSp(),
                              Rcpp::Named("suc") = cells.returnSuc());
}


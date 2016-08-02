#include <Rcpp.h>
// #include<cmath>
// #include <vector>
// #include <iterator>
// #include <algorithm>
using namespace Rcpp;

IntegerMatrix UpdateCells(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
                          NumericVector c_e, NumericVector c_l,
                          int i_max, int j_max);
void UpdateCellsTest(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
                          NumericVector c_e, NumericVector c_l,
		     int i_max, int j_max);

int Torus(int i, int ii, int i_max);

int Reflexion(int i, int ii, int i_max);

int ColonizeNeigh(int ii, int jj , IntegerMatrix mat_sp,  IntegerMatrix mat_suc,
                  NumericVector c_e, NumericVector c_l, int i_max, int j_max);

int KillCellDisturb(int ii, int jj, IntegerMatrix mat_sp,
                    double prob_distur=0.1);
int KillCellStress(int ii, int jj, IntegerMatrix mat_sp,
                   NumericVector c_e, NumericVector c_l);

int SuccCell(int ii, int jj, IntegerMatrix mat_sp, IntegerMatrix mat_suc,
             double prob_succ = 0.2);

double Prob_Law(float x, float y, float K = 10);

int Sampel_Law(int ii, int jj, int ni, int nj,
               IntegerMatrix mat_suc, NumericVector c_,
               float K = 10);


// [[Rcpp::export]]
Rcpp::List UpdateIterTest(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
                          NumericVector c_e, NumericVector c_l, int n){
int nrow = mat_sp.nrow(), ncol = mat_sp.ncol();
  for(int nn=1; nn<n+1;  nn++){
    UpdateCellsTest(mat_sp, mat_suc,  c_e, c_l, nrow, ncol);
    // mat_sp = list_mat_ss["sp"];
    // mat_suc = list_mat_ss["suc"];
   }

return Rcpp::List::create(Rcpp::Named("sp") = mat_sp,
                          Rcpp::Named("suc") = mat_suc);
}


// [[Rcpp::export]]
Rcpp::List UpdateIter(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
                         NumericVector c_e, NumericVector c_l, int n){
int nrow = mat_sp.nrow(), ncol = mat_sp.ncol();
  for(int nn=1; nn<n+1;  nn++){
    mat_sp = UpdateCells(mat_sp, mat_suc,  c_e, c_l, nrow, ncol);
    // mat_sp = list_mat_ss["sp"];
    // mat_suc = list_mat_ss["suc"];
   }

return Rcpp::List::create(Rcpp::Named("sp") = mat_sp,
                          Rcpp::Named("suc") = mat_suc);
}

IntegerMatrix UpdateCells(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
                          NumericVector c_e, NumericVector c_l,
                          int i_max, int j_max){
    for(int i=0; i<i_max; i++){
      for(int j=0; j<j_max; j++){
	mat_sp(i , j) = ColonizeNeigh(i, j, mat_sp, mat_suc,
          			     c_e, c_l,
                                     i_max, j_max);
	mat_sp(i , j) = KillCellDisturb(i, j, mat_sp);
        mat_suc(i, j) = SuccCell(i, j, mat_sp, mat_suc);

	// tt(i , j) = KillCellStress(i, j, tt, sp);
      }
    }

    return mat_sp;
// return Rcpp::List::create(Rcpp::Named("sp") = mat_sp,
//                           Rcpp::Named("suc") = mat_suc);
}


void UpdateCellsTest(IntegerMatrix mat_sp, IntegerMatrix mat_suc,
                          NumericVector c_e, NumericVector c_l,
                          int i_max, int j_max){
    for(int i=0; i<i_max; i++){
      for(int j=0; j<j_max; j++){
	mat_sp(i , j) = ColonizeNeigh(i, j, mat_sp, mat_suc,
          			     c_e, c_l,
                                     i_max, j_max);
	mat_sp(i , j) = KillCellDisturb(i, j, mat_sp);
        mat_suc(i, j) = SuccCell(i, j, mat_sp, mat_suc);

	// tt(i , j) = KillCellStress(i, j, tt, sp);
      }
    }

}


// [[Rcpp::export]]
int Torus(int i, int ii, int i_max){
 int i_n = ( i_max + ii + i ) % i_max;

 return i_n;
}
// [[Rcpp::export]]
int Reflexion(int i, int ii, int i_max){
  int i_n = ii + i;
 if(i_n < 0) i_n = 1;
 if(i_n > i_max-1) i_n = i_max -2;

 return i_n;
}


// [[Rcpp::export]]
int ColonizeNeigh(int ii, int jj , IntegerMatrix mat_sp,  IntegerMatrix mat_suc,
                  NumericVector c_e, NumericVector c_l, int i_max, int j_max) {
  int i_n, j_n, res;
    for( int i = -1; i < 2; i++ ) {
        for( int j = -1; j < 2; j++ ) {
            i_n = Torus(i, ii, i_max);
            j_n = Reflexion(j, jj, j_max);
	    switch(mat_suc(ii, jj)){
	    case 1:
              res = Sampel_Law(ii, jj, i_n, j_n, mat_sp, c_e);
            case 2:
              res = Sampel_Law(ii, jj, i_n, j_n, mat_sp, c_l);
            }
        }
    }
    return res;
}

int KillCellDisturb(int ii, int jj, IntegerMatrix mat_sp,  double prob_distur){

  return (R::runif(0,1) < prob_distur) ? 0 : mat_sp(ii , jj);
}


int KillCellStress(int ii, int jj, IntegerMatrix mat_sp,
                   NumericVector c_e, NumericVector c_l){

  return (R::runif(0,1) < (1-c_e[mat_sp(ii, jj)]) &&
              c_l[mat_sp(ii, jj)] != -100) ?
                                                           0 : mat_sp(ii , jj);
}

int SuccCell(int ii, int jj, IntegerMatrix mat_sp, IntegerMatrix mat_suc,
             double prob_succ){
  int res;
  if(mat_sp(ii , jj)>0){
  res = (R::runif(0,1) < prob_succ) ? 2 : mat_suc(ii , jj);
  } else {
  res = 1;
  }
  return res;
}

double Prob_Law(float x, float y, float K){

  return  1 / (1 + exp(-K*(x - y)));
}

int Sampel_Law(int ii, int jj, int ni, int nj,
               IntegerMatrix mat_sp, NumericVector c_,
               float K){
  bool test;
  int res;
  if(c_[mat_sp(ni, nj)] == -100) {
    res =mat_sp(ii , jj);
    }
  else {
  test = R::runif(0,1) < Prob_Law(c_[mat_sp(ni, nj)],
                                  c_[mat_sp(ii, jj)],
			 K);
  res = (test) ? mat_sp(nj , nj) : mat_sp(ii , jj);
    }
  return res;
}


// int NmatchC(NumericVector x, double v) {
//   int n = std::distance(x.begin(), std::find(x.begin(), x.end(), v))
//   return n;
// }






#include <Rcpp.h>
using namespace Rcpp;

//' Leading NA
//'
//' Returns the TMTI_infinity statistic from a pre-sorted vector of p-values
//'
//' @param pvals A NumericVector
//' @export
// [[Rcpp::export]]
List calc_TMTI_inf(NumericVector pvals) {
  List ret;
  int m = pvals.size();
  // NumericVector Y (m);
  double currentMin = 1;
  double currentY = 0;
  double previousY = 1;
  int pos = 0;
  for(int i = 0; i < m; ++i) {
    currentY = R::pbeta(pvals[i], i + 1, m + 1 - (i + 1), true, false);
    if(currentY < currentMin){
      currentMin = currentY;
      pos = i;
    }
    previousY = currentY;
  }
  // NumericVector::iterator it = std::min_element(Y.begin(), Y.end());
  ret["val"] = currentMin;
  ret["pos"] = pos + 1;
  return ret;
}

//' Leading NA
//'
//' Returns the truncated TMTI_infinity statistic from pre-sorted,
//' pre-truncated p-values.
//'
//' @param pvals A NumericVector containing the truncated sorted p-values. It
//' is important that this vector: 1) contains only the truncated p-values (i.e,
//' those that fall below the truncation point) and 2) is sorted.
//' @param m_full The total (i.e., non-truncated) number of p-values.
//' @export
// [[Rcpp::export]]
List calc_truncatedTMTI_inf(NumericVector pvals, int m_full) {
  List ret;
  int m = pvals.size();
  double currentMin = 1;
  double currentY = 0;
  double previousY = 1;
  int pos = 0;
  for(int i = 0; i < m; ++i) {
    currentY = R::pbeta(pvals[i], i + 1, m_full + 1 - (i + 1), true, false);
    if(currentY < currentMin){
      currentMin = currentY;
      pos = i;
    }
    previousY = currentY;
  }
  ret["val"] = currentMin;
  ret["pos"] = pos + 1;
  return ret;
}

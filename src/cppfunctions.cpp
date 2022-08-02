#include <Rcpp.h>
using namespace Rcpp;

//' Leading NA
//'
//' Returns the TMTI_infinity statistic from pre-sorted,
//' pre-truncated vector of p-values. If no truncation is used, set m_full = m
//'
//' @param pvals A NumericVector containing the truncated sorted p-values. It
//' is important that this vector: 1) contains only the truncated p-values (i.e,
//' those that fall below the truncation point) and 2) is sorted.
//' @param m The total (i.e., non-truncated) number of p-values.
// [[Rcpp::export]]
double MakeZ_C(NumericVector pvals, int m) {
  int m_p = pvals.size();
  double currentMin = 1;
  double currentY = 0;
  for(int i = 0; i < m_p; ++i) {
    currentY = R::pbeta(pvals[i], i + 1, m + 1 - (i + 1), true, false);
    if(currentY < currentMin){
      currentMin = currentY;
    }
  }
  return currentMin;
}

//' Leading NA
//'
//' Returns the transformed p-values (Y) from pre-sorted p-values and
//' pre-truncated p-values. If not truncation is used, set m_full = m
//'
//' @param pvals A NumericVector containing the truncated sorted p-values. It
//' is important that this vector: 1) contains only the truncated p-values (i.e,
//' those that fall below the truncation point) and 2) is sorted.
//' @param m The total (i.e., non-truncated) number of p-values.
// [[Rcpp::export]]
NumericVector MakeY_C(NumericVector pvals, int m) {
  int m_p = pvals.size();
  NumericVector Y(m_p);
  for(int i = 0; i < m_p; ++i) {
    Y[i] = R::pbeta(pvals[i], i + 1, m + 1 - (i + 1), true, false);
  }
  return Y;
}


//' Leading NA
//'
//' Returns the transformed p-values (Y) from pre-sorted p-values and
//' pre-truncated p-values when n < m - 1
//'
//' @param pvals A NumericVector containing the truncated sorted p-values. It
//' is important that this vector: 1) contains only the truncated p-values (i.e,
//' those that fall below the truncation point) and 2) is sorted.
//' @param n A positive number (or Inf) indicating which type of local minimum
//' to consider. Defaults to Infm, corresponding to the global minimum.
//' @param m The total (i.e., non-truncated) number of p-values.
// [[Rcpp::export]]
double MakeZ_C_nsmall(NumericVector pvals, int n, int m) {
  int m_p = pvals.size();
  double PreviousY = 1;
  double Y =  0;
  std::vector<double> LeadingY;
  double LeadingMin;
  double Z = -1;
  for (int i = 1; i <= n; i++) {
    LeadingY.push_back(R::pbeta(pvals[i],
                                i + 1,
                                m - i,
                                true,
                                false));
  }
  for (int i = 0; i < m_p - n; i++) {
    Y = R::pbeta(pvals[i], i + 1, m - i, true, false);
    LeadingMin = *std::min_element(LeadingY.begin(),
                                   LeadingY.end());
    if ((Y < LeadingMin) & (Y < PreviousY)) {
      Z = Y;
      break;
    }
    LeadingY.erase(LeadingY.begin());
    LeadingY.push_back(R::pbeta(pvals[i + 1 + n],
                                (i + 1 + n) + 1,
                                m - (i + 1 + n),
                                true,
                                false));
    PreviousY = Y;
  }
  if (Z == -1) {
    Z = LeadingMin;
  }
  return Z;
}


//' Leading NA
//'
//' Tests a user-specified subset in a CTP, using a user-supplied local test
//'
//'
//' @param LocalTest A function that returns a double in (0, 1).
//' @param pSub A vector with the p-values of the set to be tested.
//' @param pRest A vector containing the remaining p-values.
//' @param alpha Double indicating the significance level.
//' @param is_subset_sequence Logical indicating whether the supplied subset of
//' p_values corresponds to the pSub.size() smallest overall p-values.
//' @param EarlyStop Logical indicating whether to exit as soon as a non-significant
//' p-value is found.
//' @param verbose Logical indicating whether to print progress.
// [[Rcpp::export]]
double TestSet_C (
    Function LocalTest,
    std::vector<double> pSub,
    std::vector<double> pRest,
    double alpha,
    bool is_subset_sequence,
    bool EarlyStop,
    bool verbose
) {
  int n = pRest.size();
  int n2 = pSub.size();
  double p;
  double currentMax = 0;

  if (verbose) {
    for (int i = 0; i < n; i++) {
      Rcout << "\r" << i + 1 << " of " << n;
      auto it = pSub.begin() + n2;
      pSub.insert(it, pRest.back());
      if (not is_subset_sequence) {
        std::partial_sort(pSub.begin(), pSub.begin() + n2 + 1, pSub.end());
      }
      pRest.pop_back();
      p = *REAL(LocalTest(pSub));
      if (p > currentMax) {
        currentMax = p;
      }
      if ((p > alpha) & (EarlyStop)) {
        break;
      }
    }
  } else {
    for (int i = 0; i < n; i++) {
      auto it = pSub.begin() + n2;
      pSub.insert(it, pRest.back());
      if (not is_subset_sequence) {
        std::partial_sort(pSub.begin(), pSub.begin() + n2 + 1, pSub.end());
      }
      pRest.pop_back();
      p = *REAL(LocalTest(pSub));
      if (p > currentMax) {
        currentMax = p;
      }
      if ((p > alpha) & (EarlyStop)) {
        break;
      }
    }
  }

  return currentMax;
}


//' Leading NA
//'
//' Tests a user-specified subset in a CTP, using a user-supplied local test
//'
//'
//' @param LocalTest A function that returns a double in (0, 1).
//' @param f A function that iterates LocalTest over the relevant test tree.
//' In practice, this is called as TestSet_C.
//' @param pvals A vector of p-values.
//' @param EarlyStop Logical indicating whether to exit as soon as a non-significant
//' p-value is found.
//' @param alpha Significance level. This is only used if EarlyStop = TRUE
// [[Rcpp::export]]
std::vector<double> FullCTP_C (Function LocalTest,
                            Function f,
                            std::deque<double> pvals,
                            bool EarlyStop,
                            double alpha) {
  std::vector<double> BottomTrees;
  std::vector<double> TopTree;
  std::vector<double> out;
  double max;
  int m_max;

  int m = pvals.size();

  if (EarlyStop) {
    m_max = m - 1;
    for (int i = 0; i < m_max; i++) {
      TopTree.push_back(*REAL(LocalTest(pvals)));
      double p = pvals.front();
      pvals.pop_front();
      BottomTrees.push_back(*REAL(f(p, pvals)));
      if ((TopTree.back() > alpha) | (BottomTrees.back() > alpha)) {
        break;
      }
    }

    int mm = TopTree.size();
    for (int i = 0; i < mm; i++) {
      max = *std::max_element(TopTree.begin(), TopTree.begin() + i + 1);
      if (BottomTrees[i] > max) {
        out.push_back(BottomTrees[i]);
      } else {
        out.push_back(max);
      }
    }
  } else {
    m_max = m - 1;
    for (int i = 0; i < m_max; i++) {
      TopTree.push_back(*REAL(LocalTest(pvals)));
      double p = pvals.front();
      pvals.pop_front();
      BottomTrees.push_back(*REAL(f(p, pvals)));
    }
    TopTree.push_back(pvals.front());


    for (int i = 0; i < m; i++) {
      max = *std::max_element(TopTree.begin(), TopTree.begin() + i + 1);
      if (BottomTrees[i] > max) {
        out.push_back(BottomTrees[i]);
      } else {
        out.push_back(max);
      }
    }
  }


  return out;
}


//' Leading NA
//'
//' Computes a confidence set for the number of false hypotheses among all hypotheses
//'
//'
//' @param LocalTest A function that returns a double in (0, 1).
//' @param pvals A vector of p-values.
//' @param alpha A double indicating the significance level
// [[Rcpp::export]]
int TopDown_C (Function LocalTest,
               std::deque<double> pvals,
               double alpha) {
  int m = pvals.size();
  int t_alpha = 0;
  double p;
  for(int i = 0; i < m; i++) {
    p = *REAL(LocalTest(pvals));
    if (p > alpha) {
      t_alpha = pvals.size();
      break;
    } else {
      pvals.pop_front();
    }
  }

  return m - t_alpha;
}

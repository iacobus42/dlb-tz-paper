#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::export]]
NumericVector matchCpp(NumericVector case_ps,
                       NumericVector case_fu,
                       NumericVector control_ps,
                       NumericVector control_fu,
                       double caliper,
                       bool display_progress = true) {
  int n_cases = case_ps.size();

  NumericVector matches(n_cases);

  // trying to get rcpp to show a progress bar
  Progress p(n_cases, display_progress);
  
  for (int i = 0; i < n_cases; ++i) {
    // trying to get rcpp to show a progress bar
    p.increment(); 
    // compute the distance between the case PS and the control PS score
    // we want to take the nearest but don't care about direction, so abs()
    NumericVector distance = abs(case_ps[i] - control_ps);
    // we'll keep track of the index of the nearest minimum distance using
    // min_index and current_min. Set these to -1 for index (e.g., null match)
    // and an impossibly large distance for the current_min distance
    int min_index = -1;
    double current_min = 999999;

    for (int j = 0; j < distance.size(); ++j) {
      // check that the duration of follow-up is +/- 90 days
      if (abs(case_fu[i] - control_fu[j]) <= 90) {
        // if less than 90 days difference, check for minimum PS
        // if the distance between cases[i] and controls[j] is less than the
        // current_min distance, change min_index to reflect the new best match
        // and current_min to reflect the new shortest distance
        if (distance[j] < current_min && distance[j] < caliper) {
          min_index = j;
          current_min = distance[j];
        }
      }
    }
    // we are matching without replacement, so once we have the greedy best
    // match for cases[i], set the controls[min_index] to a very large number
    // to ensure that it will never again be matched to one of the case
    // observations
    if (min_index != -1) {
      control_ps[min_index] = 99999999;
      matches[i] = min_index + 1; // because c++ is 0 indexed and r is 1 indexed
    } else {
      matches[i] = -1;
    }
  }

  return matches;
}

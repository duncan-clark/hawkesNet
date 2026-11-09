#ifndef HAWKESNET_CRAN_RELEASE_VALIDATION_H
#define HAWKESNET_CRAN_RELEASE_VALIDATION_H

#include <Rcpp.h>
#include <cmath>
#include <climits>
#include <vector>

namespace cran_input {

inline void release_parameters(const Rcpp::NumericVector& p) {
  if (p.size() != 15) Rcpp::stop("Release evaluation requires fifteen parameters");
  for (double value : p)
    if (!R_finite(value)) Rcpp::stop("Release parameters must be finite");
  double total = std::exp(p[0]) + std::exp(p[1]) + std::exp(p[2]);
  if (!R_finite(total) || total <= 0)
    Rcpp::stop("Release proposal means must have a positive finite sum");
}

inline bool integer(double value) {
  return R_finite(value) && value == std::floor(value);
}

inline int index(double value, int size, const char* name) {
  if (!integer(value) || value < 0 || value >= size)
    Rcpp::stop("Invalid %s index", name);
  return static_cast<int>(value);
}

inline void offsets(const Rcpp::NumericVector& p, int events, int rows, const char* name) {
  if (p.size() != static_cast<R_xlen_t>(events) + 1 || p[0] != 0)
    Rcpp::stop("Invalid %s offsets", name);
  for (R_xlen_t i = 0; i < p.size(); ++i)
    if (!integer(p[i]) || p[i] < 0 || p[i] > rows || (i && p[i] < p[i-1]))
      Rcpp::stop("Invalid %s offsets", name);
}

inline void release_data(const Rcpp::List& data, const Rcpp::NumericVector& parameters,
                         int families) {
  release_parameters(parameters);
  if (families < 1 || families > INT_MAX / 6)
    Rcpp::stop("Invalid release family count");
  Rcpp::NumericVector times = data["time"], newborn = data["new_count"];
  Rcpp::NumericVector cp = data["change_ptr"], sp = data["selected_ptr"];
  Rcpp::NumericMatrix initial = data["initial"], changes = data["changes"], selected = data["selected"];
  if (times.size() > INT_MAX / 15 || initial.ncol() != 6 || changes.ncol() != 6 ||
      selected.ncol() != 3 || newborn.size() != times.size())
    Rcpp::stop("Invalid release array dimensions");
  int n = times.size();
  offsets(cp, n, changes.nrow(), "release change");
  offsets(sp, n, selected.nrow(), "release selection");
  for (int i = 0; i < n; ++i) {
    if (!R_finite(times[i]) || (i && times[i] < times[i-1]))
      Rcpp::stop("Release times must be finite and ordered");
    if (!integer(newborn[i]) || newborn[i] < 0 || newborn[i] >= std::ldexp(1.0, 63))
      Rcpp::stop("Release new counts must be finite nonnegative integers");
  }
  for (const Rcpp::NumericMatrix* rows : {&initial, &changes}) {
    for (int i = 0; i < rows->nrow(); ++i) {
      index((*rows)(i, 0), families, "release family");
      for (int j = 1; j < 6; ++j)
        if (!R_finite((*rows)(i, j))) Rcpp::stop("Release covariates must be finite");
    }
  }
  for (int i = 0; i < selected.nrow(); ++i) {
    index(selected(i, 0), families, "selected family");
    index(selected(i, 1), 2, "release channel");
    if (!integer(selected(i, 2)) || selected(i, 2) <= 0 || selected(i, 2) >= std::ldexp(1.0, 63))
      Rcpp::stop("Selected release line counts must be positive integers");
  }
  // A selected family must exist before its event's changes are applied.
  std::vector<bool> seen(families, false);
  for (int i = 0; i < initial.nrow(); ++i) seen[static_cast<int>(initial(i, 0))] = true;
  if (n && initial.nrow() == 0) Rcpp::stop("Release evaluation needs a nonempty initial family set");
  for (int event = 0; event < n; ++event) {
    for (int row = sp[event]; row < sp[event+1]; ++row)
      if (!seen[static_cast<int>(selected(row, 0))])
        Rcpp::stop("Selected release family is not yet active");
    for (int row = cp[event]; row < cp[event+1]; ++row)
      seen[static_cast<int>(changes(row, 0))] = true;
  }
}
}
#endif

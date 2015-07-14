#ifndef METEOFILL_SUB
#define METEOFILL_SUB

#include <cmath>
#include <vector>
#include <string>
#include <sstream>

#include "except/except.h"
#include "geostat.h"

using namespace std;

void station_weights(
  // Inputs
  const vector<double> x,
  const vector<double> y,
  const vector<bool> is_target,
  const double idw_power,
  const unsigned int nsectors,
  const unsigned int norigins,
  // Outputs
  vector<unsigned int> &inds_target,
  vector<unsigned int> &inds_source,
  vector<double> &weights
);

bool linearModel(
  // Inputs
  const vector<double> &z,
  const vector<double> &v,
  const double v_nodata,
  // Outputs
  double &slope,              
  double &intercept,
  double &rsquared
);

#endif

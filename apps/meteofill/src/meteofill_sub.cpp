#include "meteofill_sub.h"

////////////////////////////////////////////////////////////////////////////////
void station_weights(
  // Inputs
  const vector<double> x,
  const vector<double> y,
  const vector<bool> is_target,
  const double idw_power,
  const unsigned int nsectors,
  const unsigned int norigins,
  // Outputs
  vector<unsigned int> &inds_target,  // vector of target indices
  vector<unsigned int> &inds_source,  // vector of source indices
  vector<double> &weights             // vector of weights
) {
  // Local data
  unsigned int status;
  vector<double> x_src, y_src;
  vector<unsigned int> i_tar, i_src;
  vector<double> weights_thisTarget;

  // Check inputs
  if ((x.size() == 0) | (x.size() != y.size())) {
    except e(__PRETTY_FUNCTION__,"Arrays of x and y coordinates have different or zero length.",__FILE__,__LINE__);
    throw(e);
  }
  if ((is_target.size() != x.size())) {
    except e(__PRETTY_FUNCTION__,"Logical array defining targets and sources has wrong length.",__FILE__,__LINE__);
    throw(e);
  }
  if (nsectors < 1) {
    except e(__PRETTY_FUNCTION__,"Number of sectors must be >= 1.",__FILE__,__LINE__);
    throw(e);
  }
  if (norigins < 1) {
    except e(__PRETTY_FUNCTION__,"Number of sectors origins must be >= 1.",__FILE__,__LINE__);
    throw(e);
  }

  // Separate targets and sources but remember the indices in the original
  // (input) vector
  i_tar.clear();
  x_src.clear();
  y_src.clear();
  i_src.clear();
  for (unsigned int i=0; i<is_target.size(); i++) {
    if (is_target[i]) {
      i_tar.push_back(i);
    } else {
      x_src.push_back(x[i]);
      y_src.push_back(y[i]);
      i_src.push_back(i);
    }
  }
  if (i_tar.size() == 0) {
    except e(__PRETTY_FUNCTION__,"Number of target locations is zero.",__FILE__,__LINE__);
    throw(e);
  }
  if (i_src.size() == 0) {
    except e(__PRETTY_FUNCTION__,"Number of source locations is zero.",__FILE__,__LINE__);
    throw(e);
  }

  // Init result
  inds_target.clear();
  inds_source.clear();
  weights.clear();

  // Loop through targets
  for (unsigned int i=0; i<i_tar.size(); i++) {

    // Compute weights for all sources and current target
    // Note: weights_thisTarget.size() == x_src.size() but
    //       weights_thisTarget.size() < x.size()
    status= idweights(x[i_tar[i]], y[i_tar[i]], x_src, y_src, nsectors, norigins,
      idw_power, weights_thisTarget);
    if (status != 0) {
      except e(__PRETTY_FUNCTION__,"Failed to compute inverse distance weights.",__FILE__,__LINE__);
      throw(e);
    }

    // Add non-zero weights to result
    // Note: Indices need to be transformed back into the indices of the
    //       original (input) array
    for (unsigned int k=0; k<weights_thisTarget.size(); k++) {
      if (weights_thisTarget[k] > 0.) {
        inds_target.push_back(i_tar[i]);
        inds_source.push_back(i_src[k]);
        weights.push_back(weights_thisTarget[k]);
      }
    }

  }
}

////////////////////////////////////////////////////////////////////////////////
// Estimation of coefficients and goodness of a linear model

// Mean value of a vector
double avg(const vector<double> &x) {
  if (x.size() == 0) {
    except e(__PRETTY_FUNCTION__,"Input arrays has zero length.",__FILE__,__LINE__);
    throw(e);
  }
  double sum;
  sum= 0;
  for (unsigned int i=0; i<x.size(); i++) {
    sum= sum + x[i];
  }
  return(sum / x.size());
}

// Variance of a vector, using 1/n instead of 1/(n-1)
double var(const vector<double> &x) {
  if (x.size() == 0) {
    except e(__PRETTY_FUNCTION__,"Input arrays has zero length.",__FILE__,__LINE__);
    throw(e);
  }
  double mean;
  mean= avg(x);
  double sum;
  sum= 0;
  for (unsigned int i=0; i<x.size(); i++) {
    sum= sum + pow((x[i] - mean), 2.);
  }
  return(sum / x.size());
}

// Covariance of two vectors, using 1/n instead of 1/(n-1)
double cov(const vector<double> &x, const vector<double> &y) {
  if ((x.size() == 0) | (x.size() != y.size())) {
    except e(__PRETTY_FUNCTION__,"Input arrays have different or zero length.",__FILE__,__LINE__);
    throw(e);
  }
  double xmean;
  xmean= avg(x);
  double ymean;
  ymean= avg(y);
  double sum;
  sum= 0;
  for (unsigned int i=0; i<x.size(); i++) {
    sum= sum + (x[i] - xmean) * (y[i] - ymean);
  }
  return(sum / x.size());
}

// Model coefficients and quality
// Returns: A boolean status. FALSE indicates a singular case in model fitting.
bool linearModel(
  const vector<double> &z,     // z-values (predictor)
  const vector<double> &v,     // observed data at stations
  const double v_nodata,       // value to mark missing observations
  double &slope,              
  double &intercept,
  double &rsquared
) {
  // Constants
  const double ZERO= 0.;
  // Check input 
  if ((z.size() == 0) | (z.size() != v.size())) {
    except e(__PRETTY_FUNCTION__,"Input arrays have different or zero length.",__FILE__,__LINE__);
    throw(e);
  }
  // Filter for valid data
  vector<double> z_use(0), v_use(0);
  for (unsigned int i=0; i<z.size(); i++) {
    if (v[i] != v_nodata) {
      z_use.push_back(z[i]);
      v_use.push_back(v[i]);
    }
  }

  if (v_use.size() == 0) {        // No valid records
    except e(__PRETTY_FUNCTION__,"No valid records in input.",__FILE__,__LINE__);
    throw(e);
  } else if (v_use.size() == 1) { // Just a single valid records
    slope= ZERO;
    intercept= avg(v_use);
    rsquared= ZERO;
    return(false);  // RESULT SHOULD NOT BE USED
  } else {                        // Two or more useable records
    double var_z= var(z_use);
    double var_v= var(v_use);
    if (var_z == ZERO) {          // var_z==0: All z-values equal --> Infinite slope
      slope= ZERO;
      intercept= avg(v_use);
      rsquared= ZERO;
      return(false);  // RESULT SHOULD NOT BE USED
    } else {
      if (var_v == ZERO) {        // var_v==0: All v-values equal --> Zero slope, zero correlation
        slope= ZERO;
        intercept= avg(v_use);
        rsquared= ZERO;
        return(true);  // RESULT OK, ALTHOUGH INPUT IS STRANGE
      } else {                    // The normal well-behaved case
        slope= cov(z_use,v_use) / var_z;
        intercept= avg(v_use) - slope * avg(z_use);
        rsquared= pow(cov(z_use,v_use), 2.) / (var_z * var_v);
        return(true);  // RESULT OK!
      }
    }
  }
}

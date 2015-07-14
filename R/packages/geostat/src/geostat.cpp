
#include "geostat.h"

////////////////////////////////////////////////////////////////////////////////
double dist(
  const double x0,
  const double y0,
  const double x1,
  const double y1
){
  return(sqrt(pow(x0-x1, 2.) + pow(y0-y1, 2.)));
}

////////////////////////////////////////////////////////////////////////////////
double northangle(
  const double x_from,
  const double y_from,
  const double x_to,
  const double y_to
) {
  // Local data
  const double PI= 3.14159265;
  const double ZERO= 0.0;
  double angle;
  double dx,dy;
  // Get distances
  dx= x_to - x_from;
  dy= y_to - y_from;
  if (dy >= ZERO) {
    if (dx > ZERO) {
      angle= 90.0 - atan(dy/dx) * 180.0 / PI;
    } else if (dx < ZERO) {
      angle= 270.0 + atan(dy/abs(dx)) * 180.0 / PI;
    } else {
      angle= 0.0;
    }
  } else {
    if (dx > ZERO) {
      angle= 90.0 + atan(abs(dy)/dx) * 180.0 / PI;
    } else if (dx < ZERO) {
      angle= 270.0 - atan(abs(dy)/abs(dx)) * 180.0 / PI;
    } else {
      angle= 180.0;
    }
  }
  return(angle);
}

////////////////////////////////////////////////////////////////////////////////
unsigned int sector(
  const unsigned int nsectors,
  const double angle_of_origin,
  const double x_from,
  const double y_from,
  const double x_to,
  const double y_to
) {
  // Local data
  double angle, orig_angle;
  unsigned int nsect;
  unsigned int tmp;
  // Correct bad number of sectors
  nsect= max(static_cast<unsigned int>(1),nsectors);
  // Correct origin angles > 360°
  orig_angle= angle_of_origin - floor(angle_of_origin/360.)*360.;
  // Compute north angle and correct for angle of origin
  angle= northangle(x_from,y_from,x_to,y_to);
  angle= angle - orig_angle;
  if (angle < 0.) angle= 360. + angle;
  // Compute sector index
  tmp= floor(angle / 360. * nsect) + 1;
  if (tmp > nsect) tmp--; // When north angle is exactly 360°
  return(tmp);
}

////////////////////////////////////////////////////////////////////////////////
unsigned int neighbors(
  const double x_tar,
  const double y_tar,
  const vector<double> x_src,
  const vector<double> y_src,
  const unsigned int nsectors,
  const unsigned int norigins,
  vector<bool> &is_neighbor,
  vector<double> &distance
) {
  // Local data
  unsigned int invalid_index;
  double orig_angle;
  vector<double> all_dists;
  vector<double> all_sects;
  vector<unsigned int> nearest_in_sector;
  vector<unsigned int> selected_sources;
  double sumOfDists_current= 0.;  // not required but avoids compiler warning
  double sumOfDists_new;
  unsigned int n_occupied_sectors;
  bool save;
  
  // Clear result arrays
  is_neighbor.clear();
  distance.clear();

  // Check arguments and return with a non-zero status in case of errors
  if ((x_src.size()==0) | (y_src.size()==0)) return(1);
  if (x_src.size() != y_src.size()) return(2);
  if (nsectors<1) return(3);
  if (norigins<1) return(4);

  // Compute distances between target and sources and get maximum
  all_dists.resize(x_src.size());
  for (unsigned int isrc=0; isrc<x_src.size(); isrc++) {
    all_dists[isrc]= dist(x_tar,y_tar,x_src[isrc],y_src[isrc]);
  }

  // Allocate array to hold the nearest sources in each sector
  nearest_in_sector.resize(nsectors);

  // Define an invalid source index (used to detect unoccupied sectors)
  invalid_index= x_src.max_size();

  // Loop through the angles of origins (sector rotations)
  all_sects.resize(x_src.size());
  for (unsigned int iorig=1; iorig<=norigins; iorig++) {
    orig_angle= 360. / nsectors * (iorig-1) / norigins;
    // Compute sectors for all sources
    for (unsigned int isrc=0; isrc<x_src.size(); isrc++) {
      all_sects[isrc]= sector(nsectors,orig_angle,x_tar,y_tar,
        x_src[isrc],y_src[isrc]);
    }
    // Initialize index array of neares source in sector to an invalid index
    for (unsigned int isect=1; isect<=nsectors; isect++) {
      nearest_in_sector[isect-1]= invalid_index;
    }
    // For each sector, determine the index of the nearest source locations.
    // If a sector is empty, the index is still the invalid initial value.
    sumOfDists_new= 0.; 
    for (unsigned int isect=1; isect<=nsectors; isect++) {
      for (unsigned int isrc=0; isrc<x_src.size(); isrc++) {
        if (all_sects[isrc] == isect) {
          // Save first source in sector
          if (nearest_in_sector[isect-1] == invalid_index) {
            nearest_in_sector[isect-1]= isrc;
          // Replace by nearer source
          } else {
            if (all_dists[isrc] < all_dists[nearest_in_sector[isect-1]]) {
              nearest_in_sector[isect-1]= isrc;
            }
          }
        }
      }
      // Update sum of distances
      if (nearest_in_sector[isect-1] != invalid_index) {
        sumOfDists_new= sumOfDists_new + all_dists[nearest_in_sector[isect-1]]; 
      }
    }
    // Save result for first angle of origin
    // Overwrite, if results for another angle of origin is better
     // Overwrite if there are more occupied sectors
    if (iorig==1) {
      save= true;
    } else {
      n_occupied_sectors= 0;
      for (unsigned int isect=1; isect<=nsectors; isect++) {
        if (nearest_in_sector[isect-1] != invalid_index) n_occupied_sectors++;
      }
      if (n_occupied_sectors > selected_sources.size()) {
        save= true;
      } else if (n_occupied_sectors == selected_sources.size()) {
        if (sumOfDists_new < sumOfDists_current) {
          save= true;
        } else {
          save= false; // Same number of occupied sectors but greater dists
        }
      } else {
        save= false;   // Less occupied sectors
      }
    }
    // Save (there is at least 1 non-empty sector if the number of sources is > 0)
    if (save) {
      selected_sources.clear();
      sumOfDists_current= 0.;
      for (unsigned int isect=1; isect<=nsectors; isect++) {
        if (nearest_in_sector[isect-1] != invalid_index) {
          selected_sources.push_back(nearest_in_sector[isect-1]);
          // Save sum of distances
          sumOfDists_current= sumOfDists_current + all_dists[nearest_in_sector[isect-1]]; 
        }
      }
    }

  } // End of loop over angles of origin

  // Set result arrays and return
  is_neighbor.resize(x_src.size());
  distance.resize(x_src.size());
  for (unsigned int isrc=0; isrc<x_src.size(); isrc++) {
    is_neighbor[isrc]= false;
    distance[isrc]= all_dists[isrc];
  }
  for (unsigned int k=0; k<selected_sources.size(); k++) {
    is_neighbor[selected_sources[k]]= true;
  }

  return(0);
}

////////////////////////////////////////////////////////////////////////////////
unsigned int idweights(
  const double x_tar,
  const double y_tar,
  const vector<double> x_src,
  const vector<double> y_src,
  const unsigned int nsectors,
  const unsigned int norigins,
  const double power,
  vector<double> &weights
) {
  // Local data
  const double NEAR_ZERO= 1.0e-12;
  unsigned int status;
  vector<bool> is_neighbor;
  vector<double> distance;
  double sum_of_invdist;

  // Clear result array
  weights.clear();

  // Find neighbors and return on error
  status= neighbors(x_tar, y_tar, x_src, y_src, nsectors, norigins,
    is_neighbor, distance);
  if (status != 0) {
    return(status);
  } else {
    weights.resize(x_src.size());
    // Compute inverse distances raised to a power
    sum_of_invdist= 0.;
    for (unsigned int i=0; i<x_src.size(); i++) {
      if (is_neighbor[i]) {
        weights[i]= 1. / pow(max(distance[i],NEAR_ZERO), power);
        sum_of_invdist= sum_of_invdist + weights[i];
      } else {
        weights[i]= 0.;
      }
    }
    // Convert to final weights
    for (unsigned int i=0; i<x_src.size(); i++) {
      weights[i]= weights[i] / sum_of_invdist;
    }
  }

  return(0);
}




#include "geostat.h"

extern "C" {

  // The interfaces of the cpp-functions to be called from R go here

  //////////////////////////////////////////////////////////////////////////////
  void gs_dist(
    // In
    const double* x0,
    const double* y0,
    const double* x1,
    const double* y1,
    // Out
    double* result
  ) {
    *result= dist(*x0, *y0, *x1, *y1);
  }

  //////////////////////////////////////////////////////////////////////////////
  void gs_northangle(
    // In
    const double* x_from,
    const double* y_from,
    const double* x_to,
    const double* y_to,
    // Out
    double* result
  ) {
    *result= northangle(*x_from, *y_from, *x_to, *y_to);
  }

  //////////////////////////////////////////////////////////////////////////////
  void gs_sector(
    // In
    const int* nsectors,
    const double* angle_of_origin,
    const double* x_from,
    const double* y_from,
    const double* x_to,
    const double* y_to,
    // Out
    int* result
  ) {
    *result= static_cast<int>(sector(static_cast<unsigned int>(*nsectors),
       *angle_of_origin, *x_from, *y_from, *x_to, *y_to));
  }

  //////////////////////////////////////////////////////////////////////////////
  void gs_neighbors(
    // In
    const double* x_tar,
    const double* y_tar,
    const int* nsources,
    const double* x_src,
    const double* y_src,
    const int* nsectors,
    const int* norigins,
    // Out
    int* source_is_selected,
    double* source_distance,
    int* result
  ) {
    // Local data
    vector<double> vect_x_src;
    vector<double> vect_y_src;
    vector<bool> is_neighbor;
    vector<double> distance;

    // Convert input arrays to vectors
    vect_x_src.resize(*nsources);
    vect_y_src.resize(*nsources);
    for (unsigned int i=0; i<static_cast<unsigned int>(*nsources); i++) {
      vect_x_src[i]= x_src[i];
      vect_y_src[i]= y_src[i];
    }

    // Function call
    *result= static_cast<int>(neighbors(*x_tar, *y_tar, vect_x_src, vect_y_src,
       static_cast<unsigned int>(*nsectors),
       static_cast<unsigned int>(*norigins),
       is_neighbor,distance));

    // Set result arrays  (on error, i.e. result != 0, nothing is selected)
    if (*result == 0) {
      for (unsigned int i=0; i<static_cast<unsigned int>(*nsources); i++) {
        if (is_neighbor[i]) {
          source_is_selected[i]= 1;
        } else {
          source_is_selected[i]= 0;
        }
        source_distance[i]= distance[i]; 
      }
    } else {
      for (unsigned int i=0; i<static_cast<unsigned int>(*nsources); i++) {
        source_is_selected[i]= 0;
        source_distance[i]= -9999.;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  void gs_idweights(
    // In
    const double* x_tar,
    const double* y_tar,
    const int* nsources,
    const double* x_src,
    const double* y_src,
    const int* nsectors,
    const int* norigins,
    const double* power,
    // Out
    double* weights,
    int* result
  ) {
    // Local data
    vector<double> vect_x_src;
    vector<double> vect_y_src;
    vector<double> vect_weights;

    // Convert input arrays to vectors
    vect_x_src.resize(*nsources);
    vect_y_src.resize(*nsources);
    for (unsigned int i=0; i<static_cast<unsigned int>(*nsources); i++) {
      vect_x_src[i]= x_src[i];
      vect_y_src[i]= y_src[i];
    }

    // Function call
    *result= static_cast<int>(idweights(*x_tar, *y_tar, vect_x_src, vect_y_src,
       static_cast<unsigned int>(*nsectors),
       static_cast<unsigned int>(*norigins),
       *power,vect_weights));

    // Set result arrays  (on error, i.e. result != 0, all weights are zero)
    if (*result == 0) {
      for (unsigned int i=0; i<static_cast<unsigned int>(*nsources); i++) {
        weights[i]= vect_weights[i]; 
      }
    } else {
      for (unsigned int i=0; i<static_cast<unsigned int>(*nsources); i++) {
        weights[i]= 0.;
      }
    }
  }


} // extern "C"



#ifndef GEOSTAT_H
#define GEOSTAT_H

#include <cmath>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Horizontal distance between two points [x0,y0] and [x1,y1].
////////////////////////////////////////////////////////////////////////////////

double dist(
  const double x0,
  const double y0,
  const double x1,
  const double y1
);

////////////////////////////////////////////////////////////////////////////////
// Angle when looking from [x_from,y_from] to [x_to,y_to]. The anngle increases
// in clockwise direction with N=0° or 360°, E=90°, S=180°, W=270".
////////////////////////////////////////////////////////////////////////////////
double northangle(
  const double x_from,
  const double y_from,
  const double x_to,
  const double y_to
);

////////////////////////////////////////////////////////////////////////////////
// Sector where a point [x_to,y_to] is located when looking from point
// [x_from,y_from]. The smallest sector index is 1 and indices increase in
// clockwise direction. The number of sectors is a parameter (a value of 0 will
// silently be changed into 1). Another parameter is the angle of origin of the
// sectors. It is a value in degrees and marks the border between sector 1 and
// the sector with the highest index. It is 0° or 360° == N, 90° == E,
// 180° == S, 270" == W. Values greater than 360° are silently corrected by
// subtracting x*360 where x is chosen to yield a result in range 0...360.
////////////////////////////////////////////////////////////////////////////////
unsigned int sector(
  const unsigned int nsectors,
  const double angle_of_origin,
  const double x_from,
  const double y_from,
  const double x_to,
  const double y_to
);

////////////////////////////////////////////////////////////////////////////////
// Select source locations for a target location (neighbors for interpolation).
// A single location is selected from each sector. By trying different angles of
// origins of the sectors, one may achieve better results (more occupied sectors
// with smaller average distance between target and sources).
//
// INPUT
//
//   x_tar, y_tar        Defines the target location
//   x_src, y_src        Vectors, defining the source locations
//   nsectors            Number of sectors to split the search space into
//   norigins            Number of sector origins (rotations) to try
//
// OUTPUT
//
//   is_neighbor         Logical array (mask) with length equal to the number of
//                       source locations. TRUE for selected neighbors.
//   distance            Array of distances between target and all source locs.
//
// RETURNS
//
//   An integer error level.
//   0: Normal exit.
//   1: Vectors x_src or y_src are of zero length.
//   2: Vectors x_src and y_src are of different length.
//   3: Value of nsectors is < 1.
//   4: Value of norigins is < 1.
//
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
);

////////////////////////////////////////////////////////////////////////////////
// Compute inverse distance weights after selecting source locations.
//
// INPUT
//
//   power               The power to be applied to inverse distances.
//   other arguments     See function neighbors().
//
// OUTPUT
//
//   weight              Array of weights with length equal to the number of
//                       source locations. Weights are non-zero only for the
//                       selected source locations.
//
// RETURNS
//
//   An integer error level. See function neighbors().
//
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
);

#endif


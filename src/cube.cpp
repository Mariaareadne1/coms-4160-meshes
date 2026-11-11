#include "cube.h"

void cube(
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  // 8 corner positions of the unit cube, one row per corner (x, y, z)
  V.resize(8,3);
  V <<
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 1, 0,
    0, 0, 1,
    1, 0, 1,
    1, 1, 1,
    0, 1, 1;

  // 14 texture positions on the cross-shaped unfolding, one row per point (u, v)
  UV.resize(14,2);
  UV <<
    0.00, 0.25,
    0.25, 0.25,
    0.50, 0.25,
    0.75, 0.25,
    1.00, 0.25,
    0.00, 0.50,
    0.25, 0.50,
    0.50, 0.50,
    0.75, 0.50,
    1.00, 0.50,
    0.25, 0.75,
    0.50, 0.75,
    0.25, 0.00,
    0.50, 0.00;

  // 6 outward unit normals, one row per face, in the same order as the rows of F
  NV.resize(6,3);
  NV <<
     0,  0, -1,
    -1,  0,  0,
     0,  0,  1,
     1,  0,  0,
     0,  1,  0,
     0, -1,  0;

  // 6 quads as row numbers of V, each listed counterclockwise when seen from outside
  F.resize(6,4);
  F <<
    1, 0, 3, 2,
    0, 4, 7, 3,
    4, 5, 6, 7,
    5, 1, 2, 6,
    3, 7, 6, 2,
    1, 5, 4, 0;

  // The same 6 quads as row numbers of UV, corner by corner in the same order as F
  UF.resize(6,4);
  UF <<
     0,  1,  6,  5,
     1,  2,  7,  6,
     2,  3,  8,  7,
     3,  4,  9,  8,
     6,  7, 11, 10,
    12, 13,  2,  1;

  // Every corner of face i uses normal row i
  NF.resize(6,4);
  NF <<
    0, 0, 0, 0,
    1, 1, 1, 1,
    2, 2, 2, 2,
    3, 3, 3, 3,
    4, 4, 4, 4,
    5, 5, 5, 5;
  ////////////////////////////////////////////////////////////////////////////

  ////Hint:
  // V.resize(8,3);
  // F.resize(6,4);
  // UV.resize(14,2);
  // UF.resize(6,4);
  // NV.resize(6,3);
  // NF.resize(6,4);
}
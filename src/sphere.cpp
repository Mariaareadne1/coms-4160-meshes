#include "sphere.h"
#include <iostream>
#include <cmath>

void sphere(
  const int num_faces_u,
  const int num_faces_v,
  Eigen::MatrixXd & V,
  Eigen::MatrixXi & F,
  Eigen::MatrixXd & UV,
  Eigen::MatrixXi & UF,
  Eigen::MatrixXd & NV,
  Eigen::MatrixXi & NF)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  // The number pi, used to turn grid steps into angles
  const double pi = 3.14159265358979323846;

  // One ring of num_faces_u positions on each inner latitude line, plus the two poles
  const int num_positions = num_faces_u * (num_faces_v - 1) + 2;
  // Row of V that holds the south pole
  const int south_pole = 0;
  // Row of V that holds the north pole
  const int north_pole = num_positions - 1;

  // Make room for all 3D positions
  V.resize(num_positions,3);
  // South pole sits at the bottom of the unit sphere
  V(south_pole,0) = 0;
  V(south_pole,1) = -1;
  V(south_pole,2) = 0;
  // North pole sits at the top of the unit sphere
  V(north_pole,0) = 0;
  V(north_pole,1) = 1;
  V(north_pole,2) = 0;
  // Walk the inner latitude lines from just above the south pole to just below the north pole
  for(int j = 1; j < num_faces_v; j++)
  {
    // Latitude angle: -pi/2 at the south pole, +pi/2 at the north pole
    const double latitude = -pi/2 + pi * j / num_faces_v;
    // Walk around the ring once (no repeated point at the seam)
    for(int i = 0; i < num_faces_u; i++)
    {
      // Longitude angle: 0 up to just under 2*pi
      const double longitude = 2 * pi * i / num_faces_u;
      // Row of V for grid point (i,j)
      const int row = 1 + (j - 1) * num_faces_u + i;
      // Point on the unit sphere at this latitude and longitude (y points to the north pole)
      V(row,0) = cos(latitude) * cos(longitude);
      V(row,1) = sin(latitude);
      V(row,2) = -cos(latitude) * sin(longitude);
    }
  }

  // On a unit sphere centered at the origin, the outward unit normal equals the position
  NV = V;

  // Texture grid keeps the seam column twice (u=0 and u=1) and a full row at each pole
  UV.resize((num_faces_u + 1) * (num_faces_v + 1),2);
  // Walk every texture grid row, bottom (v=0) to top (v=1)
  for(int j = 0; j <= num_faces_v; j++)
  {
    // Walk every texture grid column, left (u=0) to right (u=1)
    for(int i = 0; i <= num_faces_u; i++)
    {
      // Row of UV for texture grid point (i,j)
      const int row = j * (num_faces_u + 1) + i;
      // u goes around with longitude
      UV(row,0) = (double)i / num_faces_u;
      // v goes up with latitude
      UV(row,1) = (double)j / num_faces_v;
    }
  }

  // One quad per grid cell, for positions, texture positions and normals
  F.resize(num_faces_u * num_faces_v,4);
  UF.resize(num_faces_u * num_faces_v,4);
  NF.resize(num_faces_u * num_faces_v,4);
  // Column offsets of the 4 corners of a cell: bottom-left, bottom-right, top-right, top-left
  const int corner_i[4] = {0, 1, 1, 0};
  // Row offsets of the same 4 corners (counterclockwise seen from outside)
  const int corner_j[4] = {0, 0, 1, 1};
  // Walk every cell row
  for(int j = 0; j < num_faces_v; j++)
  {
    // Walk every cell column
    for(int i = 0; i < num_faces_u; i++)
    {
      // Row of F, UF and NF for cell (i,j)
      const int face = j * num_faces_u + i;
      // Fill the 4 corners of this cell
      for(int k = 0; k < 4; k++)
      {
        // Texture grid column of this corner (may equal num_faces_u, the right copy of the seam)
        const int grid_i = i + corner_i[k];
        // Grid row of this corner
        const int grid_j = j + corner_j[k];
        // Ring position of this corner: column num_faces_u wraps back to column 0
        const int ring_i = grid_i % num_faces_u;
        // Row of V for this corner
        int position = 0;
        if(grid_j == 0)
        {
          // Bottom grid row is the single south pole point
          position = south_pole;
        }
        else if(grid_j == num_faces_v)
        {
          // Top grid row is the single north pole point
          position = north_pole;
        }
        else
        {
          // Inner rows index into their ring
          position = 1 + (grid_j - 1) * num_faces_u + ring_i;
        }
        // Corner position index
        F(face,k) = position;
        // Normal rows match position rows one to one
        NF(face,k) = position;
        // Texture index uses the unwrapped column, so the seam gets u=1 on its right side
        UF(face,k) = grid_j * (num_faces_u + 1) + grid_i;
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
}
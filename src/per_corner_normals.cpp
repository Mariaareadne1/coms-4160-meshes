#include "per_corner_normals.h"
#include "triangle_area_normal.h"
// Hint:
#include "vertex_triangle_adjacency.h"
#include <iostream>
#include <cmath>

void per_corner_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double corner_threshold,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(F.rows()*3,3);
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  // The number pi, used to convert degrees to radians
  const double pi = 3.14159265358979323846;
  // Smallest allowed dot product between two unit normals: cos of the threshold angle
  const double epsilon = cos(corner_threshold * pi / 180.0);

  // Area normal (direction times area) of every triangle
  Eigen::MatrixXd area_normals(F.rows(),3);
  // Unit normal (direction only) of every triangle
  Eigen::MatrixXd unit_normals(F.rows(),3);
  // Fill both for each triangle
  for(int f = 0; f < F.rows(); f++)
  {
    // Positions of the triangle's 3 corners
    const Eigen::RowVector3d a = V.row(F(f,0));
    const Eigen::RowVector3d b = V.row(F(f,1));
    const Eigen::RowVector3d c = V.row(F(f,2));
    // Area-weighted normal
    area_normals.row(f) = triangle_area_normal(a, b, c);
    // Same direction, length 1
    unit_normals.row(f) = area_normals.row(f).normalized();
  }

  // For every vertex, the list of triangles that touch it
  std::vector<std::vector<int> > VF;
  vertex_triangle_adjacency(F, V.rows(), VF);

  // Visit every triangle
  for(int f = 0; f < F.rows(); f++)
  {
    // Visit each of its 3 corners
    for(int k = 0; k < 3; k++)
    {
      // Vertex at this corner
      const int v = F(f,k);
      // Running sum of area normals of neighbors on the same side of any crease
      Eigen::RowVector3d sum(0,0,0);
      // Visit every triangle g that touches vertex v (this includes f itself)
      for(int i = 0; i < (int)VF[v].size(); i++)
      {
        // Neighbor triangle
        const int g = VF[v][i];
        // Keep g only if its normal is within the threshold angle of f's normal
        if(unit_normals.row(g).dot(unit_normals.row(f)) > epsilon)
        {
          sum += area_normals.row(g);
        }
      }
      // Corner F(f,k) is stored at row f*3+k, scaled to length 1
      N.row(f*3 + k) = sum.normalized();
    }
  }
  ////////////////////////////////////////////////////////////////////////////
}
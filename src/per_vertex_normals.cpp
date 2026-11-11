#include "per_vertex_normals.h"
#include "triangle_area_normal.h"

void per_vertex_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  N = Eigen::MatrixXd::Zero(V.rows(),3);
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  // Visit every triangle
  for(int f = 0; f < F.rows(); f++)
  {
    // Positions of the triangle's 3 corners
    const Eigen::RowVector3d a = V.row(F(f,0));
    const Eigen::RowVector3d b = V.row(F(f,1));
    const Eigen::RowVector3d c = V.row(F(f,2));
    // Area normal: unit normal times area, which is exactly the weighted term a_f * n_f
    const Eigen::RowVector3d area_normal = triangle_area_normal(a, b, c);
    // Add the weighted term to the running sum of each of the triangle's 3 vertices
    for(int k = 0; k < 3; k++)
    {
      N.row(F(f,k)) += area_normal;
    }
  }
  // Divide each sum by its length so every vertex normal has length 1
  for(int v = 0; v < N.rows(); v++)
  {
    N.row(v) = N.row(v).normalized();
  }
  ////////////////////////////////////////////////////////////////////////////
}
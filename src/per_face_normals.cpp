#include "per_face_normals.h"
#include "triangle_area_normal.h"

void per_face_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  ////////////////////////////////////////////////////////////////////////////
  // Replace with your code:
  // One normal row per triangle
  N.resize(F.rows(),3);
  // Visit every triangle
  for(int f = 0; f < F.rows(); f++)
  {
    // Positions of the triangle's 3 corners, looked up through the indices in row f of F
    const Eigen::RowVector3d a = V.row(F(f,0));
    const Eigen::RowVector3d b = V.row(F(f,1));
    const Eigen::RowVector3d c = V.row(F(f,2));
    // Area-weighted normal, scaled to length 1 so only the direction remains
    N.row(f) = triangle_area_normal(a, b, c).normalized();
  }
  ////////////////////////////////////////////////////////////////////////////
}
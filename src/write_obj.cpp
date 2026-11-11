#include "write_obj.h"
#include <fstream>
#include <cassert>
#include <iostream>

bool write_obj(
  const std::string & filename,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & UV,
  const Eigen::MatrixXi & UF,
  const Eigen::MatrixXd & NV,
  const Eigen::MatrixXi & NF)
{
  assert((F.size() == 0 || F.cols() == 3 || F.cols() == 4) && "F must have 3 or 4 columns");
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here:
  // Open the output file for writing (creates it, or overwrites an existing one)
  std::ofstream file(filename);
  // Stop with failure if the file could not be opened
  if(!file.is_open())
  {
    return false;
  }
  // One "v x y z" line per row of V (3D positions)
  for(int i = 0; i < V.rows(); i++)
  {
    file << "v " << V(i,0) << " " << V(i,1) << " " << V(i,2) << "\n";
  }
  // One "vt u v" line per row of UV (2D texture coordinates)
  for(int i = 0; i < UV.rows(); i++)
  {
    file << "vt " << UV(i,0) << " " << UV(i,1) << "\n";
  }
  // One "vn x y z" line per row of NV (3D normal vectors)
  for(int i = 0; i < NV.rows(); i++)
  {
    file << "vn " << NV(i,0) << " " << NV(i,1) << " " << NV(i,2) << "\n";
  }
  // One "f ..." line per face, with one "v/vt/vn" group per corner
  for(int i = 0; i < F.rows(); i++)
  {
    // Start the face line
    file << "f";
    // Visit each corner of face i (3 for triangles, 4 for quads)
    for(int j = 0; j < F.cols(); j++)
    {
      // Write position, UV and normal indices, shifted by 1 because .obj counts from 1
      file << " " << F(i,j) + 1 << "/" << UF(i,j) + 1 << "/" << NF(i,j) + 1;
    }
    // End the face line
    file << "\n";
  }
  // Flush everything to disk and release the file
  file.close();
  // Report success
  return true;
  ////////////////////////////////////////////////////////////////////////////
}
#include "catmull_clark.h"
#include <unordered_map>
#include <utility>
#include <functional>
#include <vector>

void catmull_clark(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int num_iters,
  Eigen::MatrixXd & SV,
  Eigen::MatrixXi & SF)
{
  ////////////////////////////////////////////////////////////////////////////
  // Replace with your code here:
  SV = V;
  SF = F;
  // Each pass replaces (SV,SF) with its subdivided version
  for(int iter = 0; iter < num_iters; iter++)
  {
    // Sizes of the mesh being subdivided in this pass
    const int num_v = SV.rows();
    const int num_f = SF.rows();

    // Face points: the average of each quad's 4 corners
    Eigen::MatrixXd face_points(num_f,3);
    for(int f = 0; f < num_f; f++)
    {
      face_points.row(f) = (SV.row(SF(f,0)) + SV.row(SF(f,1)) + SV.row(SF(f,2)) + SV.row(SF(f,3))) / 4.0;
    }

    // Lookup from an edge's key number to its edge id (0, 1, 2, ...)
    std::unordered_map<long long,int> edge_index;
    // First endpoint of each edge, by edge id
    std::vector<int> edge_a;
    // Second endpoint of each edge, by edge id
    std::vector<int> edge_b;
    // First quad that uses each edge, by edge id
    std::vector<int> edge_face1;
    // Second quad that uses each edge, by edge id (-1 until found)
    std::vector<int> edge_face2;
    // Edge id of side k of quad f, where side k runs from corner k to corner k+1
    Eigen::MatrixXi face_edges(num_f,4);
    // Visit every side of every quad
    for(int f = 0; f < num_f; f++)
    {
      for(int k = 0; k < 4; k++)
      {
        // The side's two endpoints, going around the quad
        const int a = SF(f,k);
        const int b = SF(f,(k+1)%4);
        // Smaller and larger endpoint, so both quads sharing the edge build the same key
        int small = a;
        int large = b;
        if(b < a)
        {
          small = b;
          large = a;
        }
        // One whole number that is different for every pair (small, large)
        const long long key = (long long)small * num_v + large;
        // Edge id for this side
        int id = 0;
        if(edge_index.count(key) == 0)
        {
          // First time this edge is seen: give it the next id and record its data
          id = edge_a.size();
          edge_index[key] = id;
          edge_a.push_back(a);
          edge_b.push_back(b);
          edge_face1.push_back(f);
          edge_face2.push_back(-1);
        }
        else
        {
          // Second time: this quad is the edge's other neighbor
          id = edge_index[key];
          edge_face2[id] = f;
        }
        // Remember which edge this side is
        face_edges(f,k) = id;
      }
    }
    // Number of distinct edges
    const int num_e = edge_a.size();

    // Edge points: average of the 2 endpoints and the 2 neighboring face points
    Eigen::MatrixXd edge_points(num_e,3);
    for(int e = 0; e < num_e; e++)
    {
      if(edge_face2[e] == -1)
      {
        // Edge with only one quad (open boundary): use the plain midpoint
        edge_points.row(e) = (SV.row(edge_a[e]) + SV.row(edge_b[e])) / 2.0;
      }
      else
      {
        edge_points.row(e) = (SV.row(edge_a[e]) + SV.row(edge_b[e])
          + face_points.row(edge_face1[e]) + face_points.row(edge_face2[e])) / 4.0;
      }
    }

    // Per vertex: sum of neighboring face points, and how many there are
    Eigen::MatrixXd face_point_sum = Eigen::MatrixXd::Zero(num_v,3);
    std::vector<int> face_count(num_v,0);
    for(int f = 0; f < num_f; f++)
    {
      for(int k = 0; k < 4; k++)
      {
        face_point_sum.row(SF(f,k)) += face_points.row(f);
        face_count[SF(f,k)]++;
      }
    }
    // Per vertex: sum of midpoints of its edges, and how many there are
    Eigen::MatrixXd edge_mid_sum = Eigen::MatrixXd::Zero(num_v,3);
    std::vector<int> edge_count(num_v,0);
    for(int e = 0; e < num_e; e++)
    {
      // Midpoint of edge e
      const Eigen::RowVector3d mid = (SV.row(edge_a[e]) + SV.row(edge_b[e])) / 2.0;
      // Counts toward both endpoints
      edge_mid_sum.row(edge_a[e]) += mid;
      edge_count[edge_a[e]]++;
      edge_mid_sum.row(edge_b[e]) += mid;
      edge_count[edge_b[e]]++;
    }

    // New positions: moved originals first, then face points, then edge points
    Eigen::MatrixXd new_V(num_v + num_f + num_e,3);
    for(int v = 0; v < num_v; v++)
    {
      // Number of quads around vertex v (its valence)
      const int n = face_count[v];
      if(n == 0)
      {
        // Vertex used by no quad: leave it where it is
        new_V.row(v) = SV.row(v);
      }
      else
      {
        // F: average face point around v
        const Eigen::RowVector3d avg_face = face_point_sum.row(v) / n;
        // R: average edge midpoint around v
        const Eigen::RowVector3d avg_mid = edge_mid_sum.row(v) / edge_count[v];
        // P: the old position of v
        const Eigen::RowVector3d P = SV.row(v);
        // Catmull-Clark rule: (F + 2R + (n-3)P) / n
        new_V.row(v) = (avg_face + 2.0 * avg_mid + (n - 3.0) * P) / n;
      }
    }
    // Face point of quad f lives at row num_v + f
    for(int f = 0; f < num_f; f++)
    {
      new_V.row(num_v + f) = face_points.row(f);
    }
    // Edge point of edge e lives at row num_v + num_f + e
    for(int e = 0; e < num_e; e++)
    {
      new_V.row(num_v + num_f + e) = edge_points.row(e);
    }

    // Each old quad becomes 4 new quads, one at each of its corners
    Eigen::MatrixXi new_F(num_f * 4,4);
    for(int f = 0; f < num_f; f++)
    {
      for(int k = 0; k < 4; k++)
      {
        // Side of quad f that arrives at corner k (from corner k-1)
        const int prev = (k + 3) % 4;
        // Row of new_F for the small quad at corner k of quad f
        const int row = f * 4 + k;
        // Corner k itself (its moved position)
        new_F(row,0) = SF(f,k);
        // Edge point of the side leaving corner k
        new_F(row,1) = num_v + num_f + face_edges(f,k);
        // Face point in the middle of quad f
        new_F(row,2) = num_v + f;
        // Edge point of the side arriving at corner k
        new_F(row,3) = num_v + num_f + face_edges(f,prev);
      }
    }

    // The subdivided mesh becomes the input of the next pass
    SV = new_V;
    SF = new_F;
  }
  ////////////////////////////////////////////////////////////////////////////
}
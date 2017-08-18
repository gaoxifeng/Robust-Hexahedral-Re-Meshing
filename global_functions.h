#pragma once
#include "io.h"
#include "global_types.h"
#include "Eigen/Dense"
#include <algorithm>
#include <iterator>
#include <omp.h>
#include <iostream>
#include <map>
#include <set>
#include <queue>
#include <tbb/tbb.h>
#include "igl/hausdorff.h"

using namespace Eigen;
using namespace std;

//===================================mesh connectivities===================================
void build_connectivity(Mesh &hmi);
void topology_info(Mesh &mesh, Frame &frame, Mesh_Topology & mt);
bool disk_polygon(Mesh &mesh, Frame &frame, vector<vector<uint32_t>> &fes, vector<short> &E_flag, vector<short> &V_flag, const bool &Ismesh);
bool sphere_polyhedral(Mesh &mesh, Frame &frame, vector<vector<uint32_t>> &F_nvs, vector<vector<uint32_t>> &pfs, vector<bool> &F_flag, vector<short> &E_flag, vector<short> &V_flag, const bool &Ismesh);
bool comp_topology(Mesh_Topology & mt0, Mesh_Topology & mt1);
bool redundentV_check(Mesh &meshI, Mesh &meshO);
double average_edge_length(Mesh &mesh);
void re_indexing_connectivity(Mesh &hmi, MatrixXi &H);

void extract_surface_mesh(Mesh &meshi, Mesh &mesho);
void  orient_surface_mesh(Mesh &hmi);
void  orient_triangle_mesh(Mesh &hmi);
void  orient_triangle_mesh(MatrixXd &Tri_V, MatrixXi &Tri_F);

Float	uctet(vector<Float> a, vector<Float> b, vector<Float> c, vector<Float> d);
//===================================mesh quality==========================================
bool scaled_jacobian(Mesh &hmi, Mesh_Quality &mq);
inline Float a_jacobian(Vector3d &v0, Vector3d &v1, Vector3d &v2, Vector3d &v3);
//===================================feature v tags==========================================
bool triangle_mesh_feature(Mesh_Feature &mf, Mesh &hmi);
bool initial_feature(Mesh_Feature &mf, Feature_Constraints &fc, Mesh &hmi);
bool project_surface_update_feature(Mesh_Feature &mf, Feature_Constraints &fc, MatrixXd &V, VectorXi &b, MatrixXd &bc, uint32_t Loop = 1);
bool phong_projection(vector<uint32_t> &tids, uint32_t Loop, uint32_t &tid, Vector3d &v, Vector3d &interpolP, Vector3d &interpolN, Vector3d &PreinterpolP, Vector3d &PreinterpolN);
void point_line_projection(const Vector3d &v1, const Vector3d &v2, const Vector3d &v, Vector3d &pv, double &t);
void projectPointOnQuad(const vector<Vector3d>& quad_vs, vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN);
void projectPointOnTriangle(const vector<Vector3d>& tri_vs, const vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN);
template <typename T>
T bilinear(const T& v1, const T& v2, const T& v3, const T& v4, const Vector2d& uv);
template <typename T>
T barycentric(const T& v1, const T& v2, const T& v3, const Vector2d& uv);
template <typename T>
bool num_equal(const T& v1, const T& v2, const double &precision);

Float rescale(Mesh &mesh, Float scaleI, bool inverse);
void compute_referenceMesh(MatrixXd &V, vector<Hybrid> &H, vector<uint32_t> &Hs, vector<MatrixXd> &Vout);
void hex2cuboid(MatrixXd &V, vector<uint32_t> &vs, vector<MatrixXd> &vout);
void hex2tet24(MatrixXd &V, vector<uint32_t> &vs, double & volume);
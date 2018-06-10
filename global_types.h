//    This file is part of the implementation of

//    Robust Structure Simplification for Hex Re-meshing
//    Xifeng Gao, Daniele Panozzo, Wenping Wang, Zhigang Deng, Guoning Chen
//    In ACM Transactions on Graphics (Proceedings of SIGGRAPH ASIA 2017)
// 
// Copyright (C) 2017 Xifeng Gao<gxf.xisha@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cstdlib>
#include <vector>
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

/*typedefs*/
#if defined(SINGLE_PRECISION)
typedef float Float;
#else
typedef double Float;
#endif

#define Interior_RegularE 4
#define Boundary_RegularE 2

#define Precision 1.e-7
#define Precision_Pro 1.0e-5
#define Jacobian_Bound 1.0e-4

#define PAI 3.1415926
enum Base_Set{
	SHEET=0,
	CHORD
};
typedef std::tuple<uint32_t, Base_Set, Float> Tuple_Candidate;//id, type, weight

enum Feature_V_Type {
	INTERIOR = -4,
	CORNER,
	LINE,
	REGULAR
};

const int tetra_table[4][4] =
{
	{ 0,2,1,3 },
	{ 1,0,2,3 },
	{ 2,1,0,3 },
	{ 3,0,1,2 },
};
const int tet_faces[4][3] = {
	{ 1, 0, 2 },
	{ 3, 2, 0 },
	{ 1, 2, 3 },
	{ 0, 1, 3 } 
};
const int tet_edges[6][2] = {
	{ 0, 1},
	{ 0, 2},
	{ 0, 3 },
	{ 1, 2},
	{ 1, 3 },
	{ 2, 3 }
};

const int hex_face_table[6][4] =
{
	{ 0,1,2,3 },
	{ 4,5,6,7 },
	{ 0,1,5,4 },
	{ 0,4,7,3 },
	{ 3,2,6,7 },
	{ 1,5,6,2 },
};
const int hex_tetra_table[8][4] =
{
	{ 0,3,4,1 },
	{ 1,0,5,2 },
	{ 2,1,6,3 },
	{ 3,2,7,0 },
	{ 4,7,5,0 },
	{ 5,4,6,1 },
	{ 6,5,7,2 },
	{ 7,6,4,3 },
};
//-------------------------------------------------------------------
//---For Hybrid mesh-------------------------------------------------
struct Hybrid_V
{
	uint32_t id, svid, fvid;
	vector<Float> v;
	vector<uint32_t> neighbor_vs;
	vector<uint32_t> neighbor_es;
	vector<uint32_t> neighbor_fs;
	vector<uint32_t> neighbor_hs;

	bool boundary;
};
struct Hybrid_E
{
	uint32_t id;
	vector<uint32_t> vs;
	vector<uint32_t> neighbor_fs;
	vector<uint32_t> neighbor_hs;
	
	bool boundary;

	bool hex_edge = false;
};
struct Hybrid_F
{
	uint32_t id;
	vector<uint32_t> vs;
	vector<uint32_t> es;
	vector<uint32_t> neighbor_hs;
	bool boundary;
};

struct Hybrid
{
	uint32_t id;
	vector<uint32_t> vs;
	vector<uint32_t> es;
	vector<uint32_t> fs;
	bool boundary;
};
//-------------------------------------------------------------------
//-------------------------------------------------------------------
struct Singular_V
{
	uint32_t id, hid;
	bool boundary;

	vector<uint32_t> neighbor_svs;//singular vs
	vector<uint32_t> neighbor_ses;//singular es
	bool fake;

	uint32_t which_singularity;
	uint32_t which_singularity_type;
};
struct Singular_E
{
	uint32_t id;
	std::vector<uint32_t> vs;//singular v
	vector<uint32_t> es_link;//hex e
	vector<uint32_t> vs_link;//hex v
	bool boundary;

	vector<uint32_t> neighbor_ses;//singular es
	bool circle;
};
struct Frame_V
{
	uint32_t id, hid, svid;
	uint32_t what_type;
	vector<uint32_t>  neighbor_fvs;
	vector<uint32_t>  neighbor_fes;
	vector<uint32_t>  neighbor_ffs;
	vector<uint32_t>  neighbor_fhs;
	bool boundary;
};
struct Frame_E
{
	uint32_t id;
	bool singular = false;
	std::vector<uint32_t> vs;
	bool boundary;
	vector<uint32_t>  vs_link;//v of hex_v
	vector<uint32_t>  es_link;//e of hex_e
	vector<uint32_t>  neighbor_fes;
	vector<uint32_t>  neighbor_ffs;
	vector<uint32_t>  neighbor_fhs;
};
struct Frame_F
{
	uint32_t id;
	bool boundary;
	uint32_t F_location;

	vector<uint32_t> vs;
	vector<uint32_t> es;
	vector<uint32_t>  fvs_net;
	vector<uint32_t>  ffs_net;
	vector<uint32_t>  neighbor_ffs;
	vector<uint32_t>  neighbor_fhs;
	uint32_t Color_ID;
};
struct Frame_H
{
	uint32_t id;

	std::vector<uint32_t> vs;
	std::vector<uint32_t> es;
	std::vector<uint32_t> fs;
	vector<vector<vector<uint32_t> >> vs_net;
	vector<uint32_t>  fs_net;
	vector<uint32_t>  hs_net;
	vector<uint32_t>  neighbor_fhs;//neighboring cube	
	uint32_t Color_ID;
};


enum Sheet_type {
	open = 0,
	close,
	tagent,
	intersect,
	mobius
};
struct Sheet
{
	uint32_t id;
	Sheet_type type;
	bool fake;
	std::vector<uint32_t> ns;
	std::vector<uint32_t> es;
	std::vector<uint32_t> fs;
	std::vector<uint32_t> cs;

	std::vector<uint32_t> middle_es, middle_es_b, left_es, right_es;
	std::vector<uint32_t> middle_fs, side_fs, left_fs, right_fs;

	vector<vector<uint32_t>> vs_pairs, vs_links, Vs_Group;
	VectorXi target_vs;
	MatrixXd target_coords;

	Float weight, weight_val_average, weight_val_max, weight_val_min, weight_len;
	bool valence_filter;
};
struct CHord
{
	uint32_t id;

	Sheet_type type;
	bool fake;
	bool side;//false 0, true 1 
	std::vector<uint32_t> ns;
	std::vector<uint32_t> es;
	std::vector<uint32_t> fs;
	std::vector<uint32_t> cs;

	std::vector<uint32_t> parallel_ns[4];
	std::vector<uint32_t> parallel_es[4], vertical_es[4];
	std::vector<uint32_t> parallel_fs, vertical_fs[4];;

	vector<uint32_t> tangent_vs;
	vector<uint32_t> tangent_es;
	vector<uint32_t> tangent_fs;
	vector<uint32_t> tangent_cs;


	vector<vector<uint32_t>> Vs_Group;
	VectorXi target_vs;
	MatrixXd target_coords;

	Float weight, weight_val_average, weight_val_max, weight_val_min, weight_len;
	bool valence_filter;

};

struct Singularity
{
	vector<Singular_V> SVs;
	vector<Singular_E> SEs;
};
struct Frame
{
	vector<Frame_V> FVs;
	vector<Frame_E> FEs;
	vector<Frame_F> FFs;
	vector<Frame_H> FHs;
};
enum Mesh_type {
	Tri = 0,
	Qua,
	Tet,
	Hyb,
	Hex
};
struct Mesh_Topology
{
	bool euler_problem;
	bool manifoldness_problem;

	int genus;
	int surface_euler;
	int volume_euler;
	bool surface_manifoldness;
	bool volume_manifoldness;

	bool frame_euler_problem;
	bool frame_manifoldness_problem;

	int frame_genus;
	int frame_surface_euler;
	int frame_volume_euler;
	bool frame_surface_manifoldness;
	bool frame_volume_manifoldness;
};

struct Mesh_Quality
{
	std::string Name;
	double min_Jacobian;
	double ave_Jacobian;
	double deviation_Jacobian;
	VectorXd V_Js;
	VectorXd H_Js;
	VectorXd Num_Js;

	int32_t V_num, H_num;
	int32_t BV_num, BC_num;
	int32_t RemovedSheetChord_num;
	double RemovedCuboid_ratio;
	double Hausdorff_ratio;

	double timings = -1;
};
struct Mesh
{
	short type;
	MatrixXd V;
	vector<Hybrid_V> Vs;
	vector<Hybrid_E> Es;
	vector<Hybrid_F> Fs;
	vector<Hybrid> Hs;
};
struct Mesh_Feature
{//ground-truth feature
	Mesh tri;
	vector<int> V_map, V_map_reverse;

	vector<Vector3d> Tcenters;
	double ave_length;
	double angle_threshold = 140.0 / 180.0 * PAI;

	vector<uint32_t> corners;
	vector<vector<uint32_t>> corner_curves;


	vector<vector<uint32_t>> curve_vs;
	vector<vector<uint32_t>> curve_es;
	vector<bool> circles;

	MatrixXd normal_V, normal_Tri;

	vector<int> v_types;
};
struct Feature_Constraints
{
	vector<Feature_V_Type> V_types;
	vector<int> V_ids;
	vector<bool> RV_type;

	//corner constraints
	Eigen::VectorXi ids_C;
	Eigen::MatrixXd C;
	double lamda_C = 0;
	//tagent plane constraints
	Eigen::VectorXi ids_T;
	Eigen::MatrixXd normal_T;
	Eigen::VectorXd dis_T;
	Eigen::MatrixXd V_T;
	double lamda_T = 0;
	//feature line constraints
	uint32_t num_a;
	Eigen::VectorXi ids_L;
	Eigen::MatrixXd Axa_L;
	Eigen::MatrixXd origin_L;
	double lamda_L = 0;
	//
	vector<vector<uint32_t>> curve_vs;
	vector<int> curveIds;
};
struct Tetralize_Set {
	vector<uint32_t> V_map, Reverse_V_map;
	MatrixXd V;
	MatrixXi T;
	vector<MatrixXd> RT;
	VectorXi b;
	MatrixXd bc;
	Feature_Constraints fc;

	bool projection;
	Eigen::VectorXi s;
	Eigen::MatrixXd sc;
	//global optimization
	bool global;
	double lamda_region = 1.0e+7;
	Eigen::VectorXi regionb;
	Eigen::MatrixXd regionbc;

	vector<vector<uint32_t>> Vgroups;

};
struct Collapse_Info {
	vector<vector<uint32_t>> V_Groups;	
	VectorXi target_vs;
	MatrixXd target_coords;
	vector<uint32_t> hs;

	vector<uint32_t> fs_before, before_region;
	vector<uint32_t> fs_after, after_region;
	vector<uint32_t> fs_subdivided, subd_region;

	vector<uint32_t> Hsregion;
};

extern Mesh_Feature mf;

extern vector<MatrixXi> LocalMeshTri_F;
extern vector<MatrixXd> LocalMeshTri_V;
extern vector<VectorXd> LocalMeshTri_C;
extern vector<MatrixXi> LocalMeshEs;

extern vector<Feature_Constraints> SeLocalProjection;

extern vector<MatrixXi> MeshTri_F;
extern vector<MatrixXd> MeshTri_V, MeshTri_N, MeshTri_VN;
extern vector<VectorXd> MeshTri_C;
extern vector<MatrixXi> MeshEs;

extern vector<MatrixXi> SkeletonTri_F;
extern vector<MatrixXd> SkeletonTri_V;
extern vector<MatrixXd> SkeletonTri_N;
extern vector<MatrixXd> SkeletonTri_C;


extern vector<MatrixXi> FrameTri_F;
extern vector<MatrixXd> FrameTri_V, FrameTri_N;
extern vector<VectorXd> FrameTri_C;

extern vector<MatrixXi> MeshSkeletonTri_F;
extern vector<MatrixXd> MeshSkeletonTri_V, MeshSkeletonTri_N, MeshSkeletonTri_VN;
extern vector<VectorXd> MeshSkeletonTri_C;
extern vector<VectorXd> AOSkeletons;
extern vector<VectorXd> AOs;

extern vector<MatrixXi> CollapseTri_F;
extern vector<MatrixXd> CollapseTri_V, CollapseTri_N;
extern vector<VectorXd> CollapseTri_C;

extern char path_out[300];
//parallel
extern int32_t GRAIN_SIZE;
//meshes
extern MatrixXd tenC, tenC_perm;

extern double diagonal_len;
//temporary
extern Mesh mesh_sheet, mesh_sheetS;
extern bool HEXAHEDRAL_COLLASPE;
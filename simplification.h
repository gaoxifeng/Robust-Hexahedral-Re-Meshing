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
#include "global_functions.h"
#include "base_complex.h"
#include <algorithm>
#include "igl/slim.h"
#include "metro_hausdorff.h"

class simplification
{
public:
	simplification() {
		Remove_Iteration = 0;
		SHARP_FEATURE = false;
		Slim_Iteration_base = 2;
		remove_cuboid_ratio = 1.0;
		remove_sheet_ratio = 1.0;
		Hex_Num_Threshold = 0;
		hausdorff_ratio_threshould = 0.01;
		Projection_range = 1.0;

		TOPOLOGY = true;
		width_sheet = 1;
		Slim_Iteration = 2;
		Slim_Iteration_Limit = 5;
		Projection_limit = 2;
		subdivision_project_range = 1;
		Slim_region =4;
		Slim_global_region = Slim_region;
		hausdorff_ratio = 0;

		mf.angle_threshold = 0;
		
		fc.lamda_C = 1e+3;
		fc.lamda_L = 1e+3;
		fc.lamda_T = 1e+3;
		ts.lamda_region = 1e+7;
		last_candidate_pos = 0;
	};
	~simplification() {};

	void pipeline();
	bool initialize();
	
	void set_sharp_feature(bool sharp_feature) {SHARP_FEATURE = sharp_feature;}
	void set_slim_iteration_base(uint32_t iter) {Slim_Iteration_base = iter;}
	void set_cuboid_ratio(double ratio) {remove_cuboid_ratio = ratio;}
	void set_sheet_ratio(double ratio) {remove_sheet_ratio = ratio;}
	void set_hausdorff_ratio(double ratio) { hausdorff_ratio_threshould = ratio; }
	void set_target_hex_num(uint32_t hex_num) {Hex_Num_Threshold = hex_num;}
	void set_slim_region(double ratio) {Slim_region = ratio;if (Slim_region < 0) Slim_region = 0; Slim_global_region = Slim_region * 2;}

	void extract();
	bool build_sheet_info(uint32_t sheet_id);
	CHord extract_chord(uint32_t &fid, vector<bool> &f_flag);
	bool build_chord_info(uint32_t id);

	void ranking();
	void sheet_chord_weight(Tuple_Candidate &c);
	void dihedral_angle(Float &angle, Float &k_ratio, vector<uint32_t> &cs, uint32_t eid);

	bool remove();
	bool filter_topology_feature(Tuple_Candidate &c);
	bool vs_pair_sheet(uint32_t sheet_id, vector<vector<uint32_t>> &candiate_es_links, vector<vector<uint32_t>> &v_group);
	bool target_surface_sheet(uint32_t sheet_id, vector<vector<uint32_t>> &candiate_es_links, vector<vector<uint32_t>> &v_group);

	bool vs_group_chord(uint32_t chord_id);
	bool target_surface_chord(uint32_t chord_id);

	bool topology_check();
	bool hausdorff_ratio_check(Mesh &m0, Mesh &m1);

	void optimization();
	bool direct_collapse();
	bool tetralize_mesh(Tetralize_Set &ts);
	bool tetralize_mesh_omesh(Tetralize_Set &ts, Mesh &mesh_r);
	bool tetralize_mesh_submesh(Tetralize_Set &ts, Mesh &mesh_r);
	bool grow_region(uint32_t base_num, vector<uint32_t> &frontFs, vector<uint32_t> &regionFs, vector<uint32_t> &newHs, Tetralize_Set &ts, vector<bool> &H_flag, Mesh &mesh_r, bool global);
	bool grow_region2(uint32_t base_num, vector<uint32_t> &frontFs, vector<uint32_t> &regionFs, vector<uint32_t> &newHs, Tetralize_Set &ts, vector<bool> &H_flag, Mesh &mesh_r, bool global);
	void update_feature_variable_index_newmesh(Tetralize_Set &ts, Feature_Constraints &fc,
		vector<uint32_t> &new_V_map, uint32_t new_Vsize);

	void slim_opt(Tetralize_Set &ts, const uint32_t iter);
	void localize_ts(Tetralize_Set &ts, Tetralize_Set &ts_temp, int &vN, vector<int> & mapV, vector<bool> & touchedV_flag);

	void subdivision();
	bool hex_mesh_subdivision();
	uint32_t nearest_tid(vector<uint32_t> &ts, const Vector3d &v, Vector3d &n, Vector3d &pv, double &dis);

	double cuboid_num_original;
	double sheet_num_original;

	int File_num;
	
	bool TOPOLOGY;
	bool SHARP_FEATURE;
	bool OPTIMIZATION_ONLY = false;
	
	uint32_t INVALID_V, INVALID_E;

	uint32_t Hex_Num_Threshold;
	uint32_t Slim_Iteration_base, Slim_Iteration, Slim_Iteration_Limit;
	double Slim_region, Slim_global_region;

	double remove_cuboid_ratio, remove_sheet_ratio, hausdorff_ratio_threshould;
	double hausdorff_ratio;

	int32_t width_sheet;
	int32_t Remove_Iteration;
	int32_t Projection_range, Projection_limit, subdivision_project_range;
	std::vector<Sheet> All_Sheets;
	std::vector<CHord> All_Chords;
	std::vector<Tuple_Candidate> Candidates;

	uint32_t last_candidate_pos;
public:
	h_io io;
	Singularity si;
	Frame frame;
	Mesh mesh;
	base_complex base_com;
	Feature_Constraints fc;
	Mesh_Topology mt;
	Tetralize_Set ts;
	Collapse_Info CI;

	Singularity si_;
	Frame frame_;
	Mesh mesh_;
	vector<uint32_t> V_map, RV_map;

	vector<Mesh_Quality> statistics;
	int output_file_interval = 1;
};

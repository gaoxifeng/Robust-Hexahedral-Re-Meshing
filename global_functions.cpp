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

#include "global_functions.h"
#include "global_types.h"
#include "igl/bounding_box_diagonal.h"
//===================================mesh connectivities===================================
void build_connectivity(Mesh &hmi) {
	hmi.Es.clear(); if (hmi.Hs.size()) hmi.Fs.clear();
	//either hex or tri
	if (hmi.type == Mesh_type::Tri) {
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
		temp.reserve(hmi.Fs.size() * 3);
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i) {
			for (uint32_t j = 0; j < 3; ++j) {
				uint32_t v0 = hmi.Fs[i].vs[j], v1 = hmi.Fs[i].vs[(j + 1) % 3];
				if (v0 > v1) std::swap(v0, v1);
				temp.push_back(std::make_tuple(v0, v1, i, j));
			}
			hmi.Fs[i].es.resize(3);
		}
		std::sort(temp.begin(), temp.end());
		hmi.Es.reserve(temp.size() / 2);
		uint32_t E_num = 0;
		Hybrid_E e; e.boundary = true; e.vs.resize(2);
		for (uint32_t i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				hmi.Es.push_back(e);
			}
			else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
				std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
				hmi.Es[E_num - 1].boundary = false;

			hmi.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : hmi.Vs) v.boundary = false;
		for (uint32_t i = 0; i < hmi.Es.size(); ++i)
			if (hmi.Es[i].boundary) {
				hmi.Vs[hmi.Es[i].vs[0]].boundary = hmi.Vs[hmi.Es[i].vs[1]].boundary = true;
			}
	}
	else if (hmi.type == Mesh_type::Tet) {
		std::vector<std::vector<uint32_t>> total_fs(hmi.Hs.size() * 4);
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> tempF(hmi.Hs.size() * 4);
		std::vector<uint32_t> vs(3);
		for (uint32_t i = 0; i < hmi.Hs.size(); ++i) {
			for (short j = 0; j < 4; j++) {
				for (short k = 0; k < 3; k++) vs[k] = hmi.Hs[i].vs[tet_faces[j][k]];
				uint32_t id = 4 * i + j;
				total_fs[id] = vs;
				std::sort(vs.begin(), vs.end());
				tempF[id] = std::make_tuple(vs[0], vs[1], vs[2], id, i, j);
			}
			hmi.Hs[i].fs.resize(4);
		}
		std::sort(tempF.begin(), tempF.end());
		hmi.Fs.reserve(tempF.size() / 3);
		Hybrid_F f; f.boundary = true;
		uint32_t F_num = 0;
		for (uint32_t i = 0; i < tempF.size(); ++i) {
			if (i == 0 || (i != 0 &&
				(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) || std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
					std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1])))) {
				f.id = F_num; F_num++;
				f.vs = total_fs[std::get<3>(tempF[i])];
				hmi.Fs.push_back(f);
			}
			else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) && std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
				std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1])))
				hmi.Fs[F_num - 1].boundary = false;

			hmi.Hs[std::get<4>(tempF[i])].fs[std::get<5>(tempF[i])] = F_num - 1;
		}

		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp(hmi.Fs.size() * 3);
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i) {
			for (uint32_t j = 0; j < 3; ++j) {
				uint32_t v0 = hmi.Fs[i].vs[j], v1 = hmi.Fs[i].vs[(j + 1) % 3];
				if (v0 > v1) std::swap(v0, v1);
				temp[3 * i + j] = std::make_tuple(v0, v1, i, j);
			}
			hmi.Fs[i].es.resize(3);
		}
		std::sort(temp.begin(), temp.end());
		hmi.Es.reserve(temp.size() / 2);
		uint32_t E_num = 0;
		Hybrid_E e; e.boundary = false; e.vs.resize(2);
		for (uint32_t i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				hmi.Es.push_back(e);
			}
			hmi.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : hmi.Vs) v.boundary = false;
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i)
			if (hmi.Fs[i].boundary) for (uint32_t j = 0; j < 3; ++j) {
				uint32_t eid = hmi.Fs[i].es[j];
				hmi.Es[eid].boundary = true;
				hmi.Vs[hmi.Es[eid].vs[0]].boundary = hmi.Vs[hmi.Es[eid].vs[1]].boundary = true;
			}
	}
	else if(hmi.type == Mesh_type::Hex) {

		std::vector<std::vector<uint32_t>> total_fs(hmi.Hs.size() * 6);
		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> tempF(hmi.Hs.size() * 6);
		std::vector<uint32_t> vs(4);
		for (uint32_t i = 0; i < hmi.Hs.size(); ++i) {
			for (short j = 0; j < 6; j++){
				for (short k = 0; k < 4; k++) vs[k] = hmi.Hs[i].vs[hex_face_table[j][k]];
				uint32_t id = 6 * i + j;
				total_fs[id] = vs;
				std::sort(vs.begin(), vs.end());
				tempF[id] = std::make_tuple(vs[0], vs[1], vs[2], vs[3], id, i, j);
			}
			hmi.Hs[i].fs.resize(6);
		}
		std::sort(tempF.begin(), tempF.end());
		hmi.Fs.reserve(tempF.size() / 3);
		Hybrid_F f; f.boundary = true;
		uint32_t F_num = 0;
		for (uint32_t i = 0; i < tempF.size(); ++i) {
			if (i == 0 || (i != 0 &&
				(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) || std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
					std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1]) || std::get<3>(tempF[i]) != std::get<3>(tempF[i - 1])))) {
				f.id = F_num; F_num++;
				f.vs = total_fs[std::get<4>(tempF[i])];
				hmi.Fs.push_back(f);
			}
			else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) && std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
				std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1]) && std::get<3>(tempF[i]) == std::get<3>(tempF[i - 1])))
				hmi.Fs[F_num - 1].boundary = false;

			hmi.Hs[std::get<5>(tempF[i])].fs[std::get<6>(tempF[i])] = F_num - 1;
		}

		std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp(hmi.Fs.size() * 4);
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i) {
			for (uint32_t j = 0; j < 4; ++j) {
				uint32_t v0 = hmi.Fs[i].vs[j], v1 = hmi.Fs[i].vs[(j + 1) % 4];
				if (v0 > v1) std::swap(v0, v1);
				temp[4 * i + j] = std::make_tuple(v0, v1, i, j);
			}
			hmi.Fs[i].es.resize(4);
		}
		std::sort(temp.begin(), temp.end());
		hmi.Es.reserve(temp.size() / 2);
		uint32_t E_num = 0;
		Hybrid_E e; e.boundary = false; e.vs.resize(2);
		for (uint32_t i = 0; i < temp.size(); ++i) {
			if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
				std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
				e.id = E_num; E_num++;
				e.vs[0] = std::get<0>(temp[i]);
				e.vs[1] = std::get<1>(temp[i]);
				hmi.Es.push_back(e);
			}
			hmi.Fs[std::get<2>(temp[i])].es[std::get<3>(temp[i])] = E_num - 1;
		}
		//boundary
		for (auto &v : hmi.Vs) v.boundary = false;
		for (uint32_t i = 0; i < hmi.Fs.size(); ++i)
			if (hmi.Fs[i].boundary) for (uint32_t j = 0; j < 4; ++j) {
				uint32_t eid = hmi.Fs[i].es[j];
				hmi.Es[eid].boundary = true;
				hmi.Vs[hmi.Es[eid].vs[0]].boundary = hmi.Vs[hmi.Es[eid].vs[1]].boundary = true;
			}
	}
	//f_nhs;
	for (uint32_t i = 0; i < hmi.Hs.size(); i++) {
		for (uint32_t j = 0; j < hmi.Hs[i].fs.size(); j++) hmi.Fs[hmi.Hs[i].fs[j]].neighbor_hs.push_back(i);
	}
	//e_nfs, v_nfs
	for (uint32_t i = 0; i < hmi.Fs.size(); i++) {
		for (uint32_t j = 0; j < hmi.Fs[i].es.size(); j++) hmi.Es[hmi.Fs[i].es[j]].neighbor_fs.push_back(i);
		for (uint32_t j = 0; j < hmi.Fs[i].vs.size(); j++) hmi.Vs[hmi.Fs[i].vs[j]].neighbor_fs.push_back(i);
	}
	//v_nes, v_nvs
	for (uint32_t i = 0; i < hmi.Es.size(); i++) {
		uint32_t v0 = hmi.Es[i].vs[0], v1 = hmi.Es[i].vs[1];
		hmi.Vs[v0].neighbor_es.push_back(i);
		hmi.Vs[v1].neighbor_es.push_back(i);
		hmi.Vs[v0].neighbor_vs.push_back(v1);
		hmi.Vs[v1].neighbor_vs.push_back(v0);
	}
	//e_nhs
	for (uint32_t i = 0; i < hmi.Es.size(); i++) {
		std::vector<uint32_t> nhs;
		for (uint32_t j = 0; j < hmi.Es[i].neighbor_fs.size(); j++) {
			uint32_t nfid = hmi.Es[i].neighbor_fs[j];
			nhs.insert(nhs.end(), hmi.Fs[nfid].neighbor_hs.begin(), hmi.Fs[nfid].neighbor_hs.end());
		}
		std::sort(nhs.begin(), nhs.end()); nhs.erase(std::unique(nhs.begin(), nhs.end()), nhs.end());
		hmi.Es[i].neighbor_hs = nhs;
	}
}
void topology_info(Mesh &mesh, Frame &frame, Mesh_Topology & mt) {

//==================hex-mesh==================//
	for (int i = 0; i < mesh.Vs.size(); i++) if (!mesh.Vs[i].neighbor_hs.size()) ;
	//surface_euler, volume_euler;
	uint32_t num_boundary_v = 0, num_boundary_e = 0, num_boundary_f = 0;
	for (int i = 0; i<mesh.Vs.size(); i++) if (mesh.Vs[i].boundary) num_boundary_v++;
	for (int i = 0; i<mesh.Es.size(); i++) if (mesh.Es[i].boundary) num_boundary_e++;
	for (int i = 0; i<mesh.Fs.size(); i++) if (mesh.Fs[i].boundary) num_boundary_f++;
	mt.surface_euler = num_boundary_v + num_boundary_f - num_boundary_e;
	mt.genus = (2 - mt.surface_euler) / 2;
	mt.volume_euler = mesh.Vs.size() + mesh.Fs.size() - mesh.Es.size() - mesh.Hs.size();
	mt.euler_problem = false;
	if (mt.surface_euler != 2 * mt.volume_euler) { mt.euler_problem = true; }
	//surface_manifoldness;
	mt.manifoldness_problem = false;
	mt.surface_manifoldness = true;
	vector<short> E_flag(mesh.Es.size(), 0); vector<short> V_flag(mesh.Vs.size(), 0);
	for (uint32_t i = 0; i<mesh.Vs.size(); i++){//topology disk
		if (!mesh.Vs[i].boundary) continue;
		vector<vector<uint32_t>> fes;
		for (auto fid : mesh.Vs[i].neighbor_fs) {
			if (mesh.Fs[fid].boundary) fes.push_back(mesh.Fs[fid].es);
		}
		if (!disk_polygon(mesh, frame, fes, E_flag, V_flag, true)) {
			mt.surface_manifoldness = false;
			mt.manifoldness_problem = true;
			break;
		}
	}
	//volume_manifoldness;
	mt.volume_manifoldness = true;
	vector<bool> F_flag(mesh.Fs.size(), false);
	vector<vector<uint32_t>> F_nvs;
//==================frame==================//
	for (int i = 0; i < frame.FVs.size(); i++) if (!frame.FVs[i].neighbor_fhs.size()) ;
	//surface_euler, volume_euler;
	num_boundary_v = num_boundary_e = num_boundary_f = 0;
	for (int i = 0; i<frame.FVs.size(); i++) if (frame.FVs[i].boundary) num_boundary_v++;
	for (int i = 0; i<frame.FEs.size(); i++) if (frame.FEs[i].boundary) num_boundary_e++;
	for (int i = 0; i<frame.FFs.size(); i++) if (frame.FFs[i].boundary) num_boundary_f++;
	mt.frame_surface_euler = num_boundary_v + num_boundary_f - num_boundary_e;
	mt.frame_genus = (2 - mt.frame_surface_euler) / 2;
	mt.frame_volume_euler = frame.FVs.size() + frame.FFs.size() - frame.FEs.size() - frame.FHs.size();
	mt.frame_euler_problem = false;
	if (mt.frame_surface_euler != 2 * mt.frame_volume_euler) { mt.frame_euler_problem = true; }
	//surface_manifoldness;
	mt.frame_manifoldness_problem = false;
	mt.frame_surface_manifoldness = true;
	E_flag.resize(frame.FEs.size()); V_flag.resize(frame.FVs.size());
	for (uint32_t i = 0; i<frame.FVs.size(); i++) {//topology disk
		if (!frame.FVs[i].boundary) continue;
		vector<vector<uint32_t>> fes;
		for (auto fid : frame.FVs[i].neighbor_ffs) {
			if (frame.FFs[fid].boundary) fes.push_back(frame.FFs[fid].es);
		}
		if (!disk_polygon(mesh, frame, fes, E_flag, V_flag, false)) {
			mt.frame_surface_manifoldness = false;
			mt.frame_manifoldness_problem = true;
			break;
		}
	}
	//volume_manifoldness;
	mt.frame_volume_manifoldness = true;
	F_flag.resize(frame.FFs.size()); E_flag.resize(frame.FEs.size()); V_flag.resize(frame.FVs.size());
	F_nvs.resize(frame.FVs.size());
	for (uint32_t i = 0; i<frame.FVs.size(); i++) {//topology sphere
		vector<vector<uint32_t>> pfs;
		for (auto hid : frame.FVs[i].neighbor_fhs) pfs.push_back(frame.FHs[hid].fs);

		if (!sphere_polyhedral(mesh, frame, F_nvs, pfs, F_flag, E_flag, V_flag, false)) {
			mt.frame_volume_manifoldness = false;
			mt.frame_manifoldness_problem = true;
			break;
		}
	}
}
bool disk_polygon(Mesh &mesh, Frame &frame, vector<vector<uint32_t>> &fes, vector<short> &E_flag, vector<short> &V_flag, const bool &Ismesh) {
	vector<uint32_t> pes;
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (E_flag[fes[i][j]]) E_flag[fes[i][j]] = false;
			else E_flag[fes[i][j]] = true;
	}
	for (int i = 0; i < fes.size(); i++)
		for (int j = 0; j < fes[i].size(); j++)
			if (E_flag[fes[i][j]]) { pes.push_back(fes[i][j]); E_flag[fes[i][j]] = false; }
	//test nvs for each v
	function<void(vector<uint32_t> &, const uint32_t &)> two_vs = [&](vector<uint32_t> &vs, const uint32_t & eid) {
		if (Ismesh) { vs = mesh.Es[eid].vs; }
		else { vs = frame.FEs[eid].vs; }
	};

	vector<uint32_t> vs;
	for (uint32_t i = 0; i < pes.size(); ++i) {
		two_vs(vs, pes[i]);
		V_flag[vs[0]]++; V_flag[vs[1]]++;
		if (V_flag[vs[0]] > 2 || V_flag[vs[1]] > 2) {
			for (uint32_t j = 0; j < pes.size(); ++j) {
				two_vs(vs, pes[j]);
				V_flag[vs[0]] = V_flag[vs[1]] = 0;
			}
			return false;
		}
	}
	for (uint32_t i = 0; i < pes.size(); ++i) {
		two_vs(vs, pes[i]);
		V_flag[vs[0]] = V_flag[vs[1]] = 0;
	}
	//extract the polygon	
	if (!pes.size()) return false;
	vector<uint32_t> pvs;
	pvs.reserve(pes.size());
	vector<bool> e_flag(pes.size(), false);
	two_vs(vs, pes[0]);
	pvs.insert(pvs.end(), vs.begin(),vs.end());
	e_flag[0] = true;
	uint32_t start_v = pvs[1];
	for (uint32_t i = 2; i < pes.size(); i++) {
		for (uint32_t j = 1; j < pes.size(); j++) {
			if (!e_flag[j]) {

				two_vs(vs, pes[j]);
				if (vs[0] == start_v) {
					e_flag[j] = true;
					pvs.push_back(vs[1]);
					start_v = vs[1];
					break;
				}
				else if (vs[1] == start_v) {
					e_flag[j] = true;
					pvs.push_back(vs[0]);
					start_v = vs[0];
					break;
				}
			}
		}
	}
	if (pvs.size() != pes.size()) return false;
	return true;
}
bool sphere_polyhedral(Mesh &mesh, Frame &frame, vector<vector<uint32_t>> &F_nvs, vector<vector<uint32_t>> &pfs, vector<bool> &F_flag, vector<short> &E_flag, vector<short> &V_flag, const bool &Ismesh) {

	vector<uint32_t> pf;
	for (int i = 0; i < pfs.size(); i++)
		for (int j = 0; j < pfs[i].size(); j++)
			if (F_flag[pfs[i][j]]) F_flag[pfs[i][j]] = false;
			else F_flag[pfs[i][j]] = true;
	for (int i = 0; i < pfs.size(); i++)
		for (int j = 0; j < pfs[i].size(); j++)
			if (F_flag[pfs[i][j]]) { pf.push_back(pfs[i][j]); F_flag[pfs[i][j]] = false; }
	//test each e whether non-manifold
	function<void(vector<uint32_t> &, const uint32_t &)> four_es = [&](vector<uint32_t> &es, const uint32_t & fid) {
		if (Ismesh) { es = mesh.Fs[fid].es; }
		else { es = frame.FFs[fid].es; }
	};
	function<void(vector<uint32_t> &, const uint32_t &)> four_vs = [&](vector<uint32_t> &vs, const uint32_t & fid) {
		if (Ismesh) { vs = mesh.Fs[fid].vs; }
		else { vs = frame.FFs[fid].vs; }
	};
	bool non_simple = false;
	vector<uint32_t> es;
	for (uint32_t i = 0; i < pf.size(); ++i) {
		four_es(es, pf[i]);
		for (auto eid : es) {
			E_flag[eid]++;
			if (E_flag[eid] > 2) non_simple = true;
		}
		if (non_simple) {
			for (uint32_t k = 0; k < pf.size(); ++k) {
				four_es(es, pf[k]);
				for (auto eid : es) E_flag[eid] = false;
			}
			return false;
		}
	}
	for (uint32_t k = 0; k < pf.size(); ++k){
		four_es(es, pf[k]);
		for (auto eid : es) if(E_flag[eid]!=2) non_simple = true;
	}
	for (uint32_t k = 0; k < pf.size(); ++k) {
		four_es(es, pf[k]);
		for (auto eid : es) E_flag[eid] = false;
	}
	if (non_simple) return false;
	//test each v whether non-manifold
	std::vector<uint32_t> vs_set;
	for (uint32_t i = 0; i < pf.size(); ++i) {
		std::vector<uint32_t> vs;
		four_vs(vs, pf[i]);
		for (auto vid : vs) {
			F_nvs[vid].push_back(i);
			if (!V_flag[vid]) {vs_set.push_back(vid); V_flag[vid] = true;}
		}
	}
	for (auto vid: vs_set)  V_flag[vid] = false;

	for (auto vid: vs_set) {
		sort(F_nvs[vid].begin(), F_nvs[vid].end());
		F_nvs[vid].erase(std::unique(F_nvs[vid].begin(), F_nvs[vid].end()), F_nvs[vid].end());

		vector<vector<uint32_t>> fes;
		for (auto fid : F_nvs[vid])
			if(Ismesh) fes.push_back(mesh.Fs[fid].es); else fes.push_back(frame.FFs[fid].es);

		if (!disk_polygon(mesh, frame, fes, E_flag, V_flag, Ismesh)){
			non_simple = true; break;
		}
	}
	for (auto vid : vs_set) F_nvs[vid].clear();

	if (non_simple) return false;
	return true;
}
bool comp_topology(Mesh_Topology & mt0, Mesh_Topology & mt1) {

	if (mt1.manifoldness_problem) return false;
	if (mt0.genus != mt1.genus) return false;
	if (mt0.frame_genus != mt1.frame_genus) return false;
	return true;
}
bool redundentV_check(Mesh &meshI, Mesh &meshO) {
	bool redundentV = false;
	build_connectivity(meshI);
	vector<int> V_tag(meshI.Vs.size(), -1);
	int vI = 0;
	for (auto v : meshI.Vs) {
		if (!v.neighbor_hs.size() || !v.neighbor_fs.size()) { redundentV = true; continue; }
		V_tag[v.id] = vI++;
	}
	meshO.type = meshI.type;
	meshO.V.resize(3, vI);
	for (int i = 0; i < V_tag.size();i++) {
		if (V_tag[i] == -1) continue;
		Hybrid_V vo;
		vo.id = meshO.Vs.size(); 
		vo.boundary = false;

		meshO.Vs.push_back(vo);
		meshO.V.col(V_tag[i]) = meshI.V.col(i);
	}
	meshO.Hs = meshI.Hs;
	Hybrid h;
	if (meshO.type == Mesh_type::Hex)h.vs.resize(8);
	else if (meshO.type == Mesh_type::Tet)h.vs.resize(4);

	for (auto &h:meshO.Hs) {
		for (uint32_t j = 0; j < h.vs.size(); j++) h.vs[j] = V_tag[h.vs[j]];
		for (uint32_t i = 0; i < h.vs.size(); i++) meshO.Vs[h.vs[i]].neighbor_hs.push_back(h.id);
	}
	return redundentV;
}
double average_edge_length(Mesh &mesh) {
	double ave_length = 0;
	for (uint32_t i = 0; i < mesh.Fs.size(); i++) {
		auto &vs = mesh.Fs[i].vs;
		for (uint32_t j = 0; j < 3; j++) ave_length += (mesh.V.col(vs[j]) - mesh.V.col(vs[(j + 1) % 3])).norm();
	}
	int N = mesh.Fs.size() * 3;
	ave_length /= N;
	return ave_length;
}
void re_indexing_connectivity(Mesh &hmi, MatrixXi &H){
	hmi.Vs.clear();
	hmi.Es.clear();
	hmi.Fs.clear();
	hmi.Hs.clear();
	//Vs
	hmi.Vs.resize(hmi.V.cols());
	for (uint32_t i = 0; i<hmi.V.cols(); i++)
	{
		Hybrid_V v;
		v.id = i; v.boundary = false;
		hmi.Vs[i] = v;
	}
	//Hs
	hmi.Hs.resize(H.cols());
	Hybrid h; ;
	if (hmi.type == Mesh_type::Hex)h.vs.resize(8);
	else if (hmi.type == Mesh_type::Tet)h.vs.resize(4);

	for (uint32_t i = 0; i < H.cols(); i++) {
		for (uint32_t j = 0; j < H.rows(); j++) h.vs[j] = H(j, i);
		h.id = i; hmi.Hs[h.id] = h;
		for (uint32_t i = 0; i < h.vs.size(); i++) hmi.Vs[h.vs[i]].neighbor_hs.push_back(h.id);
	}
//Es, Fs, and their connectivities
	build_connectivity(hmi);
}

void extract_surface_mesh(Mesh &meshi, Mesh &mesho) {
	mesho.Vs.clear(); mesho.Es.clear(); mesho.Fs.clear(); mesho.Hs.clear();
	mesho.type = Mesh_type::Tri;

	vector<bool> V_tag(meshi.Vs.size(), false);
	vector<int32_t> V_map(meshi.Vs.size(), -1), V_map_reverse;

	for (auto f : meshi.Fs) if (f.boundary) {
		for (auto vid : f.vs) V_tag[vid] = true;
		
		if (f.vs.size() == 3) {
			Hybrid_F hf; hf.vs = f.vs;
			mesho.Fs.push_back(hf);
		}
		else if (f.vs.size() == 4) {
			Hybrid_F hf;
			hf.vs.push_back(f.vs[0]);
			hf.vs.push_back(f.vs[1]);
			hf.vs.push_back(f.vs[2]);
			mesho.Fs.push_back(hf);
			hf.vs.clear();
			hf.vs.push_back(f.vs[2]);
			hf.vs.push_back(f.vs[3]);
			hf.vs.push_back(f.vs[0]);
			mesho.Fs.push_back(hf);
		}

	}
	//re-indexing
	uint32_t newV_id = 0;
	for (uint32_t i = 0; i < V_tag.size(); i++) if (V_tag[i]) {
		V_map[i] = newV_id++; V_map_reverse.push_back(i);
	}
	mesho.V.resize(3, newV_id);
	for (uint32_t i = 0; i < V_tag.size(); i++) if (V_tag[i]) {
		Hybrid_V v;
		v.id = mesho.Vs.size(); mesho.Vs.push_back(v);
		mesho.V.col(v.id) = meshi.V.col(i);
	}
	for (uint32_t i = 0; i < mesho.Fs.size(); i++) for (uint32_t j = 0; j < 3; j++) mesho.Fs[i].vs[j] = V_map[mesho.Fs[i].vs[j]];
	orient_surface_mesh(mesho);
}
void  orient_surface_mesh(Mesh &hmi) {

	vector<bool> flag(hmi.Fs.size(), true);
	flag[0] = false;

	std::queue<uint32_t> pf_temp; pf_temp.push(0);
	while (!pf_temp.empty()) {
		uint32_t fid = pf_temp.front(); pf_temp.pop();
		for (auto eid : hmi.Fs[fid].es) for (auto nfid : hmi.Es[eid].neighbor_fs) {
			if (!flag[nfid]) continue;
			uint32_t v0 = hmi.Es[eid].vs[0], v1 = hmi.Es[eid].vs[1];
			int32_t v0_pos = std::find(hmi.Fs[fid].vs.begin(), hmi.Fs[fid].vs.end(), v0) - hmi.Fs[fid].vs.begin();
			int32_t v1_pos = std::find(hmi.Fs[fid].vs.begin(), hmi.Fs[fid].vs.end(), v1) - hmi.Fs[fid].vs.begin();

			if ((v0_pos + 1) % hmi.Fs[fid].vs.size() != v1_pos) swap(v0, v1);

			int32_t v0_pos_ = std::find(hmi.Fs[nfid].vs.begin(), hmi.Fs[nfid].vs.end(), v0) - hmi.Fs[nfid].vs.begin();
			int32_t v1_pos_ = std::find(hmi.Fs[nfid].vs.begin(), hmi.Fs[nfid].vs.end(), v1) - hmi.Fs[nfid].vs.begin();

			if ((v0_pos_ + 1) % hmi.Fs[nfid].vs.size() == v1_pos_) std::reverse(hmi.Fs[nfid].vs.begin(), hmi.Fs[nfid].vs.end());

			pf_temp.push(nfid); flag[nfid] = false;
		}
	}
	Float res = 0;
	Vector3d ori; ori.setZero();
	for (auto f : hmi.Fs) {
		auto &fvs = f.vs;
		Vector3d center; center.setZero(); for (auto vid : fvs) center += hmi.V.col(vid); center /= fvs.size();

		for (uint32_t j = 0; j < fvs.size(); j++) {
			Vector3d x = hmi.V.col(fvs[j]) - ori, y = hmi.V.col(fvs[(j + 1) % fvs.size()]) - ori, z = center - ori;
			res += -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));
		}
	}
	if (res > 0) {
		for (uint32_t i = 0; i < hmi.Fs.size(); i++) std::reverse(hmi.Fs[i].vs.begin(), hmi.Fs[i].vs.end());
	}
}
void orient_triangle_mesh(Mesh &hmi) {
	vector<vector<Float>> p(hmi.V.cols()); vector<vector<uint32_t>> f(hmi.Fs.size());
	for (uint32_t i = 0; i<hmi.V.cols(); i++){
		vector<Float> pv;
		pv.push_back(hmi.V(0,i));
		pv.push_back(hmi.V(1,i));
		pv.push_back(hmi.V(2,i));
		p[i] = pv;
	}
	for (uint32_t i = 0; i < hmi.Fs.size(); i++) f[i] = hmi.Fs[i].vs;

	std::map<std::set<uint32_t>, vector<uint32_t>> edge_2_neb_tri;
	std::set<vector<uint32_t>> direct_edges;
	std::set<std::set<uint32_t>> no_direct_edges;
	vector<bool> whe_tri_in;

	vector<vector<uint32_t>> nf;

	nf = f;

	for (int i = 0; i<f.size(); i++)
	{
		std::set<uint32_t> xa, xb, xc;

		xa.insert(f[i][0]);
		xa.insert(f[i][1]);

		xb.insert(f[i][1]);
		xb.insert(f[i][2]);

		xc.insert(f[i][2]);
		xc.insert(f[i][0]);

		whe_tri_in.push_back(false);

		vector<uint32_t> ct;

		ct.push_back(i);

		if (edge_2_neb_tri.find(xa) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<uint32_t>, vector<uint32_t>>(xa, ct));
		}
		else
		{
			edge_2_neb_tri[xa].push_back(i);
		}

		if (edge_2_neb_tri.find(xb) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<uint32_t>, vector<uint32_t>>(xb, ct));
		}
		else
		{
			edge_2_neb_tri[xb].push_back(i);
		}

		if (edge_2_neb_tri.find(xc) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<uint32_t>, vector<uint32_t>>(xc, ct));
		}
		else
		{
			edge_2_neb_tri[xc].push_back(i);
		}
	}

	std::set<uint32_t> xa, xb, xc;
	vector<uint32_t> ya, yb, yc;

	xa.insert(f[0][0]);
	xa.insert(f[0][1]);

	xb.insert(f[0][1]);
	xb.insert(f[0][2]);

	xc.insert(f[0][2]);
	xc.insert(f[0][0]);

	ya.push_back(f[0][0]);
	ya.push_back(f[0][1]);

	yb.push_back(f[0][1]);
	yb.push_back(f[0][2]);

	yc.push_back(f[0][2]);
	yc.push_back(f[0][0]);

	no_direct_edges.insert(xa);
	no_direct_edges.insert(xb);
	no_direct_edges.insert(xc);

	direct_edges.insert(ya);
	direct_edges.insert(yb);
	direct_edges.insert(yc);

	whe_tri_in[0] = true;

	std::queue<uint32_t> queue_loop;

	for (uint32_t i = 0; i<2; i++)
	{
		if (!whe_tri_in[edge_2_neb_tri[xa][i]])
		{
			queue_loop.push(edge_2_neb_tri[xa][i]);

			whe_tri_in[edge_2_neb_tri[xa][i]] = true;
		}

		if (!whe_tri_in[edge_2_neb_tri[xb][i]])
		{
			queue_loop.push(edge_2_neb_tri[xb][i]);

			whe_tri_in[edge_2_neb_tri[xb][i]] = true;
		}

		if (!whe_tri_in[edge_2_neb_tri[xc][i]])
		{
			queue_loop.push(edge_2_neb_tri[xc][i]);

			whe_tri_in[edge_2_neb_tri[xc][i]] = true;
		}
	}

	while (!queue_loop.empty())
	{
		xa.clear();
		xb.clear();
		xc.clear();

		ya.clear();
		yb.clear();
		yc.clear();

		uint32_t c;

		c = queue_loop.front();

		xa.insert(f[c][0]);
		xa.insert(f[c][1]);

		xb.insert(f[c][1]);
		xb.insert(f[c][2]);

		xc.insert(f[c][2]);
		xc.insert(f[c][0]);

		ya.push_back(f[c][0]);
		ya.push_back(f[c][1]);

		yb.push_back(f[c][1]);
		yb.push_back(f[c][2]);

		yc.push_back(f[c][2]);
		yc.push_back(f[c][0]);

		uint32_t cnt, ct;

		cnt = 0;
		ct = 0;

		if (no_direct_edges.find(xa) != no_direct_edges.end())
		{
			cnt++;
		}

		if (no_direct_edges.find(xb) != no_direct_edges.end())
		{
			cnt++;
		}

		if (no_direct_edges.find(xc) != no_direct_edges.end())
		{
			cnt++;
		}

		if (direct_edges.find(ya) != direct_edges.end())
		{
			ct++;
		}

		if (direct_edges.find(yb) != direct_edges.end())
		{
			ct++;
		}

		if (direct_edges.find(yc) != direct_edges.end())
		{
			ct++;
		}

		if (cnt == 0)
		{
			std::cout << "Error in triangle direction solving!" << std::endl;
			exit(0);
		}

		if (ct != 0)
		{
			ya.clear();
			yb.clear();
			yc.clear();

			ya.push_back(f[c][1]);
			ya.push_back(f[c][0]);

			yb.push_back(f[c][2]);
			yb.push_back(f[c][1]);

			yc.push_back(f[c][0]);
			yc.push_back(f[c][2]);

			nf[c][0] = f[c][1];
			nf[c][1] = f[c][0];

		}
		else
		{
			// do nothing
		}

		no_direct_edges.insert(xa);
		no_direct_edges.insert(xb);
		no_direct_edges.insert(xc);

		direct_edges.insert(ya);
		direct_edges.insert(yb);
		direct_edges.insert(yc);

		for (uint32_t i = 0; i<2; i++)
		{
			if (!whe_tri_in[edge_2_neb_tri[xa][i]])
			{
				queue_loop.push(edge_2_neb_tri[xa][i]);

				whe_tri_in[edge_2_neb_tri[xa][i]] = true;
			}

			if (!whe_tri_in[edge_2_neb_tri[xb][i]])
			{
				queue_loop.push(edge_2_neb_tri[xb][i]);

				whe_tri_in[edge_2_neb_tri[xb][i]] = true;
			}

			if (!whe_tri_in[edge_2_neb_tri[xc][i]])
			{
				queue_loop.push(edge_2_neb_tri[xc][i]);

				whe_tri_in[edge_2_neb_tri[xc][i]] = true;
			}
		}

		queue_loop.pop();
	}

	Float res = 0;

	vector<Float> ori(3, 0);

	for (uint32_t i = 0; i<nf.size(); i++)
		res += uctet(ori, p[nf[i][0]], p[nf[i][1]], p[nf[i][2]]);

	if (res > 0)
	{
		uint32_t tmi;

		for (uint32_t i = 0; i<nf.size(); i++)
		{
			tmi = nf[i][0];
			nf[i][0] = nf[i][1];
			nf[i][1] = tmi;
		}
	}

	for (uint32_t i = 0; i < hmi.Fs.size(); i++)
		hmi.Fs[i].vs = nf[i];
}
void  orient_triangle_mesh(MatrixXd &Tri_V, MatrixXi &Tri_F) {

	vector<vector<Float>> p(Tri_V.cols()); vector<vector<uint32_t>> f(Tri_F.rows());
	for (uint32_t i = 0; i<Tri_V.cols(); i++) {
		vector<Float> pv;
		pv.push_back(Tri_V(0,i));
		pv.push_back(Tri_V(1,i));
		pv.push_back(Tri_V(2,i));
		p[i] = pv;
	}
	for (uint32_t i = 0; i < Tri_F.rows(); i++)
		for (uint32_t j = 0; j < 3; j++)
			f[i].push_back(Tri_F(i, j));

	std::map<std::set<uint32_t>, vector<uint32_t>> edge_2_neb_tri;
	std::set<vector<uint32_t>> direct_edges;
	std::set<std::set<uint32_t>> no_direct_edges;
	vector<bool> whe_tri_in;

	vector<vector<uint32_t>> nf;

	nf = f;

	for (int i = 0; i<f.size(); i++)
	{
		std::set<uint32_t> xa, xb, xc;

		xa.insert(f[i][0]);
		xa.insert(f[i][1]);

		xb.insert(f[i][1]);
		xb.insert(f[i][2]);

		xc.insert(f[i][2]);
		xc.insert(f[i][0]);

		whe_tri_in.push_back(false);

		vector<uint32_t> ct;

		ct.push_back(i);

		if (edge_2_neb_tri.find(xa) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<uint32_t>, vector<uint32_t>>(xa, ct));
		}
		else
		{
			edge_2_neb_tri[xa].push_back(i);
		}

		if (edge_2_neb_tri.find(xb) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<uint32_t>, vector<uint32_t>>(xb, ct));
		}
		else
		{
			edge_2_neb_tri[xb].push_back(i);
		}

		if (edge_2_neb_tri.find(xc) == edge_2_neb_tri.end())
		{
			edge_2_neb_tri.insert(std::pair<std::set<uint32_t>, vector<uint32_t>>(xc, ct));
		}
		else
		{
			edge_2_neb_tri[xc].push_back(i);
		}
	}

	std::set<uint32_t> xa, xb, xc;
	vector<uint32_t> ya, yb, yc;

	xa.insert(f[0][0]);
	xa.insert(f[0][1]);

	xb.insert(f[0][1]);
	xb.insert(f[0][2]);

	xc.insert(f[0][2]);
	xc.insert(f[0][0]);

	ya.push_back(f[0][0]);
	ya.push_back(f[0][1]);

	yb.push_back(f[0][1]);
	yb.push_back(f[0][2]);

	yc.push_back(f[0][2]);
	yc.push_back(f[0][0]);

	no_direct_edges.insert(xa);
	no_direct_edges.insert(xb);
	no_direct_edges.insert(xc);

	direct_edges.insert(ya);
	direct_edges.insert(yb);
	direct_edges.insert(yc);

	whe_tri_in[0] = true;

	std::queue<uint32_t> queue_loop;

	for (uint32_t i = 0; i<2; i++)
	{
		if (!whe_tri_in[edge_2_neb_tri[xa][i]])
		{
			queue_loop.push(edge_2_neb_tri[xa][i]);

			whe_tri_in[edge_2_neb_tri[xa][i]] = true;
		}

		if (!whe_tri_in[edge_2_neb_tri[xb][i]])
		{
			queue_loop.push(edge_2_neb_tri[xb][i]);

			whe_tri_in[edge_2_neb_tri[xb][i]] = true;
		}

		if (!whe_tri_in[edge_2_neb_tri[xc][i]])
		{
			queue_loop.push(edge_2_neb_tri[xc][i]);

			whe_tri_in[edge_2_neb_tri[xc][i]] = true;
		}
	}

	while (!queue_loop.empty())
	{
		xa.clear();
		xb.clear();
		xc.clear();

		ya.clear();
		yb.clear();
		yc.clear();

		uint32_t c;

		c = queue_loop.front();

		xa.insert(f[c][0]);
		xa.insert(f[c][1]);

		xb.insert(f[c][1]);
		xb.insert(f[c][2]);

		xc.insert(f[c][2]);
		xc.insert(f[c][0]);

		ya.push_back(f[c][0]);
		ya.push_back(f[c][1]);

		yb.push_back(f[c][1]);
		yb.push_back(f[c][2]);

		yc.push_back(f[c][2]);
		yc.push_back(f[c][0]);

		uint32_t cnt, ct;

		cnt = 0;
		ct = 0;

		if (no_direct_edges.find(xa) != no_direct_edges.end())
		{
			cnt++;
		}

		if (no_direct_edges.find(xb) != no_direct_edges.end())
		{
			cnt++;
		}

		if (no_direct_edges.find(xc) != no_direct_edges.end())
		{
			cnt++;
		}

		if (direct_edges.find(ya) != direct_edges.end())
		{
			ct++;
		}

		if (direct_edges.find(yb) != direct_edges.end())
		{
			ct++;
		}

		if (direct_edges.find(yc) != direct_edges.end())
		{
			ct++;
		}

		if (cnt == 0)
		{
			std::cout << "Error in triangle direction solving!" << std::endl;
			exit(0);
		}

		if (ct != 0)
		{
			ya.clear();
			yb.clear();
			yc.clear();

			ya.push_back(f[c][1]);
			ya.push_back(f[c][0]);

			yb.push_back(f[c][2]);
			yb.push_back(f[c][1]);

			yc.push_back(f[c][0]);
			yc.push_back(f[c][2]);

			nf[c][0] = f[c][1];
			nf[c][1] = f[c][0];

		}
		else
		{
			// do nothing
		}

		no_direct_edges.insert(xa);
		no_direct_edges.insert(xb);
		no_direct_edges.insert(xc);

		direct_edges.insert(ya);
		direct_edges.insert(yb);
		direct_edges.insert(yc);

		for (uint32_t i = 0; i<2; i++)
		{
			if (!whe_tri_in[edge_2_neb_tri[xa][i]])
			{
				queue_loop.push(edge_2_neb_tri[xa][i]);

				whe_tri_in[edge_2_neb_tri[xa][i]] = true;
			}

			if (!whe_tri_in[edge_2_neb_tri[xb][i]])
			{
				queue_loop.push(edge_2_neb_tri[xb][i]);

				whe_tri_in[edge_2_neb_tri[xb][i]] = true;
			}

			if (!whe_tri_in[edge_2_neb_tri[xc][i]])
			{
				queue_loop.push(edge_2_neb_tri[xc][i]);

				whe_tri_in[edge_2_neb_tri[xc][i]] = true;
			}
		}

		queue_loop.pop();
	}

	Float res = 0;

	vector<Float> ori(3, 0);

	for (uint32_t i = 0; i<nf.size(); i++)
		res += uctet(ori, p[nf[i][0]], p[nf[i][1]], p[nf[i][2]]);

	if (res > 0)
	{
		uint32_t tmi;

		for (uint32_t i = 0; i<nf.size(); i++)
		{
			tmi = nf[i][0];
			nf[i][0] = nf[i][1];
			nf[i][1] = tmi;
		}
	}

	for (uint32_t i = 0; i < Tri_F.rows(); i++)
		for (uint32_t j = 0; j < 3; j++)
			Tri_F(i, j) = nf[i][j];
}
Float uctet(vector<Float> a, vector<Float> b, vector<Float> c, vector<Float> d) {
	Float res = 0;
	vector<Float> x, y, z;

	for (uint32_t i = 0; i<3; i++){
		x.push_back(b[i] - a[i]);
		y.push_back(c[i] - a[i]);
		z.push_back(d[i] - a[i]);
	}
	res = -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));

	return res;
}

//===================================mesh quality==========================================
bool scaled_jacobian(Mesh &hmi, Mesh_Quality &mq)
{
	if (hmi.type != Mesh_type::Hex) return false;

	mq.ave_Jacobian = 0;
	mq.min_Jacobian = 1;
	mq.deviation_Jacobian = 0;
	mq.V_Js.resize(hmi.Hs.size() * 8); mq.V_Js.setZero();
	mq.H_Js.resize(hmi.Hs.size()); mq.H_Js.setZero();

	for (uint32_t i = 0; i<hmi.Hs.size(); i++)
	{
		double hex_minJ = 1;
		for (uint32_t j = 0; j<8; j++)
		{
			uint32_t v0, v1, v2, v3;
			v0 = hex_tetra_table[j][0]; v1 = hex_tetra_table[j][1];
			v2 = hex_tetra_table[j][2]; v3 = hex_tetra_table[j][3];

			Vector3d c0 = hmi.V.col(hmi.Hs[i].vs[v0]);
			Vector3d c1 = hmi.V.col(hmi.Hs[i].vs[v1]);
			Vector3d c2 = hmi.V.col(hmi.Hs[i].vs[v2]);
			Vector3d c3 = hmi.V.col(hmi.Hs[i].vs[v3]);

			double jacobian_value = a_jacobian(c0, c1, c2, c3);

			if (hex_minJ>jacobian_value) hex_minJ = jacobian_value;

			uint32_t id = 8 * i + j;
			mq.V_Js[id] = jacobian_value;
		}
		mq.H_Js[i] = hex_minJ;
		mq.ave_Jacobian += hex_minJ;
		if (mq.min_Jacobian > hex_minJ) mq.min_Jacobian = hex_minJ;

	}
	mq.ave_Jacobian /= hmi.Hs.size();
	for (int i = 0; i < mq.H_Js.size(); i++)
		mq.deviation_Jacobian += (mq.H_Js[i] - mq.ave_Jacobian)*(mq.H_Js[i] - mq.ave_Jacobian);
	mq.deviation_Jacobian/= hmi.Hs.size();

	return true;
}
double a_jacobian(Vector3d &v0, Vector3d &v1, Vector3d &v2, Vector3d &v3)
{
	Matrix3d Jacobian;
	
	Jacobian.col(0) = v1 - v0;
	Jacobian.col(1) = v2 - v0;
	Jacobian.col(2) = v3 - v0;


	double norm1 = Jacobian.col(0).norm();
	double norm2 = Jacobian.col(1).norm();
	double norm3 = Jacobian.col(2).norm();

	double scaled_jacobian = Jacobian.determinant();
	if (std::abs(norm1) < Precision || std::abs(norm2) < Precision || std::abs(norm3) < Precision){
		system("PAUSE");
	}
	scaled_jacobian /= norm1*norm2*norm3;
	return scaled_jacobian;
}
//===================================feature v tags==========================================
bool triangle_mesh_feature(Mesh_Feature &mf, Mesh &hmi){

	if (hmi.type != Mesh_type::Hex) return false;
	//tri-mesh
	const uint32_t INVALID_V = hmi.Vs.size();

	Mesh &mesh = mf.tri;
	mesh.Vs.clear(); mesh.Es.clear(); mesh.Fs.clear(); mesh.Hs.clear();
	mesh.type = Mesh_type::Tri;

	mf.V_map_reverse.clear();
	mf.V_map.resize(hmi.Vs.size()); std::fill(mf.V_map.begin(), mf.V_map.end(), -1);

	vector<bool> V_tag(hmi.Vs.size(), false);
	
	for (auto f : hmi.Fs) if (f.boundary) {
		for (auto vid : f.vs) V_tag[vid] = true;
		Hybrid_F t;
		t.id = mesh.Fs.size();
		t.vs.resize(3);
		t.vs[0] = f.vs[0];
		t.vs[1] = f.vs[1];
		t.vs[2] = f.vs[2];
		mesh.Fs.push_back(t);

		t.id = mesh.Fs.size();
		t.vs[0] = f.vs[2];
		t.vs[1] = f.vs[3];
		t.vs[2] = f.vs[0];
		mesh.Fs.push_back(t);
	}
	//re-indexing
	uint32_t newV_id = 0;
	for (uint32_t i = 0; i < V_tag.size(); i++) if (V_tag[i]) {
		mf.V_map[i] = newV_id++; mf.V_map_reverse.push_back(i);
	}
	mesh.V.resize(3, newV_id);
	for (uint32_t i = 0; i < V_tag.size(); i++) if (V_tag[i]) {
		Hybrid_V v;
		v.id = mesh.Vs.size(); mesh.Vs.push_back(v);
		mesh.V.col(v.id) = hmi.V.col(i);	
	}
	for (uint32_t i = 0; i < mesh.Fs.size();i++) for (uint32_t j = 0; j < 3; j++) mesh.Fs[i].vs[j] = mf.V_map[mesh.Fs[i].vs[j]];
	//orient direction
	orient_triangle_mesh(mesh);
	//connectivity
	build_connectivity(mesh);

//normals
	mf.normal_Tri.resize(3, mesh.Fs.size()); mf.normal_Tri.setZero();
	mf.normal_V.resize(3, mesh.Vs.size()); mf.normal_V.setZero();
	mf.Tcenters.clear(); mf.Tcenters.resize(mesh.Fs.size());
	for (uint32_t i = 0; i < mesh.Fs.size(); i++) {
		const auto &vs = mesh.Fs[i].vs;
		Vector3d vec0 = mesh.V.col(vs[1]) - mesh.V.col(vs[0]);
		Vector3d vec1 = mesh.V.col(vs[2]) - mesh.V.col(vs[0]);

		mf.normal_Tri.col(i) = (vec0.cross(vec1)).normalized();
		mf.normal_V.col(vs[0]) += mf.normal_Tri.col(i);
		mf.normal_V.col(vs[1]) += mf.normal_Tri.col(i);
		mf.normal_V.col(vs[2]) += mf.normal_Tri.col(i);

		mf.Tcenters[i].setZero();
		mf.Tcenters[i] += mesh.V.col(vs[0]);
		mf.Tcenters[i] += mesh.V.col(vs[1]);
		mf.Tcenters[i] += mesh.V.col(vs[2]);
		mf.Tcenters[i] /= 3;
	}
	for (uint32_t i = 0; i<mf.normal_V.cols(); ++i) 
		if (mf.normal_V.col(i) != Vector3d::Zero()) mf.normal_V.col(i) = mf.normal_V.col(i).normalized();
	mf.ave_length = 0;
	for (uint32_t i = 0; i < mesh.Es.size(); i++) {
			uint32_t v0 = mesh.Es[i].vs[0];
			uint32_t v1 = mesh.Es[i].vs[1];

			mf.ave_length += (mesh.V.col(v0) - mesh.V.col(v1)).norm();
	}
	mf.ave_length /= mesh.Es.size();
//feature edges, vs
	vector<bool> E_feature_flag(mesh.Es.size(), false);
	mf.v_types.resize(mesh.Vs.size()); fill(mf.v_types.begin(), mf.v_types.end(), 0);
	vector<Float> Dihedral_angles(mesh.Es.size());
	for (uint32_t i = 0; i < mesh.Es.size(); i++) {
		Vector3d n0 = mf.normal_Tri.col(mesh.Es[i].neighbor_fs[0]);
		Vector3d n1 = mf.normal_Tri.col(mesh.Es[i].neighbor_fs[1]);

		Dihedral_angles[i] = PAI - acos(n0.dot(n1));

		if (Dihedral_angles[i] < mf.angle_threshold) {
			E_feature_flag[i] = true;

			uint32_t v0 = mesh.Es[i].vs[0];
			uint32_t v1 = mesh.Es[i].vs[1];

			mf.v_types[v0]++;
			mf.v_types[v1]++;
		}
	}
	mf.corners.clear();
	for (uint32_t i = 0; i < mf.v_types.size(); i++) {

		if (mf.v_types[i] == 0) mf.v_types[i] = Feature_V_Type::REGULAR;
		else if (mf.v_types[i] == 1  || mf.v_types[i] > 2) {
			mf.v_types[i] = Feature_V_Type::CORNER;
			mf.corners.push_back(i);
		}
		else if (mf.v_types[i] == 2) {
			vector<Vector3d> ns; vector<int> vs;
			for (auto eid : mesh.Vs[i].neighbor_es) if (E_feature_flag[eid]) {
				uint32_t v0 = mesh.Es[eid].vs[0]; vs.push_back(v0);
				uint32_t v1 = mesh.Es[eid].vs[1]; vs.push_back(v1);
				ns.push_back((mesh.V.col(v0) - mesh.V.col(v1)).normalized());
			}

			if (vs[0] == vs[2] || vs[1] == vs[3]) ns[0] *= -1;
			double angles = PAI - acos(ns[0].dot(ns[1]));

			if (angles < mf.angle_threshold) {
				mf.v_types[i] = Feature_V_Type::CORNER;
				mf.corners.push_back(i);
			}
		}

	}
	mf.corner_curves.clear(); mf.corner_curves.resize(mf.corners.size());
//feature curves
	mf.curve_es.clear(); mf.curve_vs.clear(); mf.circles.clear();
	uint32_t INVALID_E = mesh.Es.size();
	vector<bool> E_flag(mesh.Es.size(), false);
	for (uint32_t i = 0; i < mesh.Es.size(); i++) {

		if (!E_feature_flag[i]) continue;
		if (E_flag[i]) continue;

		std::function<bool(uint32_t, uint32_t, uint32_t &)> feature_line_proceed = [&](uint32_t vid, uint32_t eid, uint32_t &neid)->bool {

			if (mf.v_types[vid] == Feature_V_Type::REGULAR || mf.v_types[vid] == Feature_V_Type::CORNER) return false;

			for (uint32_t j = 0; j < mesh.Vs[vid].neighbor_es.size(); j++) {
				uint32_t cur_e = mesh.Vs[vid].neighbor_es[j];
				if (cur_e == eid || !E_feature_flag[cur_e] || E_flag[cur_e]) continue;
				neid = cur_e;
				return true;
			}
			return false;
		};

		uint32_t v_left = mesh.Es[i].vs[0], v_right = mesh.Es[i].vs[1];
		uint32_t sv_left, sv_right;
		std::vector<uint32_t> vs_left, vs_right, es_left, es_right;

		bool is_circle = false;
		//left 
		es_left.push_back(i); vs_left.push_back(v_left);
		uint32_t cur_e = i, next_e = INVALID_E;
		while (feature_line_proceed(v_left, cur_e, next_e)) {
			cur_e = next_e;
			E_flag[next_e] = true;
			if (cur_e == i) { is_circle = true; break; }
			es_left.push_back(next_e);
			if (mesh.Es[cur_e].vs[0] == v_left) v_left = mesh.Es[cur_e].vs[1]; else v_left = mesh.Es[cur_e].vs[0];
			vs_left.push_back(v_left);
		}

		if (is_circle) {
			mf.circles.push_back(true);
			mf.curve_es.push_back(es_left); mf.curve_vs.push_back(vs_left);
			continue;
		}
		//right
		vs_right.push_back(v_right);
		cur_e = i, next_e = INVALID_E;
		while (feature_line_proceed(v_right, cur_e, next_e)) {
			cur_e = next_e;
			E_flag[next_e] = true;
			if (mesh.Es[cur_e].vs[0] == v_right) v_right = mesh.Es[cur_e].vs[1]; else v_right = mesh.Es[cur_e].vs[0];
			vs_right.push_back(v_right);
			es_right.push_back(next_e);
		}
		std::reverse(vs_left.begin(), vs_left.end());
		vector<uint32_t> vs_link;
		vs_link = vs_left; vs_link.insert(vs_link.end(), vs_right.begin(), vs_right.end());
		std::reverse(es_left.begin(), es_left.end());
		vector<uint32_t> es_link;
		es_link = es_left; es_link.insert(es_link.end(), es_right.begin(), es_right.end());
	
		mf.curve_es.push_back(es_link); mf.curve_vs.push_back(vs_link);
		mf.circles.push_back(false);
	}
//
	for (uint32_t i = 0; i < mf.curve_vs.size(); i++) for (auto vid : mf.curve_vs[i])
		if (mf.v_types[vid] != Feature_V_Type::CORNER) mf.v_types[vid] = i;
		else {
			int pos = find(mf.corners.begin(), mf.corners.end(), vid) - mf.corners.begin();
			mf.corner_curves[pos].push_back(i);
		}
}
bool initial_feature(Mesh_Feature &mf, Feature_Constraints &fc, Mesh &hmi) {

	fc.V_types.resize(hmi.Vs.size()); fill(fc.V_types.begin(), fc.V_types.end(), Feature_V_Type::INTERIOR);
	fc.V_ids.resize(hmi.Vs.size());
	fc.RV_type.resize(hmi.Vs.size());
	fill(fc.RV_type.begin(), fc.RV_type.end(), false);
	//matrices
	uint32_t num_corners = 0, num_lines = 0, num_regulars = 0;
	for (uint32_t i = 0; i < mf.v_types.size(); i++) {
		if (mf.v_types[i] == Feature_V_Type::REGULAR) {
			fc.V_types[mf.V_map_reverse[i]] = Feature_V_Type::REGULAR;
			fc.V_ids[mf.V_map_reverse[i]] = i;
			num_regulars++;
		}
		else if (mf.v_types[i] == Feature_V_Type::CORNER) {
			fc.V_types[mf.V_map_reverse[i]] = Feature_V_Type::CORNER;
			fc.V_ids[mf.V_map_reverse[i]] = i;
			num_corners++;
		}
		else if (mf.v_types[i] != Feature_V_Type::INTERIOR) {
			fc.V_types[mf.V_map_reverse[i]] = Feature_V_Type::LINE;
			fc.V_ids[mf.V_map_reverse[i]] = mf.v_types[i];
			num_lines++;
		}
	}
	fc.ids_C.resize(num_corners); fc.C.resize(num_corners, 3);
	fc.num_a = num_lines; fc.ids_L.resize(num_lines); fc.Axa_L.resize(num_lines, 3); fc.origin_L.resize(num_lines, 3);
	fc.ids_T.resize(num_regulars); fc.normal_T.resize(num_regulars, 3); fc.dis_T.resize(num_regulars); fc.V_T.resize(num_regulars, 3);
	num_corners = num_lines = num_regulars = 0;
	for (uint32_t i = 0; i < hmi.Vs.size(); i++) {
		if (fc.V_types[i] == Feature_V_Type::CORNER) {
			fc.ids_C[num_corners] = i;
			fc.C.row(num_corners++) = mf.tri.V.col(fc.V_ids[i]);
		}
		else if (fc.V_types[i] == Feature_V_Type::LINE) {
			fc.ids_L[num_lines] = i;
			fc.origin_L.row(num_lines) = hmi.V.col(i);
			uint32_t curve_id = fc.V_ids[i];
			vector<uint32_t> &curve = mf.curve_vs[curve_id];
			if (find(curve.begin(), curve.end(), mf.V_map[i]) == curve.end()) {
				cout << "ERROR in curve" << endl; system("PAUSE");
			}
			int pos = find(curve.begin(), curve.end(), mf.V_map[i]) - curve.begin();
			Vector3d tangent(0, 0, 0);
			uint32_t curve_len = curve.size();
			if (mf.circles[curve_id] || (!mf.circles[curve_id] && pos !=0 && pos != curve_len - 1)) {
				int32_t pos_0 = (pos -1 + curve_len) % curve_len, pos_1 = (pos + 1) % curve_len;
				tangent += (mf.tri.V.col(curve[pos]) - mf.tri.V.col(curve[pos_0])).normalized();
				tangent += (mf.tri.V.col(curve[pos_1]) - mf.tri.V.col(curve[pos])).normalized();
			}
			else if (!mf.circles[curve_id] && pos == 0) {
				int32_t pos_1 = (pos + 1) % curve_len;
				tangent += (mf.tri.V.col(curve[pos_1]) - mf.tri.V.col(curve[pos])).normalized();
			}
			else if (!mf.circles[curve_id] && pos == curve_len - 1) {
				int32_t pos_0 = (pos - 1 + curve_len) % curve_len;
				tangent += (mf.tri.V.col(curve[pos]) - mf.tri.V.col(curve[pos_0])).normalized();
			}
			if(tangent == Vector3d::Zero()){ 
				cout << "ERROR in curve" << endl; system("PAUSE"); 
			}
			tangent.normalize();
			fc.Axa_L.row(num_lines++) = tangent;
		}
		else if (fc.V_types[i] == Feature_V_Type::REGULAR) {
			fc.ids_T[num_regulars] = i;
			uint32_t vid = fc.V_ids[i];
			fc.normal_T.row(num_regulars) = mf.normal_V.col(vid);
			fc.V_T.row(num_regulars) = mf.tri.V.col(vid); 
			fc.dis_T[num_regulars++] = mf.normal_V.col(vid).dot(mf.tri.V.col(vid));
		}
	}

	return true;
}
bool project_surface_update_feature(Mesh_Feature &mf, Feature_Constraints &fc, MatrixXd &V, VectorXi &b, MatrixXd &bc, uint32_t Loop) {

	uint32_t bc_num = fc.ids_C.size() + fc.ids_L.size() + fc.ids_T.size();
	vector<std::tuple<Feature_V_Type, uint32_t, uint32_t, uint32_t>> CI;
	b.resize(bc_num);
	bc.resize(bc_num, 3);
	bc.setZero();
	bc_num = 0;
	uint32_t num_corners = 0, num_lines = 0, num_regulars = 0;
		
	for (uint32_t i = 0; i < fc.V_types.size(); i++) {
		if (fc.V_types[i] == Feature_V_Type::CORNER)
			CI.push_back(std::make_tuple(Feature_V_Type::CORNER, i, num_corners++, bc_num++));
		else if (fc.V_types[i] == Feature_V_Type::LINE)
			CI.push_back(std::make_tuple(Feature_V_Type::LINE, i, num_lines++, bc_num++));
		else if (fc.V_types[i] == Feature_V_Type::REGULAR)
			CI.push_back(std::make_tuple(Feature_V_Type::REGULAR, i, num_regulars++, bc_num++));
	}

	GRAIN_SIZE = 10;
#if 0
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)CI.size(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t m = range.begin(); m != range.end(); m++) {
#endif
			for (uint32_t m = 0; m < CI.size(); m++) {

			Feature_V_Type type = std::get<0>(CI[m]);
			uint32_t i = get<1>(CI[m]);
			uint32_t mi = get<2>(CI[m]);
			uint32_t bci = get<3>(CI[m]);
			Vector3d pv, v;

			if (type == Feature_V_Type::CORNER) {
				pv = mf.tri.V.col(fc.V_ids[i]);
				fc.C.row(mi) = pv;
			}
			else if (type == Feature_V_Type::LINE) {
				pv.setZero();
				v = V.row(i);
				uint32_t curve_id = fc.V_ids[i];
				vector<uint32_t> &curve = mf.curve_vs[curve_id];
				Vector3d tangent(1, 0, 0);
				uint32_t curve_len = curve.size();

				if (!mf.circles[curve_id]) curve_len--;

				vector<Vector3d> pvs, tangents;
				vector<pair<double, uint32_t>> dis_ids;

				for (uint32_t j = 0; j < curve_len; j++) {
					uint32_t pos_0 = curve[j], pos_1 = curve[(j + 1) % curve.size()];
					double t, precision_here = 1.0e1;
					point_line_projection(mf.tri.V.col(pos_0), mf.tri.V.col(pos_1), v, pv, t);
					{
						tangent = (mf.tri.V.col(pos_1) - mf.tri.V.col(pos_0)).normalized();

						dis_ids.push_back(make_pair((v - pv).norm(), pvs.size()));
						pvs.push_back(pv);
						tangents.push_back(tangent);

					}
				}
				sort(dis_ids.begin(), dis_ids.end());

				if (dis_ids.size()) {
					uint32_t cloestid = dis_ids[0].second;
					pv = pvs[cloestid];
					tangent = tangents[cloestid];
				}
				else {
					for (uint32_t j = 0; j < curve.size(); j++) {
						double dis = (mf.tri.V.col(curve[j]) - v).norm();
						dis_ids.push_back(make_pair(dis, j));
					}
					sort(dis_ids.begin(), dis_ids.end());

					int pos = dis_ids[0].second;
					pv = mf.tri.V.col(curve[pos]);

					curve_len = curve.size();
					if (mf.circles[curve_id] || (!mf.circles[curve_id] && pos != 0 && pos != curve_len - 1)) {
						uint32_t pos_0 = (pos - 1 + curve_len) % curve_len, pos_1 = (pos + 1) % curve_len;
						tangent += (mf.tri.V.col(curve[pos]) - mf.tri.V.col(curve[pos_0])).normalized();
						tangent += (mf.tri.V.col(curve[pos_1]) - mf.tri.V.col(curve[pos])).normalized();
					}
					else if (!mf.circles[curve_id] && pos == 0) {
						uint32_t pos_1 = (pos + 1) % curve_len;
						tangent += (mf.tri.V.col(curve[pos_1]) - mf.tri.V.col(curve[pos])).normalized();
					}
					else if (!mf.circles[curve_id] && pos == curve_len - 1) {
						uint32_t pos_0 = (pos - 1 + curve_len) % curve_len;
						tangent += (mf.tri.V.col(curve[pos]) - mf.tri.V.col(curve[pos_0])).normalized();
					}
					tangent.normalize();
				}
				fc.origin_L.row(mi) = pv;
				fc.Axa_L.row(mi) = tangent;
			}
			else if (type == Feature_V_Type::REGULAR) {
				//breadth-first search for the best triangle plane
				vector<uint32_t>  tids;
				uint32_t tid;
				Vector3d interpolP, interpolN;
				pv.setZero();
				v = V.row(i);

				tid = fc.V_ids[i];
				if (!fc.RV_type[i]) tids = mf.tri.Vs[tid].neighbor_fs;
				else tids.push_back(tid);

				Vector3d PreinterpolP = fc.V_T.row(mi), PreinterpolN = fc.normal_T.row(mi);
				bool found = phong_projection(tids, Loop, tid, v, interpolP, interpolN, PreinterpolP, PreinterpolN);
				
				if(!found || (v - interpolP).norm() >= mf.ave_length){
					tid = 0;
					double min_dis = (v - mf.Tcenters[0]).norm();
					for (uint32_t j = 1; j < mf.Tcenters.size(); j++) {
						if (min_dis > (v - mf.Tcenters[j]).norm()) {
							tid = j;
							min_dis = (v - mf.Tcenters[j]).norm();
						}
					}
					uint32_t tid_temp = tid;
					tids.clear();
					tids.push_back(tid_temp);
					if (phong_projection(tids, Loop, tid_temp, v, interpolP, interpolN, PreinterpolP, PreinterpolN))
						tid = tid_temp;
					else{
							interpolP = mf.Tcenters[tid_temp];
							interpolN = mf.normal_Tri.col(tid_temp);
					}
				}
				pv = interpolP;
				fc.V_ids[i] = tid;
				fc.RV_type[i] = true;
				fc.normal_T.row(mi) = interpolN;
				fc.dis_T[mi] = interpolN.dot(interpolP);
				fc.V_T.row(mi) = interpolP;
			}
			else continue;			
			b[bci] = i;
			bc.row(bci) = pv;
		}
#if 0
	}
	);
#endif
	return true;
}
bool phong_projection(vector<uint32_t> &tids, uint32_t Loop, uint32_t &tid, Vector3d &v, Vector3d &interpolP, Vector3d &interpolN, Vector3d &PreinterpolP, Vector3d &PreinterpolN) {
	vector<bool> t_flag(mf.tri.Fs.size(), false);

	vector<uint32_t> tids_;
	for (uint32_t Iter = 0; Iter < Loop; Iter++) {
		for (uint32_t j = 0; j < tids.size(); j++) {
			vector<uint32_t> &vs = mf.tri.Fs[tids[j]].vs;
			for (uint32_t k = 0; k < 3; k++) {
				for (auto ntid : mf.tri.Vs[vs[k]].neighbor_fs) {
					if (t_flag[ntid]) continue; t_flag[ntid] = true;
					tids_.push_back(ntid);
				}
			}
		}
		tids.insert(tids.end(), tids_.begin(), tids_.end()); tids_.clear();
	}
	sort(tids.begin(), tids.end()); tids.erase(unique(tids.begin(), tids.end()), tids.end());
	
	vector<uint32_t> ts;
	bool found = false;
	for (auto id : tids) { t_flag[id] = true; ts.push_back(id); }
	while (ts.size()) {
		tids.clear();
		vector<Vector3d> pvs, pns;
		vector<Vector2d> uvs; vector<pair<double, uint32_t>> dis_ids;
		for (uint32_t j = 0; j < ts.size(); j++) {
			vector<Vector3d> tri_vs(3), vs_normals(3);
			vector<uint32_t> &vs = mf.tri.Fs[ts[j]].vs;
			for (uint32_t k = 0; k < 3; k++) {
				tri_vs[k] = mf.tri.V.col(vs[k]);
				vs_normals[k] = mf.normal_V.col(vs[k]);
			}
			Vector2d uv;
			projectPointOnTriangle(tri_vs, vs_normals, v, uv, interpolP, interpolN);

			if ((uv.x() >= 0.0 || num_equal(uv.x(), 0.0, Precision_Pro)) && (uv.y() >= 0.0 || num_equal(uv.y(), 0.0, Precision_Pro)) &&
				(1 - uv.x() - uv.y() >= 0.0 || num_equal(1 - uv.x() - uv.y(), 0.0, Precision_Pro))) {

				if (PreinterpolN!=Vector3d::Zero() && interpolN.dot(PreinterpolN) <= 0)continue;

				dis_ids.push_back(make_pair((v - interpolP).norm(), tids.size()));
				pvs.push_back(interpolP); pns.push_back(interpolN); uvs.push_back(uv);
				tids.push_back(ts[j]);
			}
		}
		sort(dis_ids.begin(), dis_ids.end());

		if (dis_ids.size()) {
			found = true;
			uint32_t cloestid = dis_ids[0].second;
			interpolP = pvs[cloestid];
			tid = tids[cloestid];
			interpolN = pns[cloestid];
			return true;
		}

		vector<uint32_t> ts_;
		for (uint32_t j = 0; j < ts.size(); j++) {
			vector<uint32_t> &vs = mf.tri.Fs[ts[j]].vs;
			for (uint32_t k = 0; k < 3; k++) {
				for (auto ntid : mf.tri.Vs[vs[k]].neighbor_fs) {
					if (t_flag[ntid]) continue;
					t_flag[ntid] = true;
					ts_.push_back(ntid);
				}
			}
		}
		ts.clear();
		ts.swap(ts_);
	}
	return false;
}

void point_line_projection(const Vector3d &v1, const Vector3d &v2, const Vector3d &v, Vector3d &pv, double &t)
{
	Vector3d vv1 = v - v1, v21 = v2 - v1;
	double nv21_2 = v21.squaredNorm();
	if (nv21_2 >= Precision)
		t = vv1.dot(v21) / nv21_2;
	else
		t = 0;
	if (t >= 0.0 && t <= 1.0) pv = v1 + t * (v2 - v1);
	else if (t < 0.0)pv = v1;
	else if (t > 1.0)pv = v2;
}
void projectPointOnQuad(const vector<Vector3d>& quad_vs, vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN)
{
	Eigen::Matrix<double, 3, 2> jacobian;

	uv = Vector2d::Constant(0.5f);

	Vector3d F;

	for (int i = 0; i < 4; ++i)
	{
		interpolP = bilinear(quad_vs[0], quad_vs[1], quad_vs[2], quad_vs[3], uv);
		interpolN = bilinear(vs_normals[0], vs_normals[1], vs_normals[2], vs_normals[3], uv);

		Vector3d dPdu = (1 - uv.y()) * quad_vs[1] + uv.y() * quad_vs[2] - ((1 - uv.y()) * quad_vs[0] + uv.y() * quad_vs[3]);
		Vector3d dPdv = (1 - uv.x()) *quad_vs[3] + uv.x() * quad_vs[2] - ((1 - uv.x()) * quad_vs[0] + uv.x() * quad_vs[1]);
		Vector3d dNdu = (1 - uv.y()) * vs_normals[1] + uv.y() * vs_normals[2] - ((1 - uv.y()) * vs_normals[0] + uv.y() * vs_normals[3]);
		Vector3d dNdv = (1 - uv.x()) * vs_normals[3] + uv.x() * vs_normals[2] - ((1 - uv.x()) * vs_normals[0] + uv.x() * vs_normals[1]);

		F = (p - interpolP).cross(interpolN);
		Vector3d dFdu = (-dPdu).cross(interpolN) + (p - interpolP).cross(dNdu);
		Vector3d dFdv = (-dPdv).cross(interpolN) + (p - interpolP).cross(dNdv);

		jacobian.col(0) = dFdu;
		jacobian.col(1) = dFdv;

		Vector2d rhs = -jacobian.transpose() * F;
		auto lhs = jacobian.transpose() * jacobian;
		float norm = 1.0f / (lhs(0, 0) * lhs(1, 1) - lhs(0, 1) * lhs(1, 0));

		uv += Vector2d(lhs(1, 1) * rhs.x() - lhs(0, 1) * rhs.y(), -lhs(1, 0) * rhs.x() + lhs(0, 0) * rhs.y()) * norm;
	}

	interpolP = bilinear(quad_vs[0], quad_vs[1], quad_vs[2], quad_vs[3], uv);
	interpolN = bilinear(vs_normals[0], vs_normals[1], vs_normals[2], vs_normals[3], uv);
}
void projectPointOnTriangle(const vector<Vector3d>& tri_vs, const vector<Vector3d> & vs_normals, const Vector3d& p, Vector2d& uv, Vector3d& interpolP, Vector3d& interpolN)
{
	Eigen::Matrix<double, 3, 2> jacobian;

	uv = Vector2d::Constant(0.333f);

	Vector3d F;

	Vector3d dPdu = tri_vs[0] - tri_vs[2];
	Vector3d dPdv = tri_vs[1] - tri_vs[2];
	Vector3d dNdu = vs_normals[0] - vs_normals[2];
	Vector3d dNdv = vs_normals[1] - vs_normals[2];

	for (int i = 0; i < 20; ++i)
	{
		interpolP = barycentric(tri_vs[0], tri_vs[1], tri_vs[2], uv);
		interpolN = barycentric(vs_normals[0], vs_normals[1], vs_normals[2], uv);

		F = (p - interpolP).cross(interpolN);
		Vector3d dFdu = (-dPdu).cross(interpolN) + (p - interpolP).cross(dNdu);
		Vector3d dFdv = (-dPdv).cross(interpolN) + (p - interpolP).cross(dNdv);

		jacobian.col(0) = dFdu;
		jacobian.col(1) = dFdv;

		Vector2d rhs = -jacobian.transpose() * F;
		auto lhs = jacobian.transpose() * jacobian;
		float norm = 1.0f / (lhs(0, 0) * lhs(1, 1) - lhs(0, 1) * lhs(1, 0));

		uv += Vector2d(lhs(1, 1) * rhs.x() - lhs(0, 1) * rhs.y(), -lhs(1, 0) * rhs.x() + lhs(0, 0) * rhs.y()) * norm;
	}

	interpolP = barycentric(tri_vs[0], tri_vs[1], tri_vs[2], uv);
	interpolN = barycentric(vs_normals[0], vs_normals[1], vs_normals[2], uv);
}
template <typename T>
T bilinear(const T& v1, const T& v2, const T& v3, const T& v4, const Vector2d& uv){
	return (1 - uv.x()) * ((1 - uv.y()) * v1 + uv.y() * v4) + uv.x() * ((1 - uv.y()) * v2 + uv.y() * v3);
}
template <typename T>
T barycentric(const T& v1, const T& v2, const T& v3, const Vector2d& uv)
{
	return uv.x() * v1 + uv.y() * v2 + (1 - uv.x() - uv.y()) * v3;
}
template <typename T>
bool num_equal(const T& x, const T& y, const double &precision) {
	return std::abs(x - y) <= (std::max)(precision, precision * (std::max)(std::abs(x), std::abs(y)));
}

Float rescale(Mesh &mesh, Float scaleI, bool inverse = false) {
	if (!inverse) {
		RowVector3d c; c.setZero();
		for (uint32_t i = 0; i < mesh.V.cols(); i++)
			c += mesh.V.col(i);
		c /= mesh.V.cols();
		for (uint32_t i = 0; i < mesh.V.cols(); i++)
			mesh.V.col(i) -= c;
		Vector3d min_ = mesh.V.rowwise().minCoeff();
		Vector3d max_ = mesh.V.rowwise().maxCoeff();

		double diagonal_local = (max_ - min_).norm();
		double scale = 5.0 / diagonal_local;
		mesh.V *= scale;
		return scale;
	}
	else {
		mesh.V *= scaleI;
	}
	return 1.0;
}
void compute_referenceMesh(MatrixXd &V, vector<Hybrid> &H, vector<uint32_t> &Hs, vector<MatrixXd> &Vout) {
	vector<MatrixXd>().swap(Vout);

	for (auto hid : Hs) {
		vector<MatrixXd> vout;
		hex2cuboid(V, H[hid].vs, vout);
		Vout.insert(Vout.end(),vout.begin(), vout.end());
	}
}
void hex2cuboid(MatrixXd &V, vector<uint32_t> &vs, vector<MatrixXd> &vout) {
	double volume = 0;
	hex2tet24(V, vs, volume);

	//three types of edges
	double e0 = 0, e1 = 0, e2 = 0;
	e0 += (V.row(vs[0]) - V.row(vs[1])).norm();
	e0 += (V.row(vs[3]) - V.row(vs[2])).norm();
	e0 += (V.row(vs[4]) - V.row(vs[5])).norm();
	e0 += (V.row(vs[7]) - V.row(vs[6])).norm();
	e0 /= 4;
	e1 += (V.row(vs[0]) - V.row(vs[3])).norm();
	e1 += (V.row(vs[1]) - V.row(vs[2])).norm();
	e1 += (V.row(vs[4]) - V.row(vs[7])).norm();
	e1 += (V.row(vs[5]) - V.row(vs[6])).norm();
	e1 /= 4;
	e2 += (V.row(vs[0]) - V.row(vs[4])).norm();
	e2 += (V.row(vs[1]) - V.row(vs[5])).norm();
	e2 += (V.row(vs[2]) - V.row(vs[6])).norm();
	e2 += (V.row(vs[3]) - V.row(vs[7])).norm();
	e2 /= 4;
	double ratio = std::cbrt(volume / (e0*e1*e2));
	e0 *= ratio;
	e1 *= ratio;
	e2 *= ratio;

	//eight vertices
	MatrixXd v8(8, 3);
	v8 << 0, 0, 0,
		e0, 0, 0,
		e0, e1, 0,
		0, e1, 0,
		0, 0, e2,
		e0, 0, e2,
		e0, e1, e2,
		0, e1, e2;
	//eight tets
	for (uint32_t i = 0; i < 8; i++) {
		MatrixXd tet(4, 3);
		for (uint32_t j = 0; j < 4; j++)
			tet.row(j) = v8.row(hex_tetra_table[i][j]);
		vout.push_back(tet);
	}
}
void hex2tet24(MatrixXd &V, vector<uint32_t> &vs, double & volume) {
	//6 face center
	vector<RowVector3d> fvs(6);
	for (int i = 0; i < 6; i++) {
		fvs[i].setZero();
		for (int j = 0; j < 4; j++) {
			int vid = vs[hex_face_table[i][j]];
			fvs[i] += V.row(vid);
		}
		fvs[i] /= 4;
	}
	//hex center
	RowVector3d hv; hv.setZero();
	for (int i = 0; i < 8; i++) hv += V.row(vs[i]);
	hv /= 8;
	//tet volume
	volume = 0;
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 4; j++) {
			int vid0 = vs[hex_face_table[i][j]], vid1 = vs[hex_face_table[i][(j + 1) % 4]];

			Matrix3d Jacobian;
			Jacobian.col(0) = V.row(vid0) - hv;
			Jacobian.col(1) = V.row(vid1) - hv;
			Jacobian.col(2) = fvs[i] - hv;
			volume += std::abs(Jacobian.determinant());
		}

}

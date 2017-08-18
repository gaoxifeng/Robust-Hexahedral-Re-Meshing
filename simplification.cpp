#include "simplification.h"
#include "timer.h"
void simplification::pipeline() {

	vector<unsigned long long> timings;
	timings.push_back(0);

	unsigned long long timer0, timer1, timer2, timer3, timer4;
	timer0 = timer1 = timer2 = timer3 = timer4 = 0;
	Timer<> timer;
	timer.beginStage("looping");
	uint32_t removed_candidates = 0;
	Mesh_Quality mq;

	while (true) {

		Timer<> timer0_, timer1_, timer2_, timer3_, timer4_;

		scaled_jacobian(mesh, mq);

		if (mq.min_Jacobian < Jacobian_Bound) break;
		cout << "to remove " << removed_candidates + 1 << endl;
		//stopping criterions
		if (Remove_Iteration != 0 && removed_candidates >= Remove_Iteration) break;
		if (remove_cuboid_ratio != 0 && (double)(cuboid_num_original - frame.FHs.size()) / cuboid_num_original >= remove_cuboid_ratio) break;

		if (!remove()) {
			break; }
		timer.endStage("end removing");
		timer0 = timer.value();
		timings.push_back(timer0);

		extract();
		ranking();
		subdivision();

		double remove_cs_ratio = (double)(cuboid_num_original - frame.FHs.size()) / cuboid_num_original;
		removed_candidates++;
	}

	scaled_jacobian(mesh, mq);
	cout << "minimum scaled J: " << mq.min_Jacobian << " average scaled J: " << mq.ave_Jacobian << endl;
	std::cout << "B_V B_E B_F B_H: " << frame.FVs.size() << " " << frame.FEs.size() << " " << frame.FFs.size() << " " << frame.FHs.size() << endl;
	std::cout << "#sheets removed: " << (double)(sheet_num_original - All_Sheets.size()) << endl;
	std::cout << "Removed Component Ratio: " << (double)(cuboid_num_original - frame.FHs.size()) / cuboid_num_original << endl;

	char path[300];

	optimization();
	sprintf(path, "%s%s", path_out, "_simplified_opt.vtk");
	io.write_hybrid_mesh_VTK(mesh, path);
	cout << "Structure Simplification Finished!" << endl;
	timer.endStage("end Looping");
	std::cout << "timing: " << timer.value()<<"ms"<< endl;

}
bool simplification::initialize() {
	//topology information
	topology_info(mesh, frame, mt);
	if (mt.manifoldness_problem) {
		cout << "non-manifold input, please double-check!" << endl; return false;
	}
	else {
		cout << "surface euler: " << mt.surface_euler << endl;
		cout << "volume euler: " << mt.volume_euler << endl;
	}

	INVALID_V = (uint32_t)-1;
	INVALID_E = (uint32_t)-1;

	std::cout << "Start Structure Simplification..." << endl;

	if(SHARP_FEATURE) mf.angle_threshold = 140.0/ 180 * PAI;

	triangle_mesh_feature(mf, mesh);
	initial_feature(mf, fc, mesh);

	extract();
	ranking();

	cuboid_num_original = frame.FHs.size();
	sheet_num_original = All_Sheets.size();
	//Slim_global_region = Slim_region * 2;
	Slim_global_region = Slim_region;
	return true;
}
void simplification::extract() {
	std::vector<Sheet>().swap(All_Sheets);
	std::vector<bool> e_flag(frame.FEs.size(), false);
	while (true) {
		uint32_t eid = INVALID_E;
		for (uint32_t i = 0; i < frame.FEs.size(); i++) if (!e_flag[i]) { eid = i; break; }
		if (eid == INVALID_E) break;

		Sheet sheet; sheet.id = All_Sheets.size();

		std::queue<uint32_t> e_pool; e_pool.push(eid);
		while (!e_pool.empty()) {
			uint32_t eid = e_pool.front(); e_pool.pop();
			if (e_flag[eid]) continue; e_flag[eid] = true;
			sheet.middle_es.push_back(eid);
			for (uint32_t i = 0; i < frame.FEs[eid].neighbor_ffs.size(); i++) {
				uint32_t fid = frame.FEs[eid].neighbor_ffs[i];
				uint32_t pos = std::find(frame.FFs[fid].es.begin(), frame.FFs[fid].es.end(), eid) - frame.FFs[fid].es.begin();
				uint32_t op_eid = frame.FFs[fid].es[(pos+2)%4];
				if (e_flag[op_eid]) continue;
				e_pool.push(op_eid);
			}
		}

		All_Sheets.push_back(sheet);
	}
	for (uint32_t i = 0; i < All_Sheets.size(); i++) {
		build_sheet_info(i);
	}

	std::vector<CHord>().swap(All_Chords);
	vector<bool> f_flag(frame.FFs.size(),false);
	while (true)
	{
		uint32_t fid = -1;
		for (uint32_t i = 0; i<f_flag.size(); i++){
			if (!f_flag[i]) { fid = i; break; }
		}
		if (fid == -1) break;

		CHord cc = extract_chord(fid, f_flag);
		cc.id = All_Chords.size();
		cc.side = 0;
		All_Chords.push_back(cc);
		cc.id = All_Chords.size();
		cc.side = 1;
		All_Chords.push_back(cc);
	}
	for (uint32_t i = 0; i<All_Chords.size(); i++) {
		build_chord_info(i);
	}
}
bool simplification::build_sheet_info(uint32_t sheet_id){
	All_Sheets[sheet_id].fake = false;
	All_Sheets[sheet_id].type = Sheet_type::open;
	vector<uint32_t> &cs = All_Sheets[sheet_id].cs;
	vector<uint32_t> &middle_fs = All_Sheets[sheet_id].middle_fs;
	vector<uint32_t> &middle_es = All_Sheets[sheet_id].middle_es;
	vector<uint32_t> &middle_es_b = All_Sheets[sheet_id].middle_es_b;
	for (auto eid : middle_es) {
		cs.insert(cs.end(), frame.FEs[eid].neighbor_fhs.begin(), frame.FEs[eid].neighbor_fhs.end());
		sort(cs.begin(), cs.end());
		cs.erase(unique(cs.begin(), cs.end()), cs.end());

		middle_fs.insert(middle_fs.end(), frame.FEs[eid].neighbor_ffs.begin(), frame.FEs[eid].neighbor_ffs.end());
		sort(middle_fs.begin(), middle_fs.end());
		middle_fs.erase(unique(middle_fs.begin(), middle_fs.end()), middle_fs.end());

		if (frame.FEs[eid].boundary) middle_es_b.push_back(eid);
	}
	if(!middle_es_b.size()) All_Sheets[sheet_id].type = Sheet_type::close;
	vector<short> V_flag(frame.FVs.size(), 0);
	for (auto eid : middle_es) {
		uint32_t v1 = frame.FEs[eid].vs[0];
		uint32_t v2 = frame.FEs[eid].vs[1];
		V_flag[v1]++;
		V_flag[v2]++;

		if (V_flag[v1] > 1 || V_flag[v2] > 1) {
			All_Sheets[sheet_id].type = Sheet_type::tagent;
			break;
		}
	}
	sort(middle_fs.begin(), middle_fs.end());
	for (auto cid:cs){
		bool all_middle = true;
		vector<uint32_t> fs = frame.FHs[cid].fs;
		sort(fs.begin(), fs.end());
		vector<uint32_t> mfs;
		set_intersection(fs.begin(), fs.end(), middle_fs.begin(), middle_fs.end(), back_inserter(mfs));

		if (!mfs.size()) {
			All_Sheets[sheet_id].type = Sheet_type::intersect;
			break;
		}
	}
	vector<uint32_t> &fs = All_Sheets[sheet_id].fs;
	vector<uint32_t> &left_fs = All_Sheets[sheet_id].left_fs;
	vector<uint32_t> &right_fs = All_Sheets[sheet_id].right_fs;
	
	for (auto cid : cs) fs.insert(fs.end(), frame.FHs[cid].fs.begin(), frame.FHs[cid].fs.end());
	sort(fs.begin(), fs.end()); fs.erase(unique(fs.begin(), fs.end()), fs.end());

	vector<uint32_t> &twoside_fs = All_Sheets[sheet_id].side_fs;

	vector<bool> F_flag(frame.FFs.size(),true);
	for (auto fid : middle_fs)F_flag[fid] = false;
	for (auto fid : fs) if(F_flag[fid]) twoside_fs.push_back(fid);
	fill(F_flag.begin(), F_flag.end(), false);
	for (auto fid : twoside_fs)F_flag[fid] = true;

	if (!twoside_fs.size()) { cout << "ERROR, no side fs" << endl; system("PAUSE"); }

	queue<uint32_t> fs_pool;
	fs_pool.push(twoside_fs[0]);
	while (!fs_pool.empty()){
		auto fid = fs_pool.front(); fs_pool.pop();
		if (!F_flag[fid]) continue; F_flag[fid] = false;
		left_fs.push_back(fid);

		for (auto eid: frame.FFs[fid].es)
			for (auto nfid : frame.FEs[eid].neighbor_ffs)
				if (F_flag[nfid]) fs_pool.push(nfid);
	}
	for (auto fid : twoside_fs)if(F_flag[fid]) right_fs.push_back(fid);
	if (All_Sheets[sheet_id].type == Sheet_type::tagent || All_Sheets[sheet_id].type == Sheet_type::intersect)
		return true;
	else if (left_fs.size() != right_fs.size()) { All_Sheets[sheet_id].type = Sheet_type::mobius; return true; }
	vector<uint32_t> &left_es = All_Sheets[sheet_id].left_es;
	vector<uint32_t> &right_es = All_Sheets[sheet_id].right_es;
	
	for (auto fid : left_fs) left_es.insert(left_es.end(), frame.FFs[fid].es.begin(), frame.FFs[fid].es.end());
	sort(left_es.begin(), left_es.end()); left_es.erase(unique(left_es.begin(), left_es.end()), left_es.end());
	for (auto fid : right_fs) right_es.insert(right_es.end(), frame.FFs[fid].es.begin(), frame.FFs[fid].es.end());
	sort(right_es.begin(), right_es.end()); right_es.erase(unique(right_es.begin(), right_es.end()), right_es.end());
	bool left = false, right = false, bl = true, br = true;
	for (auto eid:left_es){
		int size = frame.FEs[eid].neighbor_fhs.size();
		if ((!frame.FEs[eid].boundary && size != Interior_RegularE) || (frame.FEs[eid].boundary  && size != Boundary_RegularE))
			left = true;
		if (!frame.FEs[eid].boundary) bl = false;
	}
	for (auto eid : right_es) {
		int size = frame.FEs[eid].neighbor_fhs.size();
		if ((!frame.FEs[eid].boundary && size != Interior_RegularE) || (frame.FEs[eid].boundary  && size != Boundary_RegularE))
			right = true;
		if (!frame.FEs[eid].boundary) br = false;
	}
	if (!(left&&right) && !bl && !br) All_Sheets[sheet_id].fake = true;

	return true;
}
CHord simplification::extract_chord(uint32_t &fid, vector<bool> &f_flag) {
	CHord cc;
	cc.type = Sheet_type::open; cc.fake = false;

	uint32_t cur_f = fid;
	if (!frame.FFs[cur_f].neighbor_fhs.size()) { cout << "isolated patch" << endl; system("PAUSE"); }
	uint32_t cur_h = frame.FFs[cur_f].neighbor_fhs[0];
	vector<uint32_t> fs, fs_right;
	fs.push_back(cur_f);
	while (true) {
		Frame_F f_ = frame.FFs[cur_f];
		Frame_H h_ = frame.FHs[cur_h];

		std::vector<uint32_t> &vs = f_.vs;
		std::sort(vs.begin(), vs.end());
		short cors_f = -1;
		for (uint32_t j = 0; j < 6; j++) {
			std::vector<uint32_t> vsj = frame.FFs[h_.fs[j]].vs;
			std::sort(vsj.begin(), vsj.end());
			std::vector<uint32_t> common_vs;
			std::set_intersection(vs.begin(), vs.end(), vsj.begin(), vsj.end(), std::back_inserter(common_vs));
			if (common_vs.size())continue; else { cors_f = j; break; }
		}
		if(cors_f == -1) {
			cout << "problematic cuboid" << endl; system("PAUSE"); 
		}
		uint32_t neighbor_f = h_.fs[cors_f];
		uint32_t neighbor_h = -1;
		Frame_F f = frame.FFs[neighbor_f];

		for (auto fhid: f.neighbor_fhs) if (fhid != cur_h) neighbor_h = fhid;

		cur_f = neighbor_f;
		cur_h = neighbor_h;

		if (find(fs.begin(),fs.end(), cur_f) != fs.end()){
			cc.type = Sheet_type::close;
			break;
		}
		fs.push_back(cur_f);
		if (cur_h == -1) break;
	}

	cur_f = fid;
	if (cc.type != Sheet_type::close && frame.FFs[cur_f].neighbor_fhs.size() > 1) {
		cur_h = frame.FFs[cur_f].neighbor_fhs[1];

		while (true) {
			Frame_F f_ = frame.FFs[cur_f];
			Frame_H h_ = frame.FHs[cur_h];

			std::vector<uint32_t> &vs = f_.vs;
			std::sort(vs.begin(), vs.end());
			short cors_f = -1;
			for (uint32_t j = 0; j < 6; j++) {
				std::vector<uint32_t> vsj = frame.FFs[h_.fs[j]].vs;
				std::sort(vsj.begin(), vsj.end());
				std::vector<uint32_t> common_vs;
				std::set_intersection(vs.begin(), vs.end(), vsj.begin(), vsj.end(), std::back_inserter(common_vs));
				if (common_vs.size())continue; else { cors_f = j; break; }
			}
			if (cors_f == -1) { 
				cout << "problematic cuboid" << endl; system("PAUSE"); 
			}
			uint32_t neighbor_f = h_.fs[cors_f];

			if (find(fs.begin(), fs.end(), neighbor_f) != fs.end()) {
				cout << "problematic cuboid" << endl; system("PAUSE"); 
			}

			uint32_t neighbor_h = -1;
			Frame_F f = frame.FFs[neighbor_f];

			for (auto fhid : f.neighbor_fhs) if (fhid != cur_h) neighbor_h = fhid;

			cur_f = neighbor_f;
			cur_h = neighbor_h;
			fs_right.push_back(cur_f);

			if (cur_h == -1) break;
		}
	}
	reverse(fs_right.begin(), fs_right.end());
	cc.parallel_fs = fs_right;
	cc.parallel_fs.insert(cc.parallel_fs.end(), fs.begin(), fs.end());
	for (auto fid : cc.parallel_fs) { f_flag[fid] = true; sort(frame.FFs[fid].neighbor_fhs.begin(), frame.FFs[fid].neighbor_fhs.end()); }

	vector<uint32_t> v_di = frame.FFs[cc.parallel_fs[0]].vs;
	for (uint32_t i = 0; i<4; i++){ cc.parallel_ns[i].push_back(v_di[i]); }
	for (uint32_t i = 1; i<cc.parallel_fs.size(); i++) {
		vector<uint32_t> sharedhs;
		set_intersection(frame.FFs[cc.parallel_fs[i - 1]].neighbor_fhs.begin(), frame.FFs[cc.parallel_fs[i - 1]].neighbor_fhs.end(),
			frame.FFs[cc.parallel_fs[i]].neighbor_fhs.begin(), frame.FFs[cc.parallel_fs[i]].neighbor_fhs.end(), 
			back_inserter(sharedhs));
		cc.cs.push_back(sharedhs[0]);
		if (cc.type == Sheet_type::close && i == cc.parallel_fs.size() - 1){
			sharedhs.clear();
			set_intersection(frame.FFs[cc.parallel_fs[0]].neighbor_fhs.begin(), frame.FFs[cc.parallel_fs[0]].neighbor_fhs.end(),
				frame.FFs[cc.parallel_fs[i]].neighbor_fhs.begin(), frame.FFs[cc.parallel_fs[i]].neighbor_fhs.end(),
				back_inserter(sharedhs));
			cc.cs.push_back(sharedhs[0]);
		}
		vector<uint32_t> &v_di_cur = frame.FFs[cc.parallel_fs[i]].vs;

		int cid = cc.cs[i - 1];
		int v0 = v_di[0], v1 = v_di[1];

		for (int j = 0; j < 4; j++) {
			int v2 = v_di_cur[j], v3 = v_di_cur[(j + 1)%4];
			vector<uint32_t> vs0; 
			vs0.push_back(v0);
			vs0.push_back(v1);
			vs0.push_back(v2);
			vs0.push_back(v3);
			sort(vs0.begin(), vs0.end());
			int whichf = -1;
			for (auto fid : frame.FHs[cid].fs) {
				vector<uint32_t> vs1 = frame.FFs[fid].vs;
				sort(vs1.begin(), vs1.end());
				if (std::equal(vs0.begin(), vs0.end(), vs1.begin())) {whichf = fid; break;}
			}
			if (whichf != -1) {
				vector<uint32_t> &nvs_pre = frame.FVs[v0].neighbor_fvs;
				if (find(nvs_pre.begin(), nvs_pre.end(), v2) != nvs_pre.end())
					for (int k = 0; k < 4; k++) cc.parallel_ns[k].push_back(v_di_cur[(j + k) % 4]);
				else for (int k = 0; k < 4; k++) cc.parallel_ns[k].push_back(v_di_cur[(j - k + 5) % 4]);

				for (int k = 0; k < 4; k++) v_di[k] = cc.parallel_ns[k][cc.parallel_ns[k].size() - 1];
				break;
			}
		}
	}
	return cc;
}
bool simplification::build_chord_info(uint32_t id) {
	CHord &cc = All_Chords[id];
	cc.fake = false;

	for (uint32_t i = 0; i < 4; i++)
		for (auto vid : cc.parallel_ns[i]) {
			sort(frame.FVs[vid].neighbor_fes.begin(), frame.FVs[vid].neighbor_fes.end());
			sort(frame.FVs[vid].neighbor_ffs.begin(), frame.FVs[vid].neighbor_ffs.end());
		}

	for (uint32_t i = 0; i<4; i++){
		for (uint32_t j = 0; j < cc.parallel_ns[i].size(); j++) {
			uint32_t v0 = cc.parallel_ns[i][j], v1 = cc.parallel_ns[(i + 1) % 4][j];
			
			int fid;
			if (cc.type == Sheet_type::close && j == cc.parallel_fs.size() - 1) 
				fid = cc.parallel_fs[0];
			else fid = cc.parallel_fs[j];
			for(uint32_t k=0;k<4;k++)
				if ((v0 == frame.FFs[fid].vs[k] && v1 == frame.FFs[fid].vs[(k + 1) % 4]) || (v1 == frame.FFs[fid].vs[k] && v0 == frame.FFs[fid].vs[(k + 1) % 4])) {
					cc.parallel_es[i].push_back(frame.FFs[fid].es[k]);
					break;
				}

			if (j > 0) {
				vector<uint32_t> vs;
				vs.push_back(cc.parallel_ns[i][j]);
				vs.push_back(cc.parallel_ns[(i + 1) % 4][j]);
				vs.push_back(cc.parallel_ns[(i + 1) % 4][j - 1]);
				vs.push_back(cc.parallel_ns[i][j - 1]);

				int hid = cc.cs[j-1]; fid = -1;
				for (uint32_t k = 0; k < 6; k++) {
					fid = frame.FHs[hid].fs[k];
					vector<uint32_t> &vs1 = frame.FFs[fid].vs;
					bool equal = true;
					for (auto vid1 : vs1) {
						bool found = false;
						for (auto vid0 : vs) { if (vid1 == vid0) { found = true; break; } }
						if (!found) { equal = false; break; }
					}
					if (equal) break;
				}
				if (fid == -1) { cout << "error in build chord info" << endl; system("PAUSE"); }
				cc.vertical_fs[i].push_back(fid);

				for (uint32_t k = 0; k<4; k++)
					if ((vs[0] == frame.FFs[fid].vs[k] && vs[3] == frame.FFs[fid].vs[(k + 1) % 4]) || (vs[3] == frame.FFs[fid].vs[k] && vs[0] == frame.FFs[fid].vs[(k + 1) % 4])) {
						cc.vertical_es[i].push_back(frame.FFs[fid].es[k]);
						break;
					}
			}
			if (cc.type == Sheet_type::close && j == cc.parallel_ns[i].size() - 1) {
				vector<uint32_t> vs;
				vs.push_back(cc.parallel_ns[i][j]);
				vs.push_back(cc.parallel_ns[(i + 1) % 4][j]);
				vs.push_back(cc.parallel_ns[(i + 1) % 4][0]);
				vs.push_back(cc.parallel_ns[i][0]);

				int hid = cc.cs[j]; fid = -1;
				for (uint32_t k = 0; k < 6; k++) {
					fid = frame.FHs[hid].fs[k];
					vector<uint32_t> &vs1 = frame.FFs[fid].vs;
					bool equal = true;
					for (auto vid1 : vs1) {
						bool found = false;
						for (auto vid0 : vs) { if (vid1 == vid0) { found = true; break; } }
						if (!found) { equal = false; break; }
					}
					if (equal) break;
				}
				if (fid == -1) { cout << "error in build chord info" << endl; system("PAUSE"); }
				cc.vertical_fs[i].push_back(fid);

				for (uint32_t k = 0; k<4; k++)
					if ((vs[0] == frame.FFs[fid].vs[k] && vs[3] == frame.FFs[fid].vs[(k + 1) % 4]) || (vs[3] == frame.FFs[fid].vs[k] && vs[0] == frame.FFs[fid].vs[(k + 1) % 4])) {
						cc.vertical_es[i].push_back(frame.FFs[fid].es[k]);
						break;
					}
			}
		}
	}
	cc.fs = cc.parallel_fs;
	for (uint32_t i = 0; i<4; i++) cc.fs.insert(cc.fs.end(),cc.vertical_fs[i].begin(), cc.vertical_fs[i].end());
	for (uint32_t i = 0; i<4; i++) cc.ns.insert(cc.ns.end(), cc.parallel_ns[i].begin(), cc.parallel_ns[i].end());
	for (uint32_t i = 0; i<4; i++) cc.es.insert(cc.es.end(), cc.parallel_es[i].begin(), cc.parallel_es[i].end());
	for (uint32_t i = 0; i<4; i++) cc.es.insert(cc.es.end(), cc.vertical_es[i].begin(), cc.vertical_es[i].end());
	vector<bool> Vs_Ind(frame.FVs.size(), false), Es_Ind(frame.FEs.size(), false), Fs_Ind(frame.FFs.size(), false), Hs_Ind(frame.FHs.size(), false);

	for (auto vid : cc.ns) {
		if (Vs_Ind[vid]) cc.tangent_vs.push_back(vid); Vs_Ind[vid] = true;
	}
	for (auto eid : cc.es) {
		if (Es_Ind[eid]) cc.tangent_es.push_back(eid); Es_Ind[eid] = true;
	}
	for (auto cid : cc.cs) {
		if (Hs_Ind[cid]) cc.tangent_cs.push_back(cid); Hs_Ind[cid] = true;
	}
	for (uint32_t i = 0; i<4; i++){
		for (auto fid : cc.vertical_fs[i]) {
			if (Fs_Ind[fid]) cc.tangent_fs.push_back(fid); Fs_Ind[fid] = true;
		}
	}
	if (cc.tangent_vs.size() || cc.tangent_es.size() || cc.tangent_fs.size() || cc.tangent_cs.size())
		cc.type = Sheet_type::tagent;

	vector<uint32_t> sheet_ids;
	for (int i = 0; i < 2; i++) {
		uint32_t feid = cc.parallel_es[i][0];
		int sheet_id = -1;
		for (auto sheet : All_Sheets)
			if (find(sheet.middle_es.begin(), sheet.middle_es.end(), feid) != sheet.middle_es.end()) {
				sheet_id = sheet.id;
				if (sheet.fake) cc.fake = true;
				sheet_ids.push_back(sheet.id); break;
			}
		if (sheet_id == -1) { cout << "doesnot find the longer sheet" << endl; system("PAUSE"); }
	}

	return true;
}

void simplification::ranking() {
	Candidates.clear();
	char path[300];

	for (uint32_t i = 0; i < All_Sheets.size(); i++) {
		Tuple_Candidate tc;
		get<0>(tc) = i; 
		get<1>(tc) = Base_Set::SHEET;
		sheet_chord_weight(tc);
		Candidates.push_back(tc);
	}
	for (uint32_t i = 0; i < All_Chords.size(); i++) {
		Tuple_Candidate tc;
		get<0>(tc) = i;
		get<1>(tc) = Base_Set::CHORD;
		sheet_chord_weight(tc);
		
		Candidates.push_back(tc);
	}

	std::function<bool(Tuple_Candidate &, Tuple_Candidate &)> rank = [&](
		Tuple_Candidate &s1, Tuple_Candidate  &s2)->bool {
		return get<2>(s1) < get<2>(s2);
	};
	std::sort(Candidates.begin(), Candidates.end(), rank);
}
void simplification::sheet_chord_weight(Tuple_Candidate &c){

	const Float LARGE_NUM = 1.e+5;

	uint32_t id = get<0>(c);
	Float &weight = get<2>(c);
	if (get<1>(c) == Base_Set::SHEET) {
		//edge_volume ratio
		float volume = 0;
		Float edge_len = 0;
		for (auto fe:All_Sheets[id].middle_es){
			Float cur_len = 0;
			for (uint32_t k = 0; k<frame.FEs[fe].es_link.size(); k++){
				uint32_t eid = frame.FEs[fe].es_link[k];
				uint32_t v0 = mesh.Es[eid].vs[0], v1 = mesh.Es[eid].vs[1];
				cur_len += (mesh.V.col(v0) - mesh.V.col(v1)).norm();
			}
			edge_len += cur_len;
		}
		edge_len /= All_Sheets[id].middle_es.size();
		All_Sheets[id].weight_len = edge_len;
		All_Sheets[id].weight = All_Sheets[id].weight_len;

		//valence 
		//Es_neighborhood
		vector<vector<uint32_t>> Es_neighborhood(frame.FEs.size()), Es_group(frame.FEs.size()), Es_group_temp;
		vector<bool> E_flag(frame.FEs.size(), false);
		for (auto fid : All_Sheets[id].middle_fs) {
			for (uint32_t i = 0; i < 2; i++) {
				uint32_t eid = frame.FFs[fid].es[i];
				if (find(All_Sheets[id].middle_es.begin(), All_Sheets[id].middle_es.end(), eid) == All_Sheets[id].middle_es.end()) {
					uint32_t eid1 = frame.FFs[fid].es[(i + 2) % 4];
					Es_neighborhood[eid].push_back(eid1);
					Es_neighborhood[eid1].push_back(eid);
				}
			}
		}
		//Es_group
		for (uint32_t i = 0; i < Es_neighborhood.size(); i++) {
			if (!Es_neighborhood[i].size() || E_flag[i]) continue;
			E_flag[i] = true; 
			Es_group[i].push_back(i);
			vector<uint32_t> pool = Es_neighborhood[i];
			while (pool.size()){
				uint32_t eid = pool[pool.size() - 1]; pool.pop_back();
				if (E_flag[eid]) continue;
				E_flag[eid] = true;
				Es_group[i].push_back(eid);
				pool.insert(pool.end(), Es_neighborhood[eid].begin(), Es_neighborhood[eid].end());
			}
		}
		for (uint32_t i = 0; i<Es_group.size(); i++){
			if (!Es_group[i].size()) continue;
			Es_group_temp.push_back(Es_group[i]);
		}
		Es_group.swap(Es_group_temp);
		if (!Es_group.size()) { std::cout << "Sheet Error" << endl; system("PAUSE"); }
		
		All_Sheets[id].weight_val_average = 0;
		All_Sheets[id].weight_val_max = 0;
		All_Sheets[id].weight_val_min = (numeric_limits<Float>::max)();
		All_Sheets[id].valence_filter = false;


		vector<bool> C_flag(frame.FHs.size(), true);
		for (auto cid : All_Sheets[id].cs) C_flag[cid] = false;
		for (auto es : Es_group){
			vector<int> comps_result, comps_after;
			for (auto eid:es) comps_after.insert(comps_after.end(), frame.FEs[eid].neighbor_fhs.begin(), frame.FEs[eid].neighbor_fhs.end());
			sort(comps_after.begin(), comps_after.end()); comps_after.erase(unique(comps_after.begin(), comps_after.end()), comps_after.end());
			for (auto cid : comps_after) if (C_flag[cid])comps_result.push_back(cid);

			int min_val = (numeric_limits<int>::max)();
			for (auto eid : es) if (min_val > frame.FEs[eid].neighbor_fhs.size()) min_val = frame.FEs[eid].neighbor_fhs.size();
			if (comps_result.size() < min_val) {
				All_Sheets[id].valence_filter = true;
			}
		}
		weight = All_Sheets[id].weight;
	}
	else if(get<1>(c) == Base_Set::CHORD) {
		CHord &cc = All_Chords[id];
		Float diagonal_len = 0;
		for (uint32_t i = 0; i < cc.parallel_ns[0].size();i++) {
			uint32_t v0, v1; Vector3d vv0, vv1;
			if (cc.side == 0) {
				v0 = frame.FVs[cc.parallel_ns[1][i]].hid;
				v1 = frame.FVs[cc.parallel_ns[3][i]].hid;
			}
			else {
				v0 = frame.FVs[cc.parallel_ns[0][i]].hid;
				v1 = frame.FVs[cc.parallel_ns[2][i]].hid;
			}
			vv0 = mesh.V.col(v0); vv1 = mesh.V.col(v1);
			diagonal_len += (vv0 - vv1).norm();
		}
		diagonal_len *= 1.0 / cc.parallel_ns[0].size();		

		All_Chords[id].weight_len = diagonal_len;
		All_Chords[id].weight = All_Chords[id].weight_len;

		vector<vector<uint32_t>> Es_group;
		//grouping
		for (uint32_t i = 0; i<cc.vertical_es[0].size(); i++){
			vector<uint32_t> a_group;
			if (cc.side == 0){
				a_group.push_back(cc.vertical_es[0][i]);
				Es_group.push_back(a_group);
				a_group.clear();
				a_group.push_back(cc.vertical_es[1][i]);
				a_group.push_back(cc.vertical_es[3][i]);
				Es_group.push_back(a_group);
				a_group.clear();
				a_group.push_back(cc.vertical_es[2][i]);
				Es_group.push_back(a_group);
			}
			else {
				a_group.push_back(cc.vertical_es[0][i]);
				Es_group.push_back(a_group);
				a_group.clear();
				a_group.push_back(cc.vertical_es[2][i]);
				a_group.push_back(cc.vertical_es[1][i]);
				Es_group.push_back(a_group);
				a_group.clear();
				a_group.push_back(cc.vertical_es[3][i]);
				Es_group.push_back(a_group);
			}
		}
		for (int i = 0; i<cc.parallel_es[0].size(); i++){
			vector<uint32_t> a_group;
			if (cc.side == 0){
				a_group.push_back(cc.parallel_es[0][i]);
				a_group.push_back(cc.parallel_es[3][i]);
				Es_group.push_back(a_group);
				a_group.clear();
				a_group.push_back(cc.parallel_es[1][i]);
				a_group.push_back(cc.parallel_es[2][i]);
				Es_group.push_back(a_group);
			}
			else{
				a_group.push_back(cc.parallel_es[0][i]);
				a_group.push_back(cc.parallel_es[1][i]);
				Es_group.push_back(a_group);
				a_group.clear();
				a_group.push_back(cc.parallel_es[2][i]);
				a_group.push_back(cc.parallel_es[3][i]);
				Es_group.push_back(a_group);
			}
		}

		All_Chords[id].weight_val_average = 0;
		All_Chords[id].weight_val_max = 0;
		All_Chords[id].weight_val_min = (numeric_limits<Float>::max)();
		All_Chords[id].valence_filter = false;
		int interior_ = 0;

		Float weight_val = 0, yita = 6; uint32_t max_valence = 0, min_valence = (numeric_limits<int32_t>::max)();
		vector<Float> weight_subs;
		vector<bool> C_flag(frame.FHs.size(), true);
		for (auto cid : All_Chords[id].cs) C_flag[cid] = false;
		for (uint32_t i = 0; i<Es_group.size(); i++){
			if (!Es_group[i].size()) continue;
			vector<uint32_t> comps_result, comps_after;
			for (uint32_t j = 0; j<Es_group[i].size(); j++)
				comps_after.insert(comps_after.end(), frame.FEs[Es_group[i][j]].neighbor_fhs.begin(), frame.FEs[Es_group[i][j]].neighbor_fhs.end());
			sort(comps_after.begin(), comps_after.end()); comps_after.erase(unique(comps_after.begin(), comps_after.end()), comps_after.end());
			for (auto cid : comps_after) if (C_flag[cid])comps_result.push_back(cid);

			if (max_valence<comps_result.size()) max_valence = comps_result.size();
			if (min_valence>comps_result.size()) min_valence = comps_result.size();

			vector<uint32_t> es_boundary;
			for (auto eid:Es_group[i]) if (frame.FEs[eid].boundary) es_boundary.push_back(eid);

			if (es_boundary.size()) {
				vector<Float> angles; Float average_angle = 0;
				for (auto eid : es_boundary) {
					Float angle = 0, k_ratio = 0;
					dihedral_angle(angle, k_ratio, frame.FEs[eid].neighbor_fhs, eid);
					angles.push_back(angle);
					average_angle += angle;
				}
				average_angle /= es_boundary.size();

				Float after = 2 * average_angle / (PAI*comps_result.size());
				Float after_val = (2 * PAI) / average_angle  * comps_result.size();
				All_Chords[id].weight_val_average += after;
			}
			else{
				interior_++;

				int min_val = (numeric_limits<int>::max)();
				for (auto eid : Es_group[i]) if (min_val > frame.FEs[eid].neighbor_fhs.size()) min_val = frame.FEs[eid].neighbor_fhs.size();
				if (comps_result.size() < min_val) {
					All_Chords[id].valence_filter = true;
				}

				Float after = 4.0 / comps_result.size();
				Float after_val = comps_result.size();
				All_Chords[id].weight_val_average += after;
				if (All_Chords[id].weight_val_max < after_val) All_Chords[id].weight_val_max = after_val;
				if (All_Chords[id].weight_val_min > after_val) All_Chords[id].weight_val_min = after_val;
			}
		}
		All_Chords[id].weight_val_average /= Es_group.size();

		weight = All_Chords[id].weight;

		if (All_Chords[id].weight_val_min < 3.0 || All_Chords[id].weight_val_max > 5.0 || !(All_Chords[id].weight_val_average >= 3.0 && All_Chords[id].weight_val_average <= 5.0))
			weight *= LARGE_NUM;
	}else{
		cout << "ERROR" << endl; system("PAUSE");
	}
}
void simplification::dihedral_angle(Float &angle, Float &k_ratio, vector<uint32_t> &cs, uint32_t eid) {
	vector<Float> angles;
	uint32_t v1 = frame.FEs[eid].vs[0], v2 = frame.FEs[eid].vs[1];
	uint32_t hv1 = frame.FVs[v1].hid;
	uint32_t hv2 = frame.FVs[v2].hid;
	vector<uint32_t> &nfs =frame.FEs[eid].neighbor_ffs;
	for (auto cid: cs){
		vector<uint32_t> &fs = frame.FHs[cid].fs, f2s;
		for (auto fid : nfs) if (find(fs.begin(), fs.end(), fid) != fs.end()) f2s.push_back(fid);

		uint32_t v3, v4;
		for (uint32_t j = 0; j<4; j++){
			if (frame.FFs[f2s[0]].vs[j] != v1&&frame.FFs[f2s[0]].vs[j] != v2) v3 = frame.FFs[f2s[0]].vs[j];
			if (frame.FFs[f2s[1]].vs[j] != v1&&frame.FFs[f2s[1]].vs[j] != v2) v4 = frame.FFs[f2s[1]].vs[j];
		}
		uint32_t hv3 = frame.FVs[v3].hid;
		uint32_t hv4 = frame.FVs[v4].hid;
		Vector3d dir1, dir2, norm1, norm2;
		dir1 = mesh.V.col(hv3) - mesh.V.col(hv1);
		dir2 = mesh.V.col(hv3) - mesh.V.col(hv2);
		norm1 = dir1.cross(dir2); norm1.normalize();
		dir1 = mesh.V.col(hv4) - mesh.V.col(hv1);
		dir2 = mesh.V.col(hv4) - mesh.V.col(hv2);
		norm2 = dir1.cross(dir2); norm2.normalize();
		angles.push_back(PAI - std::acos(norm1.dot(norm2)));
	}
	angle = 0;
	for (uint32_t i = 0; i<angles.size(); i++) angle += angles[i];
	if (angle>PAI * 2) angle = PAI * 2;
	k_ratio = 2 * angle / PAI;
}

bool simplification::remove() {
	std::vector<Tuple_Candidate> Candidates_temp;
	Candidates_temp.insert(Candidates_temp.end(), Candidates.begin() + last_candidate_pos, Candidates.end());
	Candidates_temp.insert(Candidates_temp.end(), Candidates.begin(), Candidates.begin() + last_candidate_pos);
	Candidates_temp.swap(Candidates);

	File_num = 0;
	for (uint32_t i = 0; i < Candidates.size(); i++) {
		uint32_t id = i;
		File_num = id;
		if (!filter_topology_feature(Candidates[id])) continue;
		uint32_t feid;
		if (std::get<1>(Candidates[id]) == Base_Set::SHEET) {
			uint32_t sheet_id = std::get<0>(Candidates[id]);
			feid = All_Sheets[sheet_id].middle_es[0];
		}
		else if (std::get<1>(Candidates[id]) == Base_Set::CHORD) {
			uint32_t chord_id = std::get<0>(Candidates[id]);
			feid = All_Chords[chord_id].parallel_es[0][0];
		}

		Float resolution = frame.FEs[feid].es_link.size();
		Projection_range = resolution + 1;
		Slim_Iteration = resolution * Slim_Iteration_base;
		if (Slim_Iteration > Slim_Iteration_Limit) Slim_Iteration = Slim_Iteration_Limit;

		if (Projection_range > Projection_limit)Projection_range = Projection_limit;

		width_sheet = resolution;

		if (!direct_collapse()) continue;

		if (id >= Candidates.size() - last_candidate_pos)
			last_candidate_pos = id - Candidates.size() + last_candidate_pos;
		else if (id < Candidates.size() - last_candidate_pos)
			last_candidate_pos = id + last_candidate_pos;

		if (last_candidate_pos > All_Sheets.size()*0.3 || last_candidate_pos > 20) last_candidate_pos = 0;

		return true;
	}
	cout << "cannot find candidate anymore" << endl;
	return false;
}
bool simplification::filter_topology_feature(Tuple_Candidate &c) {
	if (!TOPOLOGY && !SHARP_FEATURE) return true;
	uint32_t id = get<0>(c);
	if (get<1>(c) == Base_Set::SHEET) {
		if (All_Sheets[id].valence_filter) {
			return false;}
		if (All_Sheets[id].fake) {
			return false;}
		vector<vector<uint32_t>> candiate_es_links, v_group;
		if (!vs_pair_sheet(id, candiate_es_links, v_group)) { cout << "ERROR: no vs pairs" << endl;  system("PAUSE");}
		if (!target_surface_sheet(id, candiate_es_links, v_group)) {
			return false;
		}

		if (TOPOLOGY) {
			vector<bool> F_flag(frame.FFs.size(), false);
			for (auto fid : All_Sheets[id].middle_fs) F_flag[fid] = true;
			for (auto cid : All_Sheets[id].cs) {
				vector<uint32_t> side_fs;
				for (auto fid : frame.FHs[cid].fs) if (!F_flag[fid]) side_fs.push_back(fid);
				if (side_fs.size() == 2 && frame.FFs[side_fs[0]].boundary && frame.FFs[side_fs[1]].boundary) {
					return false;
				}
			}
		}
	}
	else if (get<1>(c) == Base_Set::CHORD) {
		if (All_Chords[id].type != Sheet_type::close && All_Chords[id].type != Sheet_type::open) return false;
		if (All_Chords[id].valence_filter || All_Chords[id].weight_val_min < 3.0 || All_Chords[id].weight_val_max > 5.0) return false;

		if (All_Chords[id].fake ) {
			return false;
		}
		if (!vs_group_chord(id)) return false;
		if (!target_surface_chord(id)) return false;
	}

	if (!topology_check()) {
		return false;
	}

	return true;
}
bool simplification::vs_pair_sheet(uint32_t sheet_id, vector<vector<uint32_t>> &candiate_es_links, vector<vector<uint32_t>> &v_group) {
	if (sheet_id >= All_Sheets.size()) return false;
	
	vector<vector<uint32_t>>es_links(All_Sheets[sheet_id].middle_es.size()), es_links_temp;
	for (uint32_t i = 0; i < All_Sheets[sheet_id].middle_es.size(); i++)
		es_links[i] = frame.FEs[All_Sheets[sheet_id].middle_es[i]].es_link;

	vector<bool> e_flag(mesh.Es.size(), false);
	vector<uint32_t> v_pair(2);
	while (true) {
		for (uint32_t i = 0; i<es_links.size(); i++) {
			if (e_flag[es_links[i][0]]) continue;

			candiate_es_links.push_back(es_links[i]);

			int v1, v2;
			if (es_links[i].size() == 1) {
				v1 = mesh.Es[es_links[i][0]].vs[0];
				v2 = mesh.Es[es_links[i][0]].vs[1];
			}
			else {
				v1 = mesh.Es[es_links[i][0]].vs[0];
				for (uint32_t j = 0; j<mesh.Vs[v1].neighbor_es.size(); j++)
					if (mesh.Vs[v1].neighbor_es[j] == es_links[i][1]) {
						v1 = mesh.Es[es_links[i][0]].vs[1]; break;
					}

				uint32_t nume = es_links[i].size() - 1;
				v2 = mesh.Es[es_links[i][nume]].vs[0];
				for (uint32_t j = 0; j<mesh.Vs[v2].neighbor_es.size(); j++)
					if (mesh.Vs[v2].neighbor_es[j] == es_links[i][nume - 1]) {
						v2 = mesh.Es[es_links[i][nume]].vs[1]; break;
					}
			}
			v_pair[0] = v1;
			v_pair[1] = v2;
			v_group.push_back(v_pair);

			for (uint32_t j = 0; j<es_links[i].size(); j++) e_flag[es_links[i][j]] = true;


			for (uint32_t j = 0; j<mesh.Es[es_links[i][0]].neighbor_fs.size(); j++) {
				uint32_t nf = mesh.Es[es_links[i][0]].neighbor_fs[j];
				uint32_t op_eid = INVALID_E;
				for (uint32_t k = 0; k<4; k++) {
					if (mesh.Fs[nf].es[k] == es_links[i][0]) {
						op_eid = mesh.Fs[nf].es[(k + 2) % 4]; break;
					}
				}
				if (op_eid == INVALID_E) { 
					cout << "error " << endl; system("PAUSE"); 
				}
				if (e_flag[op_eid]) continue;

				uint32_t pre_nf = nf;
				vector<uint32_t> es_link;
				es_link.push_back(op_eid);

				for (uint32_t k = 1; k<es_links[i].size(); k++) {
					for (uint32_t m = 0; m<mesh.Es[es_links[i][k]].neighbor_fs.size(); m++) {
						uint32_t n_f = mesh.Es[es_links[i][k]].neighbor_fs[m];
						bool havesharede = false;
						for (uint32_t p = 0; p<4; p++) {
							for (uint32_t q = 0; q<4; q++) {
								if (mesh.Fs[n_f].es[p] == mesh.Fs[pre_nf].es[q]) {
									havesharede = true; break;
								}
							}
							if (havesharede) break;
						}
						if (havesharede) {
							pre_nf = n_f;
							op_eid = INVALID_E;
							for (uint32_t p = 0; p<4; p++) {
								if (mesh.Fs[n_f].es[p] == es_links[i][k])
								{
									op_eid = mesh.Fs[n_f].es[(p + 2) % 4]; break;
								}
							}
							if (op_eid == INVALID_E) {
								cout << "error " << endl; system("PAUSE"); 
							}
							es_link.push_back(op_eid);
							break;
						}
					}
				}
				es_links_temp.push_back(es_link);
			}
		}
		if (!es_links_temp.size()) break;
		else {
			es_links.clear();
			es_links.swap(es_links_temp);
		}
	}

	return true;
}
bool simplification::target_surface_sheet(uint32_t sheet_id, vector<vector<uint32_t>> &candiate_es_links, vector<vector<uint32_t>> &v_group) {
	vector<vector<uint32_t>> vs_links_temp(candiate_es_links.size());
	for (uint32_t i = 0; i < candiate_es_links.size();i++) {
		vector<uint32_t> vs_link;
		uint32_t first_v = v_group[i][0]; vs_link.push_back(first_v);
		for (auto eid: candiate_es_links[i]){
			if (mesh.Es[eid].vs[0] == first_v) {
				first_v = mesh.Es[eid].vs[1];
				vs_link.push_back(first_v);
			}
			else if (mesh.Es[eid].vs[1] == first_v) {
				first_v = mesh.Es[eid].vs[0];
				vs_link.push_back(first_v);
			}
		}
		vs_links_temp[i] = vs_link;
	}
	vector<bool> V_flag(mesh.Vs.size(), true);
	for (auto v_link : vs_links_temp)for (uint32_t j = 1; j < v_link.size() - 1; j++) V_flag[v_link[j]] = false;

	vector<vector<uint32_t>> &vs_pairs = All_Sheets[sheet_id].vs_pairs,
							 &vs_links = All_Sheets[sheet_id].vs_links;
	vs_pairs.clear(); vs_links.clear();
	for (uint32_t i = 0; i < v_group.size();i++) {
		vector<uint32_t> &a_pair = v_group[i];
		if (!V_flag[a_pair[0]] && !V_flag[a_pair[1]]) continue;
		vs_pairs.push_back(a_pair);
		vs_links.push_back(vs_links_temp[i]);
	}
	//grouping
	vector<vector<uint32_t>> &Vs_Group = All_Sheets[sheet_id].Vs_Group;
	Vs_Group.clear();
	vector<bool> v_flag(mesh.Vs.size(), false);
	vector<vector<uint32_t>> v_nvs(mesh.Vs.size());
	for (auto a_pair: vs_pairs) {
		v_nvs[a_pair[0]].push_back(a_pair[1]);
		v_nvs[a_pair[1]].push_back(a_pair[0]);
		v_flag[a_pair[0]] = true;
		v_flag[a_pair[1]] = true;
	}
	while (true) {
		vector<uint32_t> a_set, a_pool;
		for (auto a_pair : vs_pairs) {
			if (v_flag[a_pair[0]]) { a_pool.push_back(a_pair[0]); break; }
			else if (v_flag[a_pair[1]]) { a_pool.push_back(a_pair[1]); break; }
		}
		if (!a_pool.size()) break;
		a_set = a_pool;
		v_flag[a_pool[0]] = false;
		while (a_pool.size()) {
			vector<uint32_t> a_pool_sudo;
			for (auto vid0: a_pool)
				for (auto vid0_: v_nvs[vid0])
					if (v_flag[vid0_]) {
						a_pool_sudo.push_back(vid0_); v_flag[vid0_] = false;
					}
			a_pool.clear();
			if (a_pool_sudo.size()) {
				a_pool.swap(a_pool_sudo);
				a_set.insert(a_set.end(), a_pool.begin(), a_pool.end());
			}
		}
		Vs_Group.push_back(a_set);
	}

	CI.V_Groups = Vs_Group;
	vector<uint32_t>().swap(CI.hs);
	for (auto fhid : All_Sheets[sheet_id].cs) CI.hs.insert(CI.hs.end(), frame.FHs[fhid].hs_net.begin(), frame.FHs[fhid].hs_net.end());

	//target_vs, target_coords
	All_Sheets[sheet_id].target_vs.resize(Vs_Group.size());
	All_Sheets[sheet_id].target_coords.resize(Vs_Group.size(), 3);

	vector<bool> b_tags(mesh.Vs.size(),false);
	for (auto ffid:All_Sheets[sheet_id].side_fs){
		if (!frame.FFs[ffid].boundary) continue;
		for (auto hfid : frame.FFs[ffid].ffs_net)
			for (auto hvid : mesh.Fs[hfid].vs) b_tags[hvid] = true;
	}

	bool multiple_corners = false, multiple_curves = false, corner_curve_conflict = false, multiple_boundaries = false, corner_boundary_conflict = false;

	for (uint32_t i = 0; i<Vs_Group.size(); i++){
		vector<uint32_t> &a_group = Vs_Group[i];
		uint32_t Id = a_group[0];

		bool on_boundary = false, on_feature = false;
		//if on boundary
		uint32_t b_count = 0; vector<uint32_t> bvs;
		for (auto vid: a_group) if (b_tags[vid]) {
			Id = vid;
			bvs.push_back(Id);
			on_boundary = true;
			b_count++;
		}
		if(on_boundary) All_Sheets[sheet_id].target_coords.row(i) = mesh.V.col(Id);
		if (b_count > 1) multiple_boundaries = true;
		//vs_links_group
		vector<vector<uint32_t>> vs_links_group;
		vector<uint32_t> link_ids;
		for (auto vid : a_group) {
			for (uint32_t j = 0; j < vs_pairs.size();j++)
				for (auto vvid : vs_pairs[j]) if(vvid==vid) link_ids.push_back(j);
		}
		sort(link_ids.begin(), link_ids.end()); link_ids.erase(unique(link_ids.begin(), link_ids.end()), link_ids.end());

		for (auto id : link_ids) vs_links_group.push_back(vs_links[id]);
		//if on feature line/corner
		vector<uint32_t> coner_ids, corner;
		vector<uint32_t> vs_type1, curves;
		for (auto vid : a_group) {
			if (fc.V_types[vid] == Feature_V_Type::CORNER) {
				coner_ids.push_back(vid);
				corner.push_back(fc.V_ids[vid]);
			}
			else if (fc.V_types[vid] == Feature_V_Type::LINE) {
				vs_type1.push_back(vid);
				curves.push_back(fc.V_ids[vid]);
			}
		}
		sort(curves.begin(), curves.end()); curves.erase(unique(curves.begin(), curves.end()), curves.end());
		if (coner_ids.size() || vs_type1.size()) {
			sort(coner_ids.begin(), coner_ids.end());
			coner_ids.erase(unique(coner_ids.begin(), coner_ids.end()), coner_ids.end());
			if (coner_ids.size()) {
				if (coner_ids.size()> 1) multiple_corners = true;
				Id = coner_ids[0];
				All_Sheets[sheet_id].target_coords.row(i) = mesh.V.col(Id);

				if (curves.size() > 1) multiple_curves = true;
				else if (curves.size() == 1) {
					int which_corner = std::find(mf.corners.begin(), mf.corners.end(), corner[0]) - mf.corners.begin();
					if (std::find(mf.corner_curves[which_corner].begin(), mf.corner_curves[which_corner].end(), curves[0]) == mf.corner_curves[which_corner].end())
						corner_curve_conflict = true;
				}
			}
			else  if (vs_type1.size()) {
				if (curves.size() > 1) multiple_curves = true;
				Vector3d coords(0, 0, 0);
				for (auto vid : vs_type1) coords += mesh.V.col(vid);
				coords /= vs_type1.size();
				All_Sheets[sheet_id].target_coords.row(i) = coords;
				Id = vs_type1[0];
			}
			on_feature = true;
		}
		if (on_boundary && coner_ids.size()) {
			if(bvs[0] != coner_ids[0]) corner_boundary_conflict = true;

		}
		if (!on_boundary&&(coner_ids.size() || vs_type1.size())) {
			sort(coner_ids.begin(), coner_ids.end()); 
			coner_ids.erase(unique(coner_ids.begin(), coner_ids.end()), coner_ids.end());
			if (coner_ids.size()) {
				if (coner_ids.size()> 1) multiple_corners = true;
				Id = coner_ids[0];
				All_Sheets[sheet_id].target_coords.row(i) = mesh.V.col(Id);
			}
			else  if (vs_type1.size()) {
				if (curves.size() > 1) multiple_curves = true;
				Vector3d coords(0, 0, 0);
				for (auto vid : vs_type1) coords += mesh.V.col(vid);
				coords /= vs_type1.size();
				All_Sheets[sheet_id].target_coords.row(i) = coords;
				Id = vs_type1[0];
			}
			on_feature = true;
		}

		if (!on_boundary && !on_feature){
			Vector3d coords(0, 0, 0);
			for (auto a_link : vs_links_group) {
				Id = a_link[a_link.size() / 2];
				if(a_link.size() ==2) coords +=(mesh.V.col(a_link[0]) + mesh.V.col(a_link[1])) / 2;
				else if (a_link.size() % 2 == 0) {
					coords += (mesh.V.col(a_link[a_link.size() / 2 - 1]) + mesh.V.col(Id)) * 0.5;
				}
				else if (a_link.size() % 2 == 1) {
					coords += mesh.V.col(Id);
				}
			}
			coords /= vs_links_group.size();
			All_Sheets[sheet_id].target_coords.row(i) = coords;
		}
		All_Sheets[sheet_id].target_vs[i] = Id;
	}

	CI.target_vs = All_Sheets[sheet_id].target_vs;
	CI.target_coords = All_Sheets[sheet_id].target_coords;
	if (multiple_boundaries || multiple_curves || multiple_corners || corner_curve_conflict || corner_boundary_conflict) return false;
	return true;
}
bool simplification::vs_group_chord(uint32_t chord_id) {
	vector<vector<uint32_t>> vs_chains[4];
	vector<uint32_t> links[4];
	vector<bool> v_flag(mesh.Vs.size(), true);
	//four vertical links
	for (uint32_t i = 0; i < 4; i++) for (uint32_t j = 0; j < All_Chords[chord_id].vertical_es[i].size();j++) {
		uint32_t feid = All_Chords[chord_id].vertical_es[i][j];
		vector<uint32_t> a_link = frame.FEs[feid].vs_link;
		if (a_link[0] != frame.FVs[All_Chords[chord_id].parallel_ns[i][j]].hid) reverse(a_link.begin(), a_link.end());
		for (auto vid : a_link) if (v_flag[vid]) { links[i].push_back(vid); v_flag[vid] = false; }
		}
	//four parallel chains
	uint32_t link_size = links[0].size();
	for (uint32_t i = 0; i < 4; i++) for (uint32_t j = 0; j < link_size; j++) {
		vector<uint32_t> a_chain;
		if (j == 0) {
			uint32_t feid = All_Chords[chord_id].parallel_es[i][0];
			a_chain = frame.FEs[feid].vs_link;
			if(a_chain[0] != frame.FVs[All_Chords[chord_id].parallel_ns[i][0]].hid) reverse(a_chain.begin(), a_chain.end());
		}
		else {
			a_chain.push_back(links[i][j]);
			for (uint32_t k = 1; k<vs_chains[i][j-1].size(); k++){
				uint32_t v0 = vs_chains[i][j-1][k-1], v1 = vs_chains[i][j - 1][k], v2 = a_chain[k-1], v3;
				sort(mesh.Vs[v0].neighbor_fs.begin(), mesh.Vs[v0].neighbor_fs.end());
				sort(mesh.Vs[v1].neighbor_fs.begin(), mesh.Vs[v1].neighbor_fs.end());
				sort(mesh.Vs[v2].neighbor_fs.begin(), mesh.Vs[v2].neighbor_fs.end());
				vector<uint32_t> sharedf0, sharedf1;
				set_intersection(mesh.Vs[v0].neighbor_fs.begin(), mesh.Vs[v0].neighbor_fs.end(),
					mesh.Vs[v1].neighbor_fs.begin(), mesh.Vs[v1].neighbor_fs.end(), back_inserter(sharedf0));
				set_intersection(sharedf0.begin(), sharedf0.end(), mesh.Vs[v2].neighbor_fs.begin(), mesh.Vs[v2].neighbor_fs.end(), back_inserter(sharedf1));
				for (auto vid : mesh.Fs[sharedf1[0]].vs) if (vid != v0 &&vid != v1 &&vid != v2) {v3 = vid; break;}
				a_chain.push_back(v3);
			}
		}
		vs_chains[i].push_back(a_chain);
	}
	//find longer side
	bool longer = false; int side = -1; uint32_t v_end = -1;
	uint32_t longer_num, shorter_num;
	uint32_t feid = All_Chords[chord_id].parallel_es[0][0], feid_temp = All_Chords[chord_id].parallel_es[1][0];
	longer_num = frame.FEs[feid].vs_link.size(); shorter_num = frame.FEs[feid_temp].vs_link.size();
	if (longer_num < shorter_num) {
		swap(longer_num, shorter_num); feid = feid_temp; longer = true; side = 1;
	}
	else if (longer_num > shorter_num) { longer = true; side = 0;}

	if (longer) {
		v_end = vs_chains[side][0][longer_num - 1];

		vector<vector<uint32_t>> vs_chains_temp[4];
		for (uint32_t i = 0; i < link_size; i++) {
			vector<uint32_t> a_chain;
			a_chain.insert(a_chain.end(), vs_chains[side][i].begin(), vs_chains[side][i].begin() + shorter_num);
			vs_chains_temp[side].push_back(a_chain); 

			int sub_len = longer_num - shorter_num;
			a_chain.clear();
			a_chain.insert(a_chain.end(), vs_chains[side+2][i].begin()+ sub_len, vs_chains[side + 2][i].end());
			vs_chains_temp[side+2].push_back(a_chain);
		}
		vs_chains_temp[side].swap(vs_chains[side]);
		vs_chains_temp[side + 2].swap(vs_chains[side + 2]);
	}

	vector<vector<uint32_t>> collapse_pair00, collapse_pair01, collapse_pair10, collapse_pair11;
	if (All_Chords[chord_id].side == 0) {
		for (uint32_t i = 0; i < link_size; i++) {
			reverse(vs_chains[2][i].begin(), vs_chains[2][i].end());
			reverse(vs_chains[3][i].begin(), vs_chains[3][i].end());
		}
		collapse_pair00 = vs_chains[0];
		collapse_pair01 = vs_chains[3];
		collapse_pair10 = vs_chains[1];
		collapse_pair11 = vs_chains[2];
	}
	else {
		for (uint32_t i = 0; i < link_size; i++) {
			reverse(vs_chains[1][i].begin(), vs_chains[1][i].end());
			reverse(vs_chains[3][i].begin(), vs_chains[3][i].end());
		}
		collapse_pair00 = vs_chains[0];
		collapse_pair01 = vs_chains[1];
		collapse_pair10 = vs_chains[2];
		collapse_pair11 = vs_chains[3];
	}
	vector<vector<uint32_t>> collapse_pairs;
	collapse_pairs.reserve(collapse_pair00.size() * collapse_pair00[0].size() * 2);
	for (uint32_t i = 0; i < collapse_pair00.size(); i++) {
		for (uint32_t j = 0; j < collapse_pair00[i].size(); j++) {
			vector<uint32_t> a_pair(2);
			a_pair[0] = collapse_pair00[i][j];
			a_pair[1] = collapse_pair01[i][j];
			collapse_pairs.push_back(a_pair);
			a_pair[0] = collapse_pair10[i][j];
			a_pair[1] = collapse_pair11[i][j];
			collapse_pairs.push_back(a_pair);
		}
	}
	vector<vector<uint32_t>> vs_pairs, vs_links_cut;
	if (longer) {
		uint32_t sheet_id = -1;
		for (auto sheet : All_Sheets)
			if (find(sheet.middle_es.begin(), sheet.middle_es.end(), feid) != sheet.middle_es.end()) {
				sheet_id = sheet.id; break;
			}
		if (sheet_id == -1) { cout << "doesnot find the longer sheet" << endl; system("PAUSE"); }
		if (All_Sheets[sheet_id].type != Sheet_type::close && All_Sheets[sheet_id].type != Sheet_type::open) {
			return false;
		}
		//sheet v_group 
		vector<vector<uint32_t>> candiate_es_links, v_group;
		if (!vs_pair_sheet(sheet_id, candiate_es_links, v_group)) { cout << "ERROR: no vs pairs" << endl;  system("PAUSE"); }

		std::fill(v_flag.begin(), v_flag.end(), false);
		for (auto ffid : All_Sheets[sheet_id].left_fs)
			for (auto fid : frame.FFs[ffid].ffs_net) {
				for(uint32_t i=0;i<4;i++)  v_flag[mesh.Fs[fid].vs[i]] = true;
			}

		for (auto &a_pair : v_group) { if (v_flag[a_pair[1]]) swap(a_pair[0], a_pair[1]); }

		vector<vector<uint32_t>> vs_links_temp(candiate_es_links.size());
		vs_links_cut.resize(candiate_es_links.size());
		for (uint32_t i = 0; i < candiate_es_links.size(); i++) {
			vector<uint32_t> vs_link;
			uint32_t first_v = v_group[i][0]; vs_link.push_back(first_v);
			uint32_t e0 = candiate_es_links[i][0]; 
			if (first_v != mesh.Es[e0].vs[0] && first_v != mesh.Es[e0].vs[1])
				reverse(candiate_es_links[i].begin(), candiate_es_links[i].end());
			for (auto eid : candiate_es_links[i]) {
				if (mesh.Es[eid].vs[0] == first_v) {
					first_v = mesh.Es[eid].vs[1];
					vs_link.push_back(first_v);
				}
				else if (mesh.Es[eid].vs[1] == first_v) {
					first_v = mesh.Es[eid].vs[0];
					vs_link.push_back(first_v);
				}
			}
			vs_links_temp[i] = vs_link;
		}
		//cut sheet
		uint32_t which = -1, which_side = 0;
		for (uint32_t i = 0; i < v_group.size(); i++) {
			if (v_end == v_group[i][0]) { which = longer_num - shorter_num + 1; break; }
			else if (v_end == v_group[i][1]) { which_side = 1;  which = shorter_num - 1; break; }
		}
		vs_pairs.resize(v_group.size());
		for (uint32_t i = 0; i < v_group.size(); i++) {
			if (which_side == 0) {
				vs_links_cut[i].insert(vs_links_cut[i].end(), vs_links_temp[i].begin(), vs_links_temp[i].begin() + which);
				vs_pairs[i].push_back(v_group[i][which_side]);
				vs_pairs[i].push_back(vs_links_temp[i][which - 1]);
			}
			else {
				vs_links_cut[i].insert(vs_links_cut[i].end(), vs_links_temp[i].begin() + which, vs_links_temp[i].end());
				vs_pairs[i].push_back(vs_links_temp[i][which]);
				vs_pairs[i].push_back(v_group[i][which_side]);
			}
		}
	}
	CI.V_Groups.clear();
	std::fill(v_flag.begin(), v_flag.end(), false);
	vector<vector<uint32_t>> V_neighbors(mesh.Vs.size());
	for (auto a_pair:collapse_pairs) {
		V_neighbors[a_pair[0]].push_back(a_pair[1]);
		V_neighbors[a_pair[1]].push_back(a_pair[0]);
		v_flag[a_pair[0]] = v_flag[a_pair[1]] = true;
	}
	while (true) {
		vector<uint32_t> a_set, a_pool;
		for (auto a_pair : collapse_pairs) {
			if (v_flag[a_pair[0]]) { a_pool.push_back(a_pair[0]); break; }
			else if (v_flag[a_pair[1]]) { a_pool.push_back(a_pair[1]); break; }
		}
		if (!a_pool.size()) break;
		a_set = a_pool;
		v_flag[a_pool[0]] = false;
		while (a_pool.size()) {
			vector<uint32_t> a_pool_sudo;
			for (uint32_t i = 0; i < a_pool.size(); i++) {
				for (uint32_t j = 0; j < V_neighbors[a_pool[i]].size(); j++) {
					if (v_flag[V_neighbors[a_pool[i]][j]]) a_pool_sudo.push_back(V_neighbors[a_pool[i]][j]);
				}
			}
			a_pool.clear();
			if (a_pool_sudo.size()) {
				for (auto vid : a_pool_sudo) v_flag[vid] = false;
				a_pool.swap(a_pool_sudo);
				a_set.insert(a_set.end(), a_pool.begin(), a_pool.end());
			}
		}
		sort(a_set.begin(), a_set.end());
		a_set.erase(unique(a_set.begin(), a_set.end()), a_set.end());
		CI.V_Groups.push_back(a_set);
	}
	//compose final v_groups: cutted_sheet
	vector<int32_t> v_group_tag(mesh.Vs.size(), -1);
	for (uint32_t i = 0; i < CI.V_Groups.size();i++) for (auto vid : CI.V_Groups[i]) v_group_tag[vid] = i;
	for (auto a_pair : vs_pairs) {
		int32_t which_group = -1;
		for (auto vid : a_pair) if (v_group_tag[vid] != -1) which_group = v_group_tag[vid];
		if (which_group == -1) CI.V_Groups.push_back(a_pair);
		else {
			CI.V_Groups[which_group].insert(CI.V_Groups[which_group].end(), a_pair.begin(), a_pair.end());
			sort(CI.V_Groups[which_group].begin(), CI.V_Groups[which_group].end());
			CI.V_Groups[which_group].erase(unique(CI.V_Groups[which_group].begin(), CI.V_Groups[which_group].end()), CI.V_Groups[which_group].end());
		}
	}

	CI.hs.clear();
	vector<bool> H_flag(mesh.Hs.size(), false);
	std::fill(v_flag.begin(), v_flag.end(), false);
	for (auto a_cut : vs_links_cut)for (auto vid : a_cut) v_flag[vid] = true;
	for (auto h : mesh.Hs) {
		bool removable = true;
		for (auto vid : h.vs)if (!v_flag[vid])removable = false;
		if (removable) {
			CI.hs.push_back(h.id);
			H_flag[h.id] = true;
		}
	}
	for (auto fhid : All_Chords[chord_id].cs)
		for (auto hid : frame.FHs[fhid].hs_net) if (!H_flag[hid]) CI.hs.push_back(hid);

	return true;
}
bool simplification::target_surface_chord(uint32_t chord_id) {

	//target_vs, target_coords
	CI.target_vs.resize(CI.V_Groups.size());
	CI.target_coords.resize(CI.V_Groups.size(), 3);

	vector<bool> b_tags(mesh.Vs.size(), false);
	for (uint32_t i=0;i<4;i++)
	for (auto ffid : All_Chords[chord_id].vertical_fs[i]) {
		if (!frame.FFs[ffid].boundary) continue;
		for (auto hfid : frame.FFs[ffid].ffs_net) for (auto hvid : mesh.Fs[hfid].vs) b_tags[hvid] = true;
	}
	bool multiple_corners = false, multiple_curves = false, corner_curve_conflict=false, multiple_boundaries = false, corner_boundary_conflict = false;
	for (uint32_t i = 0; i<CI.V_Groups.size(); i++) {
		vector<uint32_t> &a_group = CI.V_Groups[i];
		uint32_t Id = a_group[0];

		bool on_boundary = false, on_feature = false;
		//if on boundary
		uint32_t b_count = 0; vector<uint32_t> bvs;
		for (auto vid : a_group) if (b_tags[vid]) {
			Id = vid;
			bvs.push_back(Id);
			on_boundary = true;
			b_count++;
		}
		if (on_boundary) CI.target_coords.row(i) = mesh.V.col(Id);
		if (b_count > 1) multiple_boundaries = true;
		//if on feature line/corner
		vector<uint32_t> coner_ids, corner;
		vector<uint32_t> vs_type1, curves;
		for (auto vid : a_group) {
			if (fc.V_types[vid] == Feature_V_Type::CORNER) {
				coner_ids.push_back(vid);
				corner.push_back(fc.V_ids[vid]);
			}
			else if (fc.V_types[vid] == Feature_V_Type::LINE) {
				vs_type1.push_back(vid);
				curves.push_back(fc.V_ids[vid]);
			}
		}
		sort(curves.begin(), curves.end()); curves.erase(unique(curves.begin(), curves.end()), curves.end());
		if (coner_ids.size() || vs_type1.size()) {
			sort(coner_ids.begin(), coner_ids.end());
			coner_ids.erase(unique(coner_ids.begin(), coner_ids.end()), coner_ids.end());
			if (coner_ids.size()) {
				if (coner_ids.size()> 1) multiple_corners = true;
				Id = coner_ids[0];
				CI.target_coords.row(i) = mesh.V.col(Id);

				if (curves.size() > 1) multiple_curves = true;
				else if (curves.size() == 1) {
					int which_corner = std::find(mf.corners.begin(), mf.corners.end(), corner[0]) - mf.corners.begin();
					if (std::find(mf.corner_curves[which_corner].begin(), mf.corner_curves[which_corner].end(), curves[0]) == mf.corner_curves[which_corner].end())
						corner_curve_conflict = true;
				}
			}
			else  if (vs_type1.size()) {
				if (curves.size() > 1) multiple_curves = true;
				Vector3d coords(0, 0, 0);
				for (auto vid : vs_type1) coords += mesh.V.col(vid);
				coords /= vs_type1.size();
				CI.target_coords.row(i) = coords;
				Id = vs_type1[0];
			}
			on_feature = true;
		}
		if (on_boundary && coner_ids.size()) {
			if (bvs[0] != coner_ids[0]) corner_boundary_conflict = true;
		}
		if (!on_boundary && (coner_ids.size() || vs_type1.size())) {
			sort(coner_ids.begin(), coner_ids.end());
			coner_ids.erase(unique(coner_ids.begin(), coner_ids.end()), coner_ids.end());
			if (coner_ids.size()) {
				if (coner_ids.size()> 1) multiple_corners = true;
				Id = coner_ids[0];
				CI.target_coords.row(i) = mesh.V.col(Id);
			}
			else  if (vs_type1.size()) {
				if (curves.size() > 1) multiple_curves = true;
				Vector3d coords(0, 0, 0);
				for (auto vid : vs_type1) coords += mesh.V.col(vid);
				coords /= vs_type1.size();
				CI.target_coords.row(i) = coords;
				Id = vs_type1[0];
			}
			on_feature = true;
		}

		if (!on_boundary && !on_feature) {
			Vector3d coords(0, 0, 0);
			Id = a_group[0];
			for (auto vid : a_group) coords += mesh.V.col(vid);
			coords /= a_group.size();
			CI.target_coords.row(i) = coords;
		}
		CI.target_vs[i] = Id;
	}
	if (multiple_boundaries || multiple_corners || multiple_curves || corner_curve_conflict|| corner_boundary_conflict) return false;
	return true;
}
bool simplification::topology_check() {
	
	mesh_.type = Mesh_type::Hex;
	auto &Vs_Group = CI.V_Groups;
	vector<bool> V_flag(mesh.Vs.size(), false), H_flag(mesh.Hs.size(), false);

	for (auto hid : CI.hs) {
		H_flag[hid] = true;
		for (uint32_t i = 0; i < 8; i++) V_flag[mesh.Hs[hid].vs[i]] = true;
	}
	//V_map, RV_map
	V_map.resize(mesh.Vs.size()); fill(V_map.begin(), V_map.end(), INVALID_V);
	vector<uint32_t>().swap(RV_map);

	for (auto vs : Vs_Group) for (auto vid : vs) V_flag[vid] = true;
	for (uint32_t i = 0; i < Vs_Group.size(); i++) V_flag[CI.target_vs[i]] = false;
	uint32_t V_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++) if (!V_flag[i]) { V_map[i] = V_num++; RV_map.push_back(i); }
	for (uint32_t i = 0; i < Vs_Group.size(); i++) {
		uint32_t mid = V_map[CI.target_vs[i]];
		for (auto vid : Vs_Group[i])V_map[vid] = mid;
	}
	//Vs
	vector<Hybrid_V>().swap(mesh_.Vs);
	mesh_.V.resize(3, V_num); mesh_.V.setZero();
	mesh_.Vs.resize(V_num);
	for (uint32_t i = 0; i<V_num; i++){
		Hybrid_V v;
		v.id = i; v.boundary = false;
		mesh_.Vs[i] = v;
	}
	for (uint32_t i = 0; i < V_flag.size(); i++) if (!V_flag[i]) mesh_.V.col(V_map[i]) = mesh.V.col(i);
	for (uint32_t i = 0; i < CI.target_vs.rows(); i++) mesh_.V.col(V_map[CI.target_vs[i]]) = CI.target_coords.row(i);

	//Hs
	uint32_t H_num = 0;
	for (auto h_flag : H_flag) if (!h_flag) H_num++;

	vector<Hybrid>().swap(mesh_.Hs);

	mesh_.Hs.resize(H_num); H_num = 0;
	Hybrid h_; h_.vs.resize(8);
	for (auto h : mesh.Hs) if (!H_flag[h.id]) {
		h_.id= H_num++;
		for (uint32_t i = 0; i < 8; i++) h_.vs[i] = V_map[h.vs[i]];
		mesh_.Hs[h_.id] = h_;
		for (uint32_t i = 0; i < h_.vs.size(); i++) mesh_.Vs[h_.vs[i]].neighbor_hs.push_back(h_.id);
	}
	
	build_connectivity(mesh_);
	base_com.singularity_structure(si_, mesh_);
	if(!base_com.base_complex_extraction(si_, frame_, mesh_)) return false;

	Mesh_Topology mt_;
	topology_info(mesh_, frame_, mt_);
	return comp_topology(mt, mt_);
}

void simplification::optimization() {
	
	Mesh_Quality mq;
	//quality
	scaled_jacobian(mesh, mq);
	Slim_global_region = mesh.Hs.size();
	OPTIMIZATION_ONLY = true;
	CI.target_vs.resize(0);
	tetralize_mesh_omesh(ts, mesh);
	fc.lamda_C = 1e+3;
	fc.lamda_L = 1e+3;
	fc.lamda_T = 1e+3;
	ts.fc = fc;
	ts.global = true;
	
	Mesh_Quality mq_pre = mq;
	mesh_ = mesh;
	for (uint32_t i = 0; i < 5 * Slim_Iteration; i++) {
		ts.projection = false;

		compute_referenceMesh(ts.V, mesh.Hs, CI.Hsregion, ts.RT);

		slim_opt(ts, 1);

		mesh_.V = ts.V.transpose();

		scaled_jacobian(mesh_, mq);
		if (mq_pre.min_Jacobian > mq.min_Jacobian) break;
		if (!hausdorff_ratio_check(mf.tri, mesh_)) break;
		
		mesh.V = mesh_.V;
		mq_pre = mq;

		project_surface_update_feature(mf, ts.fc, ts.V, ts.s, ts.sc, 1);
		ts.projection = true;
	}
	scaled_jacobian(mesh, mq);
	std::cout << "after: minimum scaled J: " << mq.min_Jacobian << " average scaled J: " << mq.ave_Jacobian << endl;

	hausdorff_ratio_check(mf.tri, mesh);
}
bool simplification::direct_collapse() {
	Mesh_Quality mq;

	OPTIMIZATION_ONLY = false;

	tetralize_mesh(ts);
	ts.fc = fc;
	ts.global = true;
	ts.projection = false;
	ts.s.resize(0); ts.sc.resize(0, 3);
	for (uint32_t i = 0; i < Slim_Iteration; i++) {
		ts.projection = false;

		compute_referenceMesh(ts.V, mesh.Hs, CI.Hsregion, ts.RT);

		slim_opt(ts, 1);
		project_surface_update_feature(mf, ts.fc, ts.V, ts.s, ts.sc, Projection_range);

		ts.projection = true;
	}
	ts.Vgroups.clear();
	project_surface_update_feature(mf, ts.fc, ts.V, ts.s, ts.sc, Projection_range);

	//====================post-update====================//
	vector<bool> touchedV_flag(mesh.Vs.size(), false);
	for (uint32_t i = 0; i < ts.T.rows(); i++)for (uint32_t j = 0; j < 4; j++)touchedV_flag[ts.T(i, j)] = true;
	for (uint32_t i = 0; i < ts.b.rows(); i++) touchedV_flag[ts.b[i]] = true;
	for (uint32_t i = 0; i < touchedV_flag.size(); i++) if (touchedV_flag[i] && V_map[i] != INVALID_V) {
		mesh_.V.col(V_map[i]) = ts.V.row(i);
	}

	for (uint32_t i = 0; i < CI.V_Groups.size(); i++) {
		vector<uint32_t> &vs = CI.V_Groups[i];
		uint32_t mv_id = V_map[CI.target_vs[i]];
		mesh_.V.col(mv_id).setZero();
		for (uint32_t j = 0; j < vs.size(); j++) {
			mesh_.V.col(mv_id) += ts.V.row(vs[j]);
			if (vs[j] != CI.target_vs[i]) V_map[vs[j]] = INVALID_V;
		}
		mesh_.V.col(mv_id) /= vs.size();
	}

	scaled_jacobian(mesh_, mq);

	if (mq.min_Jacobian < Jacobian_Bound) {
		return false;
	}

	//assign fc 
	Feature_Constraints fc_temp;
	update_feature_variable_index_newmesh(ts, fc_temp, V_map, mesh_.Vs.size());
	tetralize_mesh_omesh(ts, mesh_);
	ts.fc = fc_temp;
	ts.global = true;
	ts.s.resize(0); ts.sc.resize(0, 3);
	for (uint32_t i = 0; i < Slim_Iteration; i++){
		ts.projection = false;
		//RT
		compute_referenceMesh(ts.V, mesh_.Hs, CI.Hsregion, ts.RT);

		slim_opt(ts, 1);

		project_surface_update_feature(mf, ts.fc, ts.V, ts.s, ts.sc, Projection_range);
		
		ts.projection = true;
	}
	project_surface_update_feature(mf, ts.fc, ts.V, ts.s, ts.sc, Projection_range);

	std::fill(touchedV_flag.begin(), touchedV_flag.end(),false);
	for (uint32_t i = 0; i < ts.T.rows(); i++)for (uint32_t j = 0; j < 4; j++)touchedV_flag[ts.T(i, j)] = true;
	for (uint32_t i = 0; i < ts.b.rows(); i++) touchedV_flag[ts.b[i]] = true;
	for (uint32_t i = 0; i < touchedV_flag.size(); i++) if (touchedV_flag[i]) {
		mesh_.V.col(i) = ts.V.row(i);
	}

	scaled_jacobian(mesh_, mq);
	if (mq.min_Jacobian < Jacobian_Bound) { std::cout << "double check smoothing" << endl;
	return false;  system("PAUSE"); 
	}

	if (!hausdorff_ratio_check(mf.tri, mesh_)) return false;

	fc = ts.fc;

	mesh = mesh_;
	si = si_;
	frame = frame_;
	
	return true;
}
bool simplification::tetralize_mesh(Tetralize_Set &ts) {
	//re-index
	vector<bool> V_flag(mesh.Vs.size(), false), F_flag(mesh.Fs.size(), false), H_flag(mesh.Hs.size(), false);
	CI.fs_before.clear();
	for (auto hid : CI.hs) H_flag[hid] = true;

	for (auto hid : CI.hs) {
		for (auto fid : mesh.Hs[hid].fs) if (!F_flag[fid] && !mesh.Fs[fid].boundary) {
			uint32_t h0 = mesh.Fs[fid].neighbor_hs[0], h1 = mesh.Fs[fid].neighbor_hs[1];
			if (H_flag[h0] != H_flag[h1]) {
				CI.fs_before.push_back(fid); F_flag[fid] = true;
			}
		}
	}

	vector<uint32_t> Hsregion;
	grow_region2(CI.hs.size(), CI.fs_before, CI.before_region, Hsregion, ts, H_flag, mesh, false);
	CI.Hsregion = Hsregion;
	//V
	ts.V = mesh.V.transpose();
	//T
	ts.T.resize(8 * Hsregion.size(), 4);
	Vector4i t;
	for (uint32_t k = 0; k < Hsregion.size();k++) {
		auto h = mesh.Hs[Hsregion[k]];
		for (uint32_t i = 0; i < 8; i++) {
			for (uint32_t j = 0; j < 4; j++) t[j] = h.vs[hex_tetra_table[i][j]];
			ts.T.row(8 * k + i) = t;
		}
	}

	uint32_t constraint_Num = 0;
	for (auto vs : CI.V_Groups) constraint_Num += vs.size();
	ts.b.resize(constraint_Num);
	ts.bc.resize(constraint_Num, 3); ts.bc.setZero();
	constraint_Num = 0;
	for (uint32_t i = 0; i < CI.V_Groups.size(); i++) {
		for (auto vid : CI.V_Groups[i]) {
			ts.b[constraint_Num] = vid;
			ts.bc.row(constraint_Num) = CI.target_coords.row(i);
			constraint_Num++;
		}
	}
	return true;
}
bool simplification::tetralize_mesh_omesh(Tetralize_Set &ts, Mesh &mesh_r) {

	ts.V = mesh_r.V.transpose();

	vector<bool> V_flag(mesh_r.Vs.size(), false);
	for (uint32_t i = 0; i < CI.target_vs.size();i++) V_flag[V_map[CI.target_vs[i]]] = true;
	
	CI.fs_after.clear();
	for (auto f : mesh_r.Fs) {
		bool thisone = true;
		for (auto vid:f.vs) if(!V_flag[vid]) thisone= false;
		if (thisone) CI.fs_after.push_back(f.id);
	}

	vector<uint32_t> Hsregion;
	vector<bool> H_flag(mesh_r.Hs.size(), false);
	int base_num = CI.hs.size();
	if (OPTIMIZATION_ONLY) base_num = mesh_r.Hs.size();
	grow_region2(base_num, CI.fs_after, CI.after_region, Hsregion, ts, H_flag, mesh_r, true);
	CI.Hsregion = Hsregion;
	//T_
	vector<vector<uint32_t>> T_;
	vector<uint32_t> t(4);
	for (auto hid : Hsregion) {
		auto h = mesh_r.Hs[hid];
		for (uint32_t i = 0; i < 8; i++) {
			for (uint32_t j = 0; j < 4; j++) t[j] = h.vs[hex_tetra_table[i][j]];
			T_.push_back(t);
		}
	}
	//T
	ts.T.resize(T_.size(), 4);
	for (uint32_t i = 0; i < T_.size(); i++)  for (uint32_t j = 0; j < 4; j++) ts.T(i, j) = T_[i][j];

	ts.b.resize(0);
	ts.bc.resize(0, 3); ts.bc.setZero();
	return true;
}
bool simplification::tetralize_mesh_submesh(Tetralize_Set &ts, Mesh &mesh_r){
	ts.V = mesh_r.V.transpose();

	vector<bool> V_flag(mesh_r.Vs.size(), false), F_flag(mesh_r.Fs.size(), false), H_flag(mesh_r.Hs.size(), false);
	CI.fs_subdivided.clear();
	for (auto hid : CI.hs) {
		H_flag[hid] = true;
		for (uint32_t i = 0; i < 8; i++) V_flag[mesh_r.Hs[hid].vs[i]] = true;
	}
	for (auto hid : CI.hs) {
		for (auto fid : mesh_r.Hs[hid].fs) if (!F_flag[fid] && !mesh_r.Fs[fid].boundary) {
			uint32_t h0 = mesh_r.Fs[fid].neighbor_hs[0], h1 = mesh_r.Fs[fid].neighbor_hs[1];
			if (H_flag[h0] != H_flag[h1]) {
				CI.fs_subdivided.push_back(fid); F_flag[fid] = true;
			}
		}
	}
	for (auto fid : CI.fs_subdivided) F_flag[fid] = false;//but useless

	vector<uint32_t> Hsregion;
	grow_region2(CI.hs.size(), CI.fs_subdivided, CI.subd_region, Hsregion, ts, H_flag, mesh_r, true);
	Hsregion.insert(Hsregion.end(), CI.hs.begin(), CI.hs.end());
	CI.Hsregion = Hsregion;

	//T_
	vector<vector<uint32_t>> T_;
	vector<uint32_t> t(4);
	for (auto hid : Hsregion) {
		auto h = mesh_r.Hs[hid];
		for (uint32_t i = 0; i < 8; i++) {
			for (uint32_t j = 0; j < 4; j++) t[j] = h.vs[hex_tetra_table[i][j]];
			T_.push_back(t);
		}
	}
	//T
	ts.T.resize(T_.size(), 4);
	for (uint32_t i = 0; i < T_.size(); i++)  for (uint32_t j = 0; j < 4; j++) ts.T(i, j) = T_[i][j];

	ts.b.resize(0);
	ts.bc.resize(0, 3); ts.bc.setZero();
	return true;
}
bool simplification::grow_region(uint32_t base_num, vector<uint32_t> &frontFs, vector<uint32_t> &regionFs, vector<uint32_t> &newHs, Tetralize_Set &ts, vector<bool> &H_flag, Mesh &mesh_r, bool global) {
	
	if (base_num <= 0) return false;
	int upper_limit = 0; 
	if (global) {
		upper_limit = base_num * Slim_global_region;
	}
	else {
		upper_limit = base_num * Slim_region;
	}
	if (upper_limit >= mesh_r.Hs.size()) {
		upper_limit = mesh_r.Hs.size();
		for (auto h : mesh_r.Hs) if (!H_flag[h.id]) newHs.push_back(h.id);
		ts.regionb.resize(0);
		ts.regionbc.resize(0, 3);
		regionFs.clear();
		return true;
	}
	
	vector<bool> F_flag(mesh_r.Fs.size(), false);
	newHs.clear();	 
	vector<uint32_t> FFs = frontFs;
	bool newadded = false;
	while (newHs.size() <= upper_limit) {
		newadded = false;
		for (auto fid : FFs) {
			for (auto hid : mesh_r.Fs[fid].neighbor_hs) if (!H_flag[hid]) {
				H_flag[hid] = true;
				newHs.push_back(hid);
				newadded = true;
			}
		}
		for (auto fid : FFs)F_flag[fid] = true;
		FFs.clear();
		if (!newadded) break;
		for (auto hid : newHs) {
			for (auto fid : mesh_r.Hs[hid].fs)if (!F_flag[fid] && !mesh_r.Fs[fid].boundary) {
				uint32_t h0 = mesh_r.Fs[fid].neighbor_hs[0], h1 = mesh_r.Fs[fid].neighbor_hs[1];
				if (H_flag[h0] != H_flag[h1]) {
					FFs.push_back(fid); F_flag[fid] = true;
				}
			}
		}
	}
	//region boundary constraint
	regionFs = FFs;

	vector<uint32_t> vs;
	vector<bool> V_flag(mesh_r.Vs.size(), false);
	for (auto fid : FFs) {
		for (auto vid : mesh_r.Fs[fid].vs)
			if (!V_flag[vid]) { vs.push_back(vid); V_flag[vid] = true; }
	}

	ts.regionb.resize(vs.size());
	ts.regionbc.resize(vs.size(), 3); ts.regionbc.setZero();
	for (uint32_t i = 0; i < vs.size(); i++) {
		ts.regionb[i] = vs[i];
		ts.regionbc.row(i) = mesh_r.V.col(vs[i]);
	}

	return true;
}
bool simplification::grow_region2(uint32_t base_num, vector<uint32_t> &frontFs, vector<uint32_t> &regionFs, vector<uint32_t> &newHs, Tetralize_Set &ts, vector<bool> &H_flag, Mesh &mesh_r, bool global) {

	if (base_num <= 0) return false;
	else if (base_num >= mesh_r.Hs.size()) {
		for (auto h : mesh_r.Hs) if (!H_flag[h.id]) newHs.push_back(h.id);
		ts.regionb.resize(0);
		ts.regionbc.resize(0, 3);
		regionFs.clear();
		return true;
	}

	int ringN = 0;
	if (global) ringN = width_sheet * Slim_global_region;
	else  ringN = width_sheet * Slim_region;

	vector<bool> F_flag(mesh_r.Fs.size(), false);
	newHs.clear();
	vector<uint32_t> FFs = frontFs;
	bool newadded = false;
	for (int i = 0; i < ringN;i++) {
		newadded = false;
		for (auto fid : FFs) {
			for (auto hid : mesh_r.Fs[fid].neighbor_hs) if (!H_flag[hid]) {
				H_flag[hid] = true;
				newHs.push_back(hid);
				newadded = true;
			}
		}
		for (auto fid : FFs)F_flag[fid] = true;
		FFs.clear();
		if (!newadded) break;
		for (auto hid : newHs) {
			for (auto fid : mesh_r.Hs[hid].fs)if (!F_flag[fid] && !mesh_r.Fs[fid].boundary) {
				uint32_t h0 = mesh_r.Fs[fid].neighbor_hs[0], h1 = mesh_r.Fs[fid].neighbor_hs[1];
				if (H_flag[h0] != H_flag[h1]) {
					FFs.push_back(fid); F_flag[fid] = true;
				}
			}
		}
	}
	//region boundary constraint
	regionFs = FFs;

	vector<uint32_t> vs;
	vector<bool> V_flag(mesh_r.Vs.size(), false);
	for (auto fid : FFs) {
		for (auto vid : mesh_r.Fs[fid].vs)
			if (!V_flag[vid]) { vs.push_back(vid); V_flag[vid] = true; }
	}

	ts.regionb.resize(vs.size());
	ts.regionbc.resize(vs.size(), 3); ts.regionbc.setZero();
	for (uint32_t i = 0; i < vs.size(); i++) {
		ts.regionb[i] = vs[i];
		ts.regionbc.row(i) = mesh_r.V.col(vs[i]);
	}

	return true;
}

void simplification::update_feature_variable_index_newmesh(Tetralize_Set &ts, Feature_Constraints &fcc,
	vector<uint32_t> &new_V_map, uint32_t new_Vsize) {

	vector<bool> V_flag(new_V_map.size(), false);
	for (uint32_t i = 0; i < V_flag.size(); i++) if (new_V_map[i] != INVALID_V) V_flag[i] = true;

	fcc.lamda_C = ts.fc.lamda_C;
	fcc.lamda_T = ts.fc.lamda_T;
	fcc.lamda_L = ts.fc.lamda_L;

	uint32_t C_num = 0, L_num = 0, T_num = 0;
	fcc.V_types.resize(new_Vsize);
	fcc.RV_type.resize(new_Vsize);
	fcc.V_ids.resize(new_Vsize);

	for (uint32_t i = 0; i < ts.fc.V_types.size(); i++) if (V_flag[i]) {
		uint32_t mvid = new_V_map[i];
		fcc.V_types[mvid] = ts.fc.V_types[i];
		fcc.V_ids[mvid] = ts.fc.V_ids[i];
		fcc.RV_type[mvid] = ts.fc.RV_type[i];
		if (ts.fc.V_types[i] == Feature_V_Type::CORNER) C_num++;
		else if (ts.fc.V_types[i] == Feature_V_Type::LINE) L_num++;
		else if (ts.fc.V_types[i] == Feature_V_Type::REGULAR) T_num++;
	}

	if (ts.s.size()){
		uint32_t S_num = C_num + L_num + T_num;
		MatrixXd sc(S_num, 3);
		VectorXi s(S_num); S_num = 0;
		for (uint32_t i = 0; i < ts.s.size(); i++) if (V_flag[ts.s[i]]) {
			uint32_t mvid = new_V_map[ts.s[i]];
			s[S_num] = mvid;
			sc.row(S_num++) = ts.sc.row(i);
		}
		ts.s = s; ts.sc = sc;
	}

	fcc.ids_C.resize(C_num); fcc.C.resize(C_num, 3);
	fcc.num_a = L_num; fcc.ids_L.resize(L_num); fcc.Axa_L.resize(L_num, 3); fcc.origin_L.resize(L_num, 3);
	fcc.ids_T.resize(T_num); fcc.normal_T.resize(T_num, 3); fcc.dis_T.resize(T_num); fcc.V_T.resize(T_num, 3);
	C_num = L_num = T_num = 0;
	for (uint32_t i = 0; i < ts.fc.ids_C.size(); i++) {
		uint32_t vid = ts.fc.ids_C[i];
		if (V_flag[vid]) {
			fcc.ids_C[C_num] = new_V_map[vid];
			fcc.C.row(C_num++) = ts.fc.C.row(i);
		}
	}
	for (uint32_t i = 0; i < ts.fc.ids_L.size(); i++) {
		uint32_t vid = ts.fc.ids_L[i];
		if (V_flag[vid]) {
			fcc.ids_L[L_num] = new_V_map[vid];
			fcc.Axa_L.row(L_num) = ts.fc.Axa_L.row(i);
			fcc.origin_L.row(L_num++) = ts.fc.origin_L.row(i);
		}
	}
	for (uint32_t i = 0; i < ts.fc.ids_T.size(); i++) {
		uint32_t vid = ts.fc.ids_T[i];
		if (V_flag[vid]) {
			fcc.ids_T[T_num] = new_V_map[vid];
			fcc.normal_T.row(T_num) = ts.fc.normal_T.row(i);
			fcc.V_T.row(T_num) = ts.fc.V_T.row(i);
			fcc.dis_T.row(T_num++) = ts.fc.dis_T.row(i);
		}
	}
}

void simplification::slim_opt(Tetralize_Set &ts, const uint32_t iter) {
	igl::SLIMData sData;

	int vN = 0;
	vector<int> mapV(ts.V.rows(), -1), mapRV;
	vector<bool> touchedV_flag(ts.V.rows(), false);
	for (uint32_t i = 0; i < ts.T.rows(); i++)for (uint32_t j = 0; j < 4; j++)touchedV_flag[ts.T(i, j)] = true;
	for (uint32_t i = 0; i < ts.b.rows(); i++) touchedV_flag[ts.b[i]] = true;
	for (uint32_t i = 0; i < touchedV_flag.size(); i++) if (touchedV_flag[i]) { mapV[i] = vN++; mapRV.push_back(i); }

	Tetralize_Set ts_temp;
	localize_ts(ts, ts_temp, vN, mapV, touchedV_flag);

	sData.soft_const_p = 1e5;
	sData.exp_factor = 5.0;
	sData.lamda_C = ts_temp.fc.lamda_C;
	sData.lamda_T = ts_temp.fc.lamda_T;
	sData.lamda_L = ts_temp.fc.lamda_L;
	sData.lamda_region = ts_temp.lamda_region;
	sData.Vgroups = ts_temp.Vgroups;

	igl::SLIMData::SLIM_ENERGY energy = igl::SLIMData::SYMMETRIC_DIRICHLET;

	if (ts_temp.projection)
		slim_precompute(ts_temp.V, ts_temp.T, ts_temp.V, sData, energy, ts_temp.s, ts_temp.sc,
			ts_temp.fc.ids_C, ts_temp.fc.C,
			ts_temp.fc.ids_L, ts_temp.fc.Axa_L, ts_temp.fc.origin_L,
			ts_temp.fc.ids_T, ts_temp.fc.normal_T, ts_temp.fc.dis_T, ts_temp.regionb, ts_temp.regionbc, ts_temp.projection, ts_temp.global, ts.RT);
	else
		slim_precompute(ts_temp.V, ts_temp.T, ts_temp.V, sData, energy, ts_temp.b, ts_temp.bc,
			ts_temp.fc.ids_C, ts_temp.fc.C,
			ts_temp.fc.ids_L, ts_temp.fc.Axa_L, ts_temp.fc.origin_L,
			ts_temp.fc.ids_T, ts_temp.fc.normal_T, ts_temp.fc.dis_T, ts_temp.regionb, ts_temp.regionbc, ts_temp.projection, ts_temp.global, ts.RT);

	slim_solve(sData, iter);

	for (uint32_t i = 0; i < sData.V_o.rows(); i++) ts.V.row(mapRV[i]) = sData.V_o.row(i);
}
void simplification::localize_ts(Tetralize_Set &ts, Tetralize_Set &ts_temp, int &vN, vector<int> & mapV, vector<bool> & touchedV_flag) {
	ts_temp.projection = ts.projection;
	ts_temp.global = ts.global;

	ts_temp.lamda_region = ts.lamda_region;

	//V, T
	ts_temp.V.resize(vN, 3); ts_temp.T = ts.T;
	for (uint32_t i = 0; i < touchedV_flag.size(); i++) if (touchedV_flag[i]) {
		ts_temp.V.row(mapV[i]) = ts.V.row(i);
	}
	for (uint32_t i = 0; i < ts.T.rows(); i++)for (uint32_t j = 0; j < 4; j++) ts_temp.T(i, j) = mapV[ts.T(i, j)];
	//b, bc
	ts_temp.b = ts.b; ts_temp.bc = ts.bc; vector<int> bb;
	for (uint32_t i = 0; i < ts_temp.b.size(); i++) {
		ts_temp.b[i] = mapV[ts_temp.b[i]];
		bb.push_back(ts_temp.b[i]);
	}
	//Vgroups
	ts_temp.Vgroups = ts.Vgroups;
	for (uint32_t i = 0; i < ts_temp.Vgroups.size(); i++)for (uint32_t j = 0; j < ts_temp.Vgroups[i].size(); j++) ts_temp.Vgroups[i][j] = mapV[ts_temp.Vgroups[i][j]];
	//regionb, regionbc
	ts_temp.regionb = ts.regionb; ts_temp.regionbc = ts.regionbc;
	for (uint32_t i = 0; i < ts_temp.regionb.size(); i++) ts_temp.regionb[i] = mapV[ts_temp.regionb[i]];
	//s, sc
	int sN = 0;
	for (uint32_t i = 0; i < ts.s.size(); i++) if (touchedV_flag[ts.s[i]]) sN++;
	ts_temp.s.resize(sN); ts_temp.sc.resize(sN, 3);
	sN = 0;
	for (uint32_t i = 0; i < ts.s.size(); i++)if (touchedV_flag[ts.s[i]]) {
		ts_temp.s[sN] = mapV[ts.s[i]];
		ts_temp.sc.row(sN) = ts.sc.row(i);
		sN++;  
	}

	//fc
	ts_temp.fc.lamda_C = ts.fc.lamda_C;
	ts_temp.fc.lamda_T = ts.fc.lamda_T;
	ts_temp.fc.lamda_L = ts.fc.lamda_L;

	uint32_t C_num = 0, L_num = 0, T_num = 0;

	for (uint32_t i = 0; i < ts.fc.V_types.size(); i++) if (touchedV_flag[i]) {
		if (ts.fc.V_types[i] == Feature_V_Type::CORNER) C_num++;
		else if (ts.fc.V_types[i] == Feature_V_Type::LINE) L_num++;
		else if (ts.fc.V_types[i] == Feature_V_Type::REGULAR) T_num++;
	}
	ts_temp.fc.ids_C.resize(C_num); ts_temp.fc.C.resize(C_num, 3);
	ts_temp.fc.num_a = L_num; ts_temp.fc.ids_L.resize(L_num); ts_temp.fc.Axa_L.resize(L_num, 3); ts_temp.fc.origin_L.resize(L_num, 3);
	ts_temp.fc.ids_T.resize(T_num); ts_temp.fc.normal_T.resize(T_num, 3); ts_temp.fc.dis_T.resize(T_num);
	C_num = L_num = T_num = 0;
	for (uint32_t i = 0; i < ts.fc.ids_C.size(); i++) {
		uint32_t vid = ts.fc.ids_C[i];
		if (touchedV_flag[vid]) {
			ts_temp.fc.ids_C[C_num] = mapV[vid];
			ts_temp.fc.C.row(C_num++) = ts.fc.C.row(i);
		}
	}
	for (uint32_t i = 0; i < ts.fc.ids_L.size(); i++) {
		uint32_t vid = ts.fc.ids_L[i];
		if (touchedV_flag[vid]) {
			ts_temp.fc.ids_L[L_num] = mapV[vid];
			ts_temp.fc.Axa_L.row(L_num) = ts.fc.Axa_L.row(i);
			ts_temp.fc.origin_L.row(L_num++) = ts.fc.origin_L.row(i);
		}
	}
	for (uint32_t i = 0; i < ts.fc.ids_T.size(); i++) {
		uint32_t vid = ts.fc.ids_T[i];
		if (touchedV_flag[vid]) {
			ts_temp.fc.ids_T[T_num] = mapV[vid];
			ts_temp.fc.normal_T.row(T_num) = ts.fc.normal_T.row(i);
			ts_temp.fc.dis_T.row(T_num++) = ts.fc.dis_T.row(i);
		}
	}
}

void simplification::subdivision() {
	if (Hex_Num_Threshold >= mesh.Hs.size()) {
		if (hex_mesh_subdivision()) {
			cout << "subdivided the hex-mesh " << endl;
		}

		extract();
		ranking();
	}
}
bool simplification::hex_mesh_subdivision() {
	if (Hex_Num_Threshold <= mesh.Hs.size()) return false;

	std::function<bool(pair<Float, uint32_t> &, pair<Float, uint32_t> &)> greater = [&](
		pair<Float, uint32_t> &a, pair<Float, uint32_t> &b)->bool {
		return a.first > b.first;
	};

	vector<pair<Float, uint32_t>> sids;
	for (auto wi : All_Sheets) sids.push_back(make_pair(wi.weight_len/frame.FEs[wi.middle_es[0]].es_link.size(), wi.id));
	std::sort(sids.begin(), sids.end(), greater);
	//subset of sheets
	vector<uint32_t> ids_use;
	uint32_t num = 0;
	for (auto si : sids) {
		uint32_t cur_n = 0;
		for (auto cid : All_Sheets[si.second].cs) cur_n += frame.FHs[cid].hs_net.size();
		num += cur_n;
		ids_use.push_back(si.second); 
		break;
	}
	//tag edges
	vector<bool> E_flag(mesh.Es.size(), false), H_flag(mesh.Hs.size(), false);
	for (auto id : ids_use) {
		vector<vector<uint32_t>> candiate_es_links, v_group;
		if (!vs_pair_sheet(id, candiate_es_links, v_group)) { cout << "ERROR: no vs pairs" << endl;  system("PAUSE"); }
		for (auto link : candiate_es_links) for (auto eid : link) E_flag[eid] = true;
		for (auto cid : All_Sheets[id].cs) for (auto hid : frame.FHs[cid].hs_net) H_flag[hid] = true;
	}

	//subset of meshes, subdivide
	mesh_.type = mesh.type;
	mesh_.Vs.clear();
	mesh_.Vs.resize(mesh.Vs.size());
	for (auto v : mesh.Vs) {
		Hybrid_V hv;
		hv.id = v.id;
		hv.v = v.v;
		mesh_.Vs[hv.id] = hv;
	}
	uint32_t V_id = mesh.Vs.size();

	const uint32_t INVALID_ELE = mesh.Vs.size() + mesh.Es.size() + mesh.Fs.size() + mesh.Hs.size();
	vector<uint32_t> E_map(mesh.Es.size(), INVALID_ELE), F_map(mesh.Fs.size(), INVALID_ELE), H_map(mesh.Hs.size(), INVALID_ELE);
	vector<Vector3d> V_; 
	//edge v
	for (uint32_t i = 0; i <mesh.Es.size(); i++) {
		if (!E_flag[i]) continue;

		uint32_t num = 2; Vector3d vd; vd.setZero();
		for (uint32_t j = 0; j < num; j++)
			vd += mesh.V.col(mesh.Es[i].vs[j]);
		vd /= num;
		V_.push_back(vd);
		Hybrid_V hv;
		hv.id = V_id++;
		mesh_.Vs.push_back(hv);

		E_map[i] = hv.id;
	}
	//face v
	for (uint32_t i = 0; i < mesh.Fs.size(); i++) {
		bool valid_f = true;
		for (auto eid : mesh.Fs[i].es) if (!E_flag[eid]) {valid_f = false; break;}
		if (!valid_f) continue;

		uint32_t num = 4; Vector3d vd; vd.setZero();
		for (uint32_t j = 0; j < num; j++)
			vd += mesh.V.col(mesh.Fs[i].vs[j]);
		vd /= num;
		V_.push_back(vd);
		Hybrid_V hv;
		hv.id = V_id++;
		mesh_.Vs.push_back(hv);

		F_map[i] = hv.id;
	}
	//hex v
	for (uint32_t i = 0; i < mesh.Hs.size(); i++) {
		vector<uint32_t> &es = mesh.Hs[i].es;
		if (!es.size()) {
			//std::cout << "Dont have es" << endl;
			for (auto fid : mesh.Hs[i].fs) es.insert(es.end(), mesh.Fs[fid].es.begin(), mesh.Fs[fid].es.end());
			std::sort(es.begin(), es.end()); es.erase(unique(es.begin(), es.end()), es.end());
		}
		//else std::cout << "have es" << endl;
		bool valid_h = true;
		for (auto eid : es) if (!E_flag[eid]) { valid_h = false; break; }
		if (!valid_h) continue;

		uint32_t num = 8; Vector3d vd; vd.setZero();
		for (uint32_t j = 0; j < num; j++)
			vd += mesh.V.col(mesh.Hs[i].vs[j]);
		vd /= num;
		V_.push_back(vd);
		Hybrid_V hv;
		hv.id = V_id++;
		mesh_.Vs.push_back(hv);

		H_map[i] = hv.id;
	}
	//Hs
	mesh_.Hs.clear(); CI.hs.clear();
	uint32_t H_id = 0;
	for (uint32_t i = 0; i < mesh.Hs.size(); i++) {
		if (!H_flag[i]) { 
			Hybrid h;
			h.id = H_id++;
			h.vs = mesh.Hs[i].vs;
			mesh_.Hs.push_back(h); 
			continue; 
		}
		//split three direction
		if (H_map[i] < INVALID_ELE) {

			const vector<uint32_t> &vs = mesh.Hs[i].vs;
			for (auto vid:vs) {
				std::sort(mesh.Vs[vid].neighbor_es.begin(), mesh.Vs[vid].neighbor_es.end());
				std::sort(mesh.Vs[vid].neighbor_fs.begin(), mesh.Vs[vid].neighbor_fs.end());
			}
			for (uint32_t j = 0; j < 8; j++) {
				Hybrid h;
				h.id = H_id++;
				h.vs.resize(8); h.vs[j] = vs[j];
				for (uint32_t k = 0; k < 8; k++) {
					if (k == j) continue;
					//sharedes
					vector<uint32_t> sharedes;
					std::set_intersection(mesh.Vs[vs[k]].neighbor_es.begin(), mesh.Vs[vs[k]].neighbor_es.end(),
						mesh.Vs[vs[j]].neighbor_es.begin(), mesh.Vs[vs[j]].neighbor_es.end(), back_inserter(sharedes));
					if (sharedes.size()) {
						if (sharedes.size() != 1 || E_map[sharedes[0]] == INVALID_ELE) { system("PAUSE"); }
						h.vs[k] = E_map[sharedes[0]];
						continue;
					}
					//sharedfs
					vector<uint32_t> sharedfs;
					std::set_intersection(mesh.Vs[vs[k]].neighbor_fs.begin(), mesh.Vs[vs[k]].neighbor_fs.end(),
						mesh.Vs[vs[j]].neighbor_fs.begin(), mesh.Vs[vs[j]].neighbor_fs.end(), back_inserter(sharedfs));
					if (sharedfs.size()) {
						if (sharedfs.size() != 1 || F_map[sharedfs[0]] == INVALID_ELE) { system("PAUSE"); }
						h.vs[k] = F_map[sharedfs[0]];
						continue;
					}
					//hex-diagonal v
					h.vs[k] = H_map[i];
				}
				mesh_.Hs.push_back(h);
				CI.hs.push_back(h.id);
			}
			continue;
		}
		//split two direction
		vector<uint32_t> fs;
		for (auto fid : mesh.Hs[i].fs)if (F_map[fid] < INVALID_ELE) fs.push_back(fid);
		
		if (fs.size()) {
			if(fs.size() !=2) {
			system("PAUSE"); }

			const vector<uint32_t> &hvs = mesh.Hs[i].vs;
			vector<uint32_t> &fvs = mesh.Fs[fs[0]].vs;

			for (auto vid : hvs) std::sort(mesh.Vs[vid].neighbor_es.begin(), mesh.Vs[vid].neighbor_es.end());
			
			for (uint32_t j = 0; j < 4; j++) {
				uint32_t v0_fix = fvs[j], v1_fix = INVALID_ELE;
				for (auto vid : mesh.Vs[v0_fix].neighbor_vs)
					if (std::find(mesh.Fs[fs[1]].vs.begin(), mesh.Fs[fs[1]].vs.end(), vid) != mesh.Fs[fs[1]].vs.end()) {
						v1_fix = vid; break;
					}
				if(v1_fix == INVALID_ELE) {system("PAUSE"); }

				Hybrid h;
				h.id = H_id++;
				h.vs.resize(8); 
				for (uint32_t k = 0; k < 8; k++) {
					if (hvs[k] == v0_fix) { h.vs[k] = v0_fix; continue; }
					else if(hvs[k] == v1_fix) { h.vs[k] = v1_fix; continue; }

					if (std::find(mesh.Fs[fs[0]].vs.begin(), mesh.Fs[fs[0]].vs.end(), hvs[k]) != mesh.Fs[fs[0]].vs.end()) {
						//sharedes
						vector<uint32_t> sharedes;
						std::set_intersection(mesh.Vs[v0_fix].neighbor_es.begin(), mesh.Vs[v0_fix].neighbor_es.end(),
							mesh.Vs[hvs[k]].neighbor_es.begin(), mesh.Vs[hvs[k]].neighbor_es.end(), back_inserter(sharedes));
						if (sharedes.size()) {
							if (sharedes.size() != 1 || E_map[sharedes[0]] == INVALID_ELE) {system("PAUSE"); }
							h.vs[k] = E_map[sharedes[0]];
							continue;
						}else h.vs[k] = F_map[fs[0]];

					}else if(std::find(mesh.Fs[fs[1]].vs.begin(), mesh.Fs[fs[1]].vs.end(), hvs[k]) != mesh.Fs[fs[1]].vs.end()) {
						//sharedes
						vector<uint32_t> sharedes;
						std::set_intersection(mesh.Vs[v1_fix].neighbor_es.begin(), mesh.Vs[v1_fix].neighbor_es.end(),
							mesh.Vs[hvs[k]].neighbor_es.begin(), mesh.Vs[hvs[k]].neighbor_es.end(), back_inserter(sharedes));
						if (sharedes.size()) {
							if (sharedes.size() != 1 || E_map[sharedes[0]] == INVALID_ELE) {system("PAUSE"); }
							h.vs[k] = E_map[sharedes[0]];
							continue;
						}
						else h.vs[k] = F_map[fs[1]];
					}else { system("PAUSE"); }
				}
				mesh_.Hs.push_back(h);
				CI.hs.push_back(h.id);
			}
			continue;
		}
		//split one direction
		fs.clear();
		for (auto fid : mesh.Hs[i].fs) {
			bool valid_f = true;
			for (auto eid : mesh.Fs[fid].es) if (E_flag[eid]) { valid_f = false; break;}
			if (valid_f)fs.push_back(fid);
		}
		if(fs.size() != 2) {system("PAUSE"); }
		else {
			const vector<uint32_t> &hvs = mesh.Hs[i].vs;

			for (auto vid : hvs) std::sort(mesh.Vs[vid].neighbor_es.begin(), mesh.Vs[vid].neighbor_es.end());

			for (uint32_t j = 0; j < 2; j++) {
				vector<uint32_t> &fvs = mesh.Fs[fs[j]].vs;

				Hybrid h;
				h.id = H_id++;
				h.vs.resize(8);
				for (uint32_t k = 0; k < 8; k++) {
					if (std::find(fvs.begin(), fvs.end(), hvs[k]) != fvs.end()) { h.vs[k] = hvs[k]; continue; }

					for (auto vid : fvs) {
						vector<uint32_t> sharedes;
						std::set_intersection(mesh.Vs[vid].neighbor_es.begin(), mesh.Vs[vid].neighbor_es.end(),
							mesh.Vs[hvs[k]].neighbor_es.begin(), mesh.Vs[hvs[k]].neighbor_es.end(), back_inserter(sharedes));
						if (sharedes.size()) {
							if (sharedes.size() != 1 || E_map[sharedes[0]] == INVALID_ELE) { system("PAUSE"); }
							h.vs[k] = E_map[sharedes[0]];
							break;
						}
					}
				}
				mesh_.Hs.push_back(h);
				CI.hs.push_back(h.id);
			}
		}
	}
	mesh_.V.resize(3, mesh_.Vs.size()); mesh_.V.setZero();
	mesh_.V.block(0, 0, 3, mesh.V.cols()) = mesh.V;
	for (uint32_t i = 0; i < V_.size(); i++) mesh_.V.col(mesh.V.cols() + i) = V_[i];
	for (auto h:mesh_.Hs) for (uint32_t k = 0; k < 8; k++)mesh_.Vs[h.vs[k]].neighbor_hs.push_back(h.id);
	
	Mesh mesh_temp = mesh_;
	build_connectivity(mesh_temp);
	Mesh_Quality mq;
	scaled_jacobian(mesh_temp, mq);

	if (mq.min_Jacobian < Precision_Pro) {
		return false;
	}

	//feature
	Feature_Constraints fc_;
	fc_.lamda_C = fc.lamda_C;
	fc_.lamda_T = fc.lamda_T;
	fc_.lamda_L = fc.lamda_L;

	uint32_t new_Vsize = mesh_.V.cols(), C_num = 0, L_num = 0, T_num = 0;
	fc_.V_types.resize(new_Vsize);
	fc_.RV_type.resize(new_Vsize);
	fc_.V_ids.resize(new_Vsize);

	for (uint32_t i = 0; i < fc.V_types.size(); i++) {
		fc_.V_types[i] = fc.V_types[i];
		fc_.V_ids[i] = fc.V_ids[i];
		fc_.RV_type[i] = fc.RV_type[i];
		if (fc.V_types[i] == Feature_V_Type::CORNER) C_num++;
		else if (fc.V_types[i] == Feature_V_Type::LINE) L_num++;
		else if (fc.V_types[i] == Feature_V_Type::REGULAR) T_num++;
	}
	fill(fc_.RV_type.begin() + fc.RV_type.size(), fc_.RV_type.end(), true);

	vector<Vector3d> normal_Ts, V_Ts; vector<double> dis_Ts;
	for (uint32_t i = 0; i < fc.normal_T.rows(); i++) {
		normal_Ts.push_back(fc.normal_T.row(i));
		V_Ts.push_back(fc.V_T.row(i));
		dis_Ts.push_back(fc.dis_T[i]);
	}

	vector<bool> Line_flag(mf.curve_vs.size(), false);
	for (uint32_t i = 0; i < E_map.size();i++) if (E_map[i] != INVALID_ELE) {
		
		vector<uint32_t> &vs = mesh.Es[i].vs, ts, ls;
		if (!mesh.Es[i].boundary) {
			fc_.V_types[E_map[i]] = Feature_V_Type::INTERIOR; 
			continue;
		}

		for (auto vid:vs) {
			if (fc.V_types[vid] == Feature_V_Type::CORNER) { 
				int pos = find(mf.corners.begin(), mf.corners.end(), fc.V_ids[vid]) - mf.corners.begin();
				ls.insert(ls.end(),mf.corner_curves[pos].begin(), mf.corner_curves[pos].end());
				vector<uint32_t> &fs = mf.tri.Vs[fc.V_ids[vid]].neighbor_fs;
				ts.insert(ts.end(), fs.begin(), fs.end());
			}
			else if (fc.V_types[vid] == Feature_V_Type::LINE) {
				ls.push_back(fc.V_ids[vid]);
				for (auto vid_ : mf.curve_vs[fc.V_ids[vid]]) {
					vector<uint32_t> &fs = mf.tri.Vs[vid_].neighbor_fs;
					ts.insert(ts.end(), fs.begin(), fs.end());
				}
			}
			else if (fc.V_types[vid] == Feature_V_Type::REGULAR) {
				if(fc.RV_type[vid]) ts.push_back(fc.V_ids[vid]);
				else {
					vector<uint32_t> &fs = mf.tri.Vs[fc.V_ids[vid]].neighbor_fs;
					ts.insert(ts.end(), fs.begin(), fs.end());
				}
			}
		}
		//ls
		if (ls.size()) {
			bool found = false; uint32_t line_id = INVALID_ELE;
			for (auto lid : ls) {
				if (Line_flag[lid]) { Line_flag[lid] = false; found = true; line_id = lid; break; }
				else Line_flag[lid] = true;
			}
			for (auto lid : ls) Line_flag[lid] = false;
			if (found) {
				fc_.V_types[E_map[i]] = Feature_V_Type::LINE;
				fc_.V_ids[E_map[i]] = line_id;
				L_num++;
				continue;
			}
		}
		//ts
		vector<uint32_t> tids_;

		for (uint32_t Iter = 0; Iter < subdivision_project_range; Iter++) {
			for (uint32_t j = 0; j < ts.size(); j++) {
				vector<uint32_t> &vs_ = mf.tri.Fs[ts[j]].vs;
				for (uint32_t k = 0; k < 3; k++) {
					for (auto ntid : mf.tri.Vs[vs_[k]].neighbor_fs) {
						tids_.push_back(ntid);
					}
				}
			}
			ts.insert(ts.end(), tids_.begin(), tids_.end()); tids_.clear();
			sort(ts.begin(), ts.end()); ts.erase(unique(ts.begin(), ts.end()), ts.end());
		}
		Vector3d v = mesh_.V.col(E_map[i]), interpolP, interpolN; double dis = 0;
		uint32_t tid = -1;
		Vector3d PreinterpolP, PreinterpolN;
		PreinterpolP.setZero(); PreinterpolN.setZero();
		if (phong_projection(ts, subdivision_project_range, tid, v, interpolP, interpolN, PreinterpolP, PreinterpolN)) {
			if ((v - interpolP).norm() >= mf.ave_length) {
				tid = 0;
				double min_dis = (v - mf.Tcenters[0]).norm();
				for (uint32_t j = 1; j < mf.Tcenters.size(); j++) {
					if (min_dis >(v - mf.Tcenters[j]).norm()) {
						tid = j;
						min_dis = (v - mf.Tcenters[j]).norm();
					}
				}
				uint32_t tid_temp = tid;
				ts.clear();
				ts.push_back(tid_temp);
				if (phong_projection(ts, subdivision_project_range, tid_temp, v, interpolP, interpolN, PreinterpolP, PreinterpolN))
					tid = tid_temp;
				else {
					interpolP = mf.Tcenters[tid_temp];
					interpolN = mf.normal_Tri.col(tid_temp);
				}
			}
		}
		fc_.V_types[E_map[i]] = Feature_V_Type::REGULAR;
		fc_.V_ids[E_map[i]] = tid;
		normal_Ts.push_back(interpolN);
		V_Ts.push_back(interpolP);
		dis_Ts.push_back(interpolN.dot(interpolP));
		T_num++;
	}
	for (uint32_t i = 0; i < F_map.size(); i++) if (F_map[i] != INVALID_ELE) {
		
		vector<uint32_t> &vs = mesh.Fs[i].vs, ts;
		if (!mesh.Fs[i].boundary) {
			fc_.V_types[F_map[i]] = Feature_V_Type::INTERIOR; 
			continue;
		}

		for (auto vid : vs) {
			if (fc.V_types[vid] == Feature_V_Type::CORNER) {
				vector<uint32_t> &fs = mf.tri.Vs[fc.V_ids[vid]].neighbor_fs;
				ts.insert(ts.end(), fs.begin(), fs.end());
			}
			else if (fc.V_types[vid] == Feature_V_Type::LINE) {
				for (auto vid_ : mf.curve_vs[fc.V_ids[vid]]) {
					vector<uint32_t> &fs = mf.tri.Vs[vid_].neighbor_fs;
					ts.insert(ts.end(), fs.begin(), fs.end());
				}
			}
			else if (fc.V_types[vid] == Feature_V_Type::REGULAR) {

				if (fc.RV_type[vid]) ts.push_back(fc.V_ids[vid]);
				else {
					vector<uint32_t> &fs = mf.tri.Vs[fc.V_ids[vid]].neighbor_fs;
					ts.insert(ts.end(), fs.begin(), fs.end());
				}
			}
		}
		//ts
		vector<uint32_t> tids_;
		for (uint32_t Iter = 0; Iter < subdivision_project_range; Iter++) {
			for (uint32_t j = 0; j < ts.size(); j++) {
				vector<uint32_t> &vs_ = mf.tri.Fs[ts[j]].vs;
				for (uint32_t k = 0; k < 3; k++) {
					for (auto ntid : mf.tri.Vs[vs_[k]].neighbor_fs) {
						tids_.push_back(ntid);
					}
				}
			}
			ts.insert(ts.end(), tids_.begin(), tids_.end()); tids_.clear();
			sort(ts.begin(), ts.end()); ts.erase(unique(ts.begin(), ts.end()), ts.end());
		}
		Vector3d v = mesh_.V.col(F_map[i]), interpolP, interpolN; double dis = 0;
		uint32_t tid = -1;
		Vector3d PreinterpolP, PreinterpolN;
		PreinterpolP.setZero(); PreinterpolN.setZero();
		if (phong_projection(ts, subdivision_project_range, tid, v, interpolP, interpolN, PreinterpolP, PreinterpolN)) {
			if ((v - interpolP).norm() >= mf.ave_length) {
				tid = 0;
				double min_dis = (v - mf.Tcenters[0]).norm();
				for (uint32_t j = 1; j < mf.Tcenters.size(); j++) {
					if (min_dis >(v - mf.Tcenters[j]).norm()) {
						tid = j;
						min_dis = (v - mf.Tcenters[j]).norm();
					}
				}
				uint32_t tid_temp = tid;
				ts.clear();
				ts.push_back(tid_temp);
				if (phong_projection(ts, subdivision_project_range, tid_temp, v, interpolP, interpolN, PreinterpolP, PreinterpolN))
					tid = tid_temp;
				else {
					interpolP = mf.Tcenters[tid_temp];
					interpolN = mf.normal_Tri.col(tid_temp);
				}
			}
		}
		fc_.V_types[F_map[i]] = Feature_V_Type::REGULAR;
		fc_.V_ids[F_map[i]] = tid;
		normal_Ts.push_back(interpolN);
		V_Ts.push_back(interpolP);
		dis_Ts.push_back(interpolN.dot(interpolP));
		T_num++;
	}
	for (uint32_t i = 0; i < H_map.size(); i++) if (H_map[i] != INVALID_ELE) {
		fc_.V_types[H_map[i]] = Feature_V_Type::INTERIOR;
	}
	//corner
	fc_.ids_C = fc.ids_C;
	fc_.C = fc.C;
	//
	fc_.num_a = L_num; fc_.ids_L.resize(L_num); fc_.Axa_L.resize(L_num, 3); fc_.origin_L.resize(L_num, 3);
	fc_.ids_T.resize(T_num); fc_.normal_T.resize(T_num, 3); fc_.dis_T.resize(T_num); fc_.V_T.resize(T_num, 3);
	L_num = T_num = 0;

	for (uint32_t i = 0; i < fc_.V_types.size(); i++) {
		if(fc_.V_types[i] == Feature_V_Type::LINE) fc_.ids_L[L_num++] = i;
		else if(fc_.V_types[i] == Feature_V_Type::REGULAR) fc_.ids_T[T_num++] = i;
	}
	for (uint32_t i = 0; i < fc_.normal_T.rows(); i++) {
		fc_.normal_T.row(i) = normal_Ts[i];
		fc_.V_T.row(i) = V_Ts[i];
		fc_.dis_T[i] = dis_Ts[i];
	}

	Feature_Constraints fc_temp = fc_;
	MatrixXd v_temp = mesh_.V.transpose(), sc_temp;
	VectorXi s_temp;
	project_surface_update_feature(mf, fc_temp, v_temp, ts.s, ts.sc, subdivision_project_range);	
///////////////////////////////////optimize////////////////////////////
	tetralize_mesh_submesh(ts, mesh_temp);
	ts.fc = fc_temp;
	ts.global = true;

	for (uint32_t i = 0; i < Slim_Iteration; i++) {
		ts.projection = false;

		//RT
		compute_referenceMesh(ts.V, mesh_temp.Hs, CI.Hsregion, ts.RT);

		slim_opt(ts, 1);

		project_surface_update_feature(mf, ts.fc, ts.V, ts.s, ts.sc, Projection_range);
		ts.projection = true;
	}

	vector<bool> touchedV_flag(mesh_temp.Vs.size(), false);
	for (uint32_t i = 0; i < ts.T.rows(); i++)for (uint32_t j = 0; j < 4; j++)touchedV_flag[ts.T(i, j)] = true;
	for (uint32_t i = 0; i < touchedV_flag.size(); i++) if (touchedV_flag[i]) {
		mesh_temp.V.col(i) = ts.V.row(i);
	}

	scaled_jacobian(mesh_temp, mq);
	if (mq.min_Jacobian < Jacobian_Bound) {
		return false; 
	}

	if (!hausdorff_ratio_check(mf.tri, mesh_temp)) return false;
	mesh_.V = mesh_temp.V;
	fc = ts.fc;
	mesh = mesh_;
	build_connectivity(mesh);
	base_com.singularity_structure(si, mesh);
	base_com.base_complex_extraction(si, frame, mesh);

	return true;
}
uint32_t simplification::nearest_tid(vector<uint32_t> &ts, const Vector3d &v, Vector3d &n, Vector3d &pv, double &dis) {
	uint32_t tid=-1;
	if (!ts.size()) return tid;

	Vector3d interpolP, interpolN;
	vector<uint32_t>  tids;
	vector<pair<double, uint32_t>> dis_ids;
	vector<Vector3d> pvs, pns;
	bool found = false;

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
			dis_ids.push_back(make_pair((v - interpolP).norm(), tids.size()));
			tids.push_back(ts[j]);
			pvs.push_back(interpolP); pns.push_back(interpolN);
		}
	}
	sort(dis_ids.begin(), dis_ids.end());
	if (dis_ids.size()) {
		found = true;
		tid = tids[dis_ids[0].second];
		n = pns[dis_ids[0].second];
		pv = pvs[dis_ids[0].second];
		dis = dis_ids[0].first;
	}else {
		for (uint32_t j = 0; j < ts.size(); j++) {
			vector<uint32_t> &vs = mf.tri.Fs[ts[j]].vs;
			Vector3d cv; cv.setZero();
			double dis = 0;
			for (uint32_t k = 0; k < 3; k++) cv += mf.tri.V.col(vs[k]);
			cv /= 3;
			dis = (cv - v).norm();
			dis_ids.push_back(make_pair(dis, j));
			pvs.push_back(cv); pns.push_back(mf.normal_Tri.col(ts[j]));
		}
		sort(dis_ids.begin(), dis_ids.end());
		tid = ts[dis_ids[0].second];
		n = pns[dis_ids[0].second];
		pv = pvs[dis_ids[0].second];
		dis = dis_ids[0].first;
	}
	return tid;
}
bool simplification::hausdorff_ratio_check(Mesh &m0, Mesh &m1) {

	std::function<void(Mesh &, Mesh &, vector<bool> &, int &) > re_indexing = [&](Mesh &M, Mesh &m, vector<bool> &V_flag, int & N)->void {
			m.V.resize(3, N); N = 0; vector<int> v_map(M.Vs.size(), 0);
			for (uint32_t i = 0; i < V_flag.size(); i++) if (V_flag[i]) { m.V.col(N) = M.V.col(i); v_map[i] = N++; }
			for (auto f : M.Fs) {
				if (!f.boundary) continue;
				vector<uint32_t> bvs; for (auto vid : f.vs)if (V_flag[vid]) bvs.push_back(v_map[vid]);
				if (bvs.size() == 3) {
					Hybrid_F hf; hf.vs = bvs;
					m.Fs.push_back(hf);
				}
				else if (bvs.size() == 4) {
					Hybrid_F hf;
					hf.vs.push_back(bvs[0]);
					hf.vs.push_back(bvs[1]);
					hf.vs.push_back(bvs[2]);
					m.Fs.push_back(hf);
					hf.vs.clear();
					hf.vs.push_back(bvs[2]);
					hf.vs.push_back(bvs[3]);
					hf.vs.push_back(bvs[0]);
					m.Fs.push_back(hf);
				}
			}
		};
		vector<bool> V_flag; int bvN = 0, bbvN = 0;

		Mesh Mglobal;
		V_flag.resize(m1.Vs.size()); std::fill(V_flag.begin(), V_flag.end(), false);
		bvN = 0;
		for (uint32_t i = 0; i < m1.Vs.size(); i++) if (m1.Vs[i].boundary) { V_flag[i] = true; bvN++; }
		re_indexing(m1, Mglobal, V_flag, bvN);

		if (!compute(m0, Mglobal, hausdorff_ratio, hausdorff_ratio_threshould)) {
			return false;
		}
		return true;
}